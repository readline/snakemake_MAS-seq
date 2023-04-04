from os.path import join
import os
import pandas as pd
from scripts.Load import samplesheet

snakedir = os.getcwd()
print(snakedir)
configfile: 'config.yaml'
print(config)
sampledic = samplesheet(config['samplesheet'])
workdir: config['workdir']
#sbcmd="sbatch --cpus-per-task={threads} --mem={cluster.mem}"
#sbcmd+=" --time={cluster.time} --partition={cluster.partition}"
#sbcmd+=" --out={cluster.out} {cluster.extra}"

rule all:
    input:
        lambda wildcards: ["09.seurat_prep/{}/genes_seurat/matrix.mtx".format(sample) for sample in sampledic],

rule prepare:
    output:
        genome_fa   = 'ref/genome.fasta',
        anno_gtf  = 'ref/gencode.gtf',
        polya = 'ref/polyA.list.txt',
        tssbed = 'ref/refTSS.bed',
        barcode = 'ref/barcodes.txt.gz',
        primer3 = 'ref/10x_3kit_primers.fasta',
        mas16 = 'ref/mas16.fasta'
    log:
        out = snakedir + "/logs/A0.Skera_split/prepare.o",
        err = snakedir + "/logs/A0.Skera_split/prepare.e",
    threads: 2
    resources: 
        mem = '4g',
        extra = ' --gres=lscratch:10 '
    shell:
        '''module load {config[modules][smrtanalysis]} {config[modules][sambamba]} > {log.out} 2>{log.err}
        skera --version  >>{log.out} 2>>{log.err}
        lima --version >>{log.out} 2>>{log.err}
        isoseq3 --version >>{log.out} 2>>{log.err}
        pbmm2 --version >>{log.out} 2>>{log.err}
        pigeon --version >>{log.out} 2>>{log.err}
        samtools --version >>{log.out} 2>>{log.err}
        sambamba --version >>{log.out} 2>>{log.err}
        
        wget {config[references][genome_fa]} -O {output.genome_fa} >>{log.out} 2>>{log.err}
        wget {config[references][genome_fa]}.fai -O {output.genome_fa}.fai >>{log.out} 2>>{log.err}
        wget {config[references][anno_gtf]} -O {output.anno_gtf} >>{log.out} 2>>{log.err}
        wget {config[references][anno_gtf]}.pgi -O {output.anno_gtf}.pgi >>{log.out} 2>>{log.err}
        wget {config[references][polya]} -O {output.polya} >>{log.out} 2>>{log.err}
        wget {config[references][tssbed]} -O {output.tssbed} >>{log.out} 2>>{log.err}
        wget {config[references][tssbed]}.pgi -O {output.tssbed}.pgi >>{log.out} 2>>{log.err}
        wget {config[references][barcode]} -O {output.barcode} >>{log.out} 2>>{log.err}
        wget {config[references][primer3]} -O {output.primer3} >>{log.out} 2>>{log.err}
        wget {config[references][mas16]} -O {output.mas16} >>{log.out} 2>>{log.err}
        wget {config[references][mas16]}.fai -O {output.mas16}.fai >>{log.out} 2>>{log.err}
        touch ref/ref.ready >>{log.out} 2>>{log.err}
        '''
    
rule skera_split:
    input:
        bam = lambda wildcards: sampledic[wildcards.sample],
        mas16 = 'ref/mas16.fasta',
    output:
        bam   = temp('01.skera_split/{sample}/{sample}.skera.bam'),
        bami  = temp('01.skera_split/{sample}/{sample}.skera.bam.pbi'),
        bam2  = temp('01.skera_split/{sample}/{sample}.skera.non_passing.bam'),
        bam2i = temp('01.skera_split/{sample}/{sample}.skera.non_passing.bam.pbi'),
    params:
        prefix = '01.skera_split/{sample}/{sample}.skera'
    log:
        out = snakedir + "/logs/A1.skera_split/{sample}.o",
        err = snakedir + "/logs/A1.skera_split/{sample}.e",
    threads: 16
    resources: 
        mem = '16g',
        extra = ' --gres=lscratch:10 '
    shell:
        '''module load {config[modules][smrtanalysis]} {config[modules][sambamba]} > {log.out} 2>{log.err}
    skera split -j {threads} {input.bam} {input.mas16} {output.bam} >>{log.out} 2>>{log.err}
    sambamba flagstat -t {threads} {output.bam} > {output.bam}.flagstat 2>> {log.err}
    sambamba flagstat -t {threads} {output.bam2} > {output.bam2}.flagstat 2>> {log.err}
    '''
           
rule lima_rmprimer:
    input:
        bam   = '01.skera_split/{sample}/{sample}.skera.bam',
        bami  = '01.skera_split/{sample}/{sample}.skera.bam.pbi',
        primer3 = 'ref/10x_3kit_primers.fasta',
    output:
        bam   = temp('02.lima_rmprimer/{sample}/{sample}.lima.5p--3p.bam'),
        bami  = temp('02.lima_rmprimer/{sample}/{sample}.lima.5p--3p.bam.pbi'),
    params:
        prefix = '02.lima_rmprimer/{sample}/{sample}.lima',
    log:
        out = snakedir + "/logs/A2.lima_rmprimer/{sample}.o",
        err = snakedir + "/logs/A2.lima_rmprimer/{sample}.e",
    threads: 16
    resources: 
        mem = '16g',
        extra = ' --gres=lscratch:10 '
    shell:
        '''module load {config[modules][smrtanalysis]} {config[modules][sambamba]} > {log.out} 2>{log.err}
    lima -j {threads} --per-read --isoseq {input.bam} {input.primer3} {params.prefix}.bam >>{log.out} 2>>{log.err}
    sambamba flagstat -t {threads} {output.bam} > {output.bam}.flagstat 2>>{log.err}
    '''

rule clip_umi_barcode:
    input:
        bam   = '02.lima_rmprimer/{sample}/{sample}.lima.5p--3p.bam',
        bami  = '02.lima_rmprimer/{sample}/{sample}.lima.5p--3p.bam.pbi',
    output:
        bam   = temp('03.clip_umi_barcode/{sample}/{sample}.clip_umi_barcode.flt.bam'),
        bami  = temp('03.clip_umi_barcode/{sample}/{sample}.clip_umi_barcode.flt.bam.pbi'),
    params:
        design = 'T-12U-16B',
    log:
        out = snakedir + "/logs/A3.clip_umi_barcode/{sample}.o",
        err = snakedir + "/logs/A3.clip_umi_barcode/{sample}.e",
    threads: 8
    resources: 
        mem = '16g',
        extra = ' --gres=lscratch:10 '
    shell:
        '''module load {config[modules][smrtanalysis]} {config[modules][sambamba]} > {log.out} 2>{log.err}
    isoseq3 tag {input.bam} {output.bam} --design {params.design} -j {threads} >>{log.out} 2>>{log.err}
    sambamba flagstat -t {threads} {output.bam} > {output.bam}.flagstat 2>>{log.err}
    '''

        
rule refine:
    input:
        bam   = '03.clip_umi_barcode/{sample}/{sample}.clip_umi_barcode.flt.bam',
        bami  = '03.clip_umi_barcode/{sample}/{sample}.clip_umi_barcode.flt.bam.pbi',
        primer3 = 'ref/10x_3kit_primers.fasta',
    output:
        bam   = temp('04.refine/{sample}/{sample}.refine.fltnc.bam'),
        bami  = temp('04.refine/{sample}/{sample}.refine.fltnc.bam.pbi'),
        report= '04.refine/{sample}/{sample}.refine.fltnc.report.csv'
    log:
        out = snakedir + "/logs/A4.refine/{sample}.o",
        err = snakedir + "/logs/A4.refine/{sample}.e",
    threads: 16
    resources: 
        mem = '16g',
        extra = ' --gres=lscratch:10 '
    shell:
        '''module load {config[modules][smrtanalysis]} {config[modules][sambamba]} > {log.out} 2>{log.err}
    isoseq3 refine {input.bam} {input.primer3} {output.bam} --require-polya -j {threads} >>{log.out} 2>>{log.err}
    sambamba flagstat -t {threads} {output.bam} > {output.bam}.flagstat 2>>{log.err}
    '''
    
rule correction:
    input:
        bam   = '04.refine/{sample}/{sample}.refine.fltnc.bam',
        bami  = '04.refine/{sample}/{sample}.refine.fltnc.bam.pbi',
        barcode = 'ref/barcodes.txt.gz',
    output:
        bam   = temp('05.correction/{sample}/{sample}.corrected.bam'),
        bami  = temp('05.correction/{sample}/{sample}.corrected.bam.pbi'),
    log:
        out = snakedir + "/logs/A5.correction/{sample}.o",
        err = snakedir + "/logs/A5.correction/{sample}.e",
    threads: 16
    resources: 
        mem = '32g',
        extra = ' --gres=lscratch:10 '
    shell:
        '''module load {config[modules][smrtanalysis]} {config[modules][sambamba]} > {log.out} 2>{log.err}
    isoseq3 correct -B {input.barcode} {input.bam} {output.bam} -j {threads} >>{log.out} 2>>{log.err}
    sambamba flagstat -t {threads} {output.bam} > {output.bam}.flagstat 2>>{log.err}
    '''

rule deduplication:
    input:
        bam   = '05.correction/{sample}/{sample}.corrected.bam',
        bami  = '05.correction/{sample}/{sample}.corrected.bam.pbi',
    output:
        bam   = temp('06.deduplication/{sample}/{sample}.dedup.bam'),
        bami  = temp('06.deduplication/{sample}/{sample}.dedup.bam.pbi'),
        fasta = '06.deduplication/{sample}/{sample}.dedup.fasta',
    params:
        bam   = '06.deduplication/{sample}/{sample}.sorted.bam',
    log:
        out = snakedir + "/logs/A6.deduplication/{sample}.o",
        err = snakedir + "/logs/A6.deduplication/{sample}.e",
    threads: 16
    resources: 
        mem = '32g',
        extra = ' --gres=lscratch:10 '
    shell:
        '''module load {config[modules][smrtanalysis]} {config[modules][sambamba]} > {log.out} 2>{log.err}
    samtools sort -@ {threads} -t CB {input.bam} -o {params.bam} >>{log.out} 2>>{log.err}
    isoseq3 groupdedup {params.bam} {output.bam} -j {threads} >>{log.out} 2>>{log.err}
    rm {params.bam} >>{log.out} 2>>{log.err}
    sambamba flagstat -t {threads} {output.bam} > {output.bam}.flagstat 2>>{log.err}
    '''
        
rule alignment:
    input:
        bam   = '06.deduplication/{sample}/{sample}.dedup.bam',
        bami  = '06.deduplication/{sample}/{sample}.dedup.bam.pbi',
        genome_fa   = 'ref/genome.fasta',
    output:
        bam   = '07.alignment/{sample}/{sample}.bam',
        bami  = '07.alignment/{sample}/{sample}.bam.bai',
    params:
        preset = 'ISOSEQ',
    log:
        out = snakedir + "/logs/A7.alignment/{sample}.o",
        err = snakedir + "/logs/A7.alignment/{sample}.e",
    threads: 16
    resources: 
        mem = '32g',
        extra = ' --gres=lscratch:10 '
    shell:
        '''module load {config[modules][smrtanalysis]} {config[modules][sambamba]} > {log.out} 2>{log.err}
    pbmm2 align -j {threads} --preset {params.preset} --sort {input.bam} {input.genome_fa} {output.bam} >>{log.out} 2>>{log.err}
    sambamba flagstat -t {threads} {output.bam} > {output.bam}.flagstat 2>>{log.err}
    '''
    
rule collapse:
    input:
        bam   = '07.alignment/{sample}/{sample}.bam',
        bami  = '07.alignment/{sample}/{sample}.bam.bai',
        genome_fa   = 'ref/genome.fasta',
        anno_gtf  = 'ref/gencode.gtf',
        tssbed = 'ref/refTSS.bed',
        polya = 'ref/polyA.list.txt',
    output:
        collapsed_gff = temp('08.annotation/{sample}/{sample}.collapsed.gff'),
        collapsed_group = '08.annotation/{sample}/{sample}.collapsed.group.txt',
        collapsed_sorted_gff = temp('08.annotation/{sample}/{sample}.collapsed.sorted.gff'),
        collapsed_abundance = '08.annotation/{sample}/{sample}.collapsed.abundance.txt',
        classification = '08.annotation/{sample}/{sample}_classification.txt',
        classificationflt = '08.annotation/{sample}/{sample}_classification.filtered_lite_classification.txt',
        saturation = '08.annotation/{sample}/{sample}_saturation.txt',
    params:
        prefix = '08.annotation/{sample}/{sample}',
    log:
        out = snakedir + "/logs/A8.collapse/{sample}.o",
        err = snakedir + "/logs/A8.collapse/{sample}.e",
    threads: 8
    resources: 
        mem = '32g',
        extra = ' --gres=lscratch:10 '
    shell:
        '''module load {config[modules][smrtanalysis]} {config[modules][sambamba]} > {log.out} 2>{log.err}
    isoseq3 collapse \
        {input.bam} \
        {output.collapsed_gff} >>{log.out} 2>>{log.err}
    pigeon sort \
        {output.collapsed_gff} \
        -o {output.collapsed_sorted_gff} >>{log.out} 2>>{log.err}
    pigeon classify \
        {output.collapsed_sorted_gff} \
        {input.anno_gtf} \
        {input.genome_fa} \
        --fl {output.collapsed_abundance} \
        --cage-peak {input.tssbed} \
        --poly-a {input.polya} \
        -o {params.prefix} \
        -j {threads} >>{log.out} 2>>{log.err}
    pigeon filter \
        {output.classification} \
        --isoforms {output.collapsed_sorted_gff} \
        -j {threads} >>{log.out} 2>>{log.err}
    pigeon report \
        {output.classificationflt} \
        {output.saturation} \
        -j {threads} >>{log.out} 2>>{log.err}
    '''

rule seurat_prep:
    input:
        fasta = '06.deduplication/{sample}/{sample}.dedup.fasta',
        collapsed_group = '08.annotation/{sample}/{sample}.collapsed.group.txt',
        classificationflt = '08.annotation/{sample}/{sample}_classification.filtered_lite_classification.txt',
    output:
        info = '09.seurat_prep/{sample}/genes_seurat/matrix.mtx',
    params:
        prefix = '09.seurat_prep/{sample}',
    log:
        out = snakedir + "/logs/A9.seurat_prep/{sample}.o",
        err = snakedir + "/logs/A9.seurat_prep/{sample}.e",
    threads: 16
    resources: 
        mem = '64g',
        extra = ' --gres=lscratch:10 '
    shell:
        '''module load {config[modules][smrtanalysis]} {config[modules][sambamba]} > {log.out} 2>{log.err}
        pigeon make-seurat \
            --dedup {input.fasta} \
            --group {input.collapsed_group} \
            -d {params.prefix} \
            {input.classificationflt} \
            -j {threads} >>{log.out} 2>>{log.err}
        '''
        
        


