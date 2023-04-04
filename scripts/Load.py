#!/usr/bin/env python
import sys
import pandas as pd

def samplesheet(sspath):
    df = pd.read_csv(sspath, sep='\t')
    sample = {}
    for i in df.index:
        ss = df.loc[i,'Sample']
        bam= df.loc[i,'CCS_bam']
        sample[ss] = bam
    return sample

if __name__ == '__main__':
    print(samplesheet(sys.argv[1]))
