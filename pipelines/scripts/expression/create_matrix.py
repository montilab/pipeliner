#!/usr/bin/env python3
# @Author: Anthony Federico

import pandas as pd
import sys
import os

def aggregate_counts(files, method):
    """
    Merges count files into a single pandas dataframe
    Args:
        Param #1 (list): List of paths to count files
        Param #1 (str): Quantification method
    Returns:
        Pandas dataframe
    """      
    counts = {}
    for filename in files:
        sample = filename.split('/')[-1].split('.')[0]
        with open(filename) as infile:
            for l, line in enumerate(infile):
                
                if method == 'htseq':
                    if not line.startswith('_'):
                        stripped = line.strip('\n')
                        splitted = stripped.split('\t')
                        gene = splitted[0]
                        count = splitted[1]
                        try:
                            counts[sample][gene] = count
                        except KeyError:
                            counts[sample] = {gene: count}

                elif method == 'featurecounts':
                    if not line.startswith('#') and l > 1:
                        stripped = line.strip('\n')
                        splitted = stripped.split('\t')
                        gene = splitted[0]
                        count = splitted[-1]
                        try:
                            counts[sample][gene] = count
                        except KeyError:
                            counts[sample] = {gene: count}

                else:
                    raise RuntimeError("{0} is not a valid method".format(method))

    df = pd.DataFrame(counts)
    return df

def normalize_counts(files):
    """
    Merges normalized count files into a single pandas dataframe
    Args:
        Param #1 (list): List of paths to normalized count files
    Returns:
        Pandas dataframe
    """
    fpkms, tpms = {}, {}
    for filename in files:
        sample = filename.split('/')[-1].split('.')[0]
        with open(filename) as infile:
            for l, line in enumerate(infile):
                if l > 0:
                    stripped = line.strip('\n')
                    splitted = stripped.split('\t')
                    gene = splitted[0]
                    fpkm = splitted[-2]
                    tpm = splitted[-1]
                    if gene.startswith('ENSG'):
                        try:
                            fpkms[sample][gene] = fpkm
                            tpms[sample][gene] = tpm
                        except KeyError:
                            fpkms[sample] = {gene: fpkm}
                            tpms[sample] = {gene: tpm}

    fpkms = pd.DataFrame(fpkms)
    tpms = pd.DataFrame(tpms)
    return fpkms, tpms

def reindex_samples(counts, phenotypes):
    """
    Reorders columns in count matrix based on the
    sample order found in phenotypes file
    Args:
        Param #1 (frame): Pandas dataframe of counts
    Returns:
        Pandas dataframe
    """
    df = pd.read_table(phenotypes).transpose()
    df.columns = df.iloc[0]
    df = df.drop(df.index[0])
    return counts.reindex(columns=list(df.columns))

if __name__ == '__main__':

    if sys.argv[1] == '-p':
        phenotypes = sys.argv[2]
        method     = sys.argv[3]
        files      = sys.argv[4:]
        os.popen('cp {0} phenotypes.txt'.format(phenotypes))  
        
        if method == 'htseq' or method == 'featurecounts':
            df = aggregate_counts(files, method)
            df = reindex_samples(df, phenotypes)
            df.to_csv('count_expression_matrix.txt', sep='\t')

        elif method == 'stringtie':
            fpkms, tpms = normalize_counts(files)
            fpkms = reindex_samples(fpkms, phenotypes)
            tpms = reindex_samples(tpms, phenotypes)
            fpkms.to_csv('fpkm_expression_matrix.txt', sep='\t')
            tpms.to_csv('tpm_expression_matrix.txt', sep='\t')

    else:
        method     = sys.argv[1]
        files      = sys.argv[2:]

        if method == 'htseq' or method == 'featurecounts':
            df = aggregate_counts(files, method)
            df.to_csv('count_expression_matrix.txt', sep='\t')

        elif method == 'stringtie':
            fpkms, tpms = normalize_counts(files)
            fpkms.to_csv('fpkm_expression_matrix.txt', sep='\t')
            tpms.to_csv('tpm_expression_matrix.txt', sep='\t')
