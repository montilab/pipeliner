#!/usr/bin/env python3
# @Author: Anthony Federico

import pandas as pd
import sys

def get_fields(gtf):
    """
    Extracts all additional fields in gtf

    Args:
        Param #1 (string): Path to gtf file
    Returns:
        Dictionary of features
    """      
    genes = {}
    with open(gtf) as infile:
        for line in infile:
            attributes = line.strip('\n').split('\t')[-1]
            attributes = attributes.lstrip().rstrip()
            fields = {}
            for attr in attributes.rstrip(';').split('; '):
               fields[attr.split()[0]] = attr.split()[1].strip('"')
            genes[fields['gene_id']] = fields
    return genes

def create_fdata(genes, matrix):
    """
    Creates a feature matrix and reorders rows
    to match gene order of count matrix

    Args:
        Param #1 (string): Dictionary of features
        Param #2 (string): Path to count matrix
    Returns:
        Pandas dataframe
    """    
    counts = pd.read_table(matrix, index_col=0)
    for i, gene in enumerate(counts.index):
        if i == 0:
            columns = {key: [genes[gene][key]] for key in genes[gene].keys()}
        else:
            for key in genes[gene].keys():
                columns[key].append(genes[gene][key])
    fdata = pd.DataFrame(columns)
    fdata.index = counts.index 
    del fdata['gene_id']
    return fdata

if __name__ == '__main__':

    gtf    = sys.argv[1]
    matrix = sys.argv[2]

    genes = get_fields(gtf)
    df = create_fdata(genes, matrix)

    df.to_csv('features.txt', sep='\t')