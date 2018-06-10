#!/usr/bin/env python3
# @Author: Anthony Federico

import sys
import gzip
import itertools

def file_handle(f):
    """
    Opens a file in context and returns the handle

    Args:
        Param #1 (str): Filename ending in (.fq, .fastq, .fq.gz, fastq.gz)
    Yields:
        File object
    """
    if f.endswith('.gz'):
        with gzip.open(f) as handle:
            yield handle
    else:
        with open(f) as handle:
            yield handle

def check_format(f):
    """
    Ensures every 1st and 3rd line starts with '@' and '+' respectively

    Args:
        Param #1 (str): Filename ending in (.fq, .fastq, .fq.gz, fastq.gz)
    Returns:
        Nothing
    """  
    for line in file_handle(f):
        lines = line.read().split('\n')
        for i, l in enumerate(lines):
            try:
                if (i)%4 == 0:
                    if l[0] != '@':
                        sys.exit("Error: Expected line {0}  in file '{1}' to start with '@'".format(i+1, f))
                if (i+2)%4 == 0:
                    if l[0] != '+':
                        sys.exit("Error: Expected line {0}  in file '{1}' to start with '+'".format(i+1, f))
            except IndexError:
                print("Warning: Line {0} in file '{1}' is empty".format(i+1, f))

def check_pairs(f1, f2):
    """
    Ensures paired read files have the same number of reads

    Args:
        Param #1 (str): Filename ending in (.fq, .fastq, .fq.gz, fastq.gz)
        Param #2 (str): Filename ending in (.fq, .fastq, .fq.gz, fastq.gz)
    Returns:
        Nothing
    """  
    reads1, reads2 = 0, 0
    for line in file_handle(f1):
        lines = line.read().split('\n')
        for i, l in enumerate(lines):
            try:
                if (i)%4 == 0:
                    if l[0] == '@':
                        reads1 += 1
            except IndexError:
                print("Warning: Line {0} in file '{1}' is empty".format(i+1, f1))

    for line in file_handle(f2):
        lines = line.read().split('\n')
        for i, l in enumerate(lines):
            try:
                if (i)%4 == 0:
                    if l[0] == '@':
                        reads2 += 1
            except IndexError:
                print("Warning: Line {0} in file '{1}' is empty".format(i+1, f2))

        if reads1 != reads2:
            sys.exit("Error: Paired files have different number of reads")

# ------------------------------------------------------------------------------

if __name__ == "__main__":    
    if len(sys.argv) == 2:
        file1 = sys.argv[1] 
        paired = False
    elif len(sys.argv) == 3:
        file1, file2 = sys.argv[1], sys.argv[2]
        paired = True
    else:
        sys.exit("Error: Please provide only 1 or 2 fastq files for single or paired reads")

    if paired:
        check_format(file1)
        check_format(file2)
        check_pairs(file1, file2)

    if not paired: 
        check_format(file1)