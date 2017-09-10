#!/usr/bin/env python3
# @Author: Dileep Kishore <dileep>
# @Last modified by: anfederico

import sys
import csv
from functools import partial
from operator import itemgetter
from collections import defaultdict

def sample_reader(fname, gid, fid, tid):
    fdict, tdict = dict(), dict()
    with open(fname, 'r') as file_id:
        next(file_id)
        for row in csv.reader(file_id, delimiter='\t'):
            k, fdict[k], tdict[k] = itemgetter(gid, fid, tid)(row)
    return fdict, tdict, fname.split('.')[0]

def writefiles(counts, fname):
    with open('{0}.csv'.format(fname), 'w') as counts_file:
        genes = set(counts.keys())
        samples = set(j for i in counts.values() for j in i.keys())
        sample_names = [j.rsplit('/')[-1] for j in samples]
        w = csv.writer(counts_file)
        w.writerow(["gene_id"] + sample_names)
        for gid in genes:
            row_item = [gid]
            for sid in samples:
                try:
                    row_item.append(counts[gid][sid])
                except KeyError:
                    row_item.append(str(0.0))
            w.writerow(row_item)

def main(file_list, gid, fid, tid):
    count_reader = partial(sample_reader, gid=gid, fid=fid, tid=tid)
    fpkm = defaultdict(dict)
    tpm = defaultdict(dict)
    for fdict, tdict, sname in map(count_reader, file_list):
        for k, v in fdict.items():
            fpkm[k][sname] = v
        for k, v in tdict.items():
            tpm[k][sname] = v
    writefiles(fpkm, 'fpkm')
    writefiles(tpm, 'tpm')

if __name__ == '__main__':
    FILE_LIST = sys.argv[1:]
    assert len(FILE_LIST) >= 1
    with open(FILE_LIST[0], 'r') as file_id:
        header = file_id.readline().strip('\n').split('\t')
    GID, FID, TID = (header.index(x) for x in ['Gene ID', 'FPKM', 'TPM'])
    main(FILE_LIST, GID, FID, TID)