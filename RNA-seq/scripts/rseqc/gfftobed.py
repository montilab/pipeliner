#!/usr/bin/env python
"""
Convert genome annotation data in GFF/GTF to a 12 column BED format. 
BED format typically represents the transcript models. 
Usage: python gff_to_bed.py in.gff > out.bed  
Requirement:
    gffparser.py: https://github.com/vipints/GFFtools-GX/blob/master/gffparser.py    
Copyright (C) 
    2009-2012 Friedrich Miescher Laboratory of the Max Planck Society, Tubingen, Germany.
    2012-2015 Memorial Sloan Kettering Cancer Center New York City, USA.
"""

import re
import sys
import gffparser

def limitBEDWrite(tinfo):
    """
    Write a three column BED file 
    
    @args tinfo: list of genes 
    @type tinfo: numpy object  
    """

    for contig_id, feature in tinfo.items():
        uns_line = dict()
        for tid, tloc in feature.items():
            uns_line[(int(tloc[0])-1, int(tloc[1]))]=1
        for ele in sorted(uns_line):
            pline = [contig_id,
                    str(ele[0]-1),
                    str(ele[1])]

            sys.stdout.write('\t'.join(pline)+"\n")


def writeBED(tinfo):
    """
    writing result files in bed format 
    @args tinfo: list of genes 
    @type tinfo: numpy object  
    """

    for ent1 in tinfo:
        child_flag = False  

        for idx, tid in enumerate(ent1['transcripts']):
            child_flag = True 
            exon_cnt = len(ent1['exons'][idx])
            exon_len = ''
            exon_cod = '' 
            rel_start = None 
            rel_stop = None 
            for idz, ex_cod in enumerate(ent1['exons'][idx]):#check for exons of corresponding transcript  
                exon_len += '%d,' % (ex_cod[1]-ex_cod[0]+1)
                if idz == 0: #calculate the relative start position 
                    exon_cod += '0,'
                    rel_start = int(ex_cod[0])-1 
                    rel_stop = int(ex_cod[1])
                else:
                    exon_cod += '%d,' % (ex_cod[0]-1-rel_start) ## shifting the coordinates to zero 
                    rel_stop = int(ex_cod[1])
            
            if exon_len:
                score = 0 
                score = ent1['transcript_score'][idx] if ent1['transcript_score'].any() else score ## getting the transcript score 
                out_print = [ent1['chr'],
                            str(rel_start),
                            str(rel_stop),
                            tid[0],
                            str(score), 
                            ent1['strand'], 
                            str(rel_start),
                            str(rel_stop),
                            '0',
                            str(exon_cnt),
                            exon_len,
                            exon_cod]
                sys.stdout.write('\t'.join(out_print)+"\n")
        
        if not child_flag: # file just contains only a single parent type i.e, gff3 defines only one feature type 
            score = 0 
            score = ent1['transcript_score'][0] if ent1['transcript_score'].any() else score

            out_print = [ent1['chr'], 
                        '%d' % (int(ent1['start'])-1), 
                        '%d' % int(ent1['stop']),
                        ent1['name'], 
                        str(score), 
                        ent1['strand'],
                        '%d' % int(ent1['start']), 
                        '%d' % int(ent1['stop']),
                        '0',
                        '1',
                        '%d,' % (int(ent1['stop'])-int(ent1['start'])+1), 
                        '0,']

            sys.stdout.write('\t'.join(out_print)+"\n")

    
def __main__():
    try:
        query_file = sys.argv[1]
    except:
        print __doc__
        sys.exit(-1)

    Transcriptdb = gffparser.Parse(query_file)  
    writeBED(Transcriptdb)

if __name__ == "__main__": 
    __main__() 