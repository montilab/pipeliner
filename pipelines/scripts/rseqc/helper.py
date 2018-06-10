#!/usr/bin/env python
"""
Common utility functions
"""

import os 
import re
import sys 
import bz2
import gzip 
import numpy 

def init_gene():
    """
    Initializing the gene structure 
    """

    gene_det = [('id', 'f8'), 
            ('anno_id', numpy.dtype), 
            ('confgenes_id', numpy.dtype),
            ('name', 'S25'),
            ('source', 'S25'),
            ('gene_info', numpy.dtype),
            ('alias', 'S15'),
            ('name2', numpy.dtype),
            ('strand', 'S2'), 
            ('score', 'S15'), 
            ('chr', 'S25'), 
            ('chr_num', numpy.dtype),
            ('paralogs', numpy.dtype),
            ('start', 'f8'),
            ('stop', 'f8'), 
            ('transcripts', numpy.dtype),
            ('transcript_type', numpy.dtype),
            ('transcript_info', numpy.dtype),
            ('transcript_score', numpy.dtype),
            ('transcript_status', numpy.dtype),
            ('transcript_valid', numpy.dtype),
            ('exons', numpy.dtype),
            ('exons_confirmed', numpy.dtype),
            ('cds_exons', numpy.dtype),
            ('utr5_exons', numpy.dtype),
            ('utr3_exons', numpy.dtype),
            ('tis', numpy.dtype),
            ('tis_conf', numpy.dtype),
            ('tis_info', numpy.dtype),
            ('cdsStop', numpy.dtype),
            ('cdsStop_conf', numpy.dtype),
            ('cdsStop_info', numpy.dtype),
            ('tss', numpy.dtype),
            ('tss_info', numpy.dtype),
            ('tss_conf', numpy.dtype),
            ('cleave', numpy.dtype),
            ('cleave_info', numpy.dtype),
            ('cleave_conf', numpy.dtype),
            ('polya', numpy.dtype),
            ('polya_info', numpy.dtype),
            ('polya_conf', numpy.dtype),
            ('is_alt', 'f8'), 
            ('is_alt_spliced', 'f8'), 
            ('is_valid',  numpy.dtype),
            ('transcript_complete', numpy.dtype),
            ('is_complete', numpy.dtype),
            ('is_correctly_gff3_referenced', 'S5'),
            ('splicegraph', numpy.dtype) ]

    return gene_det

def open_file(fname):
    """
    Open the file (supports .gz .bz2) and returns the handler
    @args fname: input file name for reading 
    @type fname: str
    """

    if os.path.isfile(fname):
        file_prefx, ext = os.path.splitext(fname)
    else:
        exit("error: the provided file %s is not available to read. Please check!" % fname)
    
    try:
        if ext == ".gz":
            FH = gzip.open(fname, 'rb')
        elif ext == ".bz2":
            FH = bz2.BZ2File(fname, 'rb')
        else:
            FH = open(fname, 'rU')
    except Exception as error:
        sys.exit(error)

    return FH

def add_CDS_phase(strand, cds):
    """
    Calculate CDS phase and add to the CDS exons
    @args strand: feature strand information 
    @type strand: +/- 
    @args cds: coding exon coordinates 
    @type cds: numpy array [[int, int, int]]
    """

    cds_region, cds_flag = [], 0 
    if strand == '+':
        for cdspos in cds:
            if cds_flag == 0:
                cdspos = (cdspos[0], cdspos[1], 0)
                diff = (cdspos[1]-(cdspos[0]-1))%3
            else:
                xy = 0
                if diff == 0: 
                    cdspos = (cdspos[0], cdspos[1], 0)
                elif diff == 1: 
                    cdspos = (cdspos[0], cdspos[1], 2)
                    xy = 2
                elif diff == 2: 
                    cdspos = (cdspos[0], cdspos[1], 1)
                    xy = 1
                diff = ((cdspos[1]-(cdspos[0]-1))-xy)%3
            cds_region.append(cdspos)
            cds_flag = 1 
    elif strand == '-':
        cds.reverse()
        for cdspos in cds: 
            if cds_flag == 0:
                cdspos = (cdspos[0], cdspos[1], 0)
                diff = (cdspos[1]-(cdspos[0]-1))%3
            else:  
                xy = 0 
                if diff == 0: 
                    cdspos = (cdspos[0], cdspos[1], 0)
                elif diff == 1:
                    cdspos = (cdspos[0], cdspos[1], 2)
                    xy = 2
                elif diff == 2: 
                    cdspos = (cdspos[0], cdspos[1], 1)
                    xy = 1
                diff = ((cdspos[1]-(cdspos[0]-1))-xy)%3
            cds_region.append(cdspos)
            cds_flag = 1
        cds_region.reverse()
    return cds_region

def buildUTR(cc, ec, strand):
    """
    Build UTR regions from a given set of CDS and exon coordiantes of a gene
    @args cc: coding exon coordinates 
    @type cc: numpy array [[int, int, int]]
    @args ec: exon coordinates 
    @type ec: numpy array [[int, int]]
    @args strand: feature strand information 
    @type strand: +/- 
    """

    utr5 = []
    utr3 = []
    if strand == '+':
        cds_s = cc[0][0]
        for ex in ec:
            if ex[0] <= cds_s and cds_s <= ex[1]:
                if ex[0] != cds_s:utr5.append((ex[0], cds_s-1))
                break
            else:
                utr5.append(ex)
        cds_e = cc[-1][1]
        for i in range(len(ec)):
            i += 1
            if ec[-i][0] <= cds_e and cds_e <= ec[-i][1]:
                if ec[-i][1] != cds_e:utr3.append((cds_e +1, ec[-i][1]))
                break
            else:
                utr3.append(ec[-i]) 
        utr3.reverse()
    elif strand == '-':
        cds_s = cc[-1][1]
        for i in range(len(ec)):
            i += 1
            if ec[-i][0] <= cds_s and cds_s <= ec[-i][1]:
                if ec[-i][1] != cds_s:utr5.append((cds_s+1, ec[-i][1]))
                break
            else:
                utr5.append(ec[-i])
        utr5.reverse()
        cds_e = cc[0][0] 
        for ex in ec:
            if ex[0] <= cds_e and cds_e <= ex[1]:
                if ex[0] != cds_e:utr3.append((ex[0], cds_e-1))
                break
            else:
                utr3.append(ex)
    return utr5, utr3

def make_Exon_cod(strand_p, five_p_utr, cds_cod, three_p_utr):
    """
    Create exon cordinates from UTR's and CDS region
    @args strand_p: feature strand information 
    @type strand_p: +/- 
    @args five_p_utr: five prime utr exon coordinates 
    @type five_p_utr: numpy array [[int, int]]
    @args cds_cod: coding exon coordinates 
    @type cds_cod: numpy array [[int, int, int]]
    @args three_p_utr: three prime utr exon coordinates 
    @type three_p_utr: numpy array [[int, int]]
    """

    exon_pos = []
    if strand_p == '+':        
        utr5_start, utr5_end = 0, 0
        if five_p_utr != []:
            utr5_start, utr5_end = five_p_utr[-1][0], five_p_utr[-1][1] 
        cds_5start, cds_5end = cds_cod[0][0], cds_cod[0][1]
        jun_exon = []
        if cds_5start-utr5_end == 0 or cds_5start-utr5_end == 1:
            jun_exon = [utr5_start, cds_5end]    
        if len(cds_cod) == 1:
            five_prime_flag = 0
            if jun_exon != []:
                five_p_utr = five_p_utr[:-1]
                five_prime_flag = 1
            for utr5 in five_p_utr:
                exon_pos.append(utr5)
            jun_exon = []
            utr3_start, utr3_end = 0, 0
            if three_p_utr != []: 
                utr3_start = three_p_utr[0][0]
                utr3_end = three_p_utr[0][1]
            if utr3_start-cds_5end == 0 or utr3_start-cds_5end == 1:
                jun_exon = [cds_5start, utr3_end]
            three_prime_flag = 0
            if jun_exon != []: 
                cds_cod = cds_cod[:-1]
                three_p_utr = three_p_utr[1:]
                three_prime_flag = 1
            if five_prime_flag == 1 and three_prime_flag == 1:
                exon_pos.append([utr5_start, utr3_end])
            if five_prime_flag == 1 and three_prime_flag == 0:
                exon_pos.append([utr5_start, cds_5end])
                cds_cod = cds_cod[:-1]
            if five_prime_flag == 0 and three_prime_flag == 1:
                exon_pos.append([cds_5start, utr3_end])
            for cds in cds_cod:
                exon_pos.append(cds)
            for utr3 in three_p_utr:
                exon_pos.append(utr3)
        else:    
            if jun_exon != []:
                five_p_utr = five_p_utr[:-1]
                cds_cod = cds_cod[1:]
            for utr5 in five_p_utr:
                exon_pos.append(utr5)
            exon_pos.append(jun_exon) if jun_exon != [] else ''
            jun_exon = []
            utr3_start, utr3_end = 0, 0
            if three_p_utr != []:
                utr3_start = three_p_utr[0][0]
                utr3_end = three_p_utr[0][1]
            cds_3start = cds_cod[-1][0]
            cds_3end = cds_cod[-1][1]
            if utr3_start-cds_3end == 0 or utr3_start-cds_3end == 1:       
                jun_exon = [cds_3start, utr3_end]
            if jun_exon != []:
                cds_cod = cds_cod[:-1]
                three_p_utr = three_p_utr[1:]
            for cds in cds_cod:
                exon_pos.append(cds)
            exon_pos.append(jun_exon) if jun_exon != [] else ''
            for utr3 in three_p_utr:
                exon_pos.append(utr3)
    elif strand_p == '-':
        utr3_start, utr3_end = 0, 0        
        if three_p_utr != []:
            utr3_start = three_p_utr[-1][0]
            utr3_end = three_p_utr[-1][1]
        cds_3start = cds_cod[0][0]
        cds_3end = cds_cod[0][1]
        jun_exon = []
        if cds_3start-utr3_end == 0 or cds_3start-utr3_end == 1:
            jun_exon = [utr3_start, cds_3end]  
        if len(cds_cod) == 1:    
            three_prime_flag = 0
            if jun_exon != []:
                three_p_utr = three_p_utr[:-1]
                three_prime_flag = 1
            for utr3 in three_p_utr:
                exon_pos.append(utr3)
            jun_exon = []
            (utr5_start, utr5_end) = (0, 0)
            if five_p_utr != []:
                utr5_start = five_p_utr[0][0]
                utr5_end = five_p_utr[0][1]
            if utr5_start-cds_3end == 0 or utr5_start-cds_3end == 1:
                jun_exon = [cds_3start, utr5_end]
            five_prime_flag = 0
            if jun_exon != []:
                cds_cod = cds_cod[:-1]
                five_p_utr = five_p_utr[1:]
                five_prime_flag = 1
            if three_prime_flag == 1 and five_prime_flag == 1:
                exon_pos.append([utr3_start, utr5_end])
            if three_prime_flag == 1 and five_prime_flag == 0:
                exon_pos.append([utr3_start, cds_3end])
                cds_cod = cds_cod[:-1]
            if three_prime_flag == 0 and five_prime_flag == 1:
                exon_pos.append([cds_3start, utr5_end])        
            for cds in cds_cod:
                exon_pos.append(cds)
            for utr5 in five_p_utr:
                exon_pos.append(utr5)
        else:
            if jun_exon != []:
                three_p_utr = three_p_utr[:-1]
                cds_cod = cds_cod[1:]
            for utr3 in three_p_utr:
                exon_pos.append(utr3)   
            if jun_exon != []:
                exon_pos.append(jun_exon)
            jun_exon = []
            (utr5_start, utr5_end) = (0, 0)
            if five_p_utr != []:
                utr5_start = five_p_utr[0][0]
                utr5_end = five_p_utr[0][1]    
            cds_5start = cds_cod[-1][0]
            cds_5end = cds_cod[-1][1]
            if utr5_start-cds_5end == 0 or utr5_start-cds_5end == 1:
                jun_exon = [cds_5start, utr5_end]
            if jun_exon != []:
                cds_cod = cds_cod[:-1]
                five_p_utr = five_p_utr[1:]
            for cds in cds_cod:
                exon_pos.append(cds)
            if jun_exon != []:
                exon_pos.append(jun_exon)    
            for utr5 in five_p_utr:
                exon_pos.append(utr5)
    return exon_pos