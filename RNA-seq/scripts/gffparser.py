#!/usr/bin/env python
"""
Extract genome annotation from a GFF (a tab delimited format for storing sequence features and annotations) file.
Requirements: 
    Numpy :- http://numpy.org/ 
Copyright (C)	
2009-2012 Friedrich Miescher Laboratory of the Max Planck Society, Tubingen, Germany. 
2012-2015 Memorial Sloan Kettering Cancer Center, New York City, USA.
"""

import re
import os
import sys
import urllib
import numpy as np
import helper as utils 
from collections import defaultdict

def attribute_tags(col9):
    """ 
    Split the key-value tags from the attribute column, it takes column number 9 from GTF/GFF file 
    @args col9: attribute column from GFF file 
    @type col9: str
    """
    info = defaultdict(list)
    is_gff = False
    
    if not col9:
        return is_gff, info
        
    # trim the line ending semi-colon  ucsc may have some white-space  
    col9 = col9.rstrip(';| ')
    # attributes from 9th column 
    atbs = col9.split(" ; ")
    if len(atbs) == 1:
        atbs = col9.split("; ")
        if len(atbs) == 1:
            atbs = col9.split(";")
    # check the GFF3 pattern which has key value pairs like:
    gff3_pat = re.compile("\w+=")
    # sometime GTF have: gene_id uc002zkg.1;
    gtf_pat = re.compile("\s?\w+\s")

    key_vals = []

    if gff3_pat.match(atbs[0]): # gff3 pattern 
        is_gff = True
        key_vals = [at.split('=') for at in atbs]
    elif gtf_pat.match(atbs[0]): # gtf pattern
        for at in atbs:
            key_vals.append(at.strip().split(" ",1))
    else:
        # to handle attribute column has only single value 
        key_vals.append(['ID', atbs[0]])
    # get key, val items 
    for item in key_vals:
        key, val = item
        # replace the double qoutes from feature identifier 
        val = re.sub('"', '', val)
        # replace the web formating place holders to plain text format 
        info[key].extend([urllib.unquote(v) for v in val.split(',') if v])

    return is_gff, info
                
def spec_features_keywd(gff_parts):
    """
    Specify the feature key word according to the GFF specifications
    @args gff_parts: attribute field key 
    @type gff_parts: str 
    """
    for t_id in ["transcript_id", "transcriptId", "proteinId"]:
        try:
            gff_parts["info"]["Parent"] = gff_parts["info"][t_id]
            break
        except KeyError:
            pass
    for g_id in ["gene_id", "geneid", "geneId", "name", "gene_name", "genename"]:
        try:
            gff_parts["info"]["GParent"] = gff_parts["info"][g_id]
            break
        except KeyError:
            pass
    ## TODO key words
    for flat_name in ["Transcript", "CDS"]:
        if gff_parts["info"].has_key(flat_name):
            # parents
            if gff_parts['type'] in [flat_name] or re.search(r'transcript', gff_parts['type'], re.IGNORECASE):
                if not gff_parts['id']:
                    gff_parts['id'] = gff_parts['info'][flat_name][0]
                    #gff_parts["info"]["ID"] = [gff_parts["id"]]
            # children 
            elif gff_parts["type"] in ["intron", "exon", "three_prime_UTR",
                        "coding_exon", "five_prime_UTR", "CDS", "stop_codon",
                        "start_codon"]:
                gff_parts["info"]["Parent"] = gff_parts["info"][flat_name]
            break
    return gff_parts

def Parse(ga_file):
    """
    Parsing GFF/GTF file based on feature relationship, it takes the input file.
    @args ga_file: input file name 
    @type ga_file: str 
    """
    child_map = defaultdict(list)
    parent_map = dict()

    ga_handle = utils.open_file(ga_file)

    for rec in ga_handle:
        rec = rec.strip('\n\r')
        
        # skip empty line fasta identifier and commented line
        if not rec or rec[0] in  ['#', '>']:
            continue
        # skip the genome sequence 
        if not re.search('\t', rec):
            continue

        parts = rec.split('\t')
        assert len(parts) >= 8, rec

        # process the attribute column (9th column)
        ftype, tags = attribute_tags(parts[-1])
        if not tags: # skip the line if no attribute column.
	        continue 

        # extract fields  
        if parts[1]:
            tags["source"] = parts[1]
        if parts[7]:
            tags["phase"] = parts[7]

        gff_info = dict()
        gff_info['info'] = dict(tags)
        gff_info["is_gff3"] = ftype
        gff_info['chr'] = parts[0]
        gff_info['score'] = parts[5]

        if parts[3] and parts[4]:
            gff_info['location'] = [int(parts[3]) ,
                        int(parts[4])]
            gff_info['type'] = parts[2]
            gff_info['id'] = tags.get('ID', [''])[0]
            if parts[6] in ['?', '.']:
                parts[6] = None 
            gff_info['strand'] = parts[6]

            # key word according to the GFF spec.
            # is_gff3 flag is false check this condition and get the attribute fields 
            if not ftype:
                gff_info = spec_features_keywd(gff_info)
            
            # link the feature relationships
            if gff_info['info'].has_key('Parent'): 
                for p in gff_info['info']['Parent']:
                    if p == gff_info['id']:
                        gff_info['id'] = ''
                        break
                rec_category = 'child'
            elif gff_info['id']:
                rec_category = 'parent'
            else:
                rec_category = 'record'

            # depends on the record category organize the features
            if rec_category == 'child':
                for p in gff_info['info']['Parent']:
                    # create the data structure based on source and feature id 
                    child_map[(gff_info['chr'], gff_info['info']['source'], p)].append(
                                            dict( type = gff_info['type'], 
                                            location =  gff_info['location'], 
                                            strand = gff_info['strand'], 
                                            score = gff_info['score'], 
                                            ID = gff_info['id'],
                                            gene_name = gff_info['info'].get('gene_name', [''])[0],
                                            gene_id = gff_info['info'].get('GParent', [''])[0] 
                                            ))
            elif rec_category == 'parent':
                parent_map[(gff_info['chr'], gff_info['info']['source'], gff_info['id'])] = dict( 
                                            type = gff_info['type'], 
                                            location = gff_info['location'],
                                            strand = gff_info['strand'],
                                            score = gff_info['score'], 
                                            name = tags.get('Name', [''])[0])
            elif rec_category == 'record':
                #TODO how to handle plain records?
                considering_soon = 1 
    ga_handle.close()

    # depends on file type create parent feature  
    if not ftype:
        parent_map, child_map = create_missing_feature_type(parent_map, child_map)    
    
    # connecting parent child relations  
    # essentially the parent child features are here from any type of GTF/GFF2/GFF3 file
    gene_mat = format_gene_models(parent_map, child_map) 

    return gene_mat 
    
def format_gene_models(parent_nf_map, child_nf_map): 
    """
    Genarate GeneObject based on the parsed file contents
    @args parent_nf_map: parent features with source and chromosome information 
    @type parent_nf_map: collections defaultdict
    @args child_nf_map: transctipt and exon information are encoded 
    @type child_nf_map: collections defaultdict
    """

    g_cnt = 0 
    gene = np.zeros((len(parent_nf_map),), dtype = utils.init_gene())

    for pkey, pdet in parent_nf_map.items():
        # considering only gene features 
        #if not re.search(r'gene', pdet.get('type', '')):
        #    continue 

        # infer the gene start and stop if not there in the 
        if not pdet.get('location', []):
            GNS, GNE = [], []
            # multiple number of transcripts 
            for L1 in child_nf_map[pkey]:
                GNS.append(L1.get('location', [])[0]) 
                GNE.append(L1.get('location', [])[1]) 
            GNS.sort()
            GNE.sort()
            pdet['location'] = [GNS[0], GNE[-1]]

        orient = pdet.get('strand', '')
        gene[g_cnt]['id'] = g_cnt +1 
        gene[g_cnt]['chr'] = pkey[0]
        gene[g_cnt]['source'] = pkey[1]
        gene[g_cnt]['name'] = pkey[-1]
        gene[g_cnt]['start'] = pdet.get('location', [])[0]
        gene[g_cnt]['stop'] = pdet.get('location', [])[1]
        gene[g_cnt]['strand'] = orient  
        gene[g_cnt]['score'] = pdet.get('score','')

        # default value 
        gene[g_cnt]['is_alt_spliced'] = gene[g_cnt]['is_alt'] = 0
        if len(child_nf_map[pkey]) > 1:
            gene[g_cnt]['is_alt_spliced'] = gene[g_cnt]['is_alt'] = 1

        # complete sub-feature for all transcripts 
        dim = len(child_nf_map[pkey])
        TRS = np.zeros((dim,), dtype=np.object)
        TR_TYP = np.zeros((dim,), dtype=np.object)
        EXON = np.zeros((dim,), dtype=np.object)
        UTR5 = np.zeros((dim,), dtype=np.object)
        UTR3 = np.zeros((dim,), dtype=np.object)
        CDS = np.zeros((dim,), dtype=np.object)
        TISc = np.zeros((dim,), dtype=np.object)
        TSSc = np.zeros((dim,), dtype=np.object)
        CLV = np.zeros((dim,), dtype=np.object)
        CSTOP = np.zeros((dim,), dtype=np.object)
        TSTAT = np.zeros((dim,), dtype=np.object)
        TSCORE = np.zeros((dim,), dtype=np.object)

        # fetching corresponding transcripts 
        for xq, Lv1 in enumerate(child_nf_map[pkey]):

            TID = Lv1.get('ID', '')
            TRS[xq]= np.array([TID])

            TYPE = Lv1.get('type', '')
            TR_TYP[xq] = np.array('')
            TR_TYP[xq] = np.array(TYPE) if TYPE else TR_TYP[xq]

            orient = Lv1.get('strand', '')
            tr_score = Lv1.get('score', '')

            # fetching different sub-features 
            child_feat = defaultdict(list)
            for Lv2 in child_nf_map[(pkey[0], pkey[1], TID)]:
                E_TYP = Lv2.get('type', '')
                child_feat[E_TYP].append(Lv2.get('location'))
            
            # make general ascending order of coordinates 
            if orient == '-':
                for etype, excod in child_feat.items():
                    if len(excod) > 1:
                        if excod[0][0] > excod[-1][0]:
                            excod.reverse()
                            child_feat[etype] = excod

            # make exon coordinate from cds and utr regions 
            if not child_feat.get('exon'):  
                if child_feat.get('CDS'):
                    exon_cod = utils.make_Exon_cod( orient, 
                                NonetoemptyList(child_feat.get('five_prime_UTR')), 
                                NonetoemptyList(child_feat.get('CDS')),
                                NonetoemptyList(child_feat.get('three_prime_UTR')))
                    child_feat['exon'] = exon_cod 
                else: 
                    # TODO only UTR's
                    # searching through keys to find a pattern describing exon feature 
                    ex_key_pattern = [k for k in child_feat if k.endswith("exon")]
                    if ex_key_pattern:
                        child_feat['exon'] = child_feat[ex_key_pattern[0]]

            # stop_codon are seperated from CDS, add the coordinates based on strand
            if child_feat.get('stop_codon'):
                if orient == '+':
                    if child_feat.get('stop_codon')[0][0] - child_feat.get('CDS')[-1][1] == 1:
                        child_feat['CDS'][-1] = [child_feat.get('CDS')[-1][0], child_feat.get('stop_codon')[0][1]]
                    else:
                        child_feat['CDS'].append(child_feat.get('stop_codon')[0])
                elif orient == '-':
                    if child_feat.get('CDS')[0][0] - child_feat.get('stop_codon')[0][1] == 1:
                        child_feat['CDS'][0] = [child_feat.get('stop_codon')[0][0], child_feat.get('CDS')[0][1]]
                    else:
                        child_feat['CDS'].insert(0, child_feat.get('stop_codon')[0])

            # transcript signal sites 
            TIS, cdsStop, TSS, cleave = [], [], [], []
            cds_status, exon_status, utr_status = 0, 0, 0

            if child_feat.get('exon'):
                TSS = [child_feat.get('exon')[-1][1]]
                TSS = [child_feat.get('exon')[0][0]] if orient == '+' else TSS 
                cleave = [child_feat.get('exon')[0][0]]
                cleave = [child_feat.get('exon')[-1][1]] if orient == '+' else cleave
                exon_status = 1

            if child_feat.get('CDS'):
                if orient == '+': 
                    TIS = [child_feat.get('CDS')[0][0]]
                    cdsStop = [child_feat.get('CDS')[-1][1]-3]
                else:
                    TIS = [child_feat.get('CDS')[-1][1]]
                    cdsStop = [child_feat.get('CDS')[0][0]+3]
                cds_status = 1 
                # cds phase calculation 
                child_feat['CDS'] = utils.add_CDS_phase(orient, child_feat.get('CDS'))
            
            # sub-feature status 
            if child_feat.get('three_prime_UTR') or child_feat.get('five_prime_UTR'):
                utr_status =1 
            
            if utr_status == cds_status == exon_status == 1: 
                t_status = 1
            else:
                t_status = 0
            
            # add sub-feature # make array for export to different out
            TSTAT[xq] = t_status
            EXON[xq] = np.array(child_feat.get('exon'), np.float64)
            UTR5[xq] = np.array(NonetoemptyList(child_feat.get('five_prime_UTR')))
            UTR3[xq] = np.array(NonetoemptyList(child_feat.get('three_prime_UTR')))
            CDS[xq] = np.array(NonetoemptyList(child_feat.get('CDS')))
            TISc[xq] = np.array(TIS)
            CSTOP[xq] = np.array(cdsStop)
            TSSc[xq] = np.array(TSS)
            CLV[xq] = np.array(cleave)
            TSCORE[xq] = tr_score 
            
        # add sub-features to the parent gene feature
        gene[g_cnt]['transcript_status'] = TSTAT
        gene[g_cnt]['transcripts'] = TRS 
        gene[g_cnt]['exons'] = EXON
        gene[g_cnt]['utr5_exons'] = UTR5 
        gene[g_cnt]['cds_exons'] = CDS 
        gene[g_cnt]['utr3_exons'] = UTR3 
        gene[g_cnt]['transcript_type'] = TR_TYP
        gene[g_cnt]['tis'] = TISc
        gene[g_cnt]['cdsStop'] = CSTOP
        gene[g_cnt]['tss'] = TSSc
        gene[g_cnt]['cleave'] = CLV
        gene[g_cnt]['transcript_score'] = TSCORE

        gene[g_cnt]['gene_info'] = dict( ID = pkey[-1], 
                                Name = pdet.get('name'), 
                                Source = pkey[1]) 
        # few empty fields // TODO fill this:
        gene[g_cnt]['anno_id'] = []
        gene[g_cnt]['confgenes_id'] = []
        gene[g_cnt]['alias'] = ''
        gene[g_cnt]['name2'] = []
        gene[g_cnt]['chr_num'] = []
        gene[g_cnt]['paralogs'] = []
        gene[g_cnt]['transcript_valid'] = []
        gene[g_cnt]['exons_confirmed'] = []
        gene[g_cnt]['tis_conf'] = []
        gene[g_cnt]['tis_info'] = []
        gene[g_cnt]['cdsStop_conf'] = []
        gene[g_cnt]['cdsStop_info'] = []
        gene[g_cnt]['tss_info'] = []
        gene[g_cnt]['tss_conf'] = []
        gene[g_cnt]['cleave_info'] = []
        gene[g_cnt]['cleave_conf'] = []
        gene[g_cnt]['polya_info'] = []
        gene[g_cnt]['polya_conf'] = []
        gene[g_cnt]['is_valid'] = []
        gene[g_cnt]['transcript_complete'] = []
        gene[g_cnt]['is_complete'] = []
        gene[g_cnt]['is_correctly_gff3_referenced'] = ''
        gene[g_cnt]['splicegraph'] = []
        g_cnt += 1 

    ## deleting empty gene records from the main array
    XPFLG=0
    for XP, ens in enumerate(gene):
        if ens[0]==0:
            XPFLG=1
            break
    
    if XPFLG==1:
        XQC = range(XP, len(gene)+1)
        gene = np.delete(gene, XQC)

    return gene 

def NonetoemptyList(XS):
    """
    Convert a None type to empty list 
    @args XS: None type 
    @type XS: str 
    """
    return [] if XS is None else XS 

def create_missing_feature_type(p_feat, c_feat):
    """
    GFF/GTF file defines only child features. This function tries to create 
    the parent feature from the information provided in the attribute column. 
    example: 
    chr21   hg19_knownGene  exon    9690071 9690100 0.000000        +       .       gene_id "uc002zkg.1"; transcript_id "uc002zkg.1"; 
    chr21   hg19_knownGene  exon    9692178 9692207 0.000000        +       .       gene_id "uc021wgt.1"; transcript_id "uc021wgt.1"; 
    chr21   hg19_knownGene  exon    9711935 9712038 0.000000        +       .       gene_id "uc011abu.2"; transcript_id "uc011abu.2"; 
    This function gets the parsed feature annotations. 
    
    @args p_feat: Parent feature map  
    @type p_feat: collections defaultdict
    @args c_feat: Child feature map  
    @type c_feat: collections defaultdict
    """

    child_n_map = defaultdict(list)
    for fid, det in c_feat.items():
        # get the details from grand child  
        GID = STRD = SCR = None
        SPOS, EPOS = [], [] 
        TYP = dict()
        for gchild in det:
            GID = gchild.get('gene_id', '') 
            GNAME = gchild.get('gene_name', '') 
            SPOS.append(gchild.get('location', [])[0]) 
            EPOS.append(gchild.get('location', [])[1]) 
            STRD = gchild.get('strand', '')
            SCR = gchild.get('score', '')
            if gchild.get('type', '') == "gene": ## gencode GTF file has this problem 
                continue
            TYP[gchild.get('type', '')] = 1
        SPOS.sort() 
        EPOS.sort()
        
        # infer transcript type
        transcript_type = 'transcript'
        transcript_type = 'mRNA' if TYP.get('CDS', '') or TYP.get('cds', '') else transcript_type
        
        # gene id and transcript id are same
        transcript_id = fid[-1]
        if GID == transcript_id:
            transcript_id = 'Transcript:' + str(GID)
        
        # level -1 feature type 
        p_feat[(fid[0], fid[1], GID)] = dict( type = 'gene',
                                            location = [], ## infer location based on multiple transcripts  
                                            strand = STRD,
                                            name = GNAME )
        # level -2 feature type 
        child_n_map[(fid[0], fid[1], GID)].append(
                                            dict( type = transcript_type,
                                            location =  [SPOS[0], EPOS[-1]], 
                                            strand = STRD, 
                                            score = SCR, 
                                            ID = transcript_id,
                                            gene_id = '' ))
        # reorganizing the grand child
        for gchild in det:
            child_n_map[(fid[0], fid[1], transcript_id)].append(
                                            dict( type = gchild.get('type', ''),
                                            location =  gchild.get('location'),
                                            strand = gchild.get('strand'), 
                                            ID = gchild.get('ID'),
                                            score = gchild.get('score'),
                                            gene_id = '' ))
    return p_feat, child_n_map 