#!/usr/bin/env python

if __name__ == '__main__':
    import subprocess
    import os
    
    print('Setting paths for RNA-seq pipeline')

    '''
    Copies the config file with new local indir and outdir
    '''
    path_to_pipeliner = os.path.realpath(__file__)
    path_to_pipeliner = '/'.join(path_to_pipeliner.split('/')[:-2])
    print('Updating local paths in {0}/pipelines/rnaseq.config...'.format(path_to_pipeliner))
    with open('{0}/pipelines/temp.config'.format(path_to_pipeliner), 'w') as outfile:
        with open('{0}/pipelines/rnaseq.config'.format(path_to_pipeliner), 'r') as infile:
            for line in infile:
                if line.lstrip().startswith('indir'):
                    outfile.write('  indir  = "{0}/pipelines/toy_data/rnaseq"\n'.format(path_to_pipeliner))
                elif line.lstrip().startswith('outdir'): 
                    outfile.write('  outdir = "{0}/pipelines/rnaseq_results"\n'.format(path_to_pipeliner))
                else:
                    outfile.write(line)

    subprocess.call(['rm', '-rf', '{0}/pipelines/rnaseq.config'.format(path_to_pipeliner)])
    subprocess.call(['mv', '{0}/pipelines/temp.config'.format(path_to_pipeliner),
                           '{0}/pipelines/rnaseq.config'.format(path_to_pipeliner)])

    '''
    Updates reads.csv with local paths
    '''
    print('Updating local paths in {0}/pipelines/toy_data/rnaseq/ggal_reads.csv...'.format(path_to_pipeliner))
    subprocess.call(['rm', '-rf', '{0}/pipelines/toy_data/rnaseq/ggal_reads.csv'.format(path_to_pipeliner)])
    with open('{0}/pipelines/toy_data/rnaseq/ggal_reads.csv'.format(path_to_pipeliner), 'w') as outfile:
        outfile.write('Sample_Name,Read1,Read2\n')
        for sample in ['ggal_alpha', 'ggal_theta', 'ggal_gamma']:
            outfile.write('{0},{1}/pipelines/toy_data/rnaseq/reads/{0}_1.fq.gz,{1}/pipelines/toy_data/rnaseq/reads/{0}_2.fq.gz\n'.format(sample, path_to_pipeliner))

    '''
    Updates bams.csv with local paths
    '''
    print('Updating local paths in {0}/pipelines/toy_data/rnaseq/ggal_bams.csv...'.format(path_to_pipeliner))
    subprocess.call(['rm', '-rf', '{0}/pipelines/toy_data/rnaseq/ggal_bams.csv'.format(path_to_pipeliner)])
    with open('{0}/pipelines/toy_data/rnaseq/ggal_bams.csv'.format(path_to_pipeliner), 'w') as outfile:
        outfile.write('Sample_Name,Alignment\n')
        for sample in ['ggal_alpha', 'ggal_theta', 'ggal_gamma']:
            outfile.write('{0},{1}/pipelines/toy_data/rnaseq/bams/{0}.bam\n'.format(sample, path_to_pipeliner))

    
    print('Setting paths for scRNA-seq pipeline')

    '''
    Copies the config file with new local indir and outdir
    '''
    path_to_pipeliner = os.path.realpath(__file__)
    path_to_pipeliner = '/'.join(path_to_pipeliner.split('/')[:-2])
    print('Updating local paths in {0}/pipelines/scrnaseq.config...'.format(path_to_pipeliner))
    with open('{0}/pipelines/temp.config'.format(path_to_pipeliner), 'w') as outfile:
        with open('{0}/pipelines/scrnaseq.config'.format(path_to_pipeliner), 'r') as infile:
            for line in infile:
                if line.lstrip().startswith('indir'):
                    outfile.write('  indir  = "{0}/pipelines/toy_data/scrnaseq"\n'.format(path_to_pipeliner))
                elif line.lstrip().startswith('outdir'): 
                    outfile.write('  outdir = "{0}/pipelines/scrnaseq_results"\n'.format(path_to_pipeliner))
                else:
                    outfile.write(line)

    subprocess.call(['rm', '-rf', '{0}/pipelines/scrnaseq.config'.format(path_to_pipeliner)])
    subprocess.call(['mv', '{0}/pipelines/temp.config'.format(path_to_pipeliner),
                           '{0}/pipelines/scrnaseq.config'.format(path_to_pipeliner)])

    '''
    Updates reads.csv with local paths
    '''
    print('Updating local paths in {0}/pipelines/toy_data/scrnaseq/hgmm_reads.csv...'.format(path_to_pipeliner))
    subprocess.call(['rm', '-rf', '{0}/pipelines/toy_data/scrnaseq/hgmm_reads.csv'.format(path_to_pipeliner)])
    with open('{0}/pipelines/toy_data/scrnaseq/hgmm_reads.csv'.format(path_to_pipeliner), 'w') as outfile:
        outfile.write('Sample_Name,Read1,Read2\n')
        for sample in ['hgmm_beta']:
            outfile.write('{0},{1}/pipelines/toy_data/scrnaseq/reads/{0}_1.fq.gz,{1}/pipelines/toy_data/scrnaseq/reads/{0}_2.fq.gz\n'.format(sample, path_to_pipeliner))

    '''
    Updates bams.csv with local paths
    '''
    print('Updating local paths in {0}/pipelines/toy_data/scrnaseq/hgmm_bams.csv...'.format(path_to_pipeliner))
    subprocess.call(['rm', '-rf', '{0}/pipelines/toy_data/scrnaseq/hgmm_bams.csv'.format(path_to_pipeliner)])
    with open('{0}/pipelines/toy_data/scrnaseq/hgmm_bams.csv'.format(path_to_pipeliner), 'w') as outfile:
        outfile.write('Sample_Name,Alignment\n')
        for sample in ['hgmm_beta']:
            outfile.write('{0},{1}/pipelines/toy_data/scrnaseq/bams/{0}.bam\n'.format(sample, path_to_pipeliner))

    print('Setting paths for DGE pipeline')

    '''
    Copies the config file with new local indir and outdir
    '''
    path_to_pipeliner = os.path.realpath(__file__)
    path_to_pipeliner = '/'.join(path_to_pipeliner.split('/')[:-2])
    print('Updating local paths in {0}/pipelines/dge.config...'.format(path_to_pipeliner))
    with open('{0}/pipelines/temp.config'.format(path_to_pipeliner), 'w') as outfile:
        with open('{0}/pipelines/dge.config'.format(path_to_pipeliner), 'r') as infile:
            for line in infile:
                if line.lstrip().startswith('indir'):
                    outfile.write('  indir  = "{0}/pipelines/toy_data/dge"\n'.format(path_to_pipeliner))
                elif line.lstrip().startswith('outdir'): 
                    outfile.write('  outdir = "{0}/pipelines/dge_results"\n'.format(path_to_pipeliner))
                else:
                    outfile.write(line)

    subprocess.call(['rm', '-rf', '{0}/pipelines/dge.config'.format(path_to_pipeliner)])
    subprocess.call(['mv', '{0}/pipelines/temp.config'.format(path_to_pipeliner),
                           '{0}/pipelines/dge.config'.format(path_to_pipeliner)])

    '''
    Updates bams.csv with local paths
    '''
    print('Updating local paths in {0}/pipelines/toy_data/dge/dge_bams.csv...'.format(path_to_pipeliner))
    subprocess.call(['rm', '-rf', '{0}/pipelines/toy_data/dge/dge_bams.csv'.format(path_to_pipeliner)])
    with open('{0}/pipelines/toy_data/dge/dge_bams.csv'.format(path_to_pipeliner), 'w') as outfile:
        outfile.write('Sample_Name,Alignment\n')
        for sample in ['A1_marked', 'B2_marked', 'C3_marked']:
            outfile.write('{0},{1}/pipelines/toy_data/dge/bams/{0}.bam\n'.format(sample, path_to_pipeliner))

    print('Ready to run test data!')
