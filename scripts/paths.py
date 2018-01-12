#!/usr/bin/env python

if __name__ == '__main__':
    import subprocess
    import os
    
    '''
    Copies the config file with new local indir and outdir
    '''

    path_to_pipeliner = os.path.realpath(__file__)
    path_to_pipeliner = '/'.join(path_to_pipeliner.split('/')[:-2])
    print('Updating local paths in {0}/RNA-seq/nextflow.config...'.format(path_to_pipeliner))
    with open('{0}/RNA-seq/temp.config'.format(path_to_pipeliner), 'w') as outfile:
        with open('{0}/RNA-seq/nextflow.config'.format(path_to_pipeliner), 'r') as infile:
            for line in infile:
                if line.lstrip().startswith('indir'):
                    outfile.write('  indir  = "{0}/RNA-seq/ggal_data"\n'.format(path_to_pipeliner))
                elif line.lstrip().startswith('outdir'): 
                    outfile.write('  outdir = "{0}/RNA-seq/ggal_results"\n'.format(path_to_pipeliner))
                else:
                    outfile.write(line)

    subprocess.call(['rm', '-rf', '{0}/RNA-seq/nextflow.config'.format(path_to_pipeliner)])
    subprocess.call(['mv', '{0}/RNA-seq/temp.config'.format(path_to_pipeliner),
                           '{0}/RNA-seq/nextflow.config'.format(path_to_pipeliner)])

    '''
    Updates reads.csv with local paths
    '''
    print('Updating local paths in {0}/RNA-seq/ggal_data/ggal_reads.csv...'.format(path_to_pipeliner))
    subprocess.call(['rm', '-rf', '{0}/RNA-seq/ggal_data/ggal_reads.csv'.format(path_to_pipeliner)])
    with open('{0}/RNA-seq/ggal_data/ggal_reads.csv'.format(path_to_pipeliner), 'w') as outfile:
        outfile.write('Sample_Name,Read1,Read2\n')
        for sample in ('ggal_alpha', 'ggal_theta', 'ggal_gamma'):
            outfile.write('{0},{1}/RNA-seq/ggal_data/reads/{0}_1.fq.gz,{1}/RNA-seq/ggal_data/reads/{0}_2.fq.gz\n'.format(sample, path_to_pipeliner))

    '''
    Updates alignments.csv with local paths
    '''
    print('Updating local paths in {0}/RNA-seq/ggal_data/ggal_alignments.csv...'.format(path_to_pipeliner))
    subprocess.call(['rm', '-rf', '{0}/RNA-seq/ggal_data/ggal_alignments.csv'.format(path_to_pipeliner)])
    with open('{0}/RNA-seq/ggal_data/ggal_alignments.csv'.format(path_to_pipeliner), 'w') as outfile:
        outfile.write('Sample_Name,Alignment\n')
        for sample in ('ggal_alpha', 'ggal_theta', 'ggal_gamma'):
            outfile.write('{0},{1}/RNA-seq/ggal_data/alignments/{0}.bam\n'.format(sample, path_to_pipeliner))

    print('Ready to run test data!')