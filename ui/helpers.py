import os

filexts = {'annotation-files' : ['.gtf', '.gff'],
           'reference-files'  : ['.fasta', '.fa'],
           'csv-files'        : ['.csv'],
           'alignment-files'  : ['.sam', '.bam']}

def searchFiles(path, exts):
    directory = [f for f in os.listdir(path) if not f.startswith('.')]
    return [f for f in directory if f.lower().endswith(tuple(exts))]

def searchReads(csv):
    try:
        with open(csv, 'r') as infile:
            samples = [row.split(',') for row in infile][1:]
            for sample in samples:
                sample[1] = sample[1].split('/')[-1]
                try:
                    sample[2] = sample[2].split('/')[-1]
                except IndexError:
                    sample.append('')    
            return samples
    except IndexError:
        return []

def createConfig(indir, outdir, files, settings, nfdir=''):
    project = "nextflow"
    config = 'params {{\n\tindir  = "{0}"\n\toutdir = "{1}"\n'.format(indir,outdir)
    config += '\n\tfasta = "${{params.indir}}/{0}"'.format(files['reference-files'][0])
    config += '\n\tgtf   = "${{params.indir}}/{0}"'.format(files['annotation-files'][0])
    config += '\n\treads = "${{params.indir}}/{0}"'.format(files['csv-files'][0])
    config += "\n"
    for setting in settings:
        if type(settings[setting]) == str:
            config += '\n\t{0} = "{1}"'.format(setting, settings[setting])
        if type(settings[setting]) == bool: 
            config += '\n\t{0} = {1}'.format(setting, str(settings[setting]).lower())
    config += "\n}"

    if nfdir == '':
        path = '{0}.config'.format(project)
    else:
        path = '{0}/{1}.config'.format('/'.join(nfdir.split('/')[:-1]), project)

    with open(path, 'w') as outfile:
       outfile.write(config)

def createNextflow(nfdir, pipeline, env='', resuming=False):
    with open('start.sh', 'w') as writer:
        writer.write('#!/bin/bash\n')
        if env is not '':
            writer.write('source activate {0}\n'.format(env))
        writer.write('cd ~\n')
        writer.write('cd {0}\n'.format('/'.join(nfdir.split('/')[:-1])))
        writer.write('./nextflow {0} -c nextflow.config'.format(pipeline))
        if resuming:
            writer.write(' -resume')