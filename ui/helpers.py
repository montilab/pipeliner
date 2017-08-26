import os
import csv

def createFiles():
    return {'annotation-files' : [],
            'reference-files'  : [],
            'reads-files'      : [],
            'alignment-files'  : []}

def searchFiles(path, exts):
    directory = [f for f in os.listdir(path) if not f.startswith('.')]
    return [f for f in directory if f.lower().endswith(tuple(exts))]

def parseReads(pathtocsv):
  with open(pathtocsv, 'r') as infile:
    samples = [row for row in csv.reader(infile)][1:]
    for sample in samples:
      sample[1] = sample[1].split('/')[-1]
      try:
        sample[2] = sample[2].split('/')[-1]
      except IndexError:
        sample.append('')    
    return samples

def createConfig(indir, outdir, files, settings):
    project = "config_test"
    indir = "/Users/defna/pythonvillage/pypliner3/Gallus_example/ggal_data"
    outdir = "/Users/defna/pythonvillage/pypliner3/Gallus_example/results"
    files = {'annotation-files': ['genome_annotation.gff'], 'reference-files': ['genome_reference.fa'], 'reads-files': ['ggal_alpha.csv'], 'alignment-files': []}
    settings = {'aligner': 'star', 'star_index': False, 'paired': True, 'save_reference': True}

    config = 'params {{\n\tindir  = "{0}"\n\toutdir = "{1}"\n'.format(indir,outdir)
    config += '\n\tfasta = "${{params.indir}}/{0}"'.format(files['reference-files'][0])
    config += '\n\tgtf   = "${{params.indir}}/{0}"'.format(files['annotation-files'][0])
    config += '\n\treads = "${{params.indir}}/{0}"'.format(files['reads-files'][0])

    config += "\n"
    for setting in settings:
        if type(settings[setting]) == str:
            config += '\n\t{0} = "{1}"'.format(setting, settings[setting])
        if type(settings[setting]) == bool: 
            config += '\n\t{0} = {1}'.format(setting, str(settings[setting]).lower())
    config += "\n}"
    
    print(config)

    #with open('{0}.config'.format(project), 'w') as outfile:
    #   outfile.write(config)


#createConfig(1,1,1,1)

def createNextflow(nfdir, pipeline, outfile, env=None, resuming=False):
    with open(outfile, 'w') as writer:
        writer.write('#!/bin/bash\n')
        if env is not None:
            writer.write('source activate {0}\n'.format(env))
        writer.write('cd ~\n')
        writer.write('cd {0}\n'.format(nfdir))
        writer.write('./nextflow {0}'.format(pipeline))
        if resuming:
            writer.write(' -resume')

env = 'pypliner3'
nfdir = '/home/feds/Documents/pythonvillage/pypliner3/Gallus_example'
pipeline = 'main.nf'

#createNextflow(nfdir, pipeline, outfile='start.sh')
