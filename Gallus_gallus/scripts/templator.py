#!/usr/bin/env python3
# @Author: Anthony Federico <anfederico>
# @Last modified by: anfederico

def get_params(config):
    with open(config, 'r') as config:
        params = {}
        for line in config:
            line = line.strip('\n')
            line = "".join(line.split())
            line = line.split("=")
            if len(line) > 1:
                params[line[0]] = line[1].strip('"')
        return params

def add_kwargs(command, kwargs):
    for key, val in kwargs.items():
        if val == True:
            command += (' --{0}'.format(key))
        else:
            command += (' --{0} {1}'.format(key, val))
    return command

def write_template(template, command):
    if not os.path.isdir('templates'):
        subprocess.call(['mkdir', 'templates'])
    with open('templates/{0}'.format(template), 'w') as script:
        script.write(command) 

def main():
    # Fastqc -----------------------------------------------------------------------
    command = "fastqc"
    command = add_kwargs(command, settings.fastqc)
    if params['paired'] == 'true':
        command += (' ${reads[0]} ${reads[1]}')
    else:
        command += (' ${reads[0]}')
    write_template('fastqc.sh', command)

    # Trim Galore ------------------------------------------------------------------
    command = "trim_galore"
    command = add_kwargs(command, settings.trim_galore)
    if params['paired'] == 'true':
        command += (' ${reads[0]} ${reads[1]}')
    else:
        command += (' ${reads[0]}')
    write_template('trim_galore.sh', command)

    # Star Indexing ----------------------------------------------------------------
    command = "mkdir star_index\nSTAR --runMode genomeGenerate --runThreadN ${task.cpus} --sjdbGTFfile $gtf --genomeDir star_index/ --genomeFastaFiles $fasta" 
    command = add_kwargs(command, settings.star_indexing)
    write_template('star_indexing.sh', command)

    # Star Mapping -----------------------------------------------------------------
    command = "STAR --genomeDir $index --sjdbGTFfile $gtf --readFilesIn $reads --runThreadN ${task.cpus} --outFileNamePrefix \'${sampleid}.'"
    command = add_kwargs(command, settings.star_mapping)
    write_template('star_mapping.sh', command)

if __name__ == "__main__":
    import os
    import sys
    import subprocess
    import settings

    params = get_params(sys.argv[1])

    main()