fastqc = {'outdir': '.',
          'quiet' : True}

trim_galore = {'paired': True, 
               'gzip'  : True}

star_indexing = {'sjdbOverhang': 149}

star_mapping = {'twopassMode'     : 'Basic',
                'outWigType'      : 'bedGraph',
                'outSAMtype'      : 'BAM SortedByCoordinate',
                'readFilesCommand': 'gunzip -c'}
