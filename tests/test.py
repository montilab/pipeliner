#!/usr/bin/env python

if __name__ == '__main__':
    import subprocess
    import time
    import sys
    import os

    if len(sys.argv) != 2:
        print('Please specify a pipeline to test')
        print('E.g. python3 path/to/tests/test.py rnaseq')
        sys.exit()

    pipeline = sys.argv[1]
    
    path_to_tests = os.path.realpath(__file__)
    path_to_tests = '/'.join(path_to_tests.split('/')[:-1])
    path_to_configs = '{0}/configs/{1}'.format(path_to_tests, pipeline)

    omit = []

    for file in os.listdir(path_to_configs):
        if file.endswith(".config") and file not in omit:

            path_to_file = os.path.join(path_to_configs, file)

            with open(path_to_file, 'r') as infile:
                test = infile.readline().lstrip('// ').rstrip('\n')

            print('\n')
            print('Test: {0}'.format(test))
            print('Date: {0}'.format(time.strftime("%Y/%m/%d %H:%M:%S")))
            print('Status: Running...')

            subprocess.call(['./nextflow', '{0}.nf'.format(pipeline), '-c', path_to_file])
    
    print('\nTesting complete...')
    subprocess.call("exit 1", shell=True)