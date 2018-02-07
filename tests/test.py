#!/usr/bin/env python

if __name__ == '__main__':
    import subprocess
    import time
    import os
    
    path_to_tests = os.path.realpath(__file__)
    path_to_tests = '/'.join(path_to_tests.split('/')[:-1])
    path_to_configs = '{0}/configs'.format(path_to_tests)

    for file in os.listdir(path_to_configs):
        if file.endswith(".config"):
            path_to_file = os.path.join(path_to_configs, file)
            test = file.strip('.config')
            print('_'*40)
            print('Test: {0}'.format(test))
            print('Date: {0}'.format(time.strftime("%Y/%m/%d %H:%M:%S")))
            print('Status: Running...')
            subprocess.Popen(['mkdir', '-p', '{0}/logs'.format(path_to_tests)])
            subprocess.Popen(['rm', '-f', '--', '{0}/logs/{1}.txt'.format(path_to_tests, test)])
            subprocess.Popen(['./nextflow', 'rnaseq.nf', '-c', path_to_file])
