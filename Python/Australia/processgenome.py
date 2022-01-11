'''
Run trimgalore and jellyfish
'''

import os
import subprocess
import shlex
from time import time
# Setings
# Jellyfisk k size
jfk='11'
# Jellyfish threads
threads='10'
# Jellyfish elements in hash
jfhash='8G'

def trim_galore(
    input_files, trimmed_path=None, paired=False, dontzip=False, cores=None):
    '''Runs trim galore on file. Quality cutoff 30, remaining
    settings default.
    
    Args:
    input_files - (lst) List of files to be trimmed
    trimmed_path - (str, default=None) Name of output file, if not supplied 
        uses trimgalores default naming scheme
    paired - (bol) Adds paired end quality control
    '''
    
    start_time = time()

    #Test 2 cores
    cores = 2
    trimgalore_command = './trim_galore -q 30'
    if dontzip:
        trimgalore_command += ' --dont_gzip'
    if paired:
        trimgalore_command += ' --paired'
    if cores is not None:
        trimgalore_command += f' --cores {cores}'
    if trimmed_path is not None:
        trimgalore_command += f' -o {trimmed_path}'
    if isinstance(input_files, str):
        input_files = [input_files]
    for file in input_files:
        trimgalore_command += f' {file}'

    try:
        subprocess.run(shlex.split(trimgalore_command), check=True)
    except subprocess.CalledProcessError as e:
        print('Error running trim_galore command')
        return 1
    
    print(f'Running time trim_galore: {(time()-start_time)/60} min')
    return 0

def jellyfish(input_files, output_file):
    ''' Runs Jellyfish on file 
    
    Args:
    input_files - List of zipped fasta files to be counted
    output_file - Filename to dump to
    '''
    start_time = time()

    tmp_file = 'processGenomeKmercount.tmp'

    count_command = f'jellyfish count -m {jfk} -s {jfhash}'
    count_command += f' -t {threads} -C -o {tmp_file}'
    count_command += f' <(zcat {input_files[0]}) <(zcat {input_files[1]})'
    
    print(count_command)
    try:
        subprocess.run(count_command, check=True , shell=True,
        executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        print('Error running Jellyfish count command')
        print(e)
        return 1

    try:   
        with open(output_file, 'w') as fh:
            subprocess.run(['jellyfish', 'dump', tmp_file], stdout=fh, check=True)
    except subprocess.CalledProcessError as e:
        print('Error running Jellyfish dump command')
        print(e)
        return 1
    

    print(f'Running time for jellyfish: {(time()-start_time)/60} min')
    return 0