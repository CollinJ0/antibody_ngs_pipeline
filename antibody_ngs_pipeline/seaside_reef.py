#!/usr/bin/env python
# filename: seaside_reef.py


import os
import sys


from datetime import datetime
from shutil import copyfile
import time
try:
    from abutils.utils import log
    from abutils.utils.pipeline import make_dir
    from abutils.utils.progbar import progress_bar
except ImportError:
    print('You need abutils')





#These Params last upodated 7/30/2018
ABSTAR_PARAMS       = {'assigner': ['blastn'],
                       'output_type': [['json'], 'imgt', 'minimal', 'json'],
                       'merge' : {'TRUE': True, 'FALSE': False}
                      }

#Copying From Basemount
def copy_from_basemount(basemount_directory, destination_directory):
    make_dir(os.path.abspath(destination_directory))
    fastqs = []
    for (path, dirs, files) in os.walk(basemount_directory):
        for f in files:
            if f.endswith('.fastq.gz'):
                fastqs.append(os.path.join(path, f))
    
    print('')
    print('')
    print('========================================')
    print('Copying files from BaseMount')
    print('========================================')
    print('')
    print('Found {0} FASTQ files.'.format(len(fastqs)))
    print('')
    print('Copying to {}:'.format(destination_directory))
    start = datetime.now()
    progress_bar(0, len(fastqs), start_time=start, completion_string = 'ITS DONE! good job')
    for i, fastq in enumerate(fastqs):
        dest = os.path.join(destination_directory, os.path.basename(fastq))
        copyfile(fastq, dest)
        progress_bar(i + 1, len(fastqs), start_time=start)
    print('\n')
    
    

def print_splash():
    splash1 = '''
    _          _   _ _               _
   / \   _ __ | |_(_) |__   ___   __| |_   _ 
  / _ \ | '_ \| __| | '_ \ / _ \ / _` | | | |
 / ___ \| | | | |_| | |_) | (_) | (_| | |_| |
/_/   \_\_| |_|\__|_|_.__/ \___/ \__,_|\__, |
                                       |___/ '''
    splash3 = '''
 _____  
|_____|'''
    splash4 = '''
 _____ _____ 
|_____|_____|'''
    splash5 = '''
 _____ _____ _____  
|_____|_____|_____|'''
    splash6 = '''
 _____ _____ _____ _____  
|_____|_____|_____|_____|'''
    splash7 = '''
 _____ _____ _____ _____ _____  
|_____|_____|_____|_____|_____|'''
    splash8 = '''
 _____ _____ _____ _____ _____ _____  
|_____|_____|_____|_____|_____|_____|'''
    splash9 = '''
 _____ _____ _____ _____ _____ _____ _____  
|_____|_____|_____|_____|_____|_____|_____|'''
    splash10 = '''
 _____ _____ _____ _____ _____ _____ _____ _____  
|_____|_____|_____|_____|_____|_____|_____|_____|'''
    
    splash11 = '''
 _____ _____ _____ _____ _____ _____ _____ _____ _____ 
|_____|_____|_____|_____|_____|_____|_____|_____|_____|'''

                                                                              
    splash2 = '''
 _   _  ____ ____    ____  _            _ _
| \ | |/ ___/ ___|  |  _ \(_)_ __   ___| (_)_ __   ___
|  \| | |  _\___ \  | |_) | | '_ \ / _ \ | | '_ \ / _ \\
| |\  | |_| |___) | |  __/| | |_) |  __/ | | | | |  __/
|_| \_|\____|____/  |_|   |_| .__/ \___|_|_|_| |_|\___|
                            |_|                     '''

    # print('')
    print(splash1)
    print(splash2)
    print(splash11)
    print('')
    print('')