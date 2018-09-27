#!/usr/bin/env python
# filename: horseshoe.py

import sys
import os
import warnings
import shutil
import time
from getpass import getpass
from argparse import ArgumentParser

from abstar.core.abstar import Args
from abstar.assigners.registry import ASSIGNERS
from abstar.utils import mongoimport
from abstar import run_standalone, fastqc, adapter_trim, quality_trim
from abutils.utils.pipeline import make_dir
from abutils.utils.progbar import progress_bar



from .seaside_reef import ABSTAR_PARAMS, copy_from_basemount, print_splash#, MONGO_PARAMS, S3_PARAMS


def parse_arguments(print_help=False):
    parser = ArgumentParser(prog='antibody_ngs_pipeline', description="Bulk antibody sequence preprocessing, annotation with abstar, upload to MongoDB and S3")
    parser.add_argument('-f', '--fastqc', action='store_true', dest='fastqc', default=False,
                        help="Use to generate a FASTQC quality report on raw data, if used with quality trimming \
                             '-q' a FASTQC quality report will be made for both pre and post adapter trimmed data.")
    parser.add_argument('-t', '--adapter-trim', dest='adapter_fasta', default=None,
                        help="Adapter trimming using CutAdapt, if this flag is used, must specify \
                              the location of a fasta file which contains adapter sequences for both ends.")
    parser.add_argument('-q', '--quality-trim', action='store_true', dest='quality_trim', default=False,
                        help="Quality trimming using Sickle, if this flag is used, must specify \
                              the location of a fasta file which contains adapter sequences for both ends.")
    if print_help:
        parser.print_help()
    else:
        args = parser.parse_args()
        return args

######################################################
#                      Abstar                        #
######################################################


def abstar_params(project):
    directory = os.path.join('/data', project)
    proj_dir = check_dir(directory)
    
    parameters = Args(project_dir=proj_dir, merge=True)
    print("\n========================================" \
          "\nAbstar Run Arguments\n" \
          "========================================\n")
    if print_abstar_params(parameters):
        parameters = change_abstar_params(parameters)
    return parameters

def preprocess(args, pipeline_args):
    fastqc_raw_output = os.path.join(args.project_dir, 'fastqc_raw')
    original_input = os.path.join(args.project_dir, 'input')
    adapter_output = os.path.join(args.project_dir, 'adapter_trimmed')
    quality_output = os.path.join(args.project_dir, 'quality_trimmed')
    fastqc_trimmed_output = os.path.join(args.project_dir, 'fastqc_trimmed')
    if pipeline_args.fastqc and pipeline_args.adapter_fasta == None and not pipeline_args.quality_trim:
        print("\n========================================" \
              "\nFASTQC only specified, Running FASTQC on raw data...")
        fastqc(original_input, output_directory=fastqc_raw_output)
        print("\nFastQC report on raw fastqs available in {}\n" \
              "========================================\n".format(fastqc_raw_output))
        return args
    elif pipeline_args.fastqc and pipeline_args.quality_trim and pipeline_args.adapter_fasta == None:
        print("\n========================================" \
              "\nFASTQC and quality trimming specified,\n" \
              "Running FASTQC on raw data...\n")
        fastqc(original_input, output_directory=fastqc_raw_output)
        print("\nFastQC report on raw fastqs available in {}\n".format(fastqc_raw_output))
        print("\nQuality trimming raw fastqs with Sickle...\n")
        quality_trim(original_input, output_directory=quality_output)
        print("\nRunning FASTQC on quality trimmed fastqs...\n")
        fastqc(quality_output, output_directory=fastqc_trimmed_output)
        print("\nFastQC report on quality trimmed fastqs available in {}".format(fastqc_trimmed_output))
        print("========================================\n")
        args.input = quality_output
        args.output = os.path.join(args.project_dir, 'output')
        args.temp = os.path.join(args.project_dir, 'temp')
        args.project_dir = None
        return args
    elif not pipeline_args.fastqc and not pipeline_args.quality_trim and pipeline_args.adapter_fasta != None:
        print("\n========================================" \
              "\nAdapter trimming only specified." \
              "\nTrimming adapters with CutAdapt...")
        #Try adapter trim, if a non valid file of adapters has been passed through, sys.exit()
        adapter_trim(original_input, output_directory=adapter_output, adapter_both=pipeline_args.adapter_fasta)
        print("========================================\n")
        args.input = adapter_output
        args.output = os.path.join(args.project_dir, 'output')
        args.temp = os.path.join(args.project_dir, 'temp')
        args.project_dir = None       
        return args
    elif not pipeline_args.fastqc and pipeline_args.quality_trim and pipeline_args.adapter_fasta == None:
        print("\n========================================" \
              "\nQuality trimming only specified." \
              "\nQuality trimming raw fastqs with Sickle...")
        quality_trim(original_input, output_directory=quality_output)
        args.input = quality_output
        args.output = os.path.join(args.project_dir, 'output')
        args.temp = os.path.join(args.project_dir, 'temp')
        args.project_dir = None
        return args
    elif pipeline_args.fastqc and not pipeline_args.quality_trim and pipeline_args.adapter_fasta != None:       
        print("\n========================================" \
              "\nFASTQC and adapter trimming specified,\n" \
              "Running FASTQC on raw data...\n")
        fastqc(original_input, output_directory=fastqc_raw_output)
        print("\nFastQC report on raw fastqs available in {}\n".format(fastqc_raw_output))
        print("\nAdapter trimming raw fastqs with CutAdapt...\n")
        adapter_trim(original_input, output_directory=adapter_output, adapter_both=pipeline_args.adapter_fasta)
        print("\nRunning FASTQC on adapter trimmed fastqs...\n")
        fastqc(adapter_output, output_directory=fastqc_trimmed_output)
        print("\nFastQC report on adapter trimmed fastqs available in {}".format(fastqc_trimmed_output))
        print("========================================\n")
        args.input = adapter_output
        args.output = os.path.join(args.project_dir, 'output')
        args.temp = os.path.join(args.project_dir, 'temp')
        args.project_dir = None
        return args
    elif pipeline_args.fastqc and pipeline_args.quality_trim and pipeline_args.adapter_fasta != None:       
        print("\n========================================" \
              "\nFASTQC and adapter trimming and quality trimming specified,\n" \
              "Running FASTQC on raw data...\n")
        fastqc(original_input, output_directory=fastqc_raw_output)
        print("\nFastQC report on raw fastqs available in {}\n".format(fastqc_raw_output))
        print("\nAdapter trimming raw fastqs with CutAdapt...\n")
        adapter_trim(original_input, output_directory=adapter_output, adapter_both=pipeline_args.adapter_fasta)
        print("Quality trimming raw fastqs with Sickle...")
        quality_trim(adapter_output, output_directory=quality_output)
        print("\nRunning FASTQC on adapter and quality trimmed fastqs...\n")
        fastqc(quality_output, output_directory=fastqc_trimmed_output)
        print("\nFastQC report on adapter and quality trimmed fastqs available in {}".format(fastqc_trimmed_output))
        print("========================================\n")
        args.input = quality_output
        args.output = os.path.join(args.project_dir, 'output')
        args.temp = os.path.join(args.project_dir, 'temp')
        args.project_dir = None
        return args
    elif not pipeline_args.fastqc and pipeline_args.quality_trim and pipeline_args.adapter_fasta != None:       
        print("\n========================================" \
              "\nAdapter trimming and quality trimming specified")
        print("\nAdapter trimming raw fastqs with CutAdapt...")
        adapter_trim(original_input, output_directory=adapter_output, adapter_both=pipeline_args.adapter_fasta)
        print("Quality trimming raw fastqs with Sickle...")
        quality_trim(adapter_output, output_directory=quality_output)
        print("Done")
        print("========================================\n")
        args.input = quality_output
        args.output = os.path.join(args.project_dir, 'output')
        args.temp = os.path.join(args.project_dir, 'temp')
        args.project_dir = None
        return args
    else:
        return args
    return args
    
def check_dir(directory):
    if not os.path.exists(os.path.abspath(directory)):
        return directory
    else:
        print('\n******************************************'
              '\nYour Project directory already exists, please specify new directory' \
              '\nor leave blank to keep current directory listed' \
              '\n(Warning: all files in that directory will be deleted)\n ')
        direct = input("Project Directory ({}): ".format(directory))
        if len(direct) > 0:
            return direct
        else:
            if len(directory.split('/data')[1]) > 1:
                return directory
            else:
                print('No directory specified')
                sys.exit(1)

def change_abstar_params(parameters):
    print('\n========================================' \
          '\nTo change a current abstar argument,' \
          '\ntype in new argument and press enter.' \
          '\nLeave blank and press enter if you want ' \
          '\nto keep current parameter (in parentheses).'
          '\n========================================')
    
    proj=input("Project Directory ({}):".format(parameters.project_dir))
    parameters.project_dir = check_dir(proj) if len(proj) > 0 else parameters.project_dir
    ass = input("Assigner ({}):".format(parameters.assigner))
    parameters.assigner = ass if len(ass) > 0 else parameters.assigner    
    chunk = input("Chunksize ({}):".format(parameters.chunksize))
    parameters.chunksize = chunk if len(chunk) > 0 else parameters.chunksize
    out = input("Output Type ({}):".format(parameters.output_type))
    parameters.output_type = out if len(out) > 0 else parameters.output_type
    parameters = validate_abstar_params(parameters)
    merge = input("Merge ({}):".format(parameters.merge))
    if len(merge) > 0:
        while merge.upper() not in ABSTAR_PARAMS['merge'].keys():
            merge = input('Merge input must be of type: Bool' \
                          '\nPlease re-enter!: ')
        merge_bool = ABSTAR_PARAMS['merge'][merge.upper()]
        parameters.merge = merge_bool
        
    uid = input("UID ({}):".format(parameters.uid))
    parameters.uid = uid if len(uid) > 0 else parameters.uid
    celery = input("Celery ({}):".format(parameters.cluster))
    if len(celery) > 0:
        while celery.upper() not in ABSTAR_PARAMS['merge'].keys():
            celery = input('Celery input must be of type: Bool' \
                          '\nPlease re-enter!: ')
        celery_bool = ABSTAR_PARAMS['merge'][celery.upper()]
        parameters.cluster = celery_bool    
    species = input("Species ({}): ".format(parameters.species))
    parameters.species = species if len(species) > 0 else parameters.species
    
    return parameters
        
def print_abstar_params(parameters):
        params = [
                "Project Directory: {}".format(parameters.project_dir),
                "Assigner: {}".format(parameters.assigner),
                "Chunksize: {}".format(parameters.chunksize),
                "Output Type: {}".format(parameters.output_type),
                "Merge: {}".format(parameters.merge),
                "UID: {}".format(parameters.uid),
                "Celery Cluster: {}".format(parameters.cluster),
                "Species: {}".format(parameters.species)
                ]
        for param in params:
            print(param)
        change = input('\nDo you want to change Abstar arguments? [y/N]: ')
        if change.upper() in ['Y', 'YES']:
            return True
        else:
            return False
        

def validate_abstar_params(params):
    while params.assigner not in ASSIGNERS:
        print('\nASSIGNER ERROR: Assigner is not recognized' \
              '\nYour current assigner options are:')
        print(str(ABSTAR_PARAMS['assigner']))
        params.assigner = input('Please Re-enter assigner name: ')
        validate_abstar_params(params)
    while params.output_type not in ABSTAR_PARAMS['output_type']:
        print('\nOUTPUT TYPE ERROR: Output Type is not recognized' \
              '\nYour current output options are:')
        print(str(ABSTAR_PARAMS['output_type']))
        params.output_type = input('Please Re-enter assigner name: ')
    return params

def basemount_dir(bsmnt_dir, project):
    if os.path.exists(bsmnt_dir) and not shutil.which('basemount') == None:
        base = bsmnt_dir.split('/Projects')[0]
        os.system("basemount --unmount {} > /dev/null".format(base))
        print('Restarting basemount set point')
        os.system("basemount {} > /dev/null".format(base))
        print(os.path.join(bsmnt_dir, project))
        if not os.path.exists(os.path.join(bsmnt_dir, project)):
            pro = input("ERROR: Project not found! Re-enter Project: ")
            basemount_dir(bsmnt_dir, pro)
    elif shutil.which('basemount') == None:
        print('ERROR: Basemount must be installed')
        sys.exit(2)
    elif not os.path.exists('/basemount/Projects') and not shutil.which('basemount') == None:
        bd = input('Specify Basemount Set Point Directory: ')
        bd_pro = os.path.join(bd, 'Projects')
        bsmnt_dir = basemount_dir(bd_pro, project)
    else:
        print("ERROR!")
        sys.exit(3)
    return bsmnt_dir

def run_abstar(parameters, project, pipeline_args, preprocessing=False):
    default_base_setpoint = '/basemount'
    default_base_projects = os.path.join(default_base_setpoint, 'Projects')
    input_dir = os.path.join(parameters.project_dir, 'input')
    bsmnt_dir = os.path.join(default_base_projects, project)
    if os.path.exists(bsmnt_dir):
        try:
            copy_from_basemount(bsmnt_dir, input_dir)
        except ZeroDivisionError:
            print('ERROR: No Files Found in Basemount!')
            sys.exit(4)
    else:
        bsmnt_dir = basemount_dir(default_base_projects, project)
        pro_dir = os.path.join(bsmnt_dir, project)
        try:
            copy_from_basemount(pro_dir, input_dir)
        except ZeroDivisionError:
            print('ERROR: No Files Found in Basemount!')
            sys.exit(4)
    if preprocessing:
        parameters = preprocess(parameters, pipeline_args)
    if not parameters.merge and not preprocessing:
        print("\n========================================" \
              "\nUnzipping Input Files\n" \
              "========================================\n")
        os.system('gunzip {}/input/*'.format(parameters.project_dir))
    if not parameters.merge and preprocessing:
        print("\n========================================" \
              "\nUnzipping Input Files\n" \
              "========================================\n")
        os.system('gunzip {}/*'.format(parameters.input))
    warnings.filterwarnings("ignore")
    run_standalone(parameters)

######################################################
#                      MongoDB                       #
######################################################


def mongo_params(project, abstar_args):
    print("\n========================================" \
          "\nMongo Import Run Arguments\n" \
          "========================================\n")
    
    abstar_output = os.path.join(abstar_args.project_dir, 'output') if abstar_args.project_dir != None else abstar_args.output
    logs = os.path.join(abstar_args.project_dir, 'log/mongo.log') if abstar_args.project_dir != None else os.path.join(os.path.dirname(args.input), 'log/mongo.log')
    mongo_args = mongoimport.Args(db=project, input=abstar_output, delim1='.', log=logs)
    if print_mongo_args(mongo_args):
        mongo_args = change_mongo_args(mongo_args)
    return mongo_args


def print_mongo_args(parameters):
        params = [
                "IP: {}".format(parameters.ip),
                "Port: {}".format(parameters.port),
                "User: {}".format(parameters.user),
                "Password: {}".format(parameters.password),
                "Input: {}".format(parameters.input),
                "DB: {}".format(parameters.db),
                "Delim 1: {}".format(parameters.delim1),
                "Delim 2: {}".format(parameters.delim2)
                ]
        for param in params:
            print(param)
        change = input('\nDo you want to change Mongo Import arguments? [y/N]: ')
        if change.upper() in ['Y', 'YES']:
            return True
        else:
            return False    

def change_mongo_args(parameters):
    print('\n========================================' \
          '\nTo change a current Mongo argument,' \
          '\ntype in new argument and press enter.' \
          '\nLeave blank and press enter if you want ' \
          '\nto keep current parameter (in parentheses).'
          '\n========================================')
    
    ip=input("IP ({}):".format(parameters.ip))
    parameters.ip = ip if len(ip) > 0 else parameters.ip
    port = input("Port ({}):".format(parameters.port))
    parameters.port = int(port) if len(port) > 0 else parameters.port  
    user = input("User ({}):".format(parameters.user))
    parameters.user = user if len(user) > 0 else parameters.user
    password = getpass("Password ({}):".format(parameters.password))
    parameters.password = password if len(password) > 0 else parameters.password
    Input = input("Input (don't change) ({}):".format(parameters.input))
    parameters.input = Input if len(Input) > 0 else parameters.input
    db = input("DB ({}):".format(parameters.db))
    parameters.db = db if len(db) > 0 else parameters.db
    delim1 = input("Delim 1 ({}):".format(parameters.delim1))
    parameters.delim1 = delim1 if len(delim1) > 0 else parameters.delim1
    delim2 = input("Delim 2 ({}):".format(parameters.delim2))
    parameters.delim2 = delim2 if len(delim2) > 0 else parameters.delim2
    return parameters
        

def run_mongo_import(args):
    mongoimport.run(ip=args.ip,
                    port=args.port,
                    user=args.user,
                    password=args.password,
                    input=args.input,
                    log=args.log,
                    db=args.db,
                    delim1=args.delim1,
                    delim2=args.delim2)

def print_the_splash():
    print_splash()