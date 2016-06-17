#!/usr/bin/python
# Script: .py
# Author: Daniel Desiro'
"""
Description:
    ARQ: Assessment and Comparison of RNA-seq Quantification Methods
Usage:
    
Source:
    
Reference:
    
"""
import argparse
import sys
import os
import re
import traceback
import shutil
from arq_helper import unzipFastq
from arq_helper import getTime
from arq_helper import writeFastqRef
from arq_helper import countFastqReads
from arq_helper import removeJob
from arq_helper import tophatMergeOpt
from arq_helper import errorHandling
from arq_tools import getMainDictionaries
from arq_tools import getIndexConfig
from arq_tools import getAlignConfig
from arq_tools import getQuantConfig
from arq_debug import checkParametersDebug
from arq_check import checkJobRequests
from arq_check import checkProgress
from arq_dirs import makeMainDirectory
from arq_dirs import makeToolDirectories
from arq_qsub import runQsubJobs
from arq_qsub import makeBashConfig
from arq_get import getArguments
from arq_get import getJobRequests
from arq_get import getQuantDict
from arq_get import getQuantCombinations




###############################################################################
## main function
###############################################################################

def main(project, config, prefix, dataset, fasta, gtf, fastq1, fastq2, debug, test, arq_path):
    
    
    
    
    ###############################################################################
    ## initialize data structures
    
    # parse option arguments to dictionary for improved error handling
    input_options = {'project':project, 'prefix':prefix}
    # parse file arguments to dictionary for improved error handling
    input_files = {'config':config, 'fasta_gen_file':fasta, 'gtf_index_file':gtf, 'fastq1_links':fastq1, 'fastq2_links':fastq2}
    # initialize dictionary for parameter handling
    ''' TODO: save ALL possible variables and update empty checks in code '''
    print("%s Status: Get input dictionary ..." % (getTime()))
    input_parameters = getMainDictionaries('input_parameters', debug)
    input_parameters['test_call'] = test
    input_parameters['arq_path'] = arq_path
    # tool requests with tool links 
    ''' TODO: make tool requests as input parameter and save them automatically '''
    print("%s Status: Get dependency dictionary ..." % (getTime()))
    possible_jobs = getMainDictionaries('possible_jobs', debug)
    # possible alignment - quantification combinations
    input_parameters['quant_align_combi_list'] = getMainDictionaries('quant_align', debug)
    # job requests and status
    ''' TODO: options: 'request', 'enumerate', 'overwrite', 'skip', 'skipped', 'finished' '''
    # job config data
    ''' TODO: TODO: add special config data dynamically '''
    job_config_data = {}
    # set timezone
    os.environ['TZ'] = 'EST'
    
    ''' TODO: make sam merge folder, same as fq temp folder '''
    ''' TODO: make remove option '''
    ''' TODO: remove fastq file, etc. necessity -> so it can be used only for quantification finished alignments '''
    
    
    
    
    ###############################################################################
    ## save input data and options to data structures
    
    # get check and get all provided arguments
    # TODO: get and check fastq files
    print("%s Status: Get settings ..." % (getTime()))
    input_parameters, fastq1_files, fastq2_files, fastq_ref = getArguments(input_files, input_options, input_parameters, debug)
    
    # generate project directories
    print("%s Status: Make main directories ..." % (getTime()))
    input_parameters, progress_data = makeMainDirectory(input_parameters, debug)
    # add alignment free placeholder
    progress_data.append('AFREE_align')
    
    # save requested jobs
    print("%s Status: Get job requests ..." % (getTime()))
    job_requests = getJobRequests(possible_jobs, input_parameters, debug)
    
    # check for finished jobs
    ''' check and mark dependencies '''
    ''' TODO: do better distinguishing for fastq vs quant jobs '''
    print("%s Status: Check job progress ..." % (getTime()))
    job_requests = checkJobRequests(progress_data, job_requests, input_parameters, debug)
    
    # get tool parameters
    ''' TODO: make it user dependent '''
    ''' TODO: check tool paths '''
    ''' TODO: optional global command '''
    
    # get job config data
    print("%s Status: Get index config ..." % (getTime()))
    job_config_data = getIndexConfig(job_config_data, debug)
    print("%s Status: Get alignment config ..." % (getTime()))
    job_config_data = getAlignConfig(job_config_data, debug)
    print("%s Status: Get quantification config ..." % (getTime()))
    job_config_data = getQuantConfig(job_config_data, debug)
    ''' TODO: remove individual alignment file folder '''
    
    # check all parameters in debug mode
    checkParametersDebug(input_parameters, job_requests, fastq1_files, fastq2_files, debug)
    
    
    
    
    ###############################################################################
    ## start tools - index
    ''' TODO: track time and memory consumption '''
    
    # set qsub parameters
    ''' TODO: make as input parameter -> dependent on number of threads used '''
    #input_parameters['qs_mem_free'] = '40'
    #input_parameters['qs_run_time'] = '28800'
    
    # make index tool directories
    print("%s Status: Make index directories ..." % (getTime()))
    index_jobs, input_parameters = makeToolDirectories(job_requests, input_parameters, 'index', debug)
    
    # make index bash config
    print("%s Status: Make index bash config ..." % (getTime()))
    bash_index_config = makeBashConfig(index_jobs, input_parameters, job_config_data, debug)
    
    # run index jobs
    print("%s Status: Run index jobs ..." % (getTime()))
    index_jobs = runQsubJobs(index_jobs, input_parameters, bash_index_config, {}, {}, {}, 'index', debug)
    
    
    
    
    ###############################################################################
    ## start tools - alignments
    ''' TODO: track time and memory consumption -> qacct -j <number> -> get job number and save it '''
    ''' TODO: get statistics out of qacct '''
    
    # set qsub parameters
    ''' TODO: make as input parameter -> dependent on number of threads used '''
    #input_parameters['qs_mem_free'] = '40'
    #input_parameters['qs_run_time'] = '28800'
    
    # make alignment tool directories
    print("%s Status: Make alignment directories ..." % (getTime()))
    align_jobs, input_parameters = makeToolDirectories(job_requests, input_parameters, 'align', debug)
    
    # make alignment bash config
    print("%s Status: Make alignment bash config ..." % (getTime()))
    bash_align_config = makeBashConfig(align_jobs, input_parameters, job_config_data, debug)
    
    # unzip fastq files
    ''' TODO: make optional to remove fastq at the end '''
    print("%s Status: Unzip fastq1 files ..." % (getTime()))
    fastq1_files_uz = unzipFastq(fastq1_files, input_parameters)
    print("%s Status: Unzip fastq2 files ..." % (getTime()))
    fastq2_files_uz = unzipFastq(fastq2_files, input_parameters)
    
    # check alignment progress
    # INPUT  -> align_jobs      = {'TOOL1_align':'overwrite', 'TOOL2_align':'ask', ...}
    #           fastq1_files_uz = {'fq1':'path', 'fq2':'path', ...}
    # OUTPUT -> fastq_jobs      = {'TOOL1_align':{'fq1':'overwrite', 'fq2':'request', ...},
    #                              'TOOL2_align':{'fq1':'overwrite', 'fq2':'request', ...}, ...}
    print("%s Status: Check alignment progress ..." % (getTime()))
    fastq_jobs = checkProgress(align_jobs, fastq1_files_uz, input_parameters, 'align', debug)
    
    # get tophat merge options
    input_parameters['tophat_merge_opt'] = tophatMergeOpt(input_parameters['tool_dir_toprec'])
    
    # run alignment jobs
    ''' TODO: check fastq progress, when a job breaks '''
    ''' TODO: gzip sam / bam ? '''
    print("%s Status: Run alignment jobs ..." % (getTime()))
    align_jobs = runQsubJobs(fastq_jobs, input_parameters, bash_align_config, {}, fastq1_files_uz, fastq2_files_uz, 'align', debug)
    
    # remove unzipped fastq files
    ''' TODO: make definition for removeFastq and test '''
    ''' TODO: dont remove, neede for kallisto '''
    #if input_parameters['remove_fastq'].lower() in ['remove', 'true']:
        #print("%s Status: Remove unzipped fastq1 files ..." % (getTime()))
        #removeFastq(fastq1_files_uz)
    #if input_parameters['remove_fastq'].lower() in ['remove', 'true']:
        #print("%s Status: Remove unzipped fastq2 files ..." % (getTime()))
        #removeFastq(fastq2_files_uz)
    
    
    
    
    ###############################################################################
    ## start tools - quantification
    ''' TODO: quant<->align tool combination '''
    ''' TODO: get specific number of fastq files (10M, 20M, 50M reads) -> modify fastq1_files_uz and fastq2_files_uz '''
    ''' TODO: track used fastq file names? '''
    
    # set qsub parameters
    ''' TODO: make as input parameter -> dependent on number of threads used '''
    #input_parameters['num_processors'] = '4'
    #input_parameters['qs_mem_free'] = '10'
    #input_parameters['qs_run_time'] = '28800'
    
    # make quantification tool directories
    ''' TODO: make better options -> new option for only make new jobs (skip will always skip complete atm '''
    ''' TODO:     and overwrite will always overwrite all or you have to decide for each '''
    print("%s Status: Make quantification directories ..." % (getTime()))
    quant_jobs, input_parameters = makeToolDirectories(job_requests, input_parameters, 'quant', debug)
    
    # make quantification bash config
    print("%s Status: Make quantification bash config ..." % (getTime()))
    bash_quant_config = makeBashConfig(quant_jobs, input_parameters, job_config_data, debug)
    
    # make quantification job combinations
    ''' TODO: get specific number of fastq files (10M, 20M, 50M reads) -> modify fastq1_files_uz and fastq2_files_uz '''
    print("%s Status: Count fastq reads ..." % (getTime()))
    fastq_ref = countFastqReads(fastq_ref, fastq1_files_uz, input_parameters['seq_style_fq'])
    # INPUT  -> fastq_ref       = {'fq1':['fq_name_1','num_reads'], 'fq2':['fq_name_2','num_reads'], ...}
    #           job_requests    = {'TOOL1_quant':'skip', 'TOOL2_quant':'request', 'TOOL3_align':'skip', ...}
    print("%s Status: Make alignment combinations ..." % (getTime()))
    quant_files = getQuantCombinations(job_requests, fastq_ref, input_parameters, debug)
    
    # get user defined quantification names
    print("%s Status: Get quantification names ..." % (getTime()))
    input_parameters['quant_name_list'] = getQuantDict(possible_jobs, input_parameters)
    
    # check quantification progress
    # INPUT  -> quant_jobs      = {'TOOL1_quant':'overwrite', 'TOOL2_quant':'ask', ...}
    #           quant_files     = {'ALIGNTOOL1_fq1_fq2_fq3':['fq_name_1','fq_name_2','fq_name_3'],
    #                              'ALIGNTOOL2_fq1_fq2_fq3':['fq_name_1','fq_name_2','fq_name_3'],
    #                              'ALIGNTOOL1_fq4_fq5':['fq_name_4','fq_name_5'],  ...} 
    # OUTPUT -> combi_jobs      = {'tool1_quant':{'ALIGNTOOL1_fq1_fq2_fq3':'overwrite', 'ALIGNTOOL2_fq1_fq2_fq3':'request', ...},
    #                              'tool2_quant':{'ALIGNTOOL1_fq1_fq2_fq3':'overwrite', 'ALIGNTOOL2_fq1_fq2_fq3':'request', ...}, ...} 
    print("%s Status: Check quantification progress ..." % (getTime()))
    combi_jobs = checkProgress(quant_jobs, quant_files, input_parameters, 'quant', debug)
    
    # write fastq reference file for alignment compilations
    print("%s Status: Write fastq reference file ..." % (getTime()))
    input_parameters['fastq_ref_file'] = '%s/%s_fastq.ref' % (input_parameters['quant_dir_path'], input_parameters['data_set_name'])
    writeFastqRef(fastq_ref, input_parameters)
    
    # run quantification jobs
    print("%s Status: Run quantification jobs ..." % (getTime()))
    quant_jobs = runQsubJobs(combi_jobs, input_parameters, bash_quant_config, quant_files, fastq1_files_uz, fastq2_files_uz, 'quant', debug)
    
    # remove combined alignment files
    '''TODO: manage path with input parameters '''
    if input_parameters['rem_comb_align'].lower() in ['remove', 'true']:
        print("%s Status: Remove alignment combinations ..." % (getTime()))
        removeCombAlign(quant_files, input_parameters)
    
    
    
    
    ###############################################################################
    ## clean up
    
    # remove alignment tmp folder
    '''
    print("%s Status: Removing temp alignment folders ..." % (getTime()))
    for job in fastq_jobs.keys():
        search_dir = '%s/%s' % (input_parameters['align_dir_path'], job)
        for file in os.listdir(search_dir):
            is_dir = '%s/%s' % (search_dir, file)
            if os.path.isdir(is_dir):
                removeJob(is_dir, [''], debug)
    
    # remove quantification tmp folder
    print("%s Status: Removing temp quantification folders ..." % (getTime()))
    for job in combi_jobs.keys():
        search_dir = '%s/%s' % (input_parameters['quant_dir_path'], job)
        for file in os.listdir(search_dir):
            is_dir = '%s/%s' % (search_dir, file)
            if os.path.isdir(is_dir):
                removeJob(is_dir, [''], debug)
    '''
    # done
    print("%s Status: ARQ finished ..." % (getTime()))
    




###############################################################################
## parse arguments
###############################################################################

if __name__ == '__main__':
    
    # parse input arguments
    usage_args = "%(prog)s -p PROJECT -P PREFIX -D DATASET -F FASTA -G GTF -Q1 FASTQ1 [-Q2 FASTQ2] [options]\n" + "usage: use -h or --help for more information" 
    parser = argparse.ArgumentParser(prog='ARQ.py', description="Assessment and Comparison of RNA-seq Quantification Methods.", prefix_chars='-+', epilog="Reference.", usage=usage_args)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1.1 (alpha)')
    parser.add_argument('-d', '--debug',            dest='debug',   action='store_true',    help="print debug messages")
    parser.add_argument('-p', '--project',          dest='project',                         help="name of the project; folder will be created in PROJ_DIR_PREFIX")
    parser.add_argument('-c', '--config',           dest='config',                          help="optional config file with all important parameters; this will not overwrite command line arguments")
    parser.add_argument('-P', '--PROJ_DIR_PREFIX',  dest='prefix',                          help="prefix of the project directory")
    parser.add_argument('-D', '--DATA_SET_NAME',    dest='dataset',                         help="name of the provided dataset")
    parser.add_argument('-F', '--FASTA_GEN_FILE',   dest='fasta',                           help="input fasta genome file for indexing")
    parser.add_argument('-G', '--GTF_INDEX_FILE',   dest='gtf',                             help="input gtf reference file for indexing")
    parser.add_argument('-Q1', '--FASTQ1_LINKS',    dest='fastq1',                          help="input file with fastq file paths (first end; single end reads); format: <ref_name> <fastq path> (separated by tab)")
    parser.add_argument('-Q2', '--FASTQ2_LINKS',    dest='fastq2',                          help="input file with fastq file paths (second end; paired end reads)")
    parser.add_argument('-t', '--test',             dest='test',    action='store_true',    help="use predefined test data")
    ''' TODO: make full parameter list '''
    parser._optionals.title = "arguments"
    options = parser.parse_args()
    
    
    # handle test
    arq_path = ''
    if options.test:
        arq_path = os.path.dirname(os.path.realpath(__file__))
        print("%s Warning: Test requested ... " % getTime())
        print("%s Warning: Setting parameter -p test_temp ... " % getTime())
        options.project = 'test_temp'
        print("%s Warning: Setting parameter -c %s/test_files/arq_proj_test.conf ... " % (getTime(), arq_path))
        options.config = '%s/test_files/arq_proj_test.conf' % arq_path
    
    
    # check if either a config file or all required parameters are given
    config_file_given = None not in [options.project, options.config]
    required_arguments_given = None not in [options.project, options.prefix, options.dataset, options.fasta, options.gtf, options.fastq1]
    if not (config_file_given or required_arguments_given):
        print("Error: Either a config file or all mandatory parameters are required!" % ())
        print(">>> Script terminated ... exiting ... <<<")
        sys.exit()
    
    
    # call main function and handle exceptions
    try:
        main(options.project, options.config, options.prefix, options.dataset, options.fasta, options.gtf, options.fastq1, options.fastq2, options.debug, options.test, arq_path)
    except KeyboardInterrupt:
        errorHandling('%s Warning: Script interrupted by user ... ' % getTime(), arq_path, options.test)
    except SystemExit:
        errorHandling('%s Warning: System exit ... ' % getTime(), arq_path, options.test)
    except Exception:
        errorHandling('%s Warning: Script exception ... ' % getTime(), arq_path, options.test)
        traceback.print_exc(file=sys.stdout)
    if options.test:
        print("%s Warning: This was a test ... " % getTime())
        print("%s Warning: Test successful ... " % getTime())
        #print("%s Warning: Removing test files ... " % getTime())
        #shutil.rmtree('%s/test_files/test_temp' % arq_path, True)
    sys.exit(0)

