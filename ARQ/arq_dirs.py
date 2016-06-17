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
import sys
import os
import re
import shutil
from arq_debug import printDebug
from arq_helper import userCall
from arq_helper import removeJob
from arq_helper import getTime




###############################################################################
## generate main directories
###############################################################################

def makeMainDirectory(input_parameters, debug):
    
    # initialize parameters
    project_directory = input_parameters['proj_dir_path']
    progress_data = []
    debug_list = []
    format = '{:<17}'
    if debug: debug_list.append(['def: makeMainDirectory():'])
    if debug: debug_list.append(['-'])
    
    # directory list
    dir_list = ['index', 'align', 'quant', 'qsub', 'fastq', 'comb']
    # save project directories
    input_parameters['index_dir_path'] = '%s/indexes' % (project_directory)
    input_parameters['align_dir_path'] = '%s/alignments/%s' % (project_directory, input_parameters['data_set_name'])
    input_parameters['quant_dir_path'] = '%s/quantifications/%s' % (project_directory, input_parameters['data_set_name'])
    input_parameters['qsub_dir_path'] = '%s/qsub_data' % (project_directory)
    input_parameters['fastq_dir_path'] = '%s/fastq_data' % (project_directory)
    input_parameters['comb_dir_path'] = '%s/comb_data/%s' % (project_directory, input_parameters['data_set_name'])
    
    
    # make new project folders
    new_project = False
    if not os.path.isdir(project_directory):
        # main project folder
        print("%s Status: Create \'%s\' main folder ..." % (getTime(), input_parameters['p_name_suffix']))
        os.mkdir(project_directory)
        new_project = True
    # make arq data directories
    if not new_project:
        print("%s Warning: \'%s\' already exists! Checking progress ..." % (getTime(), project_directory))
    # make data folders
    if debug: debug_list.append(['dir key', 'dir path'])
    if debug: debug_list.append(['-'])
    for data_type in dir_list:
        if not os.path.isdir(input_parameters['%s_dir_path' % data_type]):
            print("%s Status: Create \'%s\' folder ..." % (getTime(), data_type))
            if debug: debug_list.append(['%s_dir_path' % data_type, input_parameters['%s_dir_path' % data_type]])
            os.makedirs('%s' % (input_parameters['%s_dir_path' % data_type]))
        else:
            print("%s Warning: \'%s\' folder already exists ..." % (getTime(), data_type))
    # remove qsub error file if exists
    if os.path.isfile('%s/qsub.error' % (input_parameters['qsub_dir_path'])):
        print("%s Status: Remove \'qsub.error\' file ..." % (getTime()))
        os.remove('%s/qsub.error' % (input_parameters['qsub_dir_path']))
    # check tool progress
    if not new_project:
        print("%s Status: Check tool progress ..." % (getTime()))
        for mainDirs in [input_parameters['index_dir_path'],input_parameters['align_dir_path'],input_parameters['quant_dir_path']]:
            if debug: debug_list.append(['mainDirs', mainDirs])
            # list all contents and check if content is a directory
            for jobName in os.listdir(mainDirs):
                subDirs = '%s/%s' % (mainDirs, jobName)
                if debug: debug_list.append(['subDirs', subDirs])
                if os.path.isdir(subDirs):
                    # save progress into progress_data
                    print("%s Warning: Job " % getTime() + '{:<15}'.format("\'%s\'" % (jobName)) + " ... is already done!")
                    progress_data.append(jobName)
                    if debug: debug_list.append([jobName, True if jobName in progress_data else False])
    
    # print debug
    if debug: printDebug(debug_list, format)
    if debug: user = userCall("Continue? y[es] / n[o]", True)
    
    # return parameters
    return input_parameters, progress_data




###############################################################################
## make tool directories
###############################################################################

def makeToolDirectories(job_requests, input_parameters, job_type, debug):
    
    # initialize debug
    debug_list = []
    format = '{:<17}'
    if debug: debug_list.append(['def: makeToolDirectories():'])
    if debug: debug_list.append(['-'])
    if debug: debug_list.append(['job', 'setting'])
    if debug: debug_list.append(['-'])
    
    # initiate tool job requests
    job_list = {}
    for job in job_requests.keys():
        # save all tool jobs to job_list
        if re.search(r'_%s$' % job_type,job):
            job_list[job] = job_requests[job]
            if debug: debug_list.append([job, job_list[job]])
    
    # declare tool job folders
    # job options: 'request', 'overwrite', 'skip', 'skipped', 'finished', 'required', 'ask', 'new'
    if debug: debug_list.append(['-'])
    for job in job_list.keys():
        job_dir_name = '%s_dir' % job
        job_folder = '%s/%s' % (input_parameters['%s_dir_path' % job_type], job)
        input_parameters[job_dir_name] = job_folder
        if job_list[job] in ['request', 'overwrite', 'required', 'new']:
            # check if it is an alignment job and ask if all alignments should be redone
            if job_list[job] == 'overwrite' and job_type in ['align', 'quant']:
                print("%s Warning: You requested \'overwrite\' for job \'%s\'!"  % (getTime(), job))
                overwrite_all = userCall("Type: y[es] to remove all / n[o] to decide for each job part", False)
                job_list[job] = 'overwrite' if overwrite_all else 'ask'
            # check for overwrite option and remove directory
            if job_list[job] == 'overwrite':
                removeJob(job_folder, [''], debug)
            # make new folder
            if not os.path.isdir(job_folder) and job_list[job] in ['request', 'overwrite', 'new']:
                print("%s Status: Create \'%s\' folder ..." % (getTime(), job))
                os.mkdir(job_folder)
                if job_list[job] == 'overwrite': job_list[job] = 'request'
                if debug: debug_list.append([job, job_list[job]])
            elif not os.path.isdir(job_folder) and job_list[job] == 'required':
                if debug: printDebug(debug_list, format)
                print("Error: \'%s\' is required but doesn't exist!" % (job_folder))
                sys.exit()
            elif os.path.isdir(job_folder) and job_list[job] in ['required', 'ask', 'new']:
                pass
            else:
                if debug: printDebug(debug_list, format)
                print("Error: \'%s\' already exists!" % (job_folder))
                sys.exit()
        elif job_list[job] == 'skip':
            job_list[job] = 'skipped'
            if debug: debug_list.append([job, job_list[job]])
    
    # print debug
    if debug: printDebug(debug_list, format)
    if debug: user = userCall("Continue? y[es] / n[o]", True)
    
    # return index jobs
    return job_list, input_parameters



