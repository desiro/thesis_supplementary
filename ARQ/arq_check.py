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
## check progress for each tool and fastq file
###############################################################################

def checkProgress(job_list, job_parts, input_parameters, job_type, debug):
    
    # initialize debug
    debug_list = []
    format = '{:<15}'
    if debug: debug_list.append(['def: checkProgress():'])
    if debug: debug_list.append(['-'])
    if debug: debug_list.append(['job', 'fastq', 'status', 'loop'])
    if debug: debug_list.append(['-'])
    
    # initialize parameters
    todo_jobs = {}
    job_parts_keys = []
    quant_align_combi_list = input_parameters['quant_align_combi_list']
    
    # make job list for fastq jobs
    for job in job_list.keys():
        # get job part keys
        job_parts_keys = job_parts.keys()
        # make modified job_list for quantification jobs
        user_naming = ''
        if job_type == 'quant':
            job_parts_keys = []
            for job_name in job_parts.keys():
                if job in input_parameters['quant_name_list'].keys():
                    # add user defined name tag
                    user_naming = '_%s' % input_parameters['quant_name_list'][job]
                    job_name = '%s%s' % (job_name, user_naming)
                else:
                    user_naming = '_default'
                    # add 'default' name tag if no user defined name tag was given for this job
                    job_name = '%s%s' % (job_name, user_naming)
                if ('%s_%s' % (job.partition('_')[0], job_name.partition('_')[0])) in quant_align_combi_list:
                    job_parts_keys.append(job_name)
        # make dict for each job type
        job_requests = {}
        # check each job when 'ask' was requested
        if job_list[job] in ['ask', 'new']:
            # check existing tool and job part data and save it to todo_jobs dict
            job_tool_dir = '%s/%s' % (input_parameters['%s_dir_path' % job_type], job)
            for jobName in os.listdir(job_tool_dir):
                jobBaseName = os.path.basename(os.path.splitext(jobName)[0])
                if debug: debug_list.append([job, jobName, jobBaseName, 'basename'])
                # only save existing jobs which are in the job part file list
                if jobBaseName in job_parts_keys:
                    if job_type == 'quant': jobBaseName = re.sub(r'\_[a-zA-Z0-9]+$', '', jobBaseName)
                    job_requests['%s' % jobBaseName] = ''
                    if debug: debug_list.append([job, jobName, jobBaseName, 'bname in keys'])
            if debug: debug_list.append(['-'])
            # iterate through requested job part names
            for job_name in job_parts.keys():
                if (job_type in ['index', 'align']) or (job_type == 'quant' and (('%s_%s' % (job.partition('_')[0], job_name.partition('_')[0])) in quant_align_combi_list)):
                    # job string
                    job_string = '%s/%s/%s%s' % (input_parameters['%s_dir_path' % job_type], job, job_name, user_naming)
                    # check if job was already done
                    if job_name in job_requests.keys() and job_list[job] == 'ask':
                        print("%s Warning: Job part \'%s%s\' for \'%s\' already exists!" % (getTime(), job_name, user_naming, job))
                        overwrite = userCall("Type: y[es] to overwrite / n[o] to skip this job part", False)
                        if overwrite:
                            removeJob(job_string, ['','.sam','.bam','.ref'], debug)
                            job_requests[job_name] = 'overwrite'
                            print("%s Status: Create \'%s%s\' folder for \'%s\' ..." % (getTime(), job_name, user_naming, job))
                            os.mkdir(job_string)
                        else:
                            job_requests[job_name] = 'skip'
                        if debug: debug_list.append([job, job_name, job_requests[job_name], 'job exists'])
                    elif job_name not in job_requests.keys():
                        job_requests[job_name] = 'request'
                        print("%s Status: Create \'%s%s\' folder for \'%s\' ..." % (getTime(), job_name, user_naming, job))
                        os.mkdir(job_string)
                        if debug: debug_list.append([job, job_name, job_requests[job_name], 'new job'])
                    else:
                        job_requests[job_name] = 'skipped'
                        if debug: debug_list.append([job, job_name, job_requests[job_name], 'else'])
            if debug: debug_list.append(['-'])
        # save overwrite all and request jobs
        elif job_list[job] in ['overwrite', 'request']:
            for job_name in job_parts.keys():
                if (job_type in ['index', 'align']) or (job_type == 'quant' and (('%s_%s' % (job.partition('_')[0], job_name.partition('_')[0])) in quant_align_combi_list)):
                    # job string
                    job_string = '%s/%s/%s%s' % (input_parameters['%s_dir_path' % job_type], job, job_name, user_naming)
                    print("%s Status: Create \'%s%s\' folder for \'%s\' ..." % (getTime(), job_name, user_naming, job))
                    os.mkdir(job_string)
                    job_requests[job_name] = job_list[job]
                    if debug: debug_list.append([job, job_name, job_requests[job_name], 'overwrite or request'])
        # save to fastq jobs
        if not job_list[job] == 'skipped':
            todo_jobs[job] = job_requests
        if debug: debug_list.append(['-'])
    
    # add todo_jobs debug
    if debug: debug_list.append(['-'])
    if debug:
        for job in todo_jobs.keys():
            for part in todo_jobs[job]:
                debug_list.append([job, part, todo_jobs[job][part], 'todo_jobs'])
    
    # print debug
    if debug: printDebug(debug_list, format)
    if debug: user = userCall("Continue? y[es] / n[o]", True)
    
    # return fastq_jobs
    return todo_jobs




###############################################################################
## validate finished jobs
###############################################################################

def checkJobRequests(progress_data, job_requests, input_parameters, debug):
    
    # initialize parameters
    debug_list = []
    format = '{:<15}'
    if debug: debug_list.append(['def: checkJobRequests():'])
    if debug: debug_list.append(['-'])
    if debug: debug_list.append(['job', 'requests[job]', 'progress[job]', 'skip all', 'over all'])
    if debug: debug_list.append(['-'])
    
    # iterate through job requests and ask user to skip job if already done
    do_for_all = {'for_all_index':'','for_all_align':'','for_all_quant':'', 'for_all_new':''}
    # skip      -> skip all <TOOL> jobs
    # overwrite -> overwrite all <TOOL> jobs (in case of alignment/quantification, it will be asked for every alignment/quantification)
    # new       -> make only jobs which aren't done yet
    user_options = {'s':'skip', 'skip':'skip', 's all':'skip', 'skip all':'skip', 
                    'o':'overwrite', 'overwrite':'overwrite', 'o all':'overwrite', 'overwrite all':'overwrite',
                    'n':'new', 'new':'new', 'n all':'new', 'new all':'new'}
    # make config setting list for user_set_over, user_set_skip and user_set_new
    conf_settings = {}
    options = ['over', 'skip', 'new']
    for opt in options:
        for job in filter(None, input_parameters['user_set_%s' % opt].split(',')):
            conf_settings[job] = opt[0]
    # ask user
    for job in job_requests.keys():
        job_done = True if job in progress_data else False
        # call if job was already done
        if job_done and job == 'AFREE_align':
            job_requests[job] = 'skip'
        elif job_done:
            all_job = 'for_all_%s' % job.rpartition('_')[2]
            if job not in conf_settings.keys():
                user_input = do_for_all[all_job] if do_for_all[all_job] else raw_input("%s Warning: What to do with existing \'%s\' job \'%s\'? a[bort] / s[kip] / o[verwrite] / n[ew] [all]: " % (getTime(), job[-5:], job))
            else:
                user_input = conf_settings[job]
            input_valid = False
            # only accept  (case insensitive) as input
            while not input_valid:
                # handle do for all cases
                if user_input.lower() in ['s all', 'skip all', 'o all', 'overwrite all', 'n all', 'new all']: do_for_all[all_job] = user_options[user_input.lower()]
                # do input command
                if user_input.lower() in user_options.keys():
                    job_requests[job] = user_options[user_input.lower()]
                    input_valid = True
                    print("%s Status: Setting for job " % getTime() + '{:<15}'.format("\'%s\'" % (job)) + " -> \'%s\'!" % (job_requests[job]))
                    if debug: debug_list.append([job, job_requests[job], job_done, user_input, do_for_all[all_job]])
                # abort if requested
                elif user_input.lower() in ['a', 'abort']:
                    if debug: printDebug(debug_list, format)
                    print("%s Warning: Script interrupted by user ... " % getTime())
                    sys.exit()
                # repeat until input is valid
                else:
                     user_input = raw_input("%s Warning: Wrong input! Use a[bort] / s[kip] / o[verwrite] / n[ew] [all]: " % getTime())
        # job not yet done
        else:
            job_requests[job] = 'request'
            print("%s Status: Setting for job " % getTime() + '{:<15}'.format("\'%s\'" % (job)) + " -> \'%s\'!" % (job_requests[job]))
            if debug: debug_list.append([job, job_requests[job], job_done, '', ''])
    
    # print debug
    if debug: printDebug(debug_list, format)
    if debug: user = userCall("Continue? y[es] / n[o]", True)
    
    # return updated job_requests
    return job_requests




###############################################################################
## validate arguments
###############################################################################

def checkArguments(input_parameters, debug):
    
    # initialize parameters
    debug_list = []
    format = '{:<17}'
    if debug: debug_list.append(['def: checkArguments():'])
    if debug: debug_list.append(['-'])
    if debug: debug_list.append(['param', 'data'])
    if debug: debug_list.append(['-'])
    
    # check all parsed arguments in the input_parameters dictionary
    for param,data in input_parameters.iteritems():
        # save date if debug option is active
        if debug: debug_list.append([param, data])
        # check for files
        if param in ['fasta_gen_file', 'gtf_index_file', 'fastq1_links', 'fastq2_links']:
            if not os.path.isfile("%s" % (data)):
                if debug: printDebug(debug_list, format)
                print("Error: No \'%s\' file named \'%s\'!" % (param, data))
                sys.exit()
        # check for tool directories
        elif param and re.match(r'tool_dir_', param):
            if not os.path.isdir("%s" % (data)):
                if debug: printDebug(debug_list, format)
                print("Error: No \'%s\' directory named \'%s\'!" % (param, data))
                sys.exit()
        else:
            ''' TODO: check other arguments '''
            pass
    
    # print debug
    if debug: printDebug(debug_list, format)
    if debug: user = userCall("Continue? y[es] / n[o]", True)




###############################################################################
## check project argument -> DONE
###############################################################################

def checkArgsProject(project, debug):
    
    # initialize parameters
    proj_dir_prefix = ''
    debug_list = []
    format = '{:<17}'
    if debug: debug_list.append(['def: checkArgsProject():'])
    if debug: debug_list.append(['-'])
    
    # check if project has a file path
    proj_path = os.path.dirname(project)
    proj_base = os.path.basename(project)
    if debug: debug_list.append(['proj_path', proj_path])
    if debug: debug_list.append(['proj_base', proj_base])
    if not proj_base:
        print("Error: Project name not defined properly!" % ())
        sys.exit()
    elif proj_path and not os.path.isdir(proj_path):
        print("Error: \'%s\' has a path, but \'%s\' is not a directory!" % (proj_base, proj_path))
        sys.exit()
    # ask user if this file path should be used as prefix
    elif proj_path and os.path.isdir(proj_path):
        print("%s Warning: \'%s\' has the path \'%s\'!" % (getTime(), proj_base, proj_path))
        print("%s Warning: This will overwrite all other \'PROJ_DIR_PREFIX\' variables." % (getTime()))
        user = userCall("Continue? y[es] / n[o]", True)
        proj_dir_prefix = proj_path
    
    # print debug
    if debug: printDebug(debug_list, format)
    if debug: user = userCall("Continue? y[es] / n[o]", True)
    
    # return arguments
    return proj_dir_prefix, proj_base



