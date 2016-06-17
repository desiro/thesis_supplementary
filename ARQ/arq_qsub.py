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
import time
import subprocess
from arq_helper import getTime
from arq_helper import makeCombAlign
from arq_helper import writeQuantRef
from arq_tools import makeBashScript
from arq_debug import printDebug
from arq_get import userCall




###############################################################################
## prepare qsub jobs
###############################################################################

def runQsubJobs(job_list, input_parameters, bash_config, align_files, fq1_files, fq2_files, job_type, debug):
    
    # initialize debug
    debug_list = []
    format = '{:<17}'
    if debug: debug_list.append(['def: runQsubJobs():'])
    if debug: debug_list.append(['-'])
    if debug and job_type == 'align': debug_list.append(['job', 'fastq', 'mode'])
    if debug and job_type == 'index': debug_list.append(['job', 'mode'])
    if debug and job_type == 'quant': debug_list.append(['job', 'combination', 'mode'])
    if debug: debug_list.append(['-'])
    
    # initialize parameters
    file_status = {}
    job_queue = []
    queues_empty = False
    q_running = False
    bowtie_running = True if 'BOWTIE_index' in job_list.keys() else False
    TOOLS_running = []
    error_file = False
    
    # make job queues
    print("%s Status: Create queue ..." % (getTime()))
    for job in job_list.keys():
        queue = []
        # save jobs directly when index is requested
        if job_type == 'index':
            if job_list[job] == 'request':
                job_queue.append(job)
        # save fastq jobs when align is requested or combi align jobs when quant is requested
        elif job_type in ['align', 'quant']:
            for file_name in job_list[job].keys():
                if job_list[job][file_name] in ['overwrite', 'request', 'new']:
                    queue.append(file_name)
            file_status[job] = queue
            job_queue.append(job)
    
    # run jobs
    print("%s Status: Run jobs ..." % (getTime()))
    while job_queue or q_running:
        # set user naming
        user_naming = ''
        # set q_running to True, because job_queue isn't empty
        q_running = True
        # print tool queue
        print('%s Queue: %s' % (getTime(), job_queue))
        # take first job from queue
        job = job_queue.pop(0) if job_queue else ''
        # start qsub for index jobs
        if job_type == 'index' and job:
            # append TOPHAT_index to queue if bowtie isn't done yet
            if job == 'TOPHAT_index' and bowtie_running:
                job_queue.append(job)
            # start index job
            else:
                if debug: debug_list.append([job, job_list[job]])
                runQsubCommand(job, input_parameters, bash_config[job], debug)
                job_list[job] = 'running'
                if debug: debug_list.append([job, job_list[job]])
                if (input_parameters['single_tool_job'].lower() in ['single', 'true']):
                    TOOLS_running.append(job)
        elif job_type in ['align', 'quant'] and job:
            queue = file_status.pop(job) if file_status.keys() else []
            # print queue
            print('%s %s: %s' % (getTime(), re.match(r'[A-Z]+', job).group(0), queue))
            # submit quantification or alignment jobs
            while queue and not (job in TOOLS_running):
                file_name = queue.pop(0)
                if debug: debug_list.append([job, file_name, job_list[job][file_name]])
                if job_type == 'align':
                    input_parameters['fastq1_file'] = fq1_files[file_name]
                    input_parameters['fastq2_file'] = fq2_files[file_name]
                if job_type == 'quant':
                    """ TODO: get to work with combi alignments (kallisto) """
                    input_parameters['fastq1_file'] = fq1_files[align_files[file_name][0]]
                    input_parameters['fastq2_file'] = fq2_files[align_files[file_name][0]]
                    input_parameters['ALIGNER_align_dir'] = '%s/%s_align' % (input_parameters['align_dir_path'], file_name.partition('_')[0])
                    input_parameters['alignment_prefix'] = align_files[file_name][0]
                    combi_string = ''
                    for align_name in align_files[file_name]:
                        combi_string = '%s%s/%s.bam ' % (combi_string, input_parameters['ALIGNER_align_dir'], align_name)
                    combi_string = combi_string.strip()
                    input_parameters['combi_align_list'] = combi_string
                    input_parameters['comb_out_prefix'] = file_name
                    if job in input_parameters['quant_name_list'].keys():
                        user_naming = '_%s' % input_parameters['quant_name_list'][job]
                    else:
                        user_naming = '_default'
                    input_parameters['quant_ref_file'] = '%s/%s/%s%s.ref' % (input_parameters['quant_dir_path'], job, file_name, user_naming)
                input_parameters['%s_out_prefix' % job_type] = '%s%s' % (file_name, user_naming)
                input_parameters['%s_out_dir' % job_type] = '%s/%s/%s%s' % (input_parameters['%s_dir_path' % job_type], job, file_name, user_naming)
                # start qsub for alignment and quantification jobs
                if job_list[job][file_name] in ['overwrite', 'request', 'new']:
                    if job_type == 'quant':
                        if len(align_files[file_name]) > 1:
                            makeCombAlign(input_parameters)
                            input_parameters['comb_in_file'] = '%s/%s.bam' % (input_parameters['comb_dir_path'], input_parameters['comb_out_prefix'])
                        else:
                            input_parameters['comb_in_file'] = '%s/%s.bam' % (input_parameters['ALIGNER_align_dir'], input_parameters['alignment_prefix'])
                        # create quantification parameter ref file
                        writeQuantRef(input_parameters, bash_config[job])
                    qsub_name = '%s-%s' % (job, input_parameters['%s_out_prefix' % job_type])
                    runQsubCommand(qsub_name, input_parameters, bash_config[job], debug)
                    job_list[job][file_name] = 'running'
                    if debug: debug_list.append([job, file_name, job_list[job][file_name]])
                # force only one job per tool
                if (input_parameters['single_tool_job'].lower() in ['single', 'true']):
                    TOOLS_running.append(job)
            # save queue if it isn't empty (for when single_tool_job option is active)
            if queue:
                file_status[job] = queue
                job_queue.append(job)
        # print tool status
        if (input_parameters['single_tool_job'].lower() in ['single', 'true']):
            print('%s Running: %s' % (getTime(), TOOLS_running))
        print('%s Qstat: ... working ...' % getTime())
        # sleep this long when jobs are in the queue and bowtie is running
        if job_queue and bowtie_running:
            time.sleep(int(input_parameters['sleep_time_qs']))
        # sleep this long when jobs are in the queue and tools are running
        elif job_queue and TOOLS_running:
            time.sleep(int(input_parameters['sleep_time_qs']))
        # sleep this long when the queue is empty
        elif not job_queue:
            time.sleep(int(input_parameters['sleep_time_qs']))
        # get qstat job info
        (out, err) = subprocess.Popen('qstat', stdout=subprocess.PIPE).communicate()
        # set bowtie running to false if BOWTIE job is not in the queue and not running anymore
        if ('BOWTIE_index' not in job_queue) and not re.search('BOWTIE',out):
            bowtie_running = False
        # check if tool in TOOLS_running are still running
        for active in TOOLS_running:
            if not re.search(active[:8],out):
                TOOLS_running.remove(active)
        # check if any job is running
        if not out: 
            q_running = False
    
    # there was an error with qsub ... exiting
    ''' TODO: make specific error detection (make new file inside bash script) and only dismiss dependent jobs '''
    if os.path.isfile('%s/qsub.error' % (input_parameters['qsub_dir_path'])):
        print("%s Warning: There was an error with qsub! Check qsub_data/qsub.error for more information." % getTime())
        if input_parameters['ignore_qs_error'].lower() not in ['ignore', 'true']:
            user = userCall("Continue? y[es] / n[o]", True)
    
    # update job info
    if debug: debug_list.append(['-'])
    for job in job_list.keys():
        if job_type == 'index':
            if job_list[job] == 'running':
                job_list[job] = 'finished'
            if debug: debug_list.append([job, job_list[job]])
        elif job_type == 'align':
            for fq_name in job_list[job].keys():
                if job_list[job][fq_name] == 'running':
                    job_list[job][fq_name] = 'finished'
                if debug: debug_list.append([job, job_list[job]])
                if debug: debug_list.append([job, fq_name, job_list[job][fq_name]])
                
    # print debug
    if debug: printDebug(debug_list, format)
    if debug: user = userCall("Continue? y[es] / n[o]", True)
    
    return job_list




###############################################################################
## make bash tool config
###############################################################################

def makeBashConfig(job_list, input_parameters, job_config_data, debug):
    
    # initialize debug
    debug_list = []
    format = '{:<17}'
    if debug: debug_list.append(['def: makeBashConfig():'])
    if debug: debug_list.append(['-'])
    if debug: debug_list.append(['job', 'script'])
    if debug: debug_list.append(['-'])
    
    # initialize bash config
    bash_config = {}
    
    # save bash config
    for job in job_list.keys():
        if job_list[job] in ['request', 'ask', 'new']:
            # get job config for requested job
            job_config = job_config_data[job]
            # save bash pre and seq config
            config = {}
            print("%s Status: Create \'bash_pre_script\' for \'%s\' ..." % (getTime(), job))
            config['bash_pre_script'] = job_config['bash_pre_script']
            print("%s Status: Create \'bash_seq_script\' for \'%s\' ..." % (getTime(), job))
            config['bash_seq_script'] = job_config['bash_seq_script']
            # go through 'script_order' list and save script command
            print("%s Status: Create \'bash_main_script\' for \'%s\' ..." % (getTime(), job))
            script = ''
            for order_item in job_config['script_order']:
                if order_item != 'options':
                    script += job_config[order_item] + ' '
                else:
                    for option in job_config['script_options']:
                        if input_parameters[option]:
                            script += job_config[option] + ' '
            print('%s Script: %s' % (getTime(), script))
            if debug: debug_list.append([job, script])
            # save bash main script
            config['bash_main_script'] = script
            # save config for job into bash config
            bash_config[job] = config
    
    # print debug
    if debug: printDebug(debug_list, format)
    if debug: user = userCall("Continue? y[es] / n[o]", True)
    
    # return bash config
    return bash_config




###############################################################################
## run qsub
###############################################################################

def runQsubCommand(job, input_parameters, bash_config, debug):
    
    # initialize debug
    debug_list = []
    format = '{:<17}'
    if debug: debug_list.append(['def: runQsubCommand():'])
    if debug: debug_list.append(['-'])
    
    # initialize qsub parameters
    job_id_name = {}
    script_path = os.path.join(input_parameters['qsub_dir_path'], '%s.sh' % job)
    input_parameters['qs_out_name'] = os.path.join(input_parameters['qsub_dir_path'], '%s.qout' % job)
    
    # generate bash script
    BASH = makeBashScript(bash_config)
    call = BASH % input_parameters
    if debug: debug_list.append(['BASH SCRIPT:\n%s' % call])
    
    # save bash script
    with open(script_path, 'w') as tmp_script:
            tmp_script.write(call)
    
    # run qsub command and get job name
    os.system('qsub -N %s %s' % (job, script_path))
    
    # remove bash script
    if input_parameters['save_bash_call'].lower() not in ['save', 'true']:
        os.remove(script_path)
    
    # print debug
    if debug: printDebug(debug_list, format)
    if debug: user = userCall("Continue? y[es] / n[o]", True)
    
    # return job dict: name <-> id
    #return job_id_name



