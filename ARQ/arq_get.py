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
from arq_check import checkArgsProject
from arq_check import checkArguments
from arq_debug import printDebug
from arq_helper import userCall
from arq_helper import getTime




###############################################################################
## make combination of quantification jobs
###############################################################################

def getQuantCombinations(job_requests, fastq_ref, input_parameters, debug):
    
    # initialize debug
    debug_list = []
    format = '{:<15}'
    if debug: debug_list.append(['def: getQuantCombinations():'])
    if debug: debug_list.append(['-'])
    if debug: debug_list.append(['combo', 'list'])
    if debug: debug_list.append(['-'])
    
    # initialize parameters
    fastq_combos = {}
    align_fq_comb = {}
    align_jobs = []
    quant_jobs = []
    
    # get possible alignment and quantification jobs
    for job in job_requests.keys():
        if re.search(r'_align$', job): align_jobs.append(re.sub(r'_align$', '', job))
        if re.search(r'_quant$', job): quant_jobs.append(job)
    
    # get possible fastq combinations
    # transforms comb job list into
    print("%s Status: Make fastq combinations ..." % (getTime()))
    for combi_str in filter(None, input_parameters['quant_comb_list'].split(';')):
        fq_combi = ''
        combi_list = []
        for combi in sorted(filter(None, combi_str.split(','))):
            # make fq combi key
            fq_combi = '%s_%s' % (fq_combi, combi)
            # save fastq names into ordered list (according to fq_combi) from -> fastq_ref[j] = [ref_name, num_reads]
            combi_list.append(fastq_ref[combi][0])
        # get final fq_combi entry into fastq_combos -> fastq_combos['_fq1_fq2_fq3'] = ['fq_name_1','fq_name_2','fq_name_3']
        if debug: debug_list.append([fq_combi, combi_list])
        fastq_combos[fq_combi] = combi_list
    print("%s Combinations: %s" % (getTime(), fastq_combos.keys()))
    
    # get possible aligner + fastq_combo combinations
    if debug: debug_list.append(['-'])
    print("%s Status: Make aligner-fastq combinations ..." % (getTime()))
    for aligner in align_jobs:
        for fq_combi in fastq_combos.keys():
            aln_fq = '%s%s' % (aligner, fq_combi)
            # get fastq - aligner combination -> align_fq_comb['ALIGNTOOL1_fq1_fq2_fq3'] = ['fq_name_1','fq_name_2','fq_name_3']
            if debug: debug_list.append([aln_fq, fastq_combos[fq_combi]])
            align_fq_comb[aln_fq] = fastq_combos[fq_combi]
    print("%s Combinations: %s" % (getTime(), align_fq_comb.keys()))
    
    # print debug
    if debug: printDebug(debug_list, format)
    if debug: user = userCall("Continue? y[es] / n[o]", True)
    
    # return
    return align_fq_comb





###############################################################################
## get required arguments
###############################################################################

def getArguments(input_files, input_options, input_parameters, debug):
    
    # initialize parameters
    debug_list = []
    format = '{:<17}'
    if debug: debug_list.append(['def: getArguments():'])
    if debug: debug_list.append(['-'])
    if debug: debug_list.append(['key', 'value'])
    if debug: debug_list.append(['-'])
    
    # check project argument and save new project name if project had a path
    input_parameters['proj_dir_prefix'], input_parameters['p_name_suffix'] = checkArgsProject(input_options['project'], debug)
    if debug: debug_list.append(['p_name_suffix', input_parameters['p_name_suffix']])
    
    # check for command line input files and save if valid
    ''' TODO: do for input_options '''
    for type,file in input_files.iteritems():
        if file and os.path.isfile("%s" % (file)) and not type == 'config':
            input_parameters[type] = file
            if debug: debug_list.append([type, file])
        elif file and not os.path.isfile("%s" % (file)):
            print("Error: No \'%s\' file named \'%s\'!" % (type, file))
            sys.exit()
    
    # check command line prefix option and save if not already declared
    prefix = input_options['prefix']
    if prefix and not os.path.isdir(prefix):
        print("Error: \'%s\' is not a directory!" % (prefix))
        sys.exit()
    elif prefix and not input_parameters['proj_dir_prefix']:
        input_parameters['proj_dir_prefix'] = prefix
    
    # get config arguments
    print("%s Status: Save config parameters ..." % (getTime()))
    input_parameters = getArgsConfig(input_files['config'], input_parameters, debug)
    if debug: debug_list.append(['proj_dir_prefix', input_parameters['proj_dir_prefix']])
    
    # configure test data if test was requested
    if input_parameters['test_call']:
        input_parameters['proj_dir_prefix'] = '%s/test_files' % input_parameters['arq_path']
        input_parameters['fasta_gen_file'] = '%s/test_files/EF204940.fa' % input_parameters['arq_path']
        input_parameters['gtf_index_file'] = '%s/test_files/ef204940.gtf' % input_parameters['arq_path']
        input_parameters['fastq1_links'] = '%s/test_files/arq_proj_test.fq1.info' % input_parameters['arq_path']
        input_parameters['fastq2_links'] = '%s/test_files/arq_proj_test.fq2.info' % input_parameters['arq_path']
    
    # define default parameters
    ''' TODO: make definition '''
    
    # validate arguments
    checkArguments(input_parameters, debug)
    
    # define index base name without path
    if not 'index_base_name' in input_parameters.keys():
        input_parameters['index_base_name'] = os.path.basename(os.path.splitext(input_parameters['fasta_gen_file'])[0])
    if debug: debug_list.append(['index_base_name', input_parameters['index_base_name']])
    
    # define tool directories
    ''' TODO: enable global tool command '''
    
    # save and validate fastq files
    print("%s Status: Save fastq files ..." % (getTime()))
    ''' TODO: make option for automatic fastq2 file read and enable 'paired end' option as input parameter '''
    fastq1_files, fastq_ref = getArgsFastq(input_parameters['fastq1_links'], input_parameters, debug)
    if input_parameters['fastq2_links']:
        fastq2_files, fastq_ref2 = getArgsFastq(input_parameters['fastq2_links'], input_parameters, debug)
        for fq_name in fastq1_files.keys():
            if not ((re.sub(r'\_1.fq.gz$', '', fastq1_files[fq_name]) == re.sub(r'\_2.fq.gz$', '', fastq2_files[fq_name])) or
                (re.sub(r'\-1.fastq.gz$', '', fastq1_files[fq_name]) == re.sub(r'\-2.fastq.gz$', '', fastq2_files[fq_name]))):
                print("%s Warning: The links of \'%s\' are not equal after substituting the enumeration!" % (getTime(), fq_name))
                ''' TODO: DO TEST '''
                user = userCall("Continue? y[es] / n[o]", True)
    
    # declare project folder
    input_parameters['proj_dir_path'] = '%s/%s' % (input_parameters['proj_dir_prefix'], input_parameters['p_name_suffix'])
    if debug: debug_list.append(['proj_dir_path', input_parameters['proj_dir_path']])
    
    # print debug
    if debug: printDebug(debug_list, format)
    if debug: user = userCall("Continue? y[es] / n[o]", True)
    
    # return arguments
    return input_parameters, fastq1_files, fastq2_files, fastq_ref




###############################################################################
## get and validate fastq files
###############################################################################

def getArgsFastq(fastq_links, input_parameters, debug):
    
    # initialize parameters
    debug_list = []
    format = '{:<17}'
    if debug: debug_list.append(['def: getArgsFastq():'])
    if debug: debug_list.append(['-'])
    if debug: debug_list.append(['ref_name', 'fastq_path', 'fastq_num'])
    if debug: debug_list.append(['-'])
    SEQ = True if input_parameters['seq_style_fq'].lower() in ['seq', 'true'] else False
    
    # open fastq links file and save links to {<ref_name>:<fastq_path>} dictionary
    fastq_files = {}
    fastq_ref = {}
    first = True if SEQ else False
    ref_name, fastq_path, num_reads = ('',)*3
    with open(fastq_links,'r') as fastq:
        for i, line in enumerate(fastq, 1):
            # skip header
            if not first:
                line = line.strip()
                try:
                    if not SEQ:
                        ref_name, fastq_path = line.split('\t')
                    else:
                        ref_name, ref_tube, bc_name, bs_seq, lane_id, num_reads, fq_path, fq_name = line.split('\t')
                        if input_parameters['test_call']: fq_path = '%s/test_files/' % input_parameters['arq_path']
                        fastq_path = '%s%s' % (fq_path, fq_name)
                except:
                    # print debug if an error occurs
                    if debug: printDebug(debug_list, format)
                    print("Error: Fastq links file \'%s\' is invalid in line \'%s\'! Can't read \'%s\'." % (fastq_links, i, line))
                    sys.exit()
                j = ('%s' % i) if not SEQ else ('%s' % (i-1))
                if debug: debug_list.append([ref_name, fastq_path, j])
                # chech if fastq file exists
                if not os.path.isfile("%s" % (fastq_path)):
                    if debug: printDebug(debug_list, format)
                    print("Error: No fastq file named \'%s\' in \'%s\' line \'%s\'!" % (fastq_path, fastq_links, i))
                    sys.exit()
                # save ref name and fastq link to dictionary
                fastq_files[ref_name] = fastq_path
                fastq_ref[j] = [ref_name, num_reads]
            first = False
    
    # print debug
    if debug: printDebug(debug_list, format)
    if debug: user = userCall("Continue? y[es] / n[o]", True)
    
    # return dictionary with fastq links and list with fastq data
    return fastq_files, fastq_ref




###############################################################################
## handle config file
###############################################################################

def getArgsConfig(config, input_parameters, debug):
    
    # initialize parameters
    debug_list = []
    format = '{:<17}'
    if debug: debug_list.append(['def: getArgsConfig():'])
    if debug: debug_list.append(['-'])
    if debug: debug_list.append(['config param', 'vaiable name', 'config data', 'terminal data'])
    if debug: debug_list.append(['-'])
    
    # read config file and validate parameters
    with open(config,'r') as cFile:
        for i, line in enumerate(cFile, 1):
            line = line.strip()
            # skip comment lines
            if not line or re.match(r'^\#', line):
                continue
            try:
                conf_param, conf_data = line.split('\t')
            except:
                # print debug if an error occurs
                if debug: printDebug(debug_list, format)
                print("Error: Config file line \'%s\' is invalid! Can't read \'%s\'." % (i, line))
                sys.exit()
            # get mandatory parameters
            valid_conf_param = False
            for param,data in input_parameters.iteritems():
                # only save conf_data if data is empty
                if conf_param.lower() == param:
                    if debug: debug_list.append([conf_param, param, conf_data, data])
                    if not data: input_parameters[param] = conf_data
                    valid_conf_param = True
                    break
            if not valid_conf_param:
                # print debug if an error occurs
                if debug: printDebug(debug_list, format)
                print("Error: Parameter \'%s\' in your config file is not supported!" % (conf_param))
                sys.exit()
    
    # print debug
    if debug: printDebug(debug_list, format)
    if debug: user = userCall("Continue? y[es] / n[o]", True)
    
    # return parameters
    return input_parameters




###############################################################################
## make user defined quantification name dictionary
###############################################################################

def getQuantDict(possible_jobs, input_parameters):
    
    # new dict for quantification names
    quant_names_dict = {}
    
    # assign quantification names
    for name_mod in filter(None, input_parameters['quant_name_list'].split(';')):
        job, name = filter(None, name_mod.split(':'))
        if job in possible_jobs.keys():
            quant_names_dict[job] = name
        else:
            print("%s Warning: \'%s\' is not a valid job, it will be ignored!" % (getTime(), job))
    
    # return dict
    return quant_names_dict




###############################################################################
## handle job requests
###############################################################################

def getJobRequests(possible_jobs, input_parameters, debug):
    
    # initialize debug
    debug_list = []
    format = '{:<17}'
    if debug: debug_list.append(['def: getJobRequests():'])
    if debug: debug_list.append(['-'])
    
    # job dict
    job_requests = {}
    
    # top job type
    job_base = input_parameters['only_make_job'].lower() if input_parameters['only_make_job'].lower() in ['index','align','quant'] else 'quant'
    
    # get aligner and quantifier
    if debug: debug_list.append(['job combi', 'job possible'])
    if debug: debug_list.append(['-'])
    for job_request in filter(None, input_parameters['arq_job_list'].split(',')):
        job_types = ['index','align','quant']
        while job_base in job_types:
            job_type = job_types.pop(0)
            job_ext = '%s_%s' % (job_request.upper(), job_type)
            job_possible = True if job_ext in possible_jobs.keys() else False
            if debug: debug_list.append([job_ext, '%s' % job_possible])
            if job_possible:
                job_requests[job_ext] = ''
    
    # get job dependencies
    if debug: debug_list.append(['-'])
    if debug: debug_list.append(['job', 'job_requests'])
    if debug: debug_list.append(['-'])
    for job_key in job_requests.keys():
        job_requests = getDependency(possible_jobs, job_requests, job_key, [])
        if debug: debug_list.append([job_key, job_requests.keys()])
    
    # print debug
    if debug: printDebug(debug_list, format)
    if debug: user = userCall("Continue? y[es] / n[o]", True)
    
    # return all job requests
    return job_requests




###############################################################################
## get dependencies recursively
###############################################################################

def getDependency(possible_jobs, job_requests, job, requirements):
    
    # get dependencies
    for required in filter(None, possible_jobs[job].split(',')):
        requirements.append(required)
    
    # iterate through dependencies
    if requirements:
        required = requirements.pop(0)
        job_requests = getDependency(possible_jobs, job_requests, required, requirements)
    
    # save dependency
    job_requests[job] = ''
    
    # return job requests
    return job_requests



