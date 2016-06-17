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
import os
import sys
import shutil
import time
import subprocess




###############################################################################
## unzip fastq files
###############################################################################

def unzipFastq(fastq_files, input_parameters):
    
    # initialize unzip dict
    fastq_files_uz = {}
    
    # unzip fastq files and save file locations in new dictionary
    for fq_name in fastq_files.keys():
        unzip_name = os.path.basename(os.path.splitext(fastq_files[fq_name])[0])
        unzip_file = '%s/%s' % (input_parameters['fastq_dir_path'], unzip_name)
        fastq_files_uz[fq_name] = unzip_file
        # check if fastq file already exists
        if not os.path.isfile(unzip_file):
            os.system('zcat %s > %s' % (fastq_files[fq_name], unzip_file))
    
    # return unzipped fastq dict
    return fastq_files_uz




###############################################################################
## count fastq reads
###############################################################################

def countFastqReads(fastq_ref, fastq_files, seqStyle):
    
    # check SEQ
    SEQ = True if seqStyle.lower() in ['seq', 'true'] else False
    
    # initialize fastq data list 
    read_count = 0
    
    # get number of total lines
    if not SEQ:
        for fq_num in fastq_ref.keys():
            fq_name = fastq_ref[fq_num][0]
            fq_file = fastq_files[fq_name]
            (out,err) = subprocess.Popen('wc -l < %s' % fq_file, stdout=subprocess.PIPE, shell=True).communicate()
            out = out.strip()
            read_count = int(out) / 4
            fastq_ref[fq_num][1] = '%s' % read_count
    
    # return fastq ref file
    return fastq_ref




###############################################################################
## remove fastq files
###############################################################################

def removeFastq(fastq_files):
    
    # remove all fastq files
    ''' TODO: TEST '''
    for fq_name in fastq_files.keys():
        os.remove(fastq_files[fq_name])




###############################################################################
## write fastq reference file
###############################################################################

def writeFastqRef(fastq_ref, input_parameters):
    
    # write file
    if not os.path.isfile('%s' % input_parameters['fastq_ref_file']):
        with open(input_parameters['fastq_ref_file'], 'w') as ref_file: 
            for fq_num in sorted(fastq_ref.keys(), key = lambda i: int(i)):
                ref_file.write('%s\t%s\t%s\n' % (fq_num, fastq_ref[fq_num][0], fastq_ref[fq_num][1]))




###############################################################################
## write quantification parameter reference file
###############################################################################

def writeQuantRef(input_parameters, bash_config):
    
    # append to file
    with open(input_parameters['quant_ref_file'], 'w') as ref_file:
        for script in ['bash_pre_script', 'bash_main_script', 'bash_seq_script']:
            command = bash_config[script] % input_parameters
            ref_file.write('%s\t%s\n' % (script, command))




###############################################################################
## tophat merge option
###############################################################################

def tophatMergeOpt(tool_dir_toprec):
    
    # check if tophat-recondition.py is provided and valid
    if tool_dir_toprec and os.path.isfile('%stophat-recondition.py' % tool_dir_toprec):
        tophat_merge_opt = 'TRUE'
    else:
        tophat_merge_opt = 'FALSE'
    
    return tophat_merge_opt




###############################################################################
## create combined alignment files
###############################################################################

def makeCombAlign(input_parameters):
    
    # check if the alignment compilation not exists and create it
    if not os.path.isfile('%(comb_dir_path)s/%(comb_out_prefix)s.bam' % input_parameters):
        print("%s Status: Create alignment compilation \'%s.bam\' ..." % (getTime(), input_parameters['comb_out_prefix']))
        command_list = ['samtools view -H %(ALIGNER_align_dir)s/%(alignment_prefix)s.bam > %(comb_dir_path)s/%(comb_out_prefix)s-header.sam', # create header information
                        'samtools merge -h %(comb_dir_path)s/%(comb_out_prefix)s-header.sam %(comb_dir_path)s/%(comb_out_prefix)s-unsorted.bam %(combi_align_list)s', # merge alignments
                        'samtools sort -n -m 20000000000 %(comb_dir_path)s/%(comb_out_prefix)s-unsorted.bam %(comb_dir_path)s/%(comb_out_prefix)s', # sort alignments
                        'rm -f %(comb_dir_path)s/%(comb_out_prefix)s-unsorted.bam', # remove unsorted merge file
                        'rm -f %(comb_dir_path)s/%(comb_out_prefix)s-header.sam'] # remove header information
        for command in command_list:
            command = command % input_parameters
            os.system(command)




###############################################################################
## remove combined alignment files
###############################################################################

def removeCombAlign(comb_files, input_parameters):
    
    # remove combined alignment files
    ''' TODO: TEST '''
    for combi_name in comb_files.keys():
        os.remove('%s/%s.bam' % (input_parameters['comb_dir_path'], comb_files[combi_name]))




###############################################################################
## ask user if he want's to continue
###############################################################################

def userCall(question, quit):
    
    # ask user
    user_input = raw_input("%s Warning: %s: " % (getTime(), question))
    input_valid = False
    answer = False
    # only accept 'y', 'yes', 'n' or 'no' (case insensitive) as input
    while not input_valid:
        if user_input.lower() in ['y', 'yes']:
            input_valid = True
            answer = True
        elif (user_input.lower() in ['n', 'no']) and not quit:
            input_valid = True
            answer = False
        elif (user_input.lower() in ['n', 'no']) and quit:
            print("%s Warning: Script interrupted by user ... " % getTime())
            sys.exit()
        else:
            user_input = raw_input("%s Warning: Wrong input! Use y[es] or n[o]: " % getTime())
    
    # return
    return answer




###############################################################################
## remove job trees
###############################################################################

def removeJob(job_path, extensions, debug):
    
    # try to remove all possible extensions or job folder
    for ext in extensions:
        removed = False
        print("%s Warning: Trying to remove \'%s%s\' ... " % (getTime(), job_path, ext))
        try:
            # remove a folder when job_path is a directory
            if os.path.isdir(job_path):
                shutil.rmtree(job_path, True)
            # remove a file if 'jb_path.ext' is a file
            elif os.path.isfile('%s%s' % (job_path, ext)):
                os.remove('%s%s' % (job_path, ext))
            removed = True
        except:
            removed = False
        # give status about removing the folder/file
        if removed and (os.path.isdir(job_path) or os.path.isfile('%s%s' % (job_path, ext))):
            print("%s Warning: ... \'%s%s\' removed!" % (getTime(), job_path, ext))
            if debug: print("Debug: You can check if \'%s%s\' was deleted properly. " % (job_path, ext))
            if debug: raw_input("Debug: Press ENTER to continue ... ")
        elif not removed:
             print("%s Warning: Could not remove \'%s%s\'!" % (getTime(), job_path, ext))




###############################################################################
## get time
###############################################################################

def getTime():
    
    # get time stamp
    ts = time.strftime('%x %X')
    
    #return time
    return ts




###############################################################################
## handle errors
###############################################################################

def errorHandling(error_message, arq_path, is_test):
    
    # handling and removing
    if is_test:
        print("%s Warning: This was a test ... " % getTime())
        print("%s Warning: Test failed ... " % getTime())
        print("%s Warning: Removing test files ... " % getTime())
        shutil.rmtree('%s/test_files/test_temp' % arq_path, True)
    print('%s' % error_message)
    print(">>> Script terminated ... exiting ... <<<")
    sys.exit(0)



