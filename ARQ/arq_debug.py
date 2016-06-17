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
from arq_helper import userCall




###############################################################################
## print debug messages
###############################################################################

def printDebug(debug_list, format):
    
    # get console window size
    rows, columns = os.popen('stty size', 'r').read().split()
    just_size = int(columns) / 2 if int(columns) % 2 == 0 else (int(columns) - 1) / 2
    print(''.ljust(just_size-4,'#')+" Debug "+''.ljust(just_size-3,'#'))
    
    # print debug
    for i in range(len(debug_list)):
        debug_msg = "Debug: "
        if debug_list[i][0] != '-':
            for n,j in enumerate(debug_list[i]):
                debug_msg += format.format(" %s " % (j))
                if n < len(debug_list[i])-1:
                    debug_msg += "##"
        elif debug_list[i][0] == '-':
            debug_msg += ''.ljust(2*just_size-7,'-')
        print(debug_msg)
    print(''.ljust(2*just_size,'#'))




###############################################################################
## check all config data
###############################################################################

def checkParametersDebug(input_parameters, job_requests, fastq1_files, fastq2_files, debug):
    
    # initialize debug
    debug_list = []
    format = '{:<17}'
    if debug: debug_list.append(['def: main(): - parameters'])
    if debug: debug_list.append(['-'])
    if debug:
        deb_print_list = {'input_param':input_parameters, 'job_requests':job_requests, 'fastq1_files':fastq1_files, 'fastq2_files':fastq2_files}
        for j,name in enumerate(sorted(deb_print_list.keys()), 1):
            for i,param in enumerate(sorted(deb_print_list[name].keys()), 1):
                if i == 1: debug_list.append([name, '{\'%s\': \'%s\',' % (param, deb_print_list[name][param])])
                elif i == len(deb_print_list[name].keys()): debug_list.append(['', ' \'%s\': \'%s\'}' % (param, deb_print_list[name][param])])
                else: debug_list.append(['', ' \'%s\': \'%s\',' % (param, deb_print_list[name][param])])
            if j < len(deb_print_list.keys()): debug_list.append(['-'])
    if debug: printDebug(debug_list, format)
    if debug: user = userCall("Continue? y[es] / n[o]", True)



