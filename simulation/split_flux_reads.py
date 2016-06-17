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


# check command line
if len(sys.argv) < 5:
    print('Usage: split_flux_reads.py <in_folder_flux> <in_prefix> <file_nums> <out_folder_flux>')
    print('Usage: in_folder_flux    folder with flux fastq files')
    print('Usage: in_prefix         prefix file for reads, will be concatenated with file_nums')
    print('Usage: file_nums         read files file number (sep=",")')
    print('Usage: out_folder_flux   folder where output fastq files will be stored')


# get transcript data
in_folder_flux = sys.argv[1]
in_prefix = sys.argv[2]
file_nums = sys.argv[3].split(',')
out_folder_flux = sys.argv[4]
for num in file_nums:
    for sample in ['A','B']:
        # modify A files
        file_fq_0 = '%s%s_%s_%s.fastq' % (in_folder_flux, in_prefix, sample, num)
        file_fq_1 = '%s%s_%s_%s-%s.fastq' % (out_folder_flux, in_prefix, sample, num, 1)
        file_fq_2 = '%s%s_%s_%s-%s.fastq' % (out_folder_flux, in_prefix, sample, num, 2)
        with open(file_fq_0, 'r') as file_fq_0_in, \
             open(file_fq_1, 'w') as file_fq_1_out, \
             open(file_fq_2, 'w') as file_fq_2_out:
            for i,line in enumerate(file_fq_0_in, 1):
                line = line.strip()
                if i % 8 in [1,5]:
                    line = line[:-4]
                if i % 8 in [1,2,3,4]:
                    file_fq_1_out.write('%s\n' % line)
                if i % 8 in [5,6,7,0]:
                    file_fq_2_out.write('%s\n' % line)
        print('File %s ... done.' % file_fq_0)
