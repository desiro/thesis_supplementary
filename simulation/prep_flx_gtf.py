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
if len(sys.argv) < 4:
    print('Usage: prep_flx_gtf.py <in_fasta> <in_gtf> <out_gtf>')
    print('Usage: in_fasta          in fasta file')
    print('Usage: in_gtf            in sorted gtf file')
    print('Usage: out_gtf           out gtf file')


# initialize parameters
chrom_name_list = []


# get transcript data
with open(sys.argv[1], 'r') as fasta_in:
    for line in fasta_in:
        line = line.strip()
        if re.match(r'^>', line):
            chrom_name_list.append(line[1:])
    print(chrom_name_list)
with open(sys.argv[2], 'r') as gtf_in, \
     open(sys.argv[3], 'w') as gtf_out:
    for line in gtf_in:
        line = line.strip()
        chrom_id = line.split('\t')[0]
        if chrom_id in chrom_name_list:
            gtf_out.write('%s\n' % line)

