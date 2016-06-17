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
if len(sys.argv) < 3:
    print('Usage: geno2chrom.py <in_file.fasta> <out_dir>')

# open genome file and create seperate chromosome files
chrom_dir = sys.argv[2]
chrom = False
with open(sys.argv[1], 'r') as genome:
    for i,line in enumerate(genome, 1):
        line = line.strip()
        if i % 5000000 == 0: print('line - %s' % i)
        if re.match(r'>', line):
            if chrom:
                chrom.close()
            chrom_name = line.strip('>')
            chrom_file = '%s/%s.fa' % (chrom_dir, chrom_name)
            chrom = open(chrom_file, 'w')
        chrom.write('%s\n' % line)
    chrom.close()

