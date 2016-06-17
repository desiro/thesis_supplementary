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
    print('Usage: trans2gene.py <in_file.gtf> <out_file.idx>')

# open gtf file and extract gene and transcript ids
with open(sys.argv[1], 'r') as gtf:
    with open('%s.tmp' % sys.argv[2], 'w') as out:
        for i,line in enumerate(gtf, 1):
            line = line.strip()
            if not re.match(r'#', line):
                if i % 100000 == 0: print('line - %s' % i)
                split_line = line.split()
                gene_id = split_line[9]
                transcript_id = split_line[11]
                gene_id = gene_id.strip('\s;"')
                transcript_id = transcript_id.strip('\s;"')
                out.write('%s\t%s\n' % (transcript_id, gene_id))

# sort according to transcript id's
os.system('sort -k1,1 %s.tmp | uniq > %s' % (sys.argv[2],sys.argv[2]))

# remove tmp file
os.remove('%s.tmp' % sys.argv[2])