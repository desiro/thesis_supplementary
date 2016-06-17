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
    print('Usage: get_cuff_counts.py <total_reads.file> <in_folder_fpkm> <in_extension> <out_extension>')
    print('Usage: total_reads.file  file with cufflinks stdout "Normalized Map Mass" (total read count)')
    print('Usage:                   format (sep=tab): CUFFL_quant ALIGNER fq_ref_num user_naming read_count')
    print('Usage: in_folder_fpkm    folder with cufflinks fpkm transcript files (will search for files from total_reads.file.cnt')
    print('Usage:                   used columns: 1:tracking_id, 8:length, 10:FPKM')
    print('Usage: in_extension     read files file in extension (e.g. fpkm)')
    print('Usage: out_extension    read files file out extension (e.g. cnt)')


# calculate cufflinks read counts
def calculateCounts(FPKM, N, L):
    # FPKM = (10^9 * C)/(N * L)
    # C = Number of reads mapped to a gene
    # N = Total mapped reads in the experiment
    # L = gene length in base-pairs for a gene
    C = (float(FPKM) * float(N) * float(L)) / float(10**9)
    return C


# initialize parameters
total_reads_hash = {}


# get total reads
with open(sys.argv[1], 'r') as total_reads:
    for line in total_reads:
        line = line.strip()
        quantifier, aligner, fastq_ref_num, user_naming, total_count = line.split('\t')
        file_prefix = '%s_%s_%s' % (aligner, fastq_ref_num, user_naming)
        total_reads_hash[file_prefix] = total_count


# get cufflinks transcript data
for prefix in total_reads_hash.keys():
    fpkm_file = '%s/%s.%s' % (sys.argv[2], prefix, sys.argv[3])
    count_file = '%s/%s.%s' % (sys.argv[2], prefix, sys.argv[4])
    if os.path.isfile('%s' % fpkm_file):
        with open(fpkm_file, 'r') as cuff_in, \
             open(count_file, 'w') as cuff_out:
            first_line = True
            for line in cuff_in:
                if not first_line:
                    line = line.strip()
                    tracking_id, class_code, nearest_ref_id, gene_id, gene_short_name, tss_id, locus, length, coverage, FPKM, FPKM_conf_lo, FPKM_conf_hi, FPKM_status = line.split('\t')
                    counts = calculateCounts(FPKM, total_reads_hash[prefix], length) if FPKM != '0' else 0
                    cuff_out.write('%s\t%s\n' % (tracking_id, counts))
                else:
                    first_line = False

