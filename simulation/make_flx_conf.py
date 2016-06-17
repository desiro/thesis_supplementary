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
    print('Usage: make_flx_conf.py <in_folder_counts> <file_extension> <file_nums> <in_file.pro> <out_file.pro>')
    print('Usage: in_folder_counts  folder with count transcript files (will search for files with "file_extension"')
    print('Usage:                   format (sep=tab): transcript_id read_count')
    print('Usage: file_extension    read files file extension (e.g. cnt, fpkm, ...)')
    print('Usage: file_nums         read files file number (sep=","; filter reads)')
    print('Usage: in_file.pro       flux simulator .pro file expression step')
    print('Usage: out_file.pro      modified flux simulator .pro file')


# initialize parameters
transcript_hash = {}
total_trans_hash = {}
max_trans_hash = {}
occurence_hash = {}


# get transcript data
valid_numbers = sys.argv[3].split(',')
for read_file in os.listdir(sys.argv[1]):
    if os.path.splitext(read_file)[1] == '.%s' % sys.argv[2] and re.search(r'[0-9]+', read_file).group() in valid_numbers:
        print('Reading: %s' % read_file)
        with open('%s%s' % (sys.argv[1], read_file), 'r') as cuff_in:
            first_line = True
            for line in cuff_in:
                line = line.strip()
                trans_id = ''
                read_count = ''
                if sys.argv[2] == 'cnt':
                    trans_id, read_count = line.split('\t')
                elif sys.argv[2] == 'fpkm':
                    if first_line:
                        first_line = False
                        continue
                    line_items = line.split('\t')
                    trans_id = line_items[0]
                    read_count = line_items[10]
                else:
                    print('ERROR')
                    sys.exit()
                old_trans_entry = total_trans_hash.get(trans_id, '0')
                old_occur_entry = occurence_hash.get(trans_id, 0)
                old_max_entry = max_trans_hash.get(trans_id, 0)
                total_trans_hash[trans_id] = float(old_trans_entry) + float(read_count)
                occurence_hash[trans_id] = old_occur_entry + 1 #if read_count != '0' else old_occur_entry
                # save maximum transcript count
                if float(old_max_entry) < float(read_count):
                    max_trans_hash[trans_id] = float(read_count)
                else:
                    max_trans_hash[trans_id] = float(old_max_entry)


# divide total counted reads with the occurence of that read
for trans_id in total_trans_hash.keys():
    total_counts = total_trans_hash[trans_id]
    max_counts = max_trans_hash[trans_id]
    occurences = occurence_hash[trans_id] if total_counts >= 1.0 else 1
    counts = float(total_counts) / float(occurences)
    transcript_hash[trans_id] = counts


# modifies flux simulators .pro file
# get total molecule count
with open(sys.argv[4], 'r') as pro_in:
    total_mol_count = 0.0
    total_trans_count = 0.0
    total_mol_frac = 0.0
    for line in pro_in:
        line = line.strip()
        line_entries = line.split('\t')
        trans_id = line_entries[1]
        trans_len = line_entries[3]
        trans_frac = line_entries[4]
        trans_mol = line_entries[5]
        trans_count = transcript_hash.get(trans_id, 0.0)
        if sys.argv[2] == 'cnt':
            total_trans_count = total_trans_count + trans_count
            total_mol_count = total_mol_count + float(trans_mol)
        elif sys.argv[2] == 'fpkm':
            total_mol_frac = total_mol_frac + float(trans_frac)
            total_trans_count = total_trans_count + trans_count
            total_mol_count = int(round(total_mol_count + float(trans_mol)))
    if sys.argv[2] == 'cnt':
        total_trans_count = int(round(total_trans_count))
        total_mol_count = int(round(total_mol_count))
    print('Total molecule count: %s' % total_mol_count)
    print('Total transcript count: %s' % total_trans_count)
    print('Total frac count: %s' % total_mol_frac)

# validate counts
for trans_id in transcript_hash.keys():
    if sys.argv[2] == 'cnt':
        total_count = total_trans_hash.get(trans_id, 0.0)
        max_count = max_trans_hash.get(trans_id, 0.0)
        occurences = occurence_hash[trans_id] if total_counts >= 1.0 else 1
        occurences = float(occurences)
        trans_count = transcript_hash.get(trans_id, 0.0)
        # check if biggest transcript count is bigger than 1/10 of sum over all samples
        # and check if biggest transcript count is bigger than 1/1000 of total transcript count
        if (max_count > total_count/10) and (max_count > total_trans_count/1000):
            count_smoother = max_count / occurences
            new_count = trans_count - count_smoother
            total_trans_count = total_trans_count - count_smoother
            transcript_hash[trans_id] = new_count
            print('Correcting count %s of transcript %s to %s, because of max count %s ...' % (trans_count, trans_id, new_count, max_count))

# set scale factor
if sys.argv[2] == 'cnt':
    scale_factor = total_mol_count / total_trans_count
elif sys.argv[2] == 'fpkm':
    scale_factor = total_mol_frac / total_trans_count

# make new pro file
test_count = 0
test_frac = 0.0
with open(sys.argv[4], 'r') as pro_in, \
     open(sys.argv[5], 'w') as pro_out:
    for line in pro_in:
        line = line.strip()
        line_entries = line.split('\t')
        read_pos = line_entries[0]
        trans_id = line_entries[1]
        trans_type = line_entries[2]
        trans_len = line_entries[3]
        trans_count = transcript_hash.get(trans_id, 0.0)
        if trans_count > 0.0:
            if sys.argv[2] == 'cnt':
                mol_count = int(round(trans_count * scale_factor))
            elif sys.argv[2] == 'fpkm':
                mol_frac = trans_count * scale_factor
                mol_count = int(round(mol_frac * float(total_mol_count)))
                if mol_frac > 0 and mol_count < 1:
                    mol_count = 1
            mol_frac = float(mol_count) / float(total_mol_count)
            test_count = test_count + mol_count
            test_frac = test_frac + mol_frac
            mol_string = '%s' % mol_frac
            if 1e-04 > mol_frac:
                mol_string = re.sub(r'\-0', '-', mol_string.upper())
            elif 1e-03 > mol_frac:
                mol_string = re.sub(r'0.000(.)(.*)', r'\1.\2E-4', mol_string)
            pro_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (read_pos, trans_id, trans_type, trans_len, mol_string, mol_count))
        else:
            pro_out.write('%s\t%s\t%s\t%s\t0.0\t0\t0.0\t0\n' % (read_pos, trans_id, trans_type, trans_len))
    print('Test molecule count: %s' % test_count)
    print('Test frac count: %s' % test_frac)

