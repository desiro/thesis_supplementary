#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/sh
#$ -v CLASSPATH,LD_LIBRARY_PATH,PATH,PYTHONPATH
#$ -pe smp 1
#$ -o /home/desirda1/master_arq/cuffl_qsub_out/cuffl_qsub.out
#$ -l m_mem_free=10G
#$ -l h_rt=14400

# get cufflinks read mass from .out file
grep 'Normalized Map Mass' /home/desirda1/master_arq/arq_project_folder/qsub_data/CUFFL_quant-*.qout | sed 's/>.*:\s//' | sed 's/\.\..*qsub_data\///' | sed 's/\(.*\)\-\(.*\)_\(.*\)_\(.*\)\.qout:/\1\t\2\t\3\t\4\t/' | sort -k2,2 -k3,3n > cuff_read_mass.data
