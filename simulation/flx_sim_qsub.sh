#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/sh
#$ -v CLASSPATH,LD_LIBRARY_PATH,PATH,PYTHONPATH
#$ -pe smp 4
#$ -o /home/desirda1/master_arq/flx_sim/flx_sim.qout
#$ -l m_mem_free=20G
#$ -l h_rt=86400

# reduce qsub cufflinks data
cd /home/desirda1/master_arq/flx_sim

#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -x -l -s -p /home/desirda1/master_arq/flx_sim/flx_sim.conf

# set flux simulator java environment memory
export FLUX_MEM='20G'

# call flux simulator expression step
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -x --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim.conf

# modify expression values (cufflinks fpkm files have to be transformed to counts)
#python /home/desirda1/master_arq/scripts_utility/make_flx_conf.py /home/desirda1/master_arq/arq_project_folder/quantifications/SEQC/CUFFL_quant/ cnt 1,2,3,4,5,6,7,8,17,18,19,20,21,22,23,24,33,34,35,36,37,38,39,40,49,50,51,52,53,54,55,56 /home/desirda1/master_arq/flx_sim/flx_sim.pro /home/desirda1/master_arq/flx_sim/flx_sim_A.pro > flx_sim_A.ref
#python /home/desirda1/master_arq/scripts_utility/make_flx_conf.py /home/desirda1/master_arq/arq_project_folder/quantifications/SEQC/CUFFL_quant/ cnt 9,10,11,12,13,14,15,16,25,26,27,28,29,30,31,32,41,42,43,44,45,46,47,48,57,58,59,60,61,62,63,64 /home/desirda1/master_arq/flx_sim/flx_sim.pro /home/desirda1/master_arq/flx_sim/flx_sim_B.pro > flx_sim_B.ref

# make new config files
#sed -e 's/#PRO_FILE_NAME\t/PRO_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_A.pro/' -e 's/#TMP_DIR\t/TMP_DIR\t\/home\/desirda1\/master_arq\/flx_sim\/tmp_A/' /home/desirda1/master_arq/flx_sim/flx_sim.conf > /home/desirda1/master_arq/flx_sim/flx_sim_A.conf
#sed -e 's/#PRO_FILE_NAME\t/PRO_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_B.pro/' -e 's/#TMP_DIR\t/TMP_DIR\t\/home\/desirda1\/master_arq\/flx_sim\/tmp_B/' /home/desirda1/master_arq/flx_sim/flx_sim.conf > /home/desirda1/master_arq/flx_sim/flx_sim_B.conf

# call flux simulator library step (config files with modified NB_MOLECULES, READ_NUMBER and PRO_FILE_NAME)
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -l --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim_A.conf
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -l --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim_B.conf

## for 8 calls
# make directories
#mkdir tmp_A_1 tmp_A_2 tmp_A_3 tmp_A_4 tmp_A_5 tmp_A_6 tmp_A_7 tmp_A_8
#mkdir tmp_B_1 tmp_B_2 tmp_B_3 tmp_B_4 tmp_B_5 tmp_B_6 tmp_B_7 tmp_B_8

# copy par files
#cp flx_sim_A.pro flx_sim_A_1.pro
#cp flx_sim_A.pro flx_sim_A_2.pro
#cp flx_sim_A.pro flx_sim_A_3.pro
#cp flx_sim_A.pro flx_sim_A_4.pro
#cp flx_sim_A.pro flx_sim_A_5.pro
#cp flx_sim_A.pro flx_sim_A_6.pro
#cp flx_sim_A.pro flx_sim_A_7.pro
#cp flx_sim_A.pro flx_sim_A_8.pro
#cp flx_sim_B.pro flx_sim_B_1.pro
#cp flx_sim_B.pro flx_sim_B_2.pro
#cp flx_sim_B.pro flx_sim_B_3.pro
#cp flx_sim_B.pro flx_sim_B_4.pro
#cp flx_sim_B.pro flx_sim_B_5.pro
#cp flx_sim_B.pro flx_sim_B_6.pro
#cp flx_sim_B.pro flx_sim_B_7.pro
#cp flx_sim_B.pro flx_sim_B_8.pro

# modify config file
#sed -e 's/#LIB_FILE_NAME\t/LIB_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_A.lib/' -e 's/#PRO_FILE_NAME\t/PRO_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_A_1.pro/' -e 's/#TMP_DIR\t/TMP_DIR\t\/home\/desirda1\/master_arq\/flx_sim\/tmp_A_1/' /home/desirda1/master_arq/flx_sim/flx_sim.conf > /home/desirda1/master_arq/flx_sim/flx_sim_A_1.conf
#sed -e 's/#LIB_FILE_NAME\t/LIB_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_A.lib/' -e 's/#PRO_FILE_NAME\t/PRO_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_A_2.pro/' -e 's/#TMP_DIR\t/TMP_DIR\t\/home\/desirda1\/master_arq\/flx_sim\/tmp_A_2/' /home/desirda1/master_arq/flx_sim/flx_sim.conf > /home/desirda1/master_arq/flx_sim/flx_sim_A_2.conf
#sed -e 's/#LIB_FILE_NAME\t/LIB_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_A.lib/' -e 's/#PRO_FILE_NAME\t/PRO_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_A_3.pro/' -e 's/#TMP_DIR\t/TMP_DIR\t\/home\/desirda1\/master_arq\/flx_sim\/tmp_A_3/' /home/desirda1/master_arq/flx_sim/flx_sim.conf > /home/desirda1/master_arq/flx_sim/flx_sim_A_3.conf
#sed -e 's/#LIB_FILE_NAME\t/LIB_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_A.lib/' -e 's/#PRO_FILE_NAME\t/PRO_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_A_4.pro/' -e 's/#TMP_DIR\t/TMP_DIR\t\/home\/desirda1\/master_arq\/flx_sim\/tmp_A_4/' /home/desirda1/master_arq/flx_sim/flx_sim.conf > /home/desirda1/master_arq/flx_sim/flx_sim_A_4.conf
#sed -e 's/#LIB_FILE_NAME\t/LIB_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_A.lib/' -e 's/#PRO_FILE_NAME\t/PRO_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_A_5.pro/' -e 's/#TMP_DIR\t/TMP_DIR\t\/home\/desirda1\/master_arq\/flx_sim\/tmp_A_5/' /home/desirda1/master_arq/flx_sim/flx_sim.conf > /home/desirda1/master_arq/flx_sim/flx_sim_A_5.conf
#sed -e 's/#LIB_FILE_NAME\t/LIB_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_A.lib/' -e 's/#PRO_FILE_NAME\t/PRO_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_A_6.pro/' -e 's/#TMP_DIR\t/TMP_DIR\t\/home\/desirda1\/master_arq\/flx_sim\/tmp_A_6/' /home/desirda1/master_arq/flx_sim/flx_sim.conf > /home/desirda1/master_arq/flx_sim/flx_sim_A_6.conf
#sed -e 's/#LIB_FILE_NAME\t/LIB_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_A.lib/' -e 's/#PRO_FILE_NAME\t/PRO_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_A_7.pro/' -e 's/#TMP_DIR\t/TMP_DIR\t\/home\/desirda1\/master_arq\/flx_sim\/tmp_A_7/' /home/desirda1/master_arq/flx_sim/flx_sim.conf > /home/desirda1/master_arq/flx_sim/flx_sim_A_7.conf
#sed -e 's/#LIB_FILE_NAME\t/LIB_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_A.lib/' -e 's/#PRO_FILE_NAME\t/PRO_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_A_8.pro/' -e 's/#TMP_DIR\t/TMP_DIR\t\/home\/desirda1\/master_arq\/flx_sim\/tmp_A_8/' /home/desirda1/master_arq/flx_sim/flx_sim.conf > /home/desirda1/master_arq/flx_sim/flx_sim_A_8.conf
#sed -e 's/#LIB_FILE_NAME\t/LIB_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_B.lib/' -e 's/#PRO_FILE_NAME\t/PRO_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_B_1.pro/' -e 's/#TMP_DIR\t/TMP_DIR\t\/home\/desirda1\/master_arq\/flx_sim\/tmp_B_1/' /home/desirda1/master_arq/flx_sim/flx_sim.conf > /home/desirda1/master_arq/flx_sim/flx_sim_B_1.conf
#sed -e 's/#LIB_FILE_NAME\t/LIB_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_B.lib/' -e 's/#PRO_FILE_NAME\t/PRO_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_B_2.pro/' -e 's/#TMP_DIR\t/TMP_DIR\t\/home\/desirda1\/master_arq\/flx_sim\/tmp_B_2/' /home/desirda1/master_arq/flx_sim/flx_sim.conf > /home/desirda1/master_arq/flx_sim/flx_sim_B_2.conf
#sed -e 's/#LIB_FILE_NAME\t/LIB_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_B.lib/' -e 's/#PRO_FILE_NAME\t/PRO_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_B_3.pro/' -e 's/#TMP_DIR\t/TMP_DIR\t\/home\/desirda1\/master_arq\/flx_sim\/tmp_B_3/' /home/desirda1/master_arq/flx_sim/flx_sim.conf > /home/desirda1/master_arq/flx_sim/flx_sim_B_3.conf
#sed -e 's/#LIB_FILE_NAME\t/LIB_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_B.lib/' -e 's/#PRO_FILE_NAME\t/PRO_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_B_4.pro/' -e 's/#TMP_DIR\t/TMP_DIR\t\/home\/desirda1\/master_arq\/flx_sim\/tmp_B_4/' /home/desirda1/master_arq/flx_sim/flx_sim.conf > /home/desirda1/master_arq/flx_sim/flx_sim_B_4.conf
#sed -e 's/#LIB_FILE_NAME\t/LIB_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_B.lib/' -e 's/#PRO_FILE_NAME\t/PRO_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_B_5.pro/' -e 's/#TMP_DIR\t/TMP_DIR\t\/home\/desirda1\/master_arq\/flx_sim\/tmp_B_5/' /home/desirda1/master_arq/flx_sim/flx_sim.conf > /home/desirda1/master_arq/flx_sim/flx_sim_B_5.conf
#sed -e 's/#LIB_FILE_NAME\t/LIB_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_B.lib/' -e 's/#PRO_FILE_NAME\t/PRO_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_B_6.pro/' -e 's/#TMP_DIR\t/TMP_DIR\t\/home\/desirda1\/master_arq\/flx_sim\/tmp_B_6/' /home/desirda1/master_arq/flx_sim/flx_sim.conf > /home/desirda1/master_arq/flx_sim/flx_sim_B_6.conf
#sed -e 's/#LIB_FILE_NAME\t/LIB_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_B.lib/' -e 's/#PRO_FILE_NAME\t/PRO_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_B_7.pro/' -e 's/#TMP_DIR\t/TMP_DIR\t\/home\/desirda1\/master_arq\/flx_sim\/tmp_B_7/' /home/desirda1/master_arq/flx_sim/flx_sim.conf > /home/desirda1/master_arq/flx_sim/flx_sim_B_7.conf
#sed -e 's/#LIB_FILE_NAME\t/LIB_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_B.lib/' -e 's/#PRO_FILE_NAME\t/PRO_FILE_NAME\t\/home\/desirda1\/master_arq\/flx_sim\/flx_sim_B_8.pro/' -e 's/#TMP_DIR\t/TMP_DIR\t\/home\/desirda1\/master_arq\/flx_sim\/tmp_B_8/' /home/desirda1/master_arq/flx_sim/flx_sim.conf > /home/desirda1/master_arq/flx_sim/flx_sim_B_8.conf

# run sequencing step
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -s --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim_A_1.conf
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -s --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim_A_2.conf
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -s --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim_A_3.conf
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -s --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim_A_4.conf
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -s --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim_A_5.conf
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -s --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim_A_6.conf
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -s --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim_A_7.conf
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -s --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim_A_8.conf
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -s --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim_B_1.conf
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -s --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim_B_2.conf
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -s --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim_B_3.conf
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -s --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim_B_4.conf
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -s --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim_B_5.conf
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -s --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim_B_6.conf
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -s --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim_B_7.conf
#/home/desirda1/tools/flux-simulator-1.2.1/bin/flux-simulator -t simulator -s --threads 4 -p /home/desirda1/master_arq/flx_sim/flx_sim_B_8.conf

# split and zip fastq reads
#python /home/desirda1/master_arq/scripts_utility/split_flux_reads.py /home/desirda1/master_arq/flx_sim/ flx_sim 1,2 /home/desirda1/master_arq/data/flux_simulations/
#gzip /home/desirda1/master_arq/data/flux_simulations/flx_sim_*_1-*.fastq
#gzip /home/desirda1/master_arq/data/flux_simulations/flx_sim_*_2-*.fastq
#python /home/desirda1/master_arq/scripts_utility/split_flux_reads.py /home/desirda1/master_arq/flx_sim/ flx_sim 3,4 /home/desirda1/master_arq/data/flux_simulations/
#gzip /home/desirda1/master_arq/data/flux_simulations/flx_sim_*_3-*.fastq
#gzip /home/desirda1/master_arq/data/flux_simulations/flx_sim_*_4-*.fastq
#python /home/desirda1/master_arq/scripts_utility/split_flux_reads.py /home/desirda1/master_arq/flx_sim/ flx_sim 5,6 /home/desirda1/master_arq/data/flux_simulations/
#gzip /home/desirda1/master_arq/data/flux_simulations/flx_sim_*_5-*.fastq
#gzip /home/desirda1/master_arq/data/flux_simulations/flx_sim_*_6-*.fastq
#python /home/desirda1/master_arq/scripts_utility/split_flux_reads.py /home/desirda1/master_arq/flx_sim/ flx_sim 7,8 /home/desirda1/master_arq/data/flux_simulations/
#gzip /home/desirda1/master_arq/data/flux_simulations/flx_sim_*_7-*.fastq
#gzip /home/desirda1/master_arq/data/flux_simulations/flx_sim_*_8-*.fastq

