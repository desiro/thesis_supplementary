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
from arq_debug import printDebug




###############################################################################
## make config for main dictionaries
###############################################################################

def getMainDictionaries(dict_type, debug):
    
    # initialize debug
    debug_list = []
    format = '{:<17}'
    if debug: debug_list.append(['def: getMainDictionaries():'])
    if debug: debug_list.append(['-'])
    if debug: debug_list.append(['done'])
    
    # initialize return dictionary
    return_dict = {}
    
    # initialize dictionary with all possible tool parameters
    # -> only these parameters will be accepted in the config file
    # -> bash_script <-(1)- % bash_config <-(2)- % input_parameters
    ''' TODO: save ALL possible variables and update empty checks in code '''
    if dict_type == 'input_parameters':
        return_dict = { # input options
                        'proj_dir_prefix':'',   'data_set_name':'',     'fasta_gen_file':'',    'fasta_trns_file':'', 
                        'gtf_index_file':'',    'fastq1_links':'',      'fastq2_links':'',      'arq_job_list':'', 
                        'quant_comb_list':'',   'quant_name_list':'', 
                        # tool paths
                        'tool_dir_star':'',     'tool_dir_tophat':'',   'tool_dir_hisat':'',    'tool_dir_bowtie':'', 
                        'tool_dir_htseq':'',    'tool_dir_featc':'',    'tool_dir_eqpqm':'',    'tool_dir_toprec':'', 
                        'tool_dir_flxcp':'',    'tool_dir_kllst':'',    'tool_dir_cuffl':'',    'tool_dir_btseq':'', 
                        'tool_dir_rsem':'', 
                        # tool parameters: alignments
                        'num_processors':'',    'star_overhang':'',     'min_intron_len':'',    'max_intron_len':'', 
                        'max_multi_hits':'',    'mate_dist_opt':'',     'max_mismatches':'',    'star_max_edit_r':'', 
                        'star_gen_load':'',     'index_base_name':'',   'bwtrs_in_opt':'',      'bwtrs_in_strand':'',
                        'bwtrs_rsem_opt':'',
                        # tool parameters: quantifications
                        'quant_feat_type':'',   'quant_feat_id':'',     'htseq_in_strand':'',   'htseq_mode_opt':'', 
                        'featc_strand_1':'',    'featc_strand_2':'',    'featc_overlap':'',     'featc_pair_val':'', 
                        'eqpqm_in_strand':'',   'eqpqm_count_opt':'',   'eqpqm_unambig':'',     'eqpqm_no_weight':'', 
                        'cuffl_lib_type':'',    'cuffl_no_update':'',   'cuffl_quiet_mode':'',  'flxcp_tool_call':'', 
                        'flxcp_tmp_dir':'',     'btseq_parse_call':'',  'btseq_estim_call':'',  'btseq_count_call':'', 
                        'btseq_prob_out':'',    'btseq_prob_in':'',     'btseq_info_file':'',   'btseq_out_count':'', 
                        'btseq_prob_folder':'', 'btseq_thet_file':'',   'btseq_seed_opt':'',    'btseq_verb_opt':'',   
                        'btseq_strand_opt':'',  'btseq_out_type':'',    'kallst_seed_opt':'',   'kallst_out_frmt':'', 
                        # management options
                        'save_bash_call':'',    'remove_fastq':'',      'rem_comb_align':'',    'single_tool_job':'', 
                        'seq_style_fq':'',      'only_make_job':'',     'user_set_over':'',     'user_set_skip':'', 
                        'user_set_new':'',      'ignore_qs_error':'',   'sleep_time_qs':'',     'qs_mem_free':'',
                        'qs_run_time':''}
    
    # initialize dictionbary with all possible jobs and dependencies as comma separated list
    # -> only these job request will be accepted in the config file
    ''' TODO: make dynamic and make possible to add jobs via config file '''
    if dict_type == 'possible_jobs':
        return_dict = { 'STAR_index':'',                                'STAR_align':'STAR_index',
                        'TOPHAT_index':'BOWTIE_index',                  'TOPHAT_align':'TOPHAT_index',
                        'HISAT_index':'',                               'HISAT_align':'HISAT_index',
                        'BOWTIE_index':'',                              'HTSEQ_quant':'',
                        'EQPQM_index':'',                               'EQPQM_quant':'EQPQM_index',
                        'FLXCP_index':'',                               'FLXCP_quant':'FLXCP_index', 
                        'KLLST_index':'',                               'KLLST_quant':'KLLST_index',
                        'BWTRS_index':'',                               'BWTRS_align':'BWTRS_index',
                        'RSEM_index':'',                                'RSEM_quant':'RSEM_index',
                        'BTSEQ_quant':'',                               'BRSEM_align':'BWTRS_index', 
                        'CUFFL_quant':'',                               'FEATC_quant':'', 
                        'AFREE_align':''}
    
    # input_parameters that will be initialized automatically
    # placeholder: *A* = aligner name       *a* = align
    #              *Q* = quantifier name    *q* = quant
    #              t c p = temp cycling parameter (parameter will iterate through all fastq_name/*A*_dir/etc.)
    #
    #   parameter     description                                           complete directory path
    # -------------------------------------------------------------------------------------------------------------------------------------------
    #   <p_name_suffix>         command line project name                   -
    #   <proj_dir_path>         main project path                           <proj_dir_prefix>/<p_name_suffix>
    #   <index_dir_path>        path to index directory                     <proj_dir_path>/indexes/<data_set_name>
    #   <align_dir_path>        path to alignment directory                 <proj_dir_path>/alignments/<data_set_name>
    #   <quant_dir_path>        path to quantification directory            <proj_dir_path>/quantifications/<data_set_name>
    #   <comb_dir_path>         path to combined alignments directory       <proj_dir_path>/combinations/<data_set_name>
    #   <*A*_index_dir>         path to aligner index directory             <index_dir_path>/*A*_index
    #   <*A*_align_dir>         path to aligner alignment directory         <align_dir_path>/*A*_align
    #   <*A*_*Q*_quant_dir>     path to aligner-quantifier directory        <quant_dir_path>/*A*_*Q*_quant
    #   <align_out_prefix>      t c p for alignment names                   -
    #   <quant_out_prefix>      t c p for quantification names              -
    #   <align_out_dir>         t c p for fastq alignments                  <align_dir_path/*A*_align/<align_out_prefix>
    #   <quant_out_dir>         t c p for combi quant                       <quant_dir_path/*Q*_quant/<align_out_prefix>
    #   <fastq1_file>           t c p for fastq1 files                      -
    #   <fastq2_file>           t c p for fastq2 files                      -
    #   <ALIGNER_align_dir>     t c p for specific alignment dir            <align_dir_path/*A*_align>']
    #   <alignment_prefix>      t c p for specific alignment first fq name  -
    #   <combi_align_list>      list of all alignment.bam files for quant   <ALIGNER_align_dir>/<align_out_prefix>1.bam <ALIGNER_align_dir>/<align_out_prefix>2.bam ...
    
    # initialize list with all possible alignment and quantification combinations
    if dict_type == 'quant_align':
        return_dict = [ 'HTSEQ_STAR', 'HTSEQ_TOPHAT', 'HTSEQ_HISAT', 'FEATC_STAR', 'FEATC_TOPHAT', 'FEATC_HISAT', 
                        'EQPQM_STAR', 'EQPQM_TOPHAT', 'EQPQM_HISAT', 'CUFFL_STAR', 'CUFFL_TOPHAT', 'CUFFL_HISAT', 
                        'FLXCP_STAR', 'FLXCP_TOPHAT', 'FLXCP_HISAT', 'BTSEQ_BWTRS', 'RSEM_BRSEM', 'KLLST_AFREE']
    
    # return dict
    return return_dict



###############################################################################
## make config for index jobs
###############################################################################

def getIndexConfig(job_config_data, debug):
    
    # initialize debug
    debug_list = []
    format = '{:<17}'
    if debug: debug_list.append(['def: getIndexConfig():'])
    if debug: debug_list.append(['-'])
    if debug: debug_list.append(['done'])
    
    ## initialize known index data
    # initialize STAR index
    job_config_data['STAR_index'] = {       'script_order'          :['tool_dir_star', 'options', 'num_processors', 'index_base_name', 'fasta_gen_file', 'gtf_index_file', 'star_overhang'], 
                                            'script_options'        :[], 
                                            'tool_dir_star'         :'%(tool_dir_star)sSTAR --runMode genomeGenerate', 
                                            'num_processors'        :'--runThreadN %(num_processors)s', 
                                            'index_base_name'       :'--genomeDir %(STAR_index_dir)s',
                                            'fasta_gen_file'        :'--genomeFastaFiles %(fasta_gen_file)s', 
                                            'gtf_index_file'        :'--sjdbGTFfile %(gtf_index_file)s', 
                                            'star_overhang'         :'--sjdbOverhang %(star_overhang)s', 
                                            'bash_pre_script'       :'cd %(STAR_index_dir)s', 
                                            'bash_seq_script'       :''}
    # initialize TOPHAT index
    job_config_data['TOPHAT_index'] = {     'script_order'          :['tool_dir_tophat', 'options', 'gtf_index_file', 'index_base_name', 'bowtie_index_base'], 
                                            'script_options'        :[], 
                                            'tool_dir_tophat'       :'%(tool_dir_tophat)stophat2', 
                                            'bowtie_index_base'     :'%(BOWTIE_index_dir)s/%(index_base_name)s', 
                                            'gtf_index_file'        :'-G %(gtf_index_file)s', 
                                            'index_base_name'       :'--transcriptome-index %(TOPHAT_index_dir)s/%(index_base_name)s', 
                                            'bash_pre_script'       :'cd %(TOPHAT_index_dir)s', 
                                            'bash_seq_script'       :''}
    # initialize HISAT index
    job_config_data['HISAT_index'] = {      'script_order'          :['tool_dir_hisat', 'options', 'fasta_gen_file', 'index_base_name'], 
                                            'script_options'        :[], 
                                            'tool_dir_hisat'        :'%(tool_dir_hisat)shisat-build', 
                                            'index_base_name'       :'%(HISAT_index_dir)s/%(index_base_name)s', 
                                            'fasta_gen_file'        :'-f %(fasta_gen_file)s', 
                                            'bash_pre_script'       :'cd %(HISAT_index_dir)s && '
                                                                    +'python %(tool_dir_hisat)s/extract_splice_sites.py %(gtf_index_file)s > %(HISAT_index_dir)s/%(index_base_name)s.splice', 
                                            'bash_seq_script'       :''}
    # initialize BOWTIE index
    job_config_data['BOWTIE_index'] = {     'script_order'          :['tool_dir_bowtie', 'options', 'fasta_gen_file', 'index_base_name'], 
                                            'script_options'        :[], 
                                            'tool_dir_bowtie'       :'%(tool_dir_bowtie)sbowtie2-build', 
                                            'index_base_name'       :'%(BOWTIE_index_dir)s/%(index_base_name)s', 
                                            'fasta_gen_file'        :'-f %(fasta_gen_file)s', 
                                            'bash_pre_script'       :'cd %(BOWTIE_index_dir)s && '
                                                                    +'cp %(fasta_gen_file)s %(BOWTIE_index_dir)s/%(index_base_name)s.fa', 
                                            'bash_seq_script'       :''}
    # initialize BOWTIE index
    job_config_data['BWTRS_index'] = {      'script_order'          :['tool_dir_bowtie', 'options', 'fasta_trns_file', 'index_base_name'], 
                                            'script_options'        :[], 
                                            'tool_dir_bowtie'       :'%(tool_dir_bowtie)sbowtie2-build', 
                                            'index_base_name'       :'%(BWTRS_index_dir)s/%(index_base_name)s', 
                                            'fasta_trns_file'       :'-f %(fasta_trns_file)s', 
                                            'bash_pre_script'       :'cd %(BWTRS_index_dir)s && '
                                                                    +'cp %(fasta_trns_file)s %(BWTRS_index_dir)s/%(index_base_name)s.fa', 
                                            'bash_seq_script'       :''}
    # initialize BOWTIE index
    job_config_data['EQPQM_index'] = {      'script_order'          :['tool_dir_eqpqm', 'options', 'gtf_index_file', 'index_base_name'], 
                                            'script_options'        :[], 
                                            'tool_dir_eqpqm'        :'%(tool_dir_eqpqm)seqp-setup.sh', 
                                            'gtf_index_file'        :'%(gtf_index_file)s', 
                                            'index_base_name'       :'%(EQPQM_index_dir)s', 
                                            'bash_pre_script'       :'cd %(EQPQM_index_dir)s', 
                                            'bash_seq_script'       :''}
    # initialize fluxCapacitor index
    job_config_data['FLXCP_index'] = {      'script_order'          :['tool_dir_flxcp', 'options', 'flxcp_tool_call', 'gtf_index_file', 'index_base_name'], 
                                            'script_options'        :[], 
                                            'tool_dir_flxcp'        :'%(tool_dir_flxcp)sflux-capacitor', 
                                            'gtf_index_file'        :'-i %(gtf_index_file)s', 
                                            'index_base_name'       :'-o %(FLXCP_index_dir)s/%(index_base_name)s.gtf', 
                                            'flxcp_tool_call'       :'-t sortGTF --force', 
                                            'bash_pre_script'       :'cd %(FLXCP_index_dir)s', 
                                            'bash_seq_script'       :''}
    # initialize rsem index
    job_config_data['RSEM_index'] = {       'script_order'          :['tool_dir_rsem', 'options', 'fasta_trns_file', 'index_base_name'], 
                                            'script_options'        :[], 
                                            'tool_dir_rsem'         :'%(tool_dir_rsem)srsem-prepare-reference', 
                                            'fasta_trns_file'       :'%(fasta_trns_file)s', 
                                            'index_base_name'       :'%(RSEM_index_dir)s/%(index_base_name)s', 
                                            'bash_pre_script'       :'cd %(RSEM_index_dir)s', 
                                            'bash_seq_script'       :''}
    # initialize kallisto index
    job_config_data['KLLST_index'] = {      'script_order'          :['tool_dir_kllst', 'options', 'index_base_name', 'fasta_trns_file'], 
                                            'script_options'        :[], 
                                            'tool_dir_kllst'        :'%(tool_dir_kllst)skallisto index', 
                                            'index_base_name'       :'-i %(KLLST_index_dir)s/%(index_base_name)s.idx', 
                                            'fasta_trns_file'        :'%(fasta_trns_file)s', 
                                            'bash_pre_script'       :'cd %(KLLST_index_dir)s', 
                                            'bash_seq_script'       :''}
    # # initialize iReckon index
    # job_config_data['IRECK_index'] = {      'script_order'          :['tool_dir_savant', 'options', 'gtf_index_file', 'index_base_name'], 
                                            # 'script_options'        :[], 
                                            # 'tool_dir_savant'       :'%(tool_dir_savant)sFormatTool.sh', 
                                            # 'gtf_index_file'        :'%(gtf_index_file)s', 
                                            # 'index_base_name'       :'%(IRECK_index_dir)s/%(index_base_name)s.tbx', 
                                            # 'bash_pre_script'       :'', 
                                            # 'bash_seq_script'       :''}
    # # initialize eXpress index -> use bedtools to generate multi line fasta file from genome fasta and gtf file
    # job_config_data['XPRSS_index'] = {      'script_order'          :['tool_dir_savant', 'options', 'gtf_index_file', 'index_base_name'], 
                                            # 'script_options'        :[], 
                                            # 'tool_dir_savant'       :'%(tool_dir_savant)sFormatTool.sh', 
                                            # 'gtf_index_file'        :'%(gtf_index_file)s', 
                                            # 'index_base_name'       :'%(IRECK_index_dir)s/%(index_base_name)s.tbx', 
                                            # 'bash_pre_script'       :'', 
                                            # 'bash_seq_script'       :''}
    
    # print debug
    if debug: printDebug(debug_list, format)
    
    # return job config
    return job_config_data




###############################################################################
## make config for alignment jobs
###############################################################################

def getAlignConfig(job_config_data, debug):
    
    # initialize debug
    debug_list = []
    format = '{:<17}'
    if debug: debug_list.append(['def: getAlignConfig():'])
    if debug: debug_list.append(['-'])
    if debug: debug_list.append(['done'])
    
    ## initialize known alignment data
    # initialize STAR alignment
    job_config_data['STAR_align'] = {       'script_order'          :['tool_dir_star', 'options', 'keep_input_order', 'num_processors', 'align_out_option', 'index_base_name', 'fastq1_file', 'fastq2_file'], 
                                            'script_options'        :['star_gen_load', 'min_intron_len', 'max_intron_len', 'max_multi_hits', 'mate_dist_opt', 'max_mismatches', 'star_max_edit_r'], 
                                            'tool_dir_star'         :'%(tool_dir_star)sSTAR', 
                                            'num_processors'        :'--runThreadN %(num_processors)s', 
                                            'align_out_option'      :'--outFileNamePrefix %(align_out_dir)s/%(align_out_prefix)s', 
                                            'index_base_name'       :'--genomeDir %(STAR_index_dir)s', 
                                            'fastq1_file'           :'--readFilesIn %(fastq1_file)s', 
                                            'fastq2_file'           :'%(fastq2_file)s', 
                                            'star_gen_load'         :'--genomeLoad %(star_gen_load)s', 
                                            'keep_input_order'      :'--outSAMorder PairedKeepInputOrder', 
                                            'min_intron_len'        :'--alignIntronMin %(min_intron_len)s', 
                                            'max_intron_len'        :'--alignIntronMax %(max_intron_len)s', 
                                            'max_multi_hits'        :'--outFilterMultimapNmax %(max_multi_hits)s', # discard when max multihits
                                            'mate_dist_opt'         :'--alignMatesGapMax %(mate_dist_opt)s', # mate gap max length
                                            'max_mismatches'        :'--outFilterMismatchNmax %(max_mismatches)s', 
                                            'star_max_edit_r'       :'--outFilterMismatchNoverLmax %(star_max_edit_r)s', 
                                            'bash_pre_script'       :'cd %(align_out_dir)s', 
                                            'bash_seq_script'       :'samtools view -Shb %(align_out_dir)s/%(align_out_prefix)sAligned.out.sam | samtools sort -n -m 20000000000 - %(STAR_align_dir)s/%(align_out_prefix)s'}
    # initialize TOPHAT alignment
    job_config_data['TOPHAT_align'] = {     'script_order'          :['tool_dir_tophat', 'options', 'micro_exon_search', 'keep_input_order', 'num_processors', 'align_out_option', 'index_base_name', 'bowtie_index_base', 'fastq1_file', 'fastq2_file'], 
                                            'script_options'        :['min_intron_len', 'max_intron_len', 'max_multi_hits', 'mate_dist_opt', 'max_mismatches'], 
                                            'tool_dir_tophat'       :'%(tool_dir_tophat)stophat2', 
                                            'num_processors'        :'--num-threads %(num_processors)s', 
                                            'align_out_option'      :'--output-dir %(align_out_dir)s', 
                                            'index_base_name'       :'--transcriptome-index %(TOPHAT_index_dir)s/%(index_base_name)s', 
                                            'bowtie_index_base'     :'%(BOWTIE_index_dir)s/%(index_base_name)s', 
                                            'fastq1_file'           :'%(fastq1_file)s', 
                                            'fastq2_file'           :'%(fastq2_file)s', 
                                            'micro_exon_search'     :'--microexon-search', 
                                            'keep_input_order'      :'--keep-fasta-order', 
                                            'min_intron_len'        :'--min-intron-length %(min_intron_len)s', 
                                            'max_intron_len'        :'--max-intron-length %(max_intron_len)s', 
                                            'max_multi_hits'        :'--max-multihits %(max_multi_hits)s', # save this many multihits
                                            'mate_dist_opt'         :'--mate-inner-dist %(mate_dist_opt)s', # mate gap medium length
                                            'max_mismatches'        :'--read-mismatches %(max_mismatches)s '
                                                                    +'--read-gap-length %(max_mismatches)s '
                                                                    +'--read-edit-dist %(max_mismatches)s', # both need to be less or equal to --read-edit-dist
                                            'bash_pre_script'       :'cd %(align_out_dir)s', 
                                            'bash_seq_script'       :'if [ \"%(tophat_merge_opt)s\" = "TRUE" ]; '
                                                                    +'then '
                                                                    +   'python %(tool_dir_toprec)stophat-recondition.py -q %(align_out_dir)s && '
                                                                    +   'samtools view -H %(align_out_dir)s/unmapped_fixup.bam > %(align_out_dir)s/unmapped_fixup-header.sam && '
                                                                    +   'samtools merge -h %(align_out_dir)s/unmapped_fixup-header.sam %(align_out_dir)s/merged_hits-unsorted.bam %(align_out_dir)s/accepted_hits.bam %(align_out_dir)s/unmapped_fixup.bam && '
                                                                    +   'samtools sort -n -m 20000000000 %(align_out_dir)s/merged_hits-unsorted.bam %(TOPHAT_align_dir)s/%(align_out_prefix)s; '
                                                                    +'else '
                                                                    +   'samtools sort -n -m 20000000000 %(align_out_dir)s/accepted_hits.bam %(TOPHAT_align_dir)s/%(align_out_prefix)s; '
                                                                    +'fi'}
    # initialize HISAT alignment
    job_config_data['HISAT_align'] = {      'script_order'          :['tool_dir_hisat', 'options', 'use_phred_33', 'keep_input_order', 'num_processors', 'align_out_option', 'index_base_name', 'known_splice_file', 'fastq1_file', 'fastq2_file'], 
                                            'script_options'        :['min_intron_len', 'max_intron_len', 'max_multi_hits', 'mate_dist_opt'], 
                                            'tool_dir_hisat'        :'%(tool_dir_hisat)shisat', 
                                            'num_processors'        :'--threads %(num_processors)s', 
                                            'align_out_option'      :'-S %(align_out_dir)s/%(align_out_prefix)s.sam', 
                                            'index_base_name'       :'-x %(HISAT_index_dir)s/%(index_base_name)s', 
                                            'known_splice_file'     :'--known-splicesite-infile %(HISAT_index_dir)s/%(index_base_name)s.splice', 
                                            'fastq1_file'           :'-1 %(fastq1_file)s', 
                                            'fastq2_file'           :'-2 %(fastq2_file)s', 
                                            'use_phred_33'          :'--phred33', 
                                            'keep_input_order'      :'--reorder', 
                                            'min_intron_len'        :'--min-intronlen %(min_intron_len)s', 
                                            'max_intron_len'        :'--max-intronlen %(max_intron_len)s', 
                                            'max_multi_hits'        :'-k %(max_multi_hits)s', 
                                            'mate_dist_opt'         :'--maxins $MATE_TOTAL_LENGTH', 
                                            'bash_pre_script'       :'cd %(align_out_dir)s && '
                                                                    +'MATE_TOTAL_LENGTH=`expr %(mate_dist_opt)s + 200`', 
                                            'bash_seq_script'       :'samtools view -Shb %(align_out_dir)s/%(align_out_prefix)s.sam | samtools sort -n -m 20000000000 - %(HISAT_align_dir)s/%(align_out_prefix)s'}
    # initialize BOWTIE transcript alignment
    job_config_data['BWTRS_align'] = {      'script_order'          :['tool_dir_bowtie', 'options', 'num_processors', 'bwtrs_in_opt', 'index_base_name', 'fastq1_file', 'fastq2_file', 'align_out_option'], 
                                            'script_options'        :['mate_dist_opt', 'bwtrs_in_strand'], 
                                            'tool_dir_bowtie'       :'%(tool_dir_bowtie)sbowtie2', 
                                            'num_processors'        :'--threads %(num_processors)s',
                                            'align_out_option'      :'-S %(align_out_dir)s/%(align_out_prefix)s.sam',
                                            'index_base_name'       :'-x %(BWTRS_index_dir)s/%(index_base_name)s', 
                                            'fastq1_file'           :'-1 %(fastq1_file)s', 
                                            'fastq2_file'           :'-2 %(fastq2_file)s', 
                                            'bwtrs_in_opt'          :'-q',
                                            'mate_dist_opt'         :'--maxins $MATE_TOTAL_LENGTH', 
                                            'bwtrs_in_strand'       :'--%(bwtrs_in_strand)s',
                                            'bash_pre_script'       :'cd %(align_out_dir)s && '
                                                                    +'MATE_TOTAL_LENGTH=`expr %(mate_dist_opt)s + 200`',
                                            'bash_seq_script'       :'samtools view -Shb %(align_out_dir)s/%(align_out_prefix)s.sam | samtools sort -n -m 20000000000 - %(BWTRS_align_dir)s/%(align_out_prefix)s'}
    # initialize BOWTIE RSEM transcript alignment
    job_config_data['BRSEM_align'] = {      'script_order'          :['tool_dir_bowtie', 'options', 'num_processors', 'bwtrs_rsem_opt', 'bwtrs_in_opt', 'index_base_name', 'fastq1_file', 'fastq2_file', 'align_out_option'], 
                                            'script_options'        :['mate_dist_opt', 'bwtrs_in_strand'], 
                                            'tool_dir_bowtie'       :'%(tool_dir_bowtie)sbowtie2', 
                                            'num_processors'        :'--threads %(num_processors)s',
                                            'align_out_option'      :'-S %(align_out_dir)s/%(align_out_prefix)s.sam',
                                            'index_base_name'       :'-x %(BWTRS_index_dir)s/%(index_base_name)s', 
                                            'fastq1_file'           :'-1 %(fastq1_file)s', 
                                            'fastq2_file'           :'-2 %(fastq2_file)s', 
                                            'bwtrs_in_opt'          :'-q',
                                            'bwtrs_rsem_opt'        :'--no-discordant --no-mixed --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1',
                                            'mate_dist_opt'         :'--maxins $MATE_TOTAL_LENGTH', 
                                            'bwtrs_in_strand'       :'--%(bwtrs_in_strand)s',
                                            'bash_pre_script'       :'cd %(align_out_dir)s && '
                                                                    +'MATE_TOTAL_LENGTH=`expr %(mate_dist_opt)s + 200`',
                                            'bash_seq_script'       :'samtools view -Shb %(align_out_dir)s/%(align_out_prefix)s.sam | samtools sort -n -m 20000000000 - %(BRSEM_align_dir)s/%(align_out_prefix)s'}
    
    
    # option comparison
    # tophat                            star                                hisat                       bowtie2
    # ------------------------------------------------------------------------------------------------------------------------------
    # --output-dir                      --outFileNamePrefix                 -S .sam                     
    # --num-threads 1                   --runThreadN 1                      -p/--threads 1              
    # --min-intron-length 70            --alignIntronMin 21                 --min-intronlen 20          
    # --max-intron-length 500000        --alignIntronMax 589824             --max-intronlen 500000      
    # --max-multihits 20                --outFilterMultimapNmax 10          -k 5                        
    # --mate-inner-dist 50              --alignMatesGapMax 589824           --minins 0 --maxins 500     --maxins 500
    # --keep-fasta-order                --outSAMorder PairedKeepInputOrder  --reorder
    # --read-mismatches 2               --outFilterMismatchNmax 10          
    #                                   --outFilterMismatchNoverLmax 0.3    
    # --microexon-search
    #                                                                       --phred33
    #                                       
    
    # print debug
    if debug: printDebug(debug_list, format)
    
    # return job config
    return job_config_data




###############################################################################
## make config for quantification jobs
###############################################################################

def getQuantConfig(job_config_data, debug):
    
    # initialize debug
    debug_list = []
    format = '{:<17}'
    if debug: debug_list.append(['def: getQuantConfig():'])
    if debug: debug_list.append(['-'])
    if debug: debug_list.append(['done'])
    
    ## initialize known quantification data
    # initialize htseq-count quantification
    job_config_data['HTSEQ_quant'] = {      'script_order'          :['tool_dir_htseq', 'options', 'quant_sam_out', 'quant_in_order', 'quant_in_format', 'comb_align_file', 'gtf_index_file', 'quant_out_option'], 
                                            'script_options'        :['htseq_in_strand', 'quant_feat_type', 'quant_feat_id', 'htseq_mode_opt'], 
                                            'tool_dir_htseq'        :'%(tool_dir_htseq)shtseq-count', 
                                            'comb_align_file'       :'%(comb_in_file)s.tmp.sam', 
                                            'gtf_index_file'        :'%(gtf_index_file)s', 
                                            'quant_out_option'      :'> %(quant_out_dir)s/%(quant_out_prefix)s.quant',
                                            'quant_in_format'       :'--format bam',
                                            'quant_in_order'        :'--order name',
                                            'quant_feat_type'       :'--type %(quant_feat_type)s',
                                            'quant_feat_id'         :'--idattr %(quant_feat_id)s', # gene count -> gene_id ; exon count -> exon_id
                                            'htseq_in_strand'       :'--stranded %(htseq_in_strand)s', # yes/no/reverse
                                            'htseq_mode_opt'        :'--mode %(htseq_mode_opt)s', # union/intersection-strict/intersection-nonempty
                                            'quant_sam_out'         :'--samout %(quant_out_dir)s/%(quant_out_prefix)s',
                                            'bash_pre_script'       :'cd %(quant_out_dir)s && '
                                                                    +'samtools view -h %(comb_in_file)s > %(comb_in_file)s.tmp.sam',
                                            'bash_seq_script'       :'rm -f %(comb_in_file)s.tmp.sam && '
                                                                    +'cp %(quant_out_dir)s/%(quant_out_prefix)s.quant %(quant_dir_path)s/HTSEQ_quant/%(quant_out_prefix)s.cnt'}
    # initialize feature-count quantification
    job_config_data['FEATC_quant'] = {      'script_order'          :['tool_dir_featc', 'options', 'quant_in_order', 'quant_out_paired', 'gtf_index_file', 'quant_out_option', 'comb_align_file'], 
                                            'script_options'        :['quant_feat_type', 'quant_feat_id', 'featc_strand_1', 'featc_strand_2', 'featc_overlap', 'featc_pair_val'], 
                                            'tool_dir_featc'        :'%(tool_dir_featc)sfeatureCounts', 
                                            'comb_align_file'       :'%(comb_in_file)s', 
                                            'gtf_index_file'        :'-a %(gtf_index_file)s', 
                                            'quant_out_option'      :'-o %(quant_out_dir)s/%(quant_out_prefix)s.quant',
                                            'quant_out_paired'      :'-p',
                                            'quant_in_order'        :'--donotsort',
                                            'quant_feat_type'       :'-t %(quant_feat_type)s',
                                            'quant_feat_id'         :'-g %(quant_feat_id)s',
                                            'featc_strand_1'        :'-s %(featc_strand_1)s', # -s 0/1/2
                                            'featc_strand_2'        :'-S %(featc_strand_2)s', # -S ff/fr/rf -> ff
                                            'featc_overlap'         :'-O -M --fraction',
                                            'featc_pair_val'        :'-C',
                                            'bash_pre_script'       :'cd %(quant_out_dir)s',
                                            'bash_seq_script'       :'cp %(quant_out_dir)s/%(quant_out_prefix)s.quant %(quant_dir_path)s/FEATC_quant/%(quant_out_prefix)s.cnt'}
    # initialize EQP-QM quantification
    job_config_data['EQPQM_quant'] = {      'script_order'          :['tool_dir_eqpqm', 'options', 'quant_in_order', 'quant_out_string', 'index_base_name', 'quant_out_option', 'comb_align_file'], 
                                            'script_options'        :['eqpqm_in_strand', 'eqpqm_count_opt', 'eqpqm_unambig', 'eqpqm_no_weight'], 
                                            'tool_dir_eqpqm'        :'%(tool_dir_eqpqm)seqp-quantify.sh', 
                                            'comb_align_file'       :'%(comb_in_file)s', 
                                            'index_base_name'       :'-d %(EQPQM_index_dir)s', 
                                            'quant_out_option'      :'%(quant_out_dir)s',
                                            'quant_out_string'      :'-o %(quant_out_prefix)s',
                                            'quant_in_order'        :'--nosort',
                                            'eqpqm_in_strand'       :'-s %(eqpqm_in_strand)s', # forward/backward (fr/rf) -> don't use for unstranded alignment
                                            'eqpqm_count_opt'       :'%(eqpqm_count_opt)s', 
                                            'eqpqm_unambig'         :'--unambig', 
                                            'eqpqm_no_weight'       :'--unweighted', 
                                            'bash_pre_script'       :'cd %(quant_out_dir)s',
                                            'bash_seq_script'       :'cp %(quant_out_dir)s/%(quant_out_prefix)s-gene.cnt %(quant_dir_path)s/EQPQM_quant/%(quant_out_prefix)s.cnt'}
    # initialize cufflinks
    job_config_data['CUFFL_quant'] = {      'script_order'          :['tool_dir_cuffl', 'options', 'cuffl_no_update', 'cuffl_quiet_mode', 'gtf_index_file', 'quant_out_option', 'num_processors', 'comb_align_file'], 
                                            'script_options'        :['cuffl_lib_type'], 
                                            'tool_dir_cuffl'        :'%(tool_dir_cuffl)scufflinks', 
                                            'comb_align_file'       :'%(quant_out_dir)s/%(quant_out_prefix)s.bam', 
                                            'quant_out_option'      :'-o %(quant_out_dir)s', 
                                            'gtf_index_file'        :'-G %(gtf_index_file)s', 
                                            'num_processors'        :'--num-threads %(num_processors)s',
                                            'cuffl_lib_type'        :'--library-type %(cuffl_lib_type)s', # ff-unstranded
                                            'cuffl_no_update'       :'--no-update-check', 
                                            'cuffl_quiet_mode'      :'--quiet', 
                                            'min_intron_len'        :'--min-intron-length %(min_intron_len)s', 
                                            'bash_pre_script'       :'cd %(quant_out_dir)s && '
                                                                    +'samtools sort -m 20000000000 %(comb_in_file)s %(quant_out_dir)s/%(quant_out_prefix)s', 
                                            'bash_seq_script'       :'cp %(quant_out_dir)s/isoforms.fpkm_tracking %(quant_dir_path)s/CUFFL_quant/%(quant_out_prefix)s.fpkm && '
                                                                    +'rm -f %(quant_out_dir)s/%(quant_out_prefix)s.bam'}
    # initialize fluxCapacitor
    job_config_data['FLXCP_quant'] = {      'script_order'          :['tool_dir_flxcp', 'flxcp_tool_call', 'flxcp_tmp_dir', 'quant_out_option', 'comb_align_file', 'index_base_name', 'options'], 
                                            'script_options'        :['min_intron_len'], 
                                            'tool_dir_flxcp'        :'%(tool_dir_flxcp)sflux-capacitor', 
                                            'flxcp_tool_call'       :'-t capacitor --force', 
                                            'flxcp_tmp_dir'         :'--tmp-dir %(quant_out_dir)s', 
                                            'quant_out_option'      :'-o %(quant_out_dir)s/%(quant_out_prefix)s.quant', 
                                            'comb_align_file'       :'-i %(quant_out_dir)s/%(quant_out_prefix)s.bam', 
                                            'index_base_name'       :'-a %(FLXCP_index_dir)s/%(index_base_name)s.gtf', 
                                            'min_intron_len'        :'--minilen %(min_intron_len)s', 
                                            'bash_pre_script'       :'cd %(quant_out_dir)s && '
                                                                    +'samtools sort -m 20000000000 %(comb_in_file)s %(quant_out_dir)s/%(quant_out_prefix)s', 
                                            'bash_seq_script'       :'cp %(quant_out_dir)s/%(quant_out_prefix)s.quant %(quant_dir_path)s/FLXCP_quant/%(quant_out_prefix)s.cnt && '
                                                                    +'rm -f %(quant_out_dir)s/%(quant_out_prefix)s.bam'}
    # initialize BitSeq
    job_config_data['BTSEQ_quant'] = {      'script_order'          :['btseq_parse_call', 'btseq_prob_out', 'fasta_trns_file', 'btseq_info_file', 'num_processors', 'btseq_strand_opt', 'btseq_verb_opt', 'comb_align_file', 'btseq_estim_call', 'quant_out_option', 'btseq_out_type', 'num_processors', 'btseq_seed_opt', 'btseq_verb_opt', 'options', 'btseq_prob_in', 'btseq_count_call', 'btseq_out_count', 'btseq_prob_folder', 'btseq_thet_file'], 
                                            'script_options'        :[], 
                                            'btseq_parse_call'      :'%(tool_dir_btseq)sparseAlignment',
                                            'btseq_prob_out'        :'--outFile=%(quant_out_dir)s/%(quant_out_prefix)s.prob',
                                            'fasta_trns_file'       :'--trSeqFile=%(fasta_trns_file)s',
                                            'btseq_info_file'       :'--trInfoFile=%(quant_out_dir)s/%(quant_out_prefix)s.info',
                                            'comb_align_file'       :'%(comb_in_file)s',
                                            'btseq_estim_call'      :'&& %(tool_dir_btseq)sestimateExpression',
                                            'quant_out_option'      :'--outPrefix=%(quant_out_dir)s/%(quant_out_prefix)s', 
                                            'btseq_prob_in'         :'%(quant_out_dir)s/%(quant_out_prefix)s.prob',
                                            'btseq_count_call'      :'&& python %(tool_dir_btseq)sgetCounts.py',
                                            'btseq_out_count'       :'-o %(quant_out_dir)s/%(quant_out_prefix)s.count',
                                            'btseq_prob_folder'     :'-p %(quant_out_dir)s',
                                            'btseq_thet_file'       :'%(quant_out_dir)s/%(quant_out_prefix)s.thetaMeans',
                                            'btseq_seed_opt'        :'--seed=0',
                                            'btseq_verb_opt'        :'--verbose',
                                            'btseq_strand_opt'      :'--unstranded',
                                            'btseq_out_type'        :'--outType=theta',
                                            'num_processors'        :'--procN=%(num_processors)s',
                                            'bash_pre_script'       :'cd %(quant_out_dir)s', 
                                            'bash_seq_script'       :'tail -n +2 %(quant_out_dir)s/%(quant_out_prefix)s.info | cut -d" " -f2 | '
                                                                    +'paste -d" " - %(quant_out_dir)s/%(quant_out_prefix)s.count > %(quant_dir_path)s/BTSEQ_quant/%(quant_out_prefix)s.cnt'}
    # initialize RSEM
    job_config_data['RSEM_quant'] = {       'script_order'          :['tool_dir_rsem', 'options', 'rsem_tool_opt', 'comb_align_file', 'index_base_name', 'quant_out_option'], 
                                            'script_options'        :[], 
                                            'tool_dir_rsem'         :'%(tool_dir_rsem)srsem-calculate-expression',
                                            'rsem_tool_opt'         :'--seed 0 --alignments --paired-end',
                                            'comb_align_file'       :'%(comb_in_file)s',
                                            'index_base_name'       :'%(RSEM_index_dir)s/%(index_base_name)s', 
                                            'quant_out_option'      :'%(quant_out_dir)s/%(quant_out_prefix)s', 
                                            'bash_pre_script'       :'cd %(quant_out_dir)s', 
                                            'bash_seq_script'       :'cp %(quant_out_dir)s/%(quant_out_prefix)s.genes.results %(quant_dir_path)s/RSEM_quant/%(quant_out_prefix)s.cnt'}
    # initialize kallisto
    job_config_data['KLLST_quant'] = {      'script_order'          :['tool_dir_kllst', 'index_base_name', 'quant_out_option', 'kallst_seed_opt', 'kallst_out_frmt', 'num_processors', 'options', 'fastq1_file', 'fastq2_file'], 
                                            'script_options'        :[], 
                                            'tool_dir_kllst'        :'%(tool_dir_kllst)skallisto quant', 
                                            'index_base_name'       :'-i %(KLLST_index_dir)s/%(index_base_name)s.idx', 
                                            'quant_out_option'      :'-o %(quant_out_dir)s', 
                                            'kallst_seed_opt'       :'--seed=0',
                                            'kallst_out_frmt'       :'--plaintext',
                                            'num_processors'        :'--threads=%(num_processors)s',
                                            'fastq1_file'           :'%(fastq1_file)s', 
                                            'fastq2_file'           :'%(fastq2_file)s', 
                                            'bash_pre_script'       :'cd %(quant_out_dir)s', 
                                            'bash_seq_script'       :'cp %(quant_out_dir)s/abundance.tsv %(quant_dir_path)s/KLLST_quant/%(quant_out_prefix)s.cnt'}
    # TODO:
    # # initialize iReckon
    # job_config_data['IRECK_quant'] = {      'script_order'          :['tool_dir_ireck', 'comb_align_file', 'fasta_gen_file', 'index_base_name', 'fastq1_file', 'fastq2_file', 'quant_out_option', 'options'], 
                                            # 'script_options'        :['ireck_method', 'ireck_kwn_splice', 'ireck_dis_novel'], 
                                            # 'tool_dir_ireck'        :'java -Xmx15000M -jar %(tool_dir_ireck)sIReckon-1.0.8.jar', 
                                            # 'comb_align_file'       :'%(comb_in_file)s', 
                                            # 'fasta_gen_file'        :'%(fasta_gen_file)s', 
                                            # 'index_base_name'       :'%(IRECK_index_dir)s/%(index_base_name)s.tbx',
                                            # 'fastq1_file'           :'-1 %(fastq1_file)s', 
                                            # 'fastq2_file'           :'-2 %(fastq2_file)s',
                                            # 'quant_out_option'      :'-o %(quant_out_dir)s', 
                                            # 'ireck_method'          :'-m %(ireck_method)s', 
                                            # 'ireck_kwn_splice'      :'-ign', 
                                            # 'ireck_dis_novel'       :'-novel 0', 
                                            # 'bash_pre_script'       :'cd %(quant_out_dir)s', 
                                            # # create bam bai file -> samtools index -b in.bam out.bam.bai
                                            # 'bash_seq_script'       :'cp %(quant_out_dir)s/%(quant_out_prefix)s-gene.cnt %(quant_dir_path)s/IRECK_quant/%(quant_out_prefix)s.cnt'}
    # # initialize FlipFlop ## TODO
    # job_config_data['FLPFP_quant'] = {      'script_order'          :['tool_dir_flpfp', 'comb_align_file', 'fasta_gen_file', 'index_base_name', 'fastq1_file', 'fastq2_file', 'quant_out_option', 'options'], 
                                            # 'script_options'        :['ireck_method', 'ireck_kwn_splice', 'ireck_dis_novel'], 
                                            # 'tool_dir_flpfp'        :'java -Xmx15000M -jar %(tool_dir_ireck)sIReckon-1.0.8.jar', 
                                            # 'comb_align_file'       :'%(comb_in_file)s', 
                                            # 'fasta_gen_file'        :'%(fasta_gen_file)s', 
                                            # 'index_base_name'       :'%(IRECK_index_dir)s/%(index_base_name)s.tbx',
                                            # 'fastq1_file'           :'-1 %(fastq1_file)s', 
                                            # 'fastq2_file'           :'-2 %(fastq2_file)s',
                                            # 'quant_out_option'      :'-o %(quant_out_dir)s', 
                                            # 'ireck_method'          :'-m %(ireck_method)s', 
                                            # 'ireck_kwn_splice'      :'-ign', 
                                            # 'ireck_dis_novel'       :'-novel 0', 
                                            # 'bash_pre_script'       :'cd %(quant_out_dir)s', 
                                            # 'bash_seq_script'       :'cp %(quant_out_dir)s/%(quant_out_prefix)s-gene.cnt %(quant_dir_path)s/IRECK_quant/%(quant_out_prefix)s.cnt'}
    # # initialize StringTie ## TODO
    # job_config_data['STRTI_quant'] = {      'script_order'          :['tool_dir_strti', 'comb_align_file', 'quant_out_option', 'gtf_index_file', 'strti_gene_abnd', 'strti_cof_ref', 'options'], 
                                            # 'script_options'        :['strti_dis_novel'], 
                                            # 'tool_dir_strti'        :'%(tool_dir_strti)sstringtie', 
                                            # 'comb_align_file'       :'%(comb_in_file)s', 
                                            # 'quant_out_option'      :'-o %(quant_out_dir)s/%(quant_out_prefix)s_ref.gtf', 
                                            # 'gtf_index_file'        :'-G %(gtf_index_file)s', 
                                            # 'strti_gene_abnd'       :'-A %(quant_out_dir)s/%(quant_out_prefix)s_abn.tab',
                                            # 'strti_cof_ref'         :'-C %(quant_out_dir)s/%(quant_out_prefix)s_cov.gtf',
                                            # 'strti_dis_novel'       :'-e', 
                                            # 'bash_pre_script'       :'cd %(quant_out_dir)s', 
                                            # 'bash_seq_script'       :'cp %(quant_out_dir)s/%(quant_out_prefix)s-gene.cnt %(quant_dir_path)s/STRTI_quant/%(quant_out_prefix)s.cnt'}
    # # initialize eXpress ## TODO
    # job_config_data['XPRSS_quant'] = {      'script_order'          :['tool_dir_xprss', 'quant_out_option', 'xprss_no_update', 'options', 'index_base_name', 'comb_align_file'], 
                                            # 'script_options'        :['xprss_out_prob', 'xprss_out_samp'], 
                                            # 'tool_dir_xprss'        :'%(tool_dir_xprss)sexpress', 
                                            # 'index_base_name'       :'%(XPRSS_index_dir)s/%(index_base_name)s.fa', 
                                            # 'comb_align_file'       :'%(comb_in_file)s', 
                                            # 'quant_out_option'      :'-o %(quant_out_dir)s/', 
                                            # 'xprss_out_prob'        :'--output-align-prob', 
                                            # 'xprss_out_samp'        :'--output-align-samp', 
                                            # # strand options
                                            # 'xprss_no_update'       :'--no-update-check', 
                                            # 'bash_pre_script'       :'cd %(quant_out_dir)s', 
                                            # 'bash_seq_script'       :'cp %(quant_out_dir)s/%(quant_out_prefix)s-gene.cnt %(quant_dir_path)s/XPRSS_quant/%(quant_out_prefix)s.cnt'}
    
    
    
    # print debug
    if debug: printDebug(debug_list, format)
    
    # return job config
    return job_config_data




###############################################################################
## get bash script
###############################################################################

def makeBashScript(bash_config):
    
    # add bash parameters for better script overview
    bash_config['bash_param_pro'] = '%(num_processors)s'
    bash_config['bash_param_qsub'] = '%(qsub_dir_path)s'
    bash_config['bash_param_qout'] = '%(qs_out_name)s'
    bash_config['bash_param_qmem'] = '%(qs_mem_free)s'
    bash_config['bash_param_qtime'] = '%(qs_run_time)s'
    
    # uncomment to overwrite prequel, main or sequel script:
    #bash_config['bash_pre_script'] = ''
    #bash_config['bash_main_script'] = ''
    #bash_config['bash_seq_script'] = ''
    
    # define the bash script for the qsub command
    bash_script= """\
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/sh
#$ -v CLASSPATH,LD_LIBRARY_PATH,PATH,PYTHONPATH
#$ -pe smp %(bash_param_pro)s
#$ -o %(bash_param_qout)s
#$ -l m_mem_free=%(bash_param_qmem)sG
#$ -l h_rt=%(bash_param_qtime)s

############################# bash prequel script #############################
echo \"%(bash_pre_script)s\"
%(bash_pre_script)s
if [ $? -ne 0 ]; then echo "Problem with prequel script ... exiting"; echo \"%(bash_pre_script)s\" >> %(bash_param_qsub)s/qsub.error; exit 1; fi
date; echo "Prequel done."; date

############################## bash main script ###############################
echo \"%(bash_main_script)s\"
%(bash_main_script)s
if [ $? -ne 0 ]; then echo "Problem with main script ... exiting"; echo \"%(bash_main_script)s\" >> %(bash_param_qsub)s/qsub.error; exit 1; fi
date; echo "Main done."; date

############################# bash sequel script ##############################
echo \"%(bash_seq_script)s\"
%(bash_seq_script)s
if [ $? -ne 0 ]; then echo "Problem with sequel script ... exiting"; echo \"%(bash_seq_script)s\" >> %(bash_param_qsub)s/qsub.error; exit 1; fi
date; echo "Sequel done."; date
"""
    
    # make bash config data
    BASH = bash_script % bash_config
    
    # return the bash script
    return BASH



