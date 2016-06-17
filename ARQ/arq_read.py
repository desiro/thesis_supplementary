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





###############################################################################
## tool manuals
###############################################################################
    
    # htseq-count
    # --------------------------------------------------------------------------------------------------
    # command:                              -> python /home/desirda1/tools/HTSeq-0.6.1/scripts/htseq-count [options] <alignment_file> <gff_file>
    # options:  -h, --help                  -> help message
    #           -f FORMAT                   -> sam/bam (default: sam)
    #           -r ORDER                    -> sorting order: pos/name (default: name) (ignored for single-end data)
    #           -s STRANDED                 -> yes/no/reverse (default: yes) (half reads lost if not strand-specific data)
    #           -a MINAQUAL                 -> skip reads with alignment quality lower (default: 10)
    #           -t TYPE                     -> feature type (3rd column in GFF) (all others ignored) (default: exon)
    #           -i IDATTR                   -> feature ID (GFF attribute) (default: gene_id)
    #           -m MODE                     -> handle reads overlapping more than one feature (default: union)
    #           -o SAMOUT                   -> SAM alignment records into output SAM file (annotating each line with its feature assignment (optional field with tag XF)
    #           -q, --quiet                 -> suppress progress report
    # modes:                                -> for each position i in the read, set S(i) defined as set of all features overlapping position i
    #   union                               -> union of all sets S(i)
    #   intersection-strict                 -> intersection of all the sets S(i)
    #   intersection-nonempty               -> intersection of all non-empty sets S(i)
    # input:    <alignment_file>            -> aligned reads (sam/bam) (use of CIGAR field info)
    #           <gff_file>                  -> features (GFF/GTF)
    # output:   <stdout>                    -> table with counts for each feature and special counters (count reads that were not counted)
    # note:     special counters            -> __no_feature             : reads(/-pairs) not assigned (set S empty)
    #                                       -> __ambiguous              : reads(/-pairs) assigned to more than one (set S more than one element)
    #                                       -> __too_low_aQual          : reads(/-pairs) skipped due to -a option
    #                                       -> __not_aligned            : reads(/-pairs) in SAM file without alignment
    #                                       -> __alignment_not_unique   : reads(/-pairs) more than one reported alignment (NH SAM field tag)
    #                                                             (multi align reads counted multiple times, unless filtered out by -a option)
    # manual:                               -> count how many reads map to each feature
    #                                       -> a feature is an interval on a chromosome or a union of such intervals
    #                                       -> RNA-Seq: features are typically genes, each gene considered as union of all its exons
    #                                       -> (for alternative splicing: may also consider each exon as a feature)
    # counting:                             -> if S precisely one feature: read (/read pair) counted for this feature
    #                                       -> if S more than one feature: read (/read pair) counted as ambiguous (not counted for any features)
    #                                       -> if S is empty:              read (/read pair) counted as no_feature.
    
    
    # feature-count
    # --------------------------------------------------------------------------------------------------
    # command:                              -> /home/desirda1/tools/subread-1.5.0-Linux-x86_64/bin/featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ... 
    # options:  -a <string>                 -> GTF file
    #           -o <string>                 -> output file including read counts (also generated: <string>.summary)
    #           input_files                 -> list of input files (BAM/SAM format)
    #           -A <string>                 -> comma delimited file with chromosome alias names to match chromosome names (if different)
    #           -F <string>                 -> GTF/SAF/GTF (default: GTF)
    #           -t <string>                 -> specify feature type in GTF annotation (default: exon) (all others ignored)
    #           -g <string>                 -> specify attribute type in GTF annotation (default: gene_id) (meta-features used for read counting)
    #           -f                          -> read counting at feature level (counting reads for exons rather than genes)
    #           -O                          -> assign reads to all their overlapping meta-features (/features if -f specified)
    #           -s <int>                    -> strand-specific read counting (default: 0) (0: unstranded, 1: stranded, 2: reverse)
    #           -M                          -> multi-mapping reads also counted (NH tag in BAM/SAM)
    #           -R                          -> output detailed assignment result for each read (names of reads and meta-features/features reads were assigned to)
    #           --largestOverlap            -> assign reads to meta-feature/feature that has largest number of overlapping bases
    #           --minOverlap <int>          -> minimum number of overlapping bases requried between a read and a meta-feature/feature that the read is assigned to (default: 1)
    #           --read2pos <5:3>            -> reduce reads to their 5' or 3' most base (read counting performed based on single base the read is reduced to)
    #           --fraction                  -> use 1/n count, instead of 1 count, for each reported alignment of a multi-mapping read in read counting (n: total number alignments reported for the multi-mapping read) (used together with -M option)
    #           --primary                   -> count primary alignments only (identified: 0x100 bit in SAM/BAM FLAG field)
    #           --ignoreDup                 -> ignore duplicate reads in read counting (identified: Ox400 bit in BAM/SAM FLAG field) (read pair ignored if one of the reads is duplicate read for paired end data)
    #           --countSplitAlignmentsOnly  -> count split alignments only (alignments with CIGAR string containing 'N') (e.g. exon-spanning reads in RNA-seq data)
    #           -p                          -> count fragments (read pairs)
    #           -P                          -> check validity of paired-end distance when counting read pairs (-d and -D to set thresholds)
    #           -B                          -> only count read pairs with both ends successfully aligned
    #           -S <ff:fr:rf>               -> specify orientation of two reads from the same pair (default: 'fr' (=forward/reverse))
    #           -C                          -> don't count read pairs mapping to different chromosomes or mapping to same chromosome but on different strands
    #           --donotsort                 -> do not sort reads in BAM/SAM input
    # input:    <input_file1,...>           -> SAM/BAM files
    #           <annotation_file>           -> GTF annotation file including chromosomal coordinates of features
    # output:   <output_file>               -> numbers of reads assigned to features (or meta-features)
    #           <">.summary                 -> stat info: overall summarization results (number successfully assigned reads and number reads failed to be assigned)
    # manual:                               -> counts mapped reads for genomic features (genes, exons, promoter, gene bodies, genomic bins and chromosomal locations)
    # features:                             -> each entry in the provided annotation file taken as a feature (e.g. an exon)
    # meta-feature:                         -> aggregation of a set of features (e.g. a gene)
    # notes:                                -> gene_id in the GTF format to group features into meta-features
    # count:                                -> at feature level or at meta-feature level
    #                                       -> summarizing reads at meta-feature level: read counts for features in the same meta-feature will be added up to read count for meta-feature.
    # overlap:                              -> read overlap a feature if at least one read base is found to overlap the feature
    #                                       -> paired-end: a fragment overlaps a feature if any of the two reads from that fragment is found to overlap the feature
    #                                       -> (default: not count reads overlapping with more than one feature (or more than one meta-feature when meta level)
    #                                       -> -O option to count such reads (will be assigned to all their overlapping features or meta-features)
    # notes:                                -> meta-feature level counting: reads that overlap multiple features of same meta-feature: counted exactly once 
    #                                          for that meta-feature, provided there is no overlap with any other meta-feature
    #                                          (e.g.: exon-spanning read counted only once for the corresponding gene even if it overlaps with more than one exon)
    
    
    # EQP-QM (Exon quantification pipeline - quantification module)
    # --------------------------------------------------------------------------------------------------
    # command:                              -> /home/desirda1/tools/EQP-QM-master/bin/eqp-setup.sh <GTF file> <data directory>
    # manual:                               -> create annotation files (once per GTF), input: GTF file
    # command:                              -> /home/desirda1/tools/EQP-QM-master/bin/eqp-quantify.sh <options> -d <setup dir> <output dir> <SAM/BAM file>
    # options:  -d STRING                   -> STRING of directory that contains the auxilliary files (from eqp_setup.sh)
    #           -g                          -> compute gene counts [computed by default, without -g,-e,-j] 
    #           -e                          -> compute exon counts [computed by default, without -g,-e,-j]
    #           -j                          -> compute junction counts [computed by default, without -g,-e,-j]
    #           -E INT                      -> Minimal overlap of a read with an exon (default: 5)
    #           -J INT                      -> Minimal overlap of a read with both exons on a junction (default: 8)
    #           -W FLOAT                    -> Minimal weight of a read; reads with a lower weight are disregarded (default: 0.01)
    #           -s DIRECTION                -> process reads as strand-specific in direction (forward for orientation fr or backward for orientation rf)
    #           -o STRING                   -> created in the current working directory (for intermediate files)
    #                                       -> count files generated in current working directory with prefix STRING
    #           --nosort                    -> alignment file already sorted by names - do not sort it
    #           --unambig                   -> count only reads that can be assigned unambiguously to a single gene or exon when creating the gene or exon counts
    #           --unweighted                -> do not use read weights in the generation of counts
    #           -w STRING                   -> use STRING as weight file instead of computing it from the genomic alignments
    # input:    <SAM/BAM file>              -> SAM/BAM file
    #           <setup dir>                 -> setup file from eqp-setup.sh
    # output:                               -> <SAM/BAM file base>-gene.cnt, <SAM/BAM file base>-gene.cnt, <SAM/BAM file base>-junction.cnt
    #           <output dir>                -> directory used for the count files SAM/BAM file: file containing the alignments of the reads against the genome with the aligner
    #           -o STRING                   -> count files generated in current working directory with prefix STRING
    # notes:    GTF file                    -> only entries with feature type exon which contain gene_id field are used by eqp-setup.sh
    # runtime:                              -> expect ~0.5-1h for ~10M paired-end reads. EQP-QM needs at least 10GB of main memory
    
    
    # aligner -> cufflinks -> cuffmerge -> cuffquant -> cuffnorm -> = normalized expression and count tables
    # --------------------------------------------------------------------------------------------------
    # 
    # cufflinks v2.2.1
    # -----------------------------
    # Usage:                                -> cufflinks [options] <hits.sam>
    # options:  -o/--output-dir             -> output dir                                            [ default:     ./ ]    +
    #           -p/--num-threads            -> number of threads                                     [ default:      1 ]    +
    #           --seed                      -> random number seed                                    [ default:      0 ]    +
    #           -G/--GTF                    -> quantitate against reference transcript annotations                          + -> only search for genes from GTF file
    #           -g/--GTF-guide              -> use reference transcript annotation to guide assembly                        ? -> also search for novel genes
    #           -M/--mask-file              -> ignore all alignment within transcripts in this file                         ?
    #           -b/--frag-bias-correct      -> use bias correction - reference fasta required        [ default:   NULL ]    ?
    #           -u/--multi-read-correct     -> use 'rescue method' for multi-reads (more accurate)   [ default:  FALSE ]    ?
    #           --library-type              -> library prep used for input reads                     [ default:  below ]    +   -> ff-unstranded alignments?
    #                                          ff-firststrand   ff-secondstrand  ff-unstranded  fr-firststrand   fr-secondstrand  fr-unstranded (default)  transfrags
    #           --library-norm-method       -> Method used to normalize library sizes                [ default:  below ]    +   -> classic-fpkmqe
    #                                          classic-fpkmqe
    # Advanced Abundance Estimation Options:
    #           -m/--frag-len-mean          -> average fragment length (unpaired reads only)         [ default:    200 ]    -
    #           -s/--frag-len-std-dev       -> fragment length std deviation (unpaired reads only)   [ default:     80 ]    -
    #           --max-mle-iterations        -> maximum iterations allowed for MLE calculation        [ default:   5000 ]    ?
    #           --compatible-hits-norm      -> count hits compatible with reference RNAs only        [ default:  FALSE ]    ?
    #           --total-hits-norm           -> count all hits for normalization                      [ default:  TRUE  ]    ?
    #           --num-frag-count-draws      -> Number of fragment generation samples                 [ default:    100 ]    ?
    #           --num-frag-assign-draws     -> Number of fragment assignment samples per generation  [ default:     50 ]    ?
    #           --max-frag-multihits        -> Maximum number of alignments allowed per fragment     [ default: unlim  ]    ?
    #           --no-effective-length-correction    -> No effective length correction                  [ default:  FALSE ]    ?
    #           --no-length-correction      -> No length correction                                  [ default:  FALSE ]    ?
    #           -N/--upper-quartile-norm    -> Deprecated, use --library-norm-method                 [    DEPRECATED   ]    ?
    #           --raw-mapped-norm           -> Deprecated, use --library-norm-method                 [    DEPRECATED   ]    ?
    # Advanced Assembly Options:
    #           -L/--label                  -> assembled transcripts have this ID prefix             [ default:   CUFF ]    +
    #           -F/--min-isoform-fraction   -> suppress transcripts below this abundance level       [ default:   0.10 ]    ?
    #           -j/--pre-mrna-fraction      -> suppress intra-intronic transcripts below this level  [ default:   0.15 ]    ?
    #           -I/--max-intron-length      -> ignore alignments with gaps longer than this          [ default: 300000 ]    ?
    #           -a/--junc-alpha             -> alpha for junction binomial test filter               [ default:  0.001 ]    ?
    #           -A/--small-anchor-fraction  -> percent read overhang taken as 'suspiciously small'   [ default:   0.09 ]    ?
    #           --min-frags-per-transfrag   -> minimum number of fragments needed for new transfrags [ default:     10 ]    ?
    #           --overhang-tolerance        -> number of terminal exon bp to tolerate in introns     [ default:      8 ]    ?
    #           --max-bundle-length         -> maximum genomic length allowed for a given bundle     [ default:3500000 ]    ?
    #           --max-bundle-frags          -> maximum fragments allowed in a bundle before skipping [ default: 500000 ]    ?
    #           --min-intron-length         -> minimum intron size allowed in genome                 [ default:     50 ]    +   -> alignmets option
    #           --trim-3-avgcov-thresh      -> minimum avg coverage required to attempt 3' trimming  [ default:     10 ]    ?
    #           --trim-3-dropoff-frac       -> fraction of avg coverage below which to trim 3' end   [ default:    0.1 ]    ?
    #           --max-multiread-fraction    -> maximum fraction of allowed multireads per transcript [ default:   0.75 ]    ?   -> max multi hits?
    #           --overlap-radius            -> maximum gap size to fill between transfrags (in bp)   [ default:     50 ]    ?   -> max gap size align?
    # Advanced Reference Annotation Guided Assembly Options:
    #           --no-faux-reads             -> disable tiling by faux reads                          [ default:  FALSE ]    ?
    #           --3-overhang-tolerance      -> overhang allowed on 3' end when merging with reference[ default:    600 ]    ?
    #           --intron-overhang-tolerance -> overhang allowed inside reference intron when merging [ default:     30 ]    ?
    # Advanced Program Behavior Options:
    #           -v/--verbose                -> log-friendly verbose processing (no progress bar)     [ default:  FALSE ]    +   -> TRUE -> output all 
    #           -q/--quiet                  -> log-friendly quiet processing (no progress bar)       [ default:  FALSE ]    +   -> TRUE (only verbose or quiet)
    #           --no-update-check           -> do not contact server to check for update availability[ default:  FALSE ]    +   -> TRUE
    # 
    # cuffmerge v2.2.1
    # -----------------------------
    # cuffmerge takes two or more Cufflinks GTF files and merges them into a
    # single unified transcript catalog.  Optionally, you can provide the script
    # with a reference GTF, and the script will use it to attach gene names and other
    # metadata to the merged catalog.
    # Usage:                                -> cuffmerge [Options] <assembly_GTF_list.txt>
    # Options:  -h/--help                   -> Prints the help message and exits
    #           -o <output_dir>             -> Directory where merged assembly will be written       [ default: ./merged_asm  ]
    #           -g/--ref-gtf                -> An optional "reference" annotation GTF.
    #           -s/--ref-sequence <seq_dir>/<seq_fasta> -> Genomic DNA sequences for the reference.
    #           --min-isoform-fraction <0-1.0>          -> Discard isoforms with abundance below this [ default:           0.05 ]
    #           -p/--num-threads <int>      -> Use this many threads to merge assemblies.            [ default:             1  ]
    #           --keep-tmp                  -> Keep all intermediate files during merge
    #
    # cuffquant v2.2.1
    # -----------------------------
    # Usage:   cuffdiff [options] <transcripts.gtf> <sample1_hits.sam> <sample2_hits.sam> [... sampleN_hits.sam]
    # Supply replicate SAMs as comma separated lists for each condition: sample1_rep1.sam,sample1_rep2.sam,...sample1_repM.sam
    # Options:  -o/--output-dir             -> write all output files to this directory              [ default:     ./ ]
    #           -M/--mask-file              -> ignore all alignment within transcripts in this file  [ default:   NULL ]
    #           -b/--frag-bias-correct      -> use bias correction - reference fasta required        [ default:   NULL ]
    #           -u/--multi-read-correct     -> use 'rescue method' for multi-reads                   [ default:  FALSE ]
    #           -p/--num-threads            -> number of threads used during quantification          [ default:      1 ]
    #           --library-type              -> Library prep used for input reads                     [ default:  below ]
    # Advanced Options
    #           -m/--frag-len-mean          -> average fragment length (unpaired reads only)         [ default:    200 ]
    #           -s/--frag-len-std-dev       -> fragment length std deviation (unpaired reads only)   [ default:     80 ]
    #           -c/--min-alignment-count    -> minimum number of alignments in a locus for testing   [ default:   10 ]
    #           --max-mle-iterations        -> maximum iterations allowed for MLE calculation        [ default:   5000 ]
    #           -v/--verbose                -> log-friendly verbose processing (no progress bar)     [ default:  FALSE ]
    #           -q/--quiet                  -> log-friendly quiet processing (no progress bar)       [ default:  FALSE ]
    #           --seed                      -> value of random number generator seed                 [ default:      0 ]
    #           --no-update-check           -> do not contact server to check for update availability[ default:  FALSE ]
    #           --max-bundle-frags          -> maximum fragments allowed in a bundle before skipping [ default: 500000 ]
    #           --max-frag-multihits        -> Maximum number of alignments allowed per fragment     [ default: unlim  ]
    #           --no-effective-length-correction    -> No effective length correction                  [ default:  FALSE ]
    #           --no-length-correction      -> No length correction                                  [ default:  FALSE ]
    #           Debugging use only:
    #           --read-skip-fraction        -> Skip a random subset of reads this size               [ default:    0.0 ]
    #           --no-read-pairs             -> Break all read pairs                                  [ default:  FALSE ]
    #           --trim-read-length          -> Trim reads to be this long (keep 5' end)              [ default:   none ]
    #           --no-scv-correction         -> Disable SCV correction                                [ default:  FALSE ]
    #
    # cuffnorm v2.2.1
    # -----------------------------
    # Usage:   cuffnorm [options] <transcripts.gtf> <sample1_expr.cxb> <sample2_expr.cxb> [... sampleN_expr.cxb]
    # Supply replicate CXB files as comma separated lists for each condition: sample1_rep1.cxb,sample1_rep2.cxb,...sample1_repM.cxb
    # Options: -o/--output-dir              -> write all output files to this directory              [ default:     ./ ]
    #           -L/--labels                 -> comma-separated list of condition labels
    #           --norm-standards-file       -> Housekeeping/spike genes to normalize libraries       [ default:   NULL ]
    #           -p/--num-threads            -> number of threads used during quantification          [ default:      1 ]
    #           --library-type              -> Library prep used for input reads                     [ default:  below ]
    #           --library-norm-method       -> Method used to normalize library sizes                [ default:  below ]
    #                                          classic-fpkm  geometric (default)  geometric  quartile
    #           --output-format             -> Format for output tables                              [ default:  below ]
    #                                          cuffdiff  simple-table (default)
    # Advanced Options:
    #           --compatible-hits-norm      -> count hits compatible with reference RNAs only        [ default:   TRUE ]
    #           --total-hits-norm           -> count all hits for normalization                      [ default:  FALSE ]
    #           -v/--verbose                -> log-friendly verbose processing (no progress bar)     [ default:  FALSE ]
    #           -q/--quiet                  -> log-friendly quiet processing (no progress bar)       [ default:  FALSE ]
    #           --seed                      -> value of random number generator seed                 [ default:      0 ]
    #           --no-update-check           -> do not contact server to check for update availability[ default:  FALSE ]
    
    
    # Flux-Capacitor v1.6.1 (Flux Library: 1.29)
    # --------------------------------------------------------------------------------------------------
    #
    # Tools:    -t <toolname>               -> sortBED      -   Sort a BED file. If no output file is specified, result is printed to standard out
    #                                          capacitor    -   The Flux Capacitor
    #                                          sortGTF      -   Sort a GTF file. If no output file is specified, result is printed to standard out
    #                                          diffexp      -   Differential Expression between two samples
    #                                          subsetter    -   Extract a random subset of lines from a file
    # 
    # Capacitor:
    # Options:
    #   --profile-output <PROFILE_OUTPUT>   -> The file for outputting profiles                                                                                 +
    #   --ignore-sam-flags                  -> Ignore SAM flags when scanning the mapping file                                                                  -
    #   --ignore-sam-pairing-information    -> Ignore SAM pairing information in the quantification                                                             -
    #   --sam-validation-stringency <SAM_V> -> Set SAMtools validation stringency for validating records. One of [STRICT|LENIENT|SILENT]                        ?
    #   -r|--sort-in-ram                    -> Sort reads in RAM memory, not on disk. This is the default behaviour. You can force sorting on disk              -
    #                                          with the --sort-on-disk option                                                                                   
    #   --count-elements <COUNT_ELEMENTS>   -> Count specified elements. Possible elements are: [SPLICE_JUNCTIONS,INTRONS]                                      ?
    #   --disable-multimap-weighting        -> Disable weighted counts for multi-maps                                                                           ?
    #   --stats-file <STATS_FILE>           -> The file to which the run characteristics are written                                                            +
    #   --tmp-dir <TMP_DIR>                 -> The temporary directory                                                                                          +
    #   --profile-file <PROFILE_FILE>       -> The profile file                                                                                                 +
    #   --min-obs <MIN_OBS>                 -> Minimum number of reads to be left in a locus by the linear solver. It can be expressed as an absolute           -
    #                                          number of reads (>1) or as fraction of the observed read count [0..1]. With MIN_OBS > 0 the lp-solver            
    #                                          will be forced to provide a deconvolution in all loci with > 0 annotation-mapped reads or pairs. The             
    #                                          derived quantifications may be variable in loci with less observations than changes that are necessary           
    #                                          to perform the deconvolution.                                                                                    
    #   -o|--output <STDOUT_FILE>           -> The file for default output                                                                                      +
    #   --coverage-file <COVERAGE_FILE>     -> Calculate coverage profile write it to the specified file                                                        ?
    #   --profile-exclude <PROFILE_EXCLUDE> -> Annotation sources (GTF field 2) that are disregarded during profiling                                           -
    #   --read-strand <READ_STRAND>         -> Information about read strandedness                                                                              ?
    #   --sam-unique-only                   -> Only use unique alignments for quantification                                                                    ?
    #   --disable-deconvolution             -> Disable the deconvolution step                                                                                   ?
    #   --minilen <MIN_ILEN>                -> Minimum length of introns of the annotation that are considered real and not indels/gaps, "introns" with         -
    #                                          a length < MIN_ILEN are removed by joining the flanking exons                                                    
    #   -i|--input <MAPPING_FILE>           -> The mapping file                                                                                                 +
    #   -q|--min-score <MIN_SCORE>          -> Minimum mapping score. Mappings with score < min_score are discarded (mapq for BAM, score for BED)               -
    #   -a|--annotation   <ANNOTATION_FILE> -> The annotation file                                                                                              +
    #   --stderr-file <STDERR_FILE>         -> The file for log messages                                                                                        +
    #   --profile-include <PROFILE_INCLUDE> -> Annotation sources (GTF field 2) that are considered during profiling                                            ?
    #   --keep-sorted <KEEP_SORTED>         -> Keeps input files (ANNOTATION_FILE, MAPPING_FILE)                                                                ?
    #   --sort-on-disk                      -> Sort reads on disk                                                                                               -
    #   --disable-file-check                -> Disable scanning of input files before the run                                                                   -
    #   --sam-primary-only                  -> Only use primary alignments for quantification                                                                   -
    #   -d|--read-descriptor <READ_DESCRIP> -> Expression how to parse the read IDs, or one of the shorthand names (ANTISENSE,STRAND_MATE,BARNA,MATE2_SENSE,    ?
    #                                          MATE1_SENSE,MATE_STRAND_CSHL,SENSE,SIMULATOR,PAIRED,SIMPLE,CASAVA18)                                             
    #   --insert-file <INSERT_FILE>         -> The file for output of inserts                                                                                   +
    #   -m|--annotation-mapping <ANNO_MAP>  -> Information from the read descriptor that will be used for annotation mapping                                    ?
    #   --profile                           -> do the profiling                                                                                                 ?
    #   -p|--parameter <file>               -> specify parameter file (PAR file)                                                                                -
    #   --printParameters                   -> Print default parameters                                                                                         -
    #
    # General Flux Options: 
    #   -t|--tool   <tool>                  -> Select a tool (default: capacitor)                                                                               +
    #   -h|--help                           -> Show help                                                                                                        -
    #   --list-tools                        -> List available tools                                                                                             -
    #   --threads <threads>                 -> Maximum number of threads to use. Default 2 (default: 2)                                                         -
    #   --log <level>                       -> Log level (NONE|INFO|ERROR|DEBUG) (default: INFO)                                                                -
    #   --force                             -> Disable interactivity. No questions will be asked                                                                +
    #   -v|--version                        -> Show version information                                                                                         -
    #
    # sortGTF:
    # Options:
    #   -t sortGTF                          -> start tool
    #   -i|--input <gtf>                    -> GTF input file
    #   -o|--output <output>                -> GTF output file. Sorts to stdout if no file is given
    #   -c|--check                          -> Check if the file is sorted before sorting
    #   --force                             -> Disable interactivity. No questions will be asked
    
    
    # BitSeq
    # --------------------------------------------------------------------------------------------------
    # gtf2bed.pl/py
    # -----------------------------
    # gtf file to bed file
    # usage:    gtf2bed.pl/py <in_file> > <out_file>
    #
    # bedtools getfasta
    # -----------------------------
    # bedtools getfasta extracts sequences from a FASTA file for each of the intervals defined in a BED/GFF/VCF file.
    # Usage and option summary
    # Usage:                                -> bedtools getfasta [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF> -fo <output FASTA>
    # Options:  -name                       -> Use the “name” column in the BED file for the FASTA headers in the output FASTA file.
    #           -tab                        -> Report extract sequences in a tab-delimited format instead of in FASTA format.
    #           -s                          -> Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented. 
    #                                          Default: strand information is ignored.
    #           -split                      -> Given BED12 input, extract and concatenate the sequences from the BED “blocks” (e.g., exons)
    #
    # BitSeq
    # -----------------------------
    # getExpression              package:BitSeq              R Documentation
    # Estimate transcript expression
    # Description:                          -> Estimate expression of transcripts. Starting from alignment and reference files function function 
    #                                          handles the entire process of expression analysis resulting in transcript expression means and
    #                                          standard deviation together with file containing all the expression samples.
    # Usage:                                -> getExpression(alignFile, trSeqFile, outPrefix=NULL, uniform=TRUE, type="RPKM", 
    #                                          log=FALSE, limitA=NULL, seed=NULL, pretend=FALSE, ... )
    # Options:  alignFile                   -> File containing read alignments.
    #           trSeqFile                   -> File containing transcript sequence in FASTA format.
    #           outPrefix                   -> Prefix for the output files. Otherwise program creates temporary files, which are only valid for current R session.
    #           uniform                     -> Use uniform read distribution.
    #           type                        -> Output type, possible values: 'theta', 'RPKM', 'counts', 'tau'.
    #           log                         -> Report mean and expression of logged expression samples.
    #           limitA                      -> Limit maximum number of alignments per read. Reads with more alignments than limit will be discarded.
    #           seed                        -> Sets the initial random seed for repeatable experiments.
    #           pretend                     -> Do not execute, only print out command line calls for the C++ version of the program.
    #           ...                         -> Other arguments are passed to 'estimateExpression', please see 'estimateExpression' for more details
    # Details:                              -> This function uses 'parseAlignment' function to compute alignment probabilities and the function 
    #                                          'estimateExpression' to produce the expression samples.
    #                                       -> In case of non-uniform read distribution, it first produces approximate estimates of expression using 
    #                                          uniform read distribution with VB inference and subsequently uses these estimates to compute read 
    #                                          distribution bias-corrected alignment probabilities, which are used in the 'estimateExpression' function 
    #                                          to produce expression estimates.
    #                                       -> The order of transcripts in the results is always the same as in the alignment file. The transcripts can 
    #                                          be identified by names stored in the 'trInfo' part of the result.
    # Value:    'list' with items:
    #           exp                         -> 'DataFrame' with transcript expression mean and standard deviation
    #           fn                          -> name of the file containing all the expression samples
    #           counts                      -> vector of estimated read counts per transcript
    #           trInfo                      -> 'DataFrame' with transcript information, contains: transcript name, possibly gene name, transcript length, 
    #                                          and adjusted transcript length
    #
    # BitSeq c++
    # -----------------------------
      
    
    # parseAlignment -o <outFileName> -s <trSeqFileName>  [OPTIONS] [alignment file]
    #   Pre-computes probabilities of (observed) reads' alignments.
    #   [alignment file] should be in either SAM or BAM format.
    #
    # Options:
    # --help                                                Show this help information.
    # --distributionFile=<distributionFileName>             Name of file to which read-distribution should be saved.
    # --excludeSingletons                                   Exclude single mate alignments for paired-end reads. (default: Off)
    # -e <expFileName> ,   --expressionFile=<expFileName>   Transcript relative expression estimates --- for better non-uniform read distribution estimation.
    # --failed=<failed>                                     File name where to save names of reads that failed to align.
    # -f <format> ,   --format=<format>                     Input format: either SAM, BAM.
    # --lenMu=<lenMu>                                       Set mean of log fragment length distribution. (l_frag ~ LogNormal(mu,sigma^2))
    # --lenSigma=<lenSigma>                                 Set sigma^2 (or variance) of log fragment length distribution. (l_frag ~ LogNormal(mu,sigma^2))
    # --mateNamesDiffer                                     Mates from paired-end reads have different names. (default: Off)
    # -l <maxAlignments> ,   --limitA=<maxAlignments>       Limit maximum number of alignments per read. (Reads with more alignments are skipped.)
    # --noiseMismatches=<numNoiseMismatches>                Number of mismatches to be considered as noise. (default: 6)
    # -o <outFileName> ,   --outFile=<outFileName>          Name of the output file.
    # -P <procN> ,   --procN=<procN>                        Maximum number of threads to be used. This provides speedup mostly when using non-uniform read distribution model (i.e. no --uniform flag). (default: 4)
    # -N <readsN> ,   --readsN=<readsN>                     Total number of reads. This is not necessary if [SB]AM contains also reads with no valid alignments.
    # --show1warning                                        Show first alignments that are considered wrong (TID unknown, TID mismatch, wrong strand). (default: Off)
    # -t <trInfoFileName> ,   --trInfoFile=<trInfoFileName> File to save transcript information extracted from [BS]AM file and reference.
    # -s <trSeqFileName> ,   --trSeqFile=<trSeqFileName>    Transcript sequence in FASTA format --- for non-uniform read distribution estimation.
    # --trSeqHeader=<trSeqHeader>                           Transcript sequence header format enables gene name extraction (standard/gencode). (default: standard)
    # --uniform                                             Use uniform read distribution. (default: Off)
    # --unstranded                                          Paired read are not strand specific. (default: Off)
    # -v ,   --verbose                                      Verbose output. (default: Off)
    # -V ,   --veryVerbose                                  Very verbose output. (default: Off)
    #
    # -----------------------------
    # estimateExpression -o <outFilePrefix>  [OPTIONS] [prob file]
    #   Estimates expression given precomputed probabilities of (observed) reads' alignments.
    #   Uses MCMC sampling algorithm to produce relative abundance or RPKM.
    #
    # Options:
    # --help                                                Show this help information.
    # --MCMC_burnIn=<MCMC_burnIn>                           Length of sampler's burn in period. (default: 1000)
    # --MCMC_chainsN=<MCMC_chainsN>                         Number of parallel chains used. At least two chains will be used. (default: 4)
    # --MCMC_dirAlpha=<MCMC_dirAlpha>                       Alpha parameter for the Dirichlet distribution. (default: 1)
    # --MCMC_samplesDOmax                                   Produce maximum number of samples (samplesNmax) in second iteration and quit. (default: Off)
    # --MCMC_samplesN=<MCMC_samplesN>                       Initial number of samples produced. Doubles after every iteration. (default: 1000)
    # --MCMC_samplesNmax=<MCMC_samplesNmax>                 Maximum number of samples produced in one iteration. After producing samplesNmax samples sampler finishes. (default: 50000)
    # --MCMC_samplesSave=<MCMC_samplesSave>                 Number of samples recorder in total. (default: 1000)
    # --MCMC_scaleReduction=<MCMC_scaleReduction>           Target scale reduction, sampler finishes after this value is met. (default: 1.2)
    # -G ,   --gibbs                                        Use Gibbs sampling instead of collapsed Gibbs sampling. (default: Off)
    # -o <outFilePrefix> ,   --outPrefix=<outFilePrefix>    Prefix for the output files.
    # -O <outputType> ,   --outType=<outputType>            Output type (theta, RPKM, counts, tau). (default: theta)
    # -p <parFileName> ,   --parFile=<parFileName>          File containing parameters for the sampler, which can be otherwise specified by --MCMC* options. As the file is checked after every MCMC iteration, the parameters can be adjusted while running.
    # -P <procN> ,   --procN=<procN>                        Limit the maximum number of threads to be used. (Default is the number of MCMC chains.)
    # --scaleReduction                                      Use scale reduction as stopping criterion, instead of computing effective sample size. (default: Off)
    # -s <seed> ,   --seed=<seed>                           Random initialization seed.
    # --thetaActFile=<thetaActFileName>                     File for logging noise parameter theta^{act}.
    # -t <trInfoFileName> ,   --trInfoFile=<trInfoFileName> File containing transcript information. (Necessary for RPKM)
    # -v ,   --verbose                                      Verbose output. (default: Off)
    #
    # -----------------------------
    # getVariance -o <outFileName>  [OPTIONS] [sampleFiles]
    #   Estimates variance of MCMC samples from 1 or multiple replicates.
    #   [sample Files] should contain transposed MCMC samples from replicates.
    #
    # Options:
    # --help                                                Show this help information.
    # -l ,   --log                                          Use logged values. (default: Off)
    # --norm=<normalization>                                Normalization constants for each input file provided as comma separated list of doubles (e.g. 1.0017,1.0,0.9999 ).
    # -o <outFileName> ,   --outFile=<outFileName>          Name of the output file.
    # -t <type> ,   --type=<type>                           Type of variance, possible values: [sample,sqDif] for sample variance or squared difference. (default: sample)
    # -v ,   --verbose                                      Verbose output. (default: Off)
    
    
    
    # RSEM
    # --------------------------------------------------------------------------------------------------
    # To take advantage of RSEM’s built-in support for the Bowtie/Bowtie 2/STAR alignment program, you must have Bowtie/Bowtie 2/STAR installed.
    # 
    # I. Preparing Reference Sequences      -> rsem-prepare-reference --gtf Homo_sapiens.GRCh38.83.gtf --bowtie Homo_sapiens.GRCh38.dna.primary_assembly.fa ref/human_ensembl
    #                                       -> RSEM can extract reference transcripts from a genome if you provide it with gene annotations in a GTF/GFF3 file. 
    #                                          Alternatively, you can provide RSEM with transcript sequences directly.
    #                                          Please note that GTF files generated from the UCSC Table Browser do not contain isoform-gene relationship information. 
    #                                          However, if you use the UCSC Genes annotation track, this information can be recovered by downloading the knownIsoforms.txt 
    #                                          file for the appropriate genome.
    # rsem-prepare-reference --help         -> Build RSEM references using RefSeq, Ensembl, or GENCODE annotations
    #                                          RefSeq and Ensembl are two frequently used annotations. For human and mouse, GENCODE annotaions are also available. 
    #                                          In this section, we show how to build RSEM references using these annotations. Note that it is important to pair the 
    #                                          genome with the annotation file for each annotation source. In addition, we recommend users to use the primary assemblies 
    #                                          of genomes. Without loss of generality, we use human genome as an example and in addition build Bowtie indices.
    # For Ensembl:                          -> genome and annotation files can be found at Ensembl FTP.
    #                                          ftp://ftp.ensembl.org/pub/release-83/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    #                                          ftp://ftp.ensembl.org/pub/release-83/gtf/homo_sapiens/Homo_sapiens.GRCh38.83.gtf.gz
    # 
    # II. Calculating Expression Values     -> 
    #                                       -> To calculate expression values, you should run the rsem-calculate-expression program. Run
    # rsem-calculate-expression --help      
    # Expression values from single-end     -> For single-end models, users have the option of providing a fragment length distribution via the --fragment-length-mean 
    #                                          and --fragment-length-sd options. The specification of an accurate fragment length distribution is important for the accuracy 
    #                                          of expression level estimates from single-end data. If the fragment length mean and sd are not provided, RSEM will not take a 
    #                                          fragment length distribution into consideration.
    # Using an alternative aligner          -> By default, RSEM automates the alignment of reads to reference transcripts using the Bowtie aligner. Turn on --bowtie2 for 
    #                                          rsem-prepare-reference and rsem-calculate-expression will allow RSEM to use the Bowtie 2 alignment program instead. 
    #                                          Please note that indel alignments, local alignments and discordant alignments are disallowed when RSEM uses Bowtie 2 since 
    #                                          RSEM currently cannot handle them. See the description of --bowtie2 option in rsem-calculate-expression for more details. 
    #                                          Similarly, turn on --star will allow RSEM to use the STAR aligner. To use an alternative alignment program, align the input 
    #                                          reads against the file reference_name.idx.fa generated by rsem-prepare-reference, and format the alignment output in SAM/BAM/CRAM 
    #                                          format. Then, instead of providing reads to rsem-calculate-expression, specify the --alignments option and provide the SAM/BAM/CRAM 
    #                                          file as an argument.
    #                                       -> RSEM requires the alignments of a read to be adjacent. For paired-end reads, RSEM also requires the two mates of any alignment 
    #                                          be adjacent. To check if your SAM/BAM/CRAM file satisfy the requirements, run
    #                                       -> rsem-sam-validator <input.sam/input.bam/input.cram>
    #                                          If your file does not satisfy the requirements, you can use convert-sam-for-rsem to convert it into a BAM file which RSEM can process. Run
    #                                           convert-sam-for-rsem --help
    #                                           to get usage information or visit the convert-sam-for-rsem documentation page.
    #                                           Note that RSEM does ** not ** support gapped alignments. So make sure that your aligner does not produce alignments with intersions/deletions. 
    #                                           In addition, you should make sure that you use reference_name.idx.fa, which is generated by RSEM, to build your aligner’s indices.
    #
    # convert-sam-for-rsem
    # -----------------------------
    # Make a RSEM compatible BAM file.
    # Usage:                                -> convert-sam-for-rsem [options] <input.sam/input.bam/input.cram> output_file_name
    # input.sam/input.bam/input.cram        -> The SAM/BAM/CRAM file generated by user's aligner. We require this file contains the header section.
    # output_file_name                      -> The output name for the converted file. 'convert-sam-for-rsem' will output a BAM with the name 'output_file_name.bam'.
    # Options:  -p/--num-threads <int>      -> Set the number of threads to be used for converting. (Default: 1)
    #           --memory-per-thread <string>-> Set the maximum allowable memory per thread. <string> represents the memory and accepts suffices 'K/M/G'. (Default: 1G)
    #           -h/--help                   -> Show help information.
    # Description:                          -> This program converts the SAM/BAM/CRAM file generated by user's aligner into a BAM file which RSEM can process. However, 
    #                                          users should make sure their aligners use 'reference_name.idx.fa' generated by 'rsem-prepare-reference' as their references 
    #                                          and output header sections. After the conversion, this program will call 'rsem-sam-validator' to validate the resulting BAM file.
    #                                          Note: You do not need to run this script if 'rsem-sam-validator' reports that your SAM/BAM/CRAM file is valid.
    #
    # rsem-prepare-reference
    # -----------------------------
    # Prepare transcript references for RSEM and optionally build BOWTIE/BOWTIE2/STAR indices.
    # Usage:                                -> rsem-prepare-reference [options] reference_fasta_file(s) reference_name
    # reference_fasta_file(s)               -> Either a comma-separated list of Multi-FASTA formatted files OR a directory name. If a directory name is specified, RSEM will 
    #                                          read all files with suffix ".fa" or ".fasta" in this directory. The files should contain either the sequences of transcripts 
    #                                          or an entire genome, depending on whether the '--gtf' option is used.
    # reference name                        -> The name of the reference used. RSEM will generate several reference-related files that are prefixed by this name. This name 
    #                                          can contain path information (e.g. '/ref/mm9').
    # Options:  --gtf <file>                -> If this option is on, RSEM assumes that 'reference_fasta_file(s)' contains the sequence of a genome, and will extract transcript 
    #                                          reference sequences using the gene annotations specified in <file>, which should be in GTF format.
    #                                          If this and '--gff3' options are off, RSEM will assume 'reference_fasta_file(s)' contains the reference transcripts. In this case, 
    #                                          RSEM assumes that name of each sequence in the Multi-FASTA files is its transcript_id. (Default: off)
    #           --gff3 <file>               -> The annotation file is in GFF3 format instead of GTF format. RSEM will first convert it to GTF format with the file name 
    #                                          'reference_name.gtf'. Please make sure that 'reference_name.gtf' does not exist. (Default: off)
    #           --gff3-RNA-patterns         -> <pattern> is a comma-separated list of transcript categories, e.g. "mRNA,rRNA". Only transcripts that match the <pattern> will 
    #           <pattern>                      be extracted. (Default: "mRNA")
    #           --trusted-sources <sources> -> <sources> is a comma-separated list of trusted sources, e.g. "ENSEMBL,HAVANA". Only transcripts coming from these sources will be 
    #                                          extracted. If this option is off, all sources are accepted. (Default: off)
    #           --transcript-to-gene-map    -> Use information from <file> to map   from transcript (isoform) ids to gene ids.  Each line of <file> should be of the form
    #           <file>                         "gene_id transcript_id" with the two fields separated by a tab character.
    #                                          If you are using a GTF file for the "UCSC Genes" gene set from the UCSC Genome Browser, then the "knownIsoforms.txt" file 
    #                                          (obtained from the "Downloads" section of the UCSC Genome Browser site) is of this format. If this option is off, then the 
    #                                          mapping of isoforms to genes depends on whether the '--gtf' option is specified.  If '--gtf' is specified, then RSEM uses the 
    #                                          "gene_id" and "transcript_id" attributes in the GTF file.  Otherwise, RSEM assumes that each sequence in the reference sequence 
    #                                          files is a separate gene. (Default: off)
    #           --allele-to-gene-map <file> -> Use information from <file> to provide gene_id and transcript_id information for each allele-specific transcript. Each line of 
    #                                          <file> should be of the form "gene_id transcript_id allele_id" with the fields separated by a tab character. This option is 
    #                                          designed for quantifying allele-specific expression. It is only valid if '--gtf' option is not specified. allele_id should be 
    #                                          the sequence names presented in the Multi-FASTA-formatted files. (Default: off)
    #           --polyA                     -> Add poly(A) tails to the end of all reference isoforms. The length of poly(A) tail added is specified by '--polyA-length' option. 
    #                                          STAR aligner users may not want to use this option. (Default: do not add poly(A) tail to any of the isoforms)
    #           --polyA-length <int>        -> The length of the poly(A) tails to be added. (Default: 125)
    #           --no-polyA-subset <file>    -> Only meaningful if '--polyA' is specified. Do not add poly(A) tails to those transcripts listed in <file>. <file> is a file 
    #                                          containing a list of transcript_ids. (Default: off)
    #           --bowtie                    -> Build Bowtie indices. (Default: off)
    #           --bowtie-path <path>        -> The path to the Bowtie executables. (Default: the path to Bowtie executables is assumed to be in the user's PATH environment variable)
    #           --bowtie2                   -> Build Bowtie 2 indices. (Default: off)
    #           --bowtie2-path              -> The path to the Bowtie 2 executables. (Default: the path to Bowtie 2 executables is assumed to be in the user's PATH environment variable)
    #           --star                      -> Build STAR indices. (Default: off)
    #           --star-path <path>          -> The path to STAR's executable. (Default: the path to STAR executable is assumed to be in user's PATH environment varaible)
    #           --star-sjdboverhang <int>   -> Length of the genomic sequence around annotated junction. It is only used for STAT to build splice junctions database and not 
    #                                          needed for Bowtie or Bowtie2. It will be passed as the value of 100 will work as well as the ideal value. (Default: 100)
    #           -p/--num-threads <int>      -> Number of threads to use for building STAR's genome indices. (Default: 1)
    #           -q/--quiet                  -> Suppress the output of logging information. (Default: off)
    #           -h/--help                   -> Show help information.
    # Description:                          -> This program extracts/preprocesses the reference sequences for RSEM. It can optionally build Bowtie indices (with '--bowtie' option) 
    #                                          and/or Bowtie 2 indices (with '--bowtie2' option) using their default parameters. It can also optionally build STAR indices (with 
    #                                          '--star' option) using parameters from ENCODE3's STAR-RSEM pipeline. If an alternative aligner is to be used, indices for that
    #                                          particular aligner can be built from either 'reference_name.idx.fa' or 'reference_name.n2g.idx.fa' (see OUTPUT for details). This 
    #                                          program is used in conjunction with the 'rsem-calculate-expression' program.
    # Output:                               -> This program will generate 'reference_name.grp', 'reference_name.ti', 'reference_name.transcripts.fa', 'reference_name.seq', 
    #                                          'reference_name.chrlist' (if '--gtf' is on), 'reference_name.idx.fa', 'reference_name.n2g.idx.fa', optional Bowtie/Bowtie 2 index 
    #                                          files, and optional STAR index files.
    #                                          'reference_name.grp', 'reference_name.ti', 'reference_name.seq', and 'reference_name.chrlist' are used by RSEM internally. 
    #                                          'reference_name.transcripts.fa' contains the extracted reference transcripts in Multi-FASTA format. Poly(A) tails are not added 
    #                                          and it may contain lower case bases in its sequences if the corresponding genomic regions are soft-masked.
    #                                          'reference_name.idx.fa' and 'reference_name.n2g.idx.fa' are used by aligners to build their own indices. In these two files, 
    #                                          all sequence bases are converted into upper case. In addition, poly(A) tails are added if '--polyA' option is set. The only 
    #                                          difference between 'reference_name.idx.fa' and 'reference_name.n2g.idx.fa' is that 'reference_name.n2g.idx.fa' in addition converts
    #                                          all 'N' characters to 'G' characters. This conversion is in particular desired for aligners (e.g. Bowtie) that do not allow reads to 
    #                                          overlap with 'N' characters in the reference sequences. Otherwise, 'reference_name.idx.fa' should be used to build the aligner's 
    #                                          index files. RSEM uses 'reference_name.idx.fa' to build Bowtie 2 indices and 'reference_name.n2g.idx.fa' to build Bowtie indices. 
    #                                          For visualizing the transcript-coordinate-based BAM files generated by RSEM in IGV, 'reference_name.idx.fa' should be imported as a 
    #                                          "genome" (see Visualization section in README.md for details).
    #
    # rsem-calculate-expression
    # -----------------------------
    # Estimate gene and isoform expression from RNA-Seq data.
    # Usage:                                -> rsem-calculate-expression [options] upstream_read_file(s) reference_name sample_name
    #                                       -> rsem-calculate-expression [options] --paired-end upstream_read_file(s) downstream_read_file(s) reference_name sample_name
    #                                       -> rsem-calculate-expression [options] --alignments [--paired-end] input reference_name sample_name
    # upstream_read_files(s)                -> Comma-separated list of files containing single-end reads or upstream reads for paired-end data. By default, these files are 
    #                                          assumed to be in FASTQ format. If the --no-qualities option is specified, then FASTA format is expected.
    # downstream_read_file(s)               -> Comma-separated list of files containing downstream reads which are paired with the upstream reads. By default, these files are 
    #                                          assumed to be in FASTQ format. If the --no-qualities option is specified, then FASTA format is expected.
    # input                                 -> SAM/BAM/CRAM formatted input file. If "-" is specified for the filename, the input is instead assumed to come from standard 
    #                                          input. RSEM requires all alignments of the same read group together. For paired-end reads, RSEM also requires the two mates of 
    #                                          any alignment be adjacent. In addition, RSEM does not allow the SEQ and QUAL fields to be empty. See Description section for 
    #                                          how to make input file obey RSEM's requirements.
    # reference_name                        -> The name of the reference used.  The user must have run 'rsem-prepare-reference' with this reference_name before running this program.
    # sample_name                           -> The name of the sample analyzed. All output files are prefixed by this name (e.g., sample_name.genes.results)
    # Options:  --paired-end                -> Input reads are paired-end reads. (Default: off)
    #           --no-qualities              -> Input reads do not contain quality scores. (Default: off)
    #           --strand-specific           -> The RNA-Seq protocol used to generate the reads is strand specific, i.e., all (upstream) reads are derived from the forward strand. 
    #                                          This option is equivalent to --forward-prob=1.0. With this option set, if RSEM runs the Bowtie/Bowtie 2 aligner, the '--norc' 
    #                                          Bowtie/Bowtie 2 option will be used, which disables alignment to the reverse strand of transcripts.  (Default: off)
    #           -p/--num-threads <int>      -> Number of threads to use. Both Bowtie/Bowtie2, expression estimation and 'samtools sort' will use this many threads. (Default: 1)
    #           --alignments                -> Input file contains alignments in SAM/BAM/CRAM format. The exact file format will be determined automatically. (Default: off)
    #           --fai <file>                -> RSEM reads header information from input by default. If this option is on, header information is read from the specified file. 
    #                                          For the format of the file, please see SAM official website. (Default: off)
    #           --bowtie2                   -> Use Bowtie 2 instead of Bowtie to align reads. Since currently RSEM does not handle indel, local and discordant alignments, the 
    #                                          Bowtie2 parameters are set in a way to avoid those alignments. In particular, we use options '--sensitive --dpad 0 --gbar 99999999 
    #                                       -> --mp 1,1 --np 1 --score-min L,0,-0.1' by default. The last parameter of '--score-min', '-0.1', is the negative of maximum mismatch 
    #                                          rate. This rate can be set by option '--bowtie2-mismatch-rate'. If reads are paired-end, we additionally use options '--no-mixed' 
    #                                          and '--no-discordant'. (Default: off)
    #           --star                      -> Use STAR to align reads. Alignment parameters are from ENCODE3's STAR-RSEM pipeline. To save computational time and memory resources, 
    #                                          STAR's Output BAM file is unsorted. It is stored in RSEM's temporary directory with name as 'sample_name.bam'. Each STAR job will 
    #                                          have its own private copy of the genome in memory. (Default: off)
    #           --append-names              -> If gene_name/transcript_name is available, append it to the end of gene_id/transcript_id (separated by '_') in files 
    #                                          'sample_name.isoforms.results' and 'sample_name.genes.results'. (Default: off)
    #           --seed <uint32>             -> Set the seed for the random number generators used in calculating posterior mean estimates and credibility intervals. The seed must 
    #                                          be a non-negative 32 bit integer. (Default: off)
    #           --single-cell-prior         -> By default, RSEM uses Dirichlet(1) as the prior to calculate posterior mean estimates and credibility intervals. However, much less 
    #                                          genes are expressed in single cell RNA-Seq data. Thus, if you want to compute posterior mean estimates and/or credibility intervals 
    #                                          and you have single-cell RNA-Seq data, you are recommended to turn on this option. Then RSEM will use Dirichlet(0.1) as the prior 
    #                                          which encourage the sparsity of the expression levels. (Default: off)
    #           --calc-pme                  -> Run RSEM's collapsed Gibbs sampler to calculate posterior mean estimates. (Default: off)
    #           --calc-ci                   -> Calculate 95% credibility intervals and posterior mean estimates. The credibility level can be changed by setting '--ci-credibility-level'. 
    #                                          (Default: off)
    #           -q/--quiet                  -> Suppress the output of logging information. (Default: off)
    #           -h/--help                   -> Show help information.
    #           --version                   -> Show version information.
    # Output Options:
    #           --sort-bam-by-read-name     -> Sort BAM file aligned under transcript coordidate by read name. Setting this option on will produce deterministic maximum likelihood 
    #                                          estimations from independent runs. Note that sorting will take long time and lots of memory. (Default: off)
    #           --no-bam-output             -> Do not output any BAM file. (Default: off)
    #           --sampling-for-bam          -> When RSEM generates a BAM file, instead of outputting all alignments a read has with their posterior probabilities, one alignment is 
    #                                          sampled according to the posterior probabilities. The sampling procedure includes the alignment to the "noise" transcript, which does 
    #                                          not appear in the BAM file. Only the sampled alignment has a weight of 1. All other alignments have weight 0. If the "noise" transcript 
    #                                          is sampled, all alignments appeared in the BAM file should have weight 0. (Default: off)
    #           --output-genome-bam         -> Generate a BAM file, 'sample_name.genome.bam', with alignments mapped to genomic coordinates and annotated with their posterior 
    #                                          probabilities. In addition, RSEM will call samtools (included in RSEM package) to sort and index the bam file. 'sample_name.genome.sorted.bam' 
    #                                          and 'sample_name.genome.sorted.bam.bai' will be generated. (Default: off)
    #           --sort-bam-by-coordinate    -> Sort RSEM generated transcript and genome BAM files by coordinates and build associated indices. (Default: off)
    #           --sort-bam-memory-per-thread-> Set the maximum memory per thread that can be used by 'samtools sort'. <string> represents the memory and accepts suffices 'K/M/G'. RSEM 
    #           <string>                       will pass <string> to the '-m' option of 'samtools sort'. Note that the default used here is different from the default used by samtools. 
    #                                          (Default: 1G)
    # Aligner Options:
    #           --seed-length <int>         -> Seed length used by the read aligner. Providing the correct value is important for RSEM. If RSEM runs Bowtie, it uses this value for 
    #                                          Bowtie's seed length parameter. Any read with its or at least one of its mates' (for paired-end reads) length less than this value will 
    #                                          be ignored. If the references are not added poly(A) tails, the minimum allowed value is 5, otherwise, the minimum allowed value is 25. 
    #                                          Note that this script will only check if the value >= 5 and give a warning message if the value < 25 but >= 5. (Default: 25)
    #           --phred33-quals             -> Input quality scores are encoded as Phred+33. (Default: on)
    #           --phred64-quals             -> Input quality scores are encoded as Phred+64 (default for GA Pipeline ver. >= 1.3). (Default: off)
    #           --solexa-quals              -> Input quality scores are solexa encoded (from GA Pipeline ver. < 1.3). (Default: off)
    #           --bowtie-path <path>        -> The path to the Bowtie executables. (Default: the path to the Bowtie executables is assumed to be in the user's PATH environment variable)
    #           --bowtie-n <int>            -> (Bowtie parameter) max # of mismatches in the seed. (Range: 0-3, Default: 2)
    #           --bowtie-e <int>            -> (Bowtie parameter) max sum of mismatch quality scores across the alignment. (Default: 99999999)
    #           --bowtie-m <int>            -> (Bowtie parameter) suppress all alignments for a read if > <int> valid alignments exist. (Default: 200)
    #           --bowtie-chunkmbs <int>     -> (Bowtie parameter) memory allocated for best first alignment calculation (Default: 0 - use Bowtie's default)
    #           --bowtie2-path <path>       -> (Bowtie 2 parameter) The path to the Bowtie 2 executables. (Default: the path to the Bowtie 2 executables is assumed to be in the user's 
    #                                          PATH environment variable)
    #           --bowtie2-mismatch-rate     -> (Bowtie 2 parameter) The maximum mismatch rate allowed. (Default: 0.1)
    #           <double>
    #           --bowtie2-k <int>           -> (Bowtie 2 parameter) Find up to <int> alignments per read. (Default: 200)
    #           --bowtie2-sensitivity-level -> (Bowtie 2 parameter) Set Bowtie 2's preset options in --end-to-end mode. This option controls how hard Bowtie 2 tries to find alignments. 
    #            <string>                      <string> must be one of "very_fast", "fast", "sensitive" and "very_sensitive". The four candidates correspond to Bowtie 2's "--very-fast", 
    #                                          "--fast", "--sensitive" and "--very-sensitive" options. (Default: "sensitive" - use Bowtie 2's default)
    #           --star-path <path>          -> The path to STAR's executable. (Default: the path to STAR executable is assumed to be in user's PATH environment variable)
    #           --star-gzipped-read-file    -> (STAR parameter) Input read file(s) is compressed by gzip. (Default: off)
    #           --star-bzipped-read-file    -> (STAR parameter) Input read file(s) is compressed by bzip2. (Default: off)
    #           --star-output-genome-bam    -> (STAR parameter) Save the BAM file from STAR alignment under genomic coordinate to 'sample_name.STAR.genome.bam'. This file is NOT sorted 
    #                                          by genomic coordinate. In this file, according to STAR's manual, 'paired ends of an alignment are always adjacent, and multiple alignments 
    #                                          of a read are adjacent as well'. (Default: off)
    # Advanced Options:
    #           --tag <string>              -> The name of the optional field used in the SAM input for identifying a read with too many valid alignments. The field should have the 
    #                                          format <tagName>:i:<value>, where a <value> bigger than 0 indicates a read with too many alignments. (Default: "")
    #           --forward-prob <double>     -> Probability of generating a read from the forward strand of a transcript. Set to 1 for a strand-specific protocol where all (upstream) 
    #                                          reads are derived from the forward strand, 0 for a strand-specific protocol where all (upstream) read are derived from the reverse strand, 
    #                                          or 0.5 for a non-strand-specific protocol. (Default: 0.5)
    #           --fragment-length-min <int> -> Minimum read/insert length allowed. This is also the value for the Bowtie/Bowtie2 -I option. (Default: 1)
    #           --fragment-length-max <int> -> Maximum read/insert length allowed. This is also the value for the Bowtie/Bowtie 2 -X option. (Default: 1000)
    #           --fragment-length-mean      -> (single-end data only) The mean of the fragment length distribution, which is assumed to be a Gaussian. (Default: -1, which disables use 
    #           <double>                       of the fragment length distribution)
    #           --fragment-length-sd        -> (single-end data only) The standard deviation of the fragment length distribution, which is assumed to be a Gaussian.  (Default: 0, 
    #           <double>                       which assumes that all fragments are of the same length, given by the rounded value of --fragment-length-mean)
    #           --estimate-rspd             -> Set this option if you want to estimate the read start position distribution (RSPD) from data. Otherwise, RSEM will use a uniform RSPD. 
    #                                          (Default: off)
    #           --num-rspd-bins <int>       -> Number of bins in the RSPD. Only relevant when '--estimate-rspd' is specified.  Use of the default setting is recommended. (Default: 20)
    #           --gibbs-burnin <int>        -> The number of burn-in rounds for RSEM's Gibbs sampler. Each round passes over the entire data set once. If RSEM can use multiple threads, 
    #                                          multiple Gibbs samplers will start at the same time and all samplers share the same burn-in number. (Default: 200)
    #           --gibbs-number-of-samples   -> The total number of count vectors RSEM will collect from its Gibbs samplers. (Default: 1000)
    #           <int>
    #           --gibbs-sampling-gap <int>  -> The number of rounds between two succinct count vectors RSEM collects. If the count vector after round N is collected, the count vector 
    #                                          after round N + <int> will also be collected. (Default: 1)
    #           --ci-credibility-level      -> The credibility level for credibility intervals. (Default: 0.95)
    #           <double>
    #           --ci-memory <int>           -> Maximum size (in memory, MB) of the auxiliary buffer used for computing credibility intervals (CI). (Default: 1024)
    #           --ci-number-of-samples-per- -> The number of read generating probability vectors sampled per sampled count vector. The crebility intervals are calculated by 
    #           count-vector <int>          -> first sampling P(C | D) and then sampling P(Theta | C) for each sampled count vector. This option controls how many Theta vectors 
    #                                          are sampled per sampled count vector. (Default: 50)
    #           --keep-intermediate-files   -> Keep temporary files generated by RSEM.  RSEM creates a temporary directory, 'sample_name.temp', into which it puts all intermediate 
    #                                          output files. If this directory already exists, RSEM overwrites all files generated by previous RSEM runs inside of it. By default, 
    #                                          after RSEM finishes, the temporary directory is deleted.  Set this option to prevent the deletion of this directory and the 
    #                                          intermediate files inside of it. (Default: off)
    #           --temporary-folder <string> -> Set where to put the temporary files generated by RSEM. If the folder specified does not exist, RSEM will try to create it. 
    #                                          (Default: sample_name.temp)
    #           --time                      -> Output time consumed by each step of RSEM to 'sample_name.time'. (Default: off)
    # Description:                          -> In its default mode, this program aligns input reads against a reference transcriptome with Bowtie and calculates expression values 
    #                                          using the alignments.  RSEM assumes the data are single-end reads with quality scores, unless the '--paired-end' or '--no-qualities' 
    #                                          options are specified. Alternatively, users can use STAR to align reads using the '--star' option. RSEM has provided options in 
    #                                          'rsem-prepare-reference' to prepare STAR's genome indices. Users may use an alternative aligner by specifying '--alignments', and 
    #                                          providing an alignment file in SAM/BAM/CRAM format. However, users should make sure that they align against the indices generated 
    #                                          by 'rsem-prepare-reference' and the alignment file satisfies the requirements mentioned in ARGUMENTS section. One simple way to 
    #                                          make the alignment file satisfying RSEM's requirements is to use the 'convert-sam-for-rsem' script. This script accepts SAM/BAM/CRAM 
    #                                          files as input and outputs a BAM file. For example, type the following command to convert a SAM file, 'input.sam', to a ready-for-use 
    #                                          BAM file, 'input_for_rsem.bam': convert-sam-for-rsem input.sam input_for_rsem For details, please refer to 'convert-sam-for-rsem's 
    #                                          documentation page.
    # Notes:                                -> 1. Users must run 'rsem-prepare-reference' with the appropriate reference before using this program.
    #                                       -> 2. For single-end data, it is strongly recommended that the user provide the fragment length distribution parameters 
    #                                             (--fragment-length-mean and --fragment-length-sd).  For paired-end data, RSEM will automatically learn a fragment length 
    #                                             distribution from the data.
    #                                       -> 3. Some aligner parameters have default values different from their original settings.
    #                                       -> 4. With the '--calc-pme' option, posterior mean estimates will be calculated in addition to maximum likelihood estimates.
    #                                       -> 5. With the '--calc-ci' option, 95% credibility intervals and posterior mean estimates will be calculated in addition to 
    #                                             maximum likelihood estimates.
    #                                       -> 6. The temporary directory and all intermediate files will be removed when RSEM finishes unless '--keep-intermediate-files' is specified.
    # Output:
    # sample_name.isoforms.results          -> File containing isoform level expression estimates. The first line contains column names separated by the tab character. The 
    #                                          format of each line in the rest of this file is:
    #                                            transcript_id gene_id length effective_length expected_count TPM FPKM IsoPct [posterior_mean_count posterior_standard_deviation_of_count 
    #                                            pme_TPM pme_FPKM IsoPct_from_pme_TPM TPM_ci_lower_bound TPM_ci_upper_bound TPM_coefficient_of_quartile_variation FPKM_ci_lower_bound 
    #                                            FPKM_ci_upper_bound FPKM_coefficient_of_quartile_variation]
    #                                          Fields are separated by the tab character. Fields within "[]" are optional. They will not be presented if neither '--calc-pme' nor 
    #                                          '--calc-ci' is set.
    #                                          - 'transcript_id' is the transcript name of this transcript. 
    #                                          - 'gene_id' is the gene name of the gene which this transcript belongs to (denote this gene as its parent gene). If no gene information is 
    #                                            provided, 'gene_id' and 'transcript_id' are the same.
    #                                          - 'length' is this transcript's sequence length (poly(A) tail is not counted). 
    #                                          - 'effective_length' counts only the positions that can generate a valid fragment. If no poly(A) tail is added,
    #                                          - 'effective_length' is equal to transcript length - mean fragment length + 1. If one transcript's effective length is less than 1, 
    #                                            this transcript's both effective length and abundance estimates are set to 0.
    #                                          - 'expected_count' is the sum of the posterior probability of each read comes from this transcript over all reads. Because 1) each read 
    #                                            aligning to this transcript has a probability of being generated from background noise; 2) RSEM may filter some alignable low quality 
    #                                            reads, the sum of expected counts for all transcript are generally less than the total number of reads aligned.
    #                                          - 'TPM' stands for Transcripts Per Million. It is a relative measure of transcript abundance. The sum of all transcripts' TPM is 1 million. 
    #                                          - 'FPKM' stands for Fragments Per Kilobase of transcript per Million mapped reads. It is another relative measure of transcript abundance. 
    #                                            If we define l_bar be the mean transcript length in a sample, which can be calculated as l_bar = \sum_i TPM_i / 10^6 * effective_length_i 
    #                                            (i goes through every transcript), the following equation is hold: FPKM_i = 10^3 / l_bar * TPM_i. We can see that the sum of FPKM is not 
    #                                            a constant across samples.
    #                                          - 'IsoPct' stands for isoform percentage. It is the percentage of this transcript's abandunce over its parent gene's abandunce. If its 
    #                                            parent gene has only one isoform or the gene information is not provided, this field will be set to 100.
    #                                          - 'posterior_mean_count', 'pme_TPM', 'pme_FPKM' are posterior mean estimates calculated by RSEM's Gibbs sampler. 
    #                                            'posterior_standard_deviation_of_count' is the posterior standard deviation of counts. 'IsoPct_from_pme_TPM' is the isoform percentage 
    #                                            calculated from 'pme_TPM' values.
    #                                          - 'TPM_ci_lower_bound', 'TPM_ci_upper_bound', 'FPKM_ci_lower_bound' and 'FPKM_ci_upper_bound' are lower(l) and upper(u) bounds of 95% 
    #                                            credibility intervals for TPM and FPKM values. The bounds are inclusive (i.e. [l, u]).
    #                                          - 'TPM_coefficient_of_quartile_variation' and 'FPKM_coefficient_of_quartile_variation' are coefficients of quartile variation (CQV) for 
    #                                            TPM and FPKM values. CQV is a robust way of measuring the ratio between the standard deviation and the mean. It is defined as
    #                                            CQV := (Q3 - Q1) / (Q3 + Q1), where Q1 and Q3 are the first and third quartiles.
    # sample_name.genes.results             -> File containing gene level expression estimates. The first line contains column names separated by the tab character. The format of each 
    #                                          line in the rest of this file is:
    #                                            gene_id transcript_id(s) length effective_length expected_count TPM FPKM [posterior_mean_count posterior_standard_deviation_of_count 
    #                                            pme_TPM pme_FPKM TPM_ci_lower_bound TPM_ci_upper_bound TPM_coefficient_of_quartile_variation FPKM_ci_lower_bound FPKM_ci_upper_bound 
    #                                            FPKM_coefficient_of_quartile_variation]
    #                                          Fields are separated by the tab character. Fields within "[]" are optional. They will not be presented if neither '--calc-pme' nor 
    #                                          '--calc-ci' is set. 'transcript_id(s)' is a comma-separated list of transcript_ids belonging to this gene. If no gene information is 
    #                                          provided, 'gene_id' and 'transcript_id(s)' are identical (the 'transcript_id').
    #                                          A gene's 'length' and 'effective_length' are defined as the weighted average of its transcripts' lengths and effective lengths 
    #                                          (weighted by 'IsoPct'). A gene's abundance estimates are just the sum of its transcripts' abundance estimates.
    # sample_name.alleles.results           -> Only generated when the RSEM references are built with allele-specific transcripts. This file contains allele level expression estimates 
    #                                          for allele-specific expression calculation. The first line contains column names separated by the tab character. The format of each line
    #                                          in the rest of this file is:
    #                                            allele_id transcript_id gene_id length effective_length expected_count TPM FPKM AlleleIsoPct AlleleGenePct [posterior_mean_count 
    #                                            posterior_standard_deviation_of_count pme_TPM pme_FPKM AlleleIsoPct_from_pme_TPM AlleleGenePct_from_pme_TPM TPM_ci_lower_bound 
    #                                            TPM_ci_upper_bound TPM_coefficient_of_quartile_variation FPKM_ci_lower_bound FPKM_ci_upper_bound FPKM_coefficient_of_quartile_variation]
    #                                          Fields are separated by the tab character. Fields within "[]" are optional. They will not be presented if neither '--calc-pme' nor 
    #                                          '--calc-ci' is set.
    #                                          - 'allele_id' is the allele-specific name of this allele-specific transcript. 
    #                                          - 'AlleleIsoPct' stands for allele-specific percentage on isoform level. It is the percentage of this allele-specific transcript's 
    #                                            abundance over its parent transcript's abundance. If its parent transcript has only one allele variant form, this field will be set to 100.
    #                                          - 'AlleleGenePct' stands for allele-specific percentage on gene level. It is the percentage of this allele-specific transcript's 
    #                                            abundance over its parent gene's abundance.
    #                                          - 'AlleleIsoPct_from_pme_TPM' and 'AlleleGenePct_from_pme_TPM' have similar meanings. They are calculated based on posterior mean estimates.
    #                                            Please note that if this file is present, the fields 'length' and 'effective_length' in 'sample_name.isoforms.results' should be 
    #                                            interpreted similarly as the corresponding definitions in 'sample_name.genes.results'.
    # sample_name.transcript.bam            -> Only generated when --no-bam-output is not specified. 'sample_name.transcript.bam' is a BAM-formatted file of read alignments in 
    #                                          transcript coordinates. The MAPQ field of each alignment is set to min(100, floor(-10 * log10(1.0 - w) + 0.5)), where w is the posterior 
    #                                          probability of that alignment being the true mapping of a read.  In addition, RSEM pads a new tag ZW:f:value, where value is a single 
    #                                          precision floating number representing the posterior probability. Because this file contains all alignment lines produced by bowtie 
    #                                          or user-specified aligners, it can also be used as a replacement of the aligner generated BAM/SAM file.
    # sample_name.transcript.sorted.bam and -> Only generated when --no-bam-output is not specified and --sort-bam-by-coordinate is specified. 'sample_name.transcript.sorted.bam' and
    # sample_name.transcript.sorted.bam.bai    'sample_name.transcript.sorted.bam.bai' are the sorted BAM file and indices generated by samtools (included in RSEM package).
    # sample_name.genome.bam                -> Only generated when --no-bam-output is not specified and --output-genome-bam is specified. 'sample_name.genome.bam' is a BAM-formatted 
    #                                          file of read alignments in genomic coordinates. Alignments of reads that have identical genomic coordinates (i.e., alignments to different
    #                                          isoforms that share the same genomic region) are collapsed into one alignment.  The MAPQ field of each alignment is set to min(100, 
    #                                          floor(-10 * log10(1.0 - w) + 0.5)), where w is the posterior probability of that alignment being the true mapping of a read.  In 
    #                                          addition, RSEM pads a new tag ZW:f:value, where value is a single precision floating number representing the posterior probability. 
    #                                          If an alignment is spliced, a XS:A:value tag is also added, where value is either '+' or '-' indicating the strand of the transcript it 
    #                                          aligns to.
    # sample_name.genome.sorted.bam and     -> Only generated when --no-bam-output is not specified, and --sort-bam-by-coordinate and --output-genome-bam are specified. 
    # sample_name.genome.sorted.bam.bai        'sample_name.genome.sorted.bam' and 'sample_name.genome.sorted.bam.bai' are the sorted BAM file and indices generated by samtools 
    #                                          (included in RSEM package).
    # sample_name.time                      -> Only generated when --time is specified. It contains time (in seconds) consumed by aligning reads, estimating expression levels and 
    #                                          calculating credibility intervals. 
    # sample_name.stat                      -> This is a folder instead of a file. All model related statistics are stored in this folder. Use 'rsem-plot-model' can generate plots 
    #                                          using this folder. 'sample_name.stat/sample_name.cnt' contains alignment statistics. The format and meanings of each field are described 
    #                                          in 'cnt_file_description.txt' under RSEM directory. 'sample_name.stat/sample_name.model' stores RNA-Seq model parameters learned 
    #                                          from the data. The format and meanings of each filed of this file are described in 'model_file_description.txt' under RSEM directory.
    
    # kallisto
    # --------------------------------------------------------------------------------------------------
    # Usage:                                -> kallisto <CMD> [arguments] ..
    #           <CMD>                       -> index         Builds a kallisto index 
    #                                       -> quant         Runs the quantification algorithm 
    #                                       -> h5dump        Converts HDF5-formatted results to plaintext
    #                                       -> version       Prints version information
    # 
    # kallisto index - Builds a kallisto index
    # -----------------------------
    # Usage:                                -> kallisto index [arguments] FASTA-files
    # Options:  -i, --index=STRING          -> Filename for the kallisto index to be constructed 
    #           -k, --kmer-size=INT         -> k-mer (odd) length (default: 31, max value: 31)
    #           --make-unique               -> Replace repeated target names with unique names
    #
    # kallisto quant - Computes equivalence classes for reads and quantifies abundances
    # -----------------------------
    # Usage: kallisto quant [arguments] FASTQ-files
    # Options:  -i, --index=STRING          -> Filename for the kallisto index to be used for quantification
    #           -o, --output-dir=STRING     -> Directory to write output to
    #           --bias                      -> Perform sequence based bias correction
    #           -b, --bootstrap-samples=INT -> Number of bootstrap samples (default: 0)
    #           --seed=INT                  -> Seed for the bootstrap sampling (default: 42)
    #           --plaintext                 -> Output plaintext instead of HDF5
    #           --single                    -> Quantify single-end reads
    #           -l,--fragment-length=DOUBLE -> Estimated average fragment length
    #           -s, --sd=DOUBLE             -> Estimated standard deviation of fragment length      (default: value is estimated from the input data)
    #           -t, --threads=INT           -> Number of threads to use (default: 1)
    #           --pseudobam                 -> Output pseudoalignments in SAM format to stdout
    #
    # kallisto h5dump - Converts HDF5-formatted results to plaintext
    # -----------------------------
    # Usage:                                -> kallisto h5dump [arguments] abundance.h5
    # Options:  -o, --output-dir=STRING     -> Directory to write output to

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # iReckon version 1.0.8
    # --------------------------------------------------------------------------------------------------
    # Requirements:                         -> java-1.6, latest version of BWA and added to PATH
    # Memory:                               -> 60M read pairs: 16G
    # Time:                                 -> 60M rp: 24 h
    # Space:                                -> 60M rp: 100G
    # Usage:                                -> java -Xmx15000M -jar iReckon-1.0.7.jar <bam_file> <reference_genome> annotations -1 reads_file [options]
    # bam_file:                             -> Alignment file in bam format. Need to be indexed (bam.bai). (Should contain split reads alignments (Like TopHat output))
    # reference_genome:                     -> In fasta format. It must be the same reference used for creating the bam_file.
    # annotations                           -> Known transcripts annotations +/ Known start sites/end sites can be preprocessed with FormatTools or Savant
    #                                          Genome Browser (FormatTools is available with Savant) preprocessing can also be automatically done by iReckon 
    #                                          (functionality tested with UCSC tab-delimited refGene files but should work with 12-columns BED format as well). 
    #                                          GTF and other formats are not accepted       
    # Output:                               -> IReckon output is a description of the found isoforms and their quantities contained in the file output_dir/result.gtf
    # Options:  -1 reads_file               -> The RNA-Seq reads file in fasta or fastq format.
    #           -2 reads_file2              -> The read mates file. If -2 is not used then the reads in reads_file1 have to be listed pair by pair. 
    #           -o output_dir               -> Select the output repository . Anticipate enough memory.
    #           -m method                   -> The tool used to align reads to isoforms. 0 for Shrimp , 2 for BWA. By default =2.
    #           -n nb_thread                -> Determine the number of threads.
    #           -b bias                     -> Select the magnitude for the gene-borders bias correction. Bias correction:  0=none, 1=weak, 2=strong. Default=0.
    #                                          Similar to Effective Isoform Length Correction (Cufflinks). Works only for isoforms longer than 200 nt.
    #           -d                          -> Use simultaneous duplicates removal.
    #           -q                          -> Quick re_run. Useful if you have recently run iReckon on the same data, in the same repository(different options).
    #           -ign                        -> Ignore splicing that do not start/end at a known splicing border site. Useless with -q. Default=false
    #           -novel                      -> Enable/ disable novel isoforms discovery. 0=Disabled, 1=Enabled. Useless with -q. Default=Enabled        
    #           -nbi nbIso_Max              -> Maximum number of correctly constructed isoforms per gene. If this number is exceeded heuristics will be applied to
    #                                          remove rare junctions in order to reduce the complexity. Useless with -q. Default=100.
    #           -minrec nbrec               -> Minimum number of records aligned to a gene to allow study. Useless with -q. Default=10.
    #           -maxref refsize             -> Maximum size of a reference file to be indexed by bwa (GO). If transcriptome is bigger, it will be split. It is 
    #                                         advisable to use higher values when enough RAM is available. (The RAM usage is determined by BWA indexing 
    #                                         algorithm: linear with reference size). Useless with -q. Default=6
    #           -chr chr-Name               -> Give a comma-separated list of the chromosomes to be investigated by iReckon. Be careful that the FPKM computed
    #                                          are normalized by the number of reads mapped there. Useless with -q.
    #           -start int                  -> Used with -chr and -end to specify the location studied by iReckon. Useless unless exactly one chromosome is         
    #                                          specified. Useless with -q. Default=1.
    #           -end int                    -> Used with -chr and -start to specify the location studied. Useless unless exactly one chromosome is specified (-chr)               
    #                                          Useless with -q. Default=End_of_Chromosome.
    # Example:                              -> java -Xmx15000M  -jar IReckon-1.0.8.jar alignment.bam reference.fa.savant hg19.refGene.gz -1 reads.fastq  -o /disk/ireckon_output/ -d -n 8 -b 2 -nbi 100 > logs.txt
    
    
    # FlipFlop
    # --------------------------------------------------------------------------------------------------
    # Then in R, you can load the package with library(flipflop) . The command example(flipflop) will show you a working example.
    # Options:  data.file [input]           -> Input alignment file in SAM format. The SAM file must be sorted according to chromosome name and starting position. 
    #                                          If you use a multi-sample strategy, please give a single SAM file that results from the "samtools merge" command with 
    #                                          the "-r" option (i.e attach RG tag infered from file name). 
    #           out.file [output]           -> Output gtf file storing the structure of the transcripts which are found to be expressed together with their abundances
    #                                          (in FPKM and expected count).
    #           output.type [output]        -> Type of output when using several samples simultaneously. When equal to "gtf" the output corresponds to a gtf file per 
    #                                          sample with specific abundances. When equal to "table" the output corresponds to a single gtf file storing the structure
    #                                          of the transcripts, and an associated table with the abundances of all samples (transcripts in rows and samples in columns). Default "gtf".
    #           annot.file [input]          -> Optional annotation file in BED12 format. If given, exon boundaries will be taken from the annotations. The BED file should 
    #                                          be sorted according to chromosome name and starting position of transcripts.
    #           samples.file [multi-samples]-> Optional samples file (one line per sample name, as shown here. ). The names should be the one present in the RG 
    #                                          tag of the provided SAM file. 
    #           mergerefit [multi-samples]  -> If TRUE use a simple refit strategy instead of the group-lasso. Default FALSE.
    #           paired [paired]             -> Boolean for paired-end reads. If FALSE your reads will be considered as single-end reads. Default FALSE.
    #           frag [paired]               -> Mean fragment size. Only used if paired is set to TRUE. Default 400. NOTE: the fragment size is equal to insert + 2*read_length. 
    #                                          You can estimate it using PicardTools or Wei Li's script.
    #           std [paired]                -> Standard deviation of fragment size. Only used if paired is set to TRUE. Default 20.
    #           OnlyPreprocess [pre-pro]    -> Boolean for performing only the pre-processing step. Output is two files: one file '.instance' and one other file 
    #                                          '.totalnumread'. Default FALSE. NOTE: The pre-processing is mainly borrowed from the IsoLasso software maintained by Wei Li.
    #           preprocess.instance [pre-pro]-> Give directly the pre-processed '.instance' input file created when using the OnlyPreprocess option. If non empty, the 
    #                                          data.file and annot.file fields are ignored.
    #           minReadNum [pre-pro]        -> The minimum number of clustered reads to output. Default 40. If you give an annotation file it will be the minimum number of 
    #                                          mapped reads to process a gene. You might want to increase this value if you want to filter genes with only few reads.
    #           minFragNum [pre-pro]        -> The minimum number of mapped read pairs to process a gene. Only used if paired is TRUE. Default 20. You might want to increase 
    #                                          this value if you want to filter genes with only few read pairs.
    #           minCvgCut [pre-pro]         -> The fraction for coverage cutoff, should be between 0-1. A higher value will be more sensitive to coverage discrepancies in 
    #                                          one gene. Default 0.05.
    #           minJuncCout [pre-pro]       -> The minimum number of reads to consider a junction as valid. Default 1. You might want to increase this value, in particular 
    #                                          when you use many samples or high coverage data.
    #           sliceCount [paralleliz]     -> Number of slices in which the pre-processing '.instance' file is cut. It creates several instance files with the extension 
    #                                          '_jj.instance' where jj is the number of the slice. If you set OnlyPreprocess to TRUE, it will create those slices and you 
    #                                           can run FlipFlop independently afterwards on each one of the slice. Default 1.
    #           mc.cores [paralleliz]       -> Number of cores. If you give sliceCount>1 with OnlyPreprocess set to FALSE, it will distribute the slices on several cores. 
    #                                          If you give a preprocess.instance file as input (which might be a slice of an original instance file), it will use several 
    #                                          cores when you are using a multi-samples strategy. Default 1.
    #           verbose [verbosity]         -> Verbosity. Default 0 (little verbosity). Put 1 for more verbosity.
    #           verbosepath [verbosity]     -> Verbosity of the optimization part. Default 0 (little verbosity). Put 1 for more verbosity. Mainly used for de-bugging.
    #           cutoff [parameter]          -> Do not report isoforms whose expression level is less than cutoff percent of the most expressed transcripts. Default 1.
    #           BICcst [parameter]          -> Constant used for model selection with the BIC criterion. Default 50.
    #           delta [parameter]           -> Loss parameter, Poisson offset. Default 1e-7.
    #           use_TSSPAS                  -> Do we restrict the candidate TSS and PAS sites. 1 is yes and 0 is no. Default 0 i.e each exon can possibly starts or ends an isoform.
    #           max_isoforms                -> Maximum number of isoforms given during regularization path. Default 10.
    # Usage:                                -> run flipflop on a SAM file from 1 sample (named in.sam), and output a GTF out.gtf: 
    #                                          - res=flipflop(data.file='in.sam', out.file='out.gtf')
    #                                          first run the pre-processing step only (it writes 2 files, named in.instance and in.totalnumread), then run the prediction step: 
    #                                          - flipflop(data.file='in.sam', OnlyPreprocess=TRUE) 
    #                                          - res=flipflop(preprocess.instance='in.instance', out.file='out.gtf')
    #                                          run flipflop on a merged SAM file from a few sample (named merge.sam), with the sample names stored in the file mysamples: 
    #                                          - res=flipflop(data.file='merge.sam', samples.file='mysamples', out.file='out.gtf') 
    #                                          This should output 3 GTF, named out_sample1.gtf, out_sample2.gtf and out_sample3.gtf.
    # different output type:                -> - res=flipflop(data.file='merge.sam', samples.file='mysamples', out.file='out.gtf', output.type='table') 
    #                                          This should output 1 GTF out.gtf storing the structure of the transcripts, and 2 tables out_count and out_fpkm storing the 
    #                                          abundances of the transcripts in each sample, in both count value and FPKM.
    # NOTE: if you are using strand-specific libraries, you should separate reads from forward or reverse strands, and run flipflop separately on two SAM files. 
    # For instance, you can separate such reads after mapping with Tophat with grep XS:A:+ in.sam > forward.sam and with grep XS:A:- in.sam > reverse.sam .
    
    
    # Running StringTie
    # --------------------------------------------------------------------------------------------------
    # Usage:                                -> stringtie <aligned_reads.bam> [options]*
    #                                          The main input of the program is a BAM file with RNA-Seq read mappings which must be sorted by their genomic location 
    #                                          (for example the accepted_hits.bam file produced by TopHat or the output of HISAT2 after sorting and converting it using 
    #                                          samtools as explained below).
    # Options:  -h/--help                   -> Prints help message and exits.
    #           -v                          -> Turns on verbose mode, printing bundle processing details.
    #           -o [<path/>]<out.gtf>       -> Sets the name of the output GTF file where StringTie will write the assembled transcripts. This can be specified as a full path, 
    #                                          in which case directories will be created as needed. By default StringTie writes the GTF at standard output.
    #           -p <int>                    -> Specify the number of processing threads (CPUs) to use for transcript assembly. The default is 1.
    #           -G <ref_ann.gff>            -> Use the reference annotation file (in GTF or GFF3 format) to guide the assembly process. The output will include expressed reference 
    #                                          transcripts as well as any novel transcripts that are assembled. This option is required by options -B, -b, -e, -C (see below).
    #           -l <label>                  -> Sets <label> as the prefix for the name of the output transcripts. Default: STRG
    #           -f <0.0-1.0>                -> Sets the minimum isoform abundance of the predicted transcripts as a fraction of the most abundant transcript assembled at a given 
    #                                          locus. Lower abundance transcripts are often artifacts of incompletely spliced precursors of processed transcripts. Default: 0.1
    #           -m <int>                    -> Sets the minimum length allowed for the predicted transcripts. Default: 200
    #           -A <gene_abund.tab>         -> Gene abundances will be reported (tab delimited format) in the output file with the given name.
    #           -C <cov_refs.gtf>           -> StringTie outputs a file with the given name with all transcripts in the provided reference file that are fully covered by reads (requires -G).
    #           -a <int>                    -> Junctions that don't have spliced reads that align across them with at least this amount of bases on both sides are filtered out. Default: 10
    #           -j <float>                  -> There should be at least this many spliced reads that align across a junction (i.e. junction coverage). This number can be fractional, 
    #                                          since some reads align in more than one place. A read that aligns in n places will contribute 1/n to the junction coverage. Default: 1
    #           -t                          -> This parameter disables trimming at the ends of the assembled transcripts. By default StringTie adjusts the predicted transcript's start 
    #                                          and/or stop coordinates based on sudden drops in coverage of the assembled transcript.
    #           -c <float>                  -> Sets the minimum read coverage allowed for the predicted transcripts. A transcript with a lower coverage than this value is not 
    #                                          shown in the output. Default: 2.5
    #           -g <int>                    -> Minimum locus gap separation value. Reads that are mapped closer than this distance are merged together in the same processing bundle. 
    #                                          Default: 50 (bp)
    #           -B                          -> This switch enables the output of Ballgown input table files (*.ctab) containing coverage data for the reference transcripts given with 
    #                                          the -G option. (See the Ballgown documentation for a description of these files.) With this option StringTie can be used as a direct 
    #                                          replacement of the tablemaker program included with the Ballgown distribution. 
    #                                          If the option -o is given as a full path to the output transcript file, StringTie will write the *.ctab files in the same directory 
    #                                          as the output GTF.
    #           -b <path>                   -> Just like -B this option enables the output of *.ctab files for Ballgown, but these files will be created in the provided directory 
    #                                          <path> instead of the directory specified by the -o option. Note: adding the -e option is recommended with the -B/-b options, unless 
    #                                          novel transcripts are still wanted in the StringTie GTF output.
    #           -e                          -> Limits the processing of read alignments to only estimate and output the assembled transcripts matching the reference transcripts 
    #                                          given with the -G option (requires -G, recommended for -B/-b). With this option, read bundles with no reference transcripts will 
    #                                          be entirely skipped, which may provide a considerable speed boost when the given set of reference transcripts is limited to a set 
    #                                          of target genes, for example.
    #           -M <0.0-1.0>                -> Sets the maximum fraction of muliple-location-mapped reads that are allowed to be present at a given locus. Default: 0.95.
    #           -x <seqid_list>             -> Ignore all read alignments (and thus do not attempt to perform transcript assembly) on the specified reference sequences. Parameter 
    #                                          <seqid_list> can be a single reference sequence name (e.g. -x chrM) or a comma-delimited list of sequence names (e.g. -x 'chrM,chrX,chrY'). 
    #                                          This can speed up StringTie especially in the case of excluding the mitochondrial genome, whose genes may have very high coverage in some 
    #                                          cases, even though they may be of no interest for a particular RNA-Seq analysis. The reference sequence names are case sensitive, they must 
    #                                          match identically the names of chromosomes/contigs of the target genome against which the RNA-Seq reads were aligned in the first place.
    #           --merge                     -> Transcript merge mode. This is a special usage mode of StringTie, distinct from the assembly usage mode described above. In the merge mode, 
    #                                          StringTie takes as input a list of GTF/GFF files and merges/assembles these transcripts into a non-redundant set of transcripts. This mode 
    #                                          is used in the new differential analysis pipeline to generate a global, unified set of transcripts (isoforms) across multiple RNA-Seq samples. 
    #
    # Transfrags:                           -> If the -G option (reference annotation) is provided, StringTie will assemble the transfrags from the input GTF files with the reference transcripts. 
    # Options:  -G <guide_gff>              -> reference annotation to include in the merging (GTF/GFF3)
    #           -o <out_gtf>                -> output file name for the merged transcripts GTF (default: stdout)
    #           -m <min_len>                -> minimum input transcript length to include in the merge (default: 50)
    #           -c <min_cov>                -> minimum input transcript coverage to include in the merge (default: 0)
    #           -F <min_fpkm>               -> minimum input transcript FPKM to include in the merge (default: 0)
    #           -T <min_tpm>                -> minimum input transcript TPM to include in the merge (default: 0)
    #           -f <min_iso>                -> minimum isoform fraction (default: 0.01)
    #           -i                          -> keep merged transcripts with retained introns (default: these are not kept unless there is strong evidence for them)
    #           -l <label>                  -> name prefix for output transcripts (default: MSTRG)
    # Input files:
    #                                       -> StringTie takes as input a binary SAM (BAM) file sorted by reference position. This file contains spliced read alignments and can be 
    #                                          produced directly by programs such as TopHat or it can be obtained by converting and sorting the output of HISAT2. We recommend using 
    #                                          HISAT2 as it is a fast and accurate alignment program. A text file in SAM format which was produced by HISAT2 must be sorted and 
    #                                          converted to BAM format using the samtools program:
    #                                          samtools view -Su alns.sam | samtools sort - alns.sorted
    #                                       -> The file resulted from the above command (alns.sorted.bam) can be used as input for StringTie.
    #                                          Every spliced read alignment (i.e. an alignment across at least one junction) in the input SAM file must contain the tag XS to indicate 
    #                                          the genomic strand that produced the RNA from which the read was sequenced. Alignments produced by TopHat and HISAT2 (when ran with 
    #                                          --dta option) already include this tag, but if you use a different read mapper you should check that this XS tag is included for spliced 
    #                                          alignments.
    #                                       -> Optionally, a reference annotation file in GTF/GFF3 format can be provided to StringTie. In this case, StringTie will check to see if the 
    #                                          reference transcripts are expressed in the RNA-Seq data, and for the ones that are expressed it will compute coverage, TPM and FPKM values. 
    #                                          Note that if option -e is not used the reference transcripts need to be fully covered by reads in order to be included in StringTie's output. 
    #                                          In that case, other transcripts assembled from the data by StringTie and not present in the reference file will be printed as well.
    # Evaluating transcript assemblies:     -> A simple way of getting more information about the transcripts assembled by StringTie (summary of gene and transcript counts, novel vs. 
    #                                          known etc.), or even performing basic tracking of assembled isoforms across multiple RNA-Seq experiments, is to use the gffcompare program. 
    #                                          Basic usage information and download options for this program can be found on the GFF utilities page.
    # Differential expression analysis:     -> StringTie can be used in a RNA-Seq analysis protocol similar to the one described originally for Cufflinks. Instead of Cuffdiff, StringTie 
    #                                          can use Ballgown as the final analysis step for estimating differential expression across multiple RNA-Seq samples and generating plots 
    #                                          and differential expression tables. DE workflow
    # Recommended workflow:                 -> for each RNA-Seq sample, map the reads to the genome with HISAT2 using the --dta option. It is highly recommended to use the reference 
    #                                          annotation information when mapping the reads, which can be either embedded in the genome index (built with the --ss and --exon options, 
    #                                          see HISAT2 manual), or provided separately at run time (using the --known-splicesite-infile option of HISAT2). The SAM output of each HISAT2 
    #                                          run must be sorted and converted to BAM using samtools as explained above.
    #                                       -> for each RNA-Seq sample, run StringTie to assemble the read alignments obtained in the previous step; it is recommended to run StringTie 
    #                                          with the -G option if the reference annotation is available.
    #                                       -> run StringTie with --merge in order to generate a non-redundant set of transcripts observed in all the RNA-Seq samples assembled previously. 
    #                                          The stringtie --merge mode takes as input a list of all the assembled transcripts files (in GTF format) previously obtained for each sample, 
    #                                          as well as a reference annotation file (-G option) if available.
    #                                       -> for each RNA-Seq sample, run StringTie using the -B/-b and -e options in order to estimate transcript abundances and generate read coverage 
    #                                          tables for Ballgown. The -e option is not required but recommended for this run in order to produce more accurate abundance estimations of 
    #                                          the input transcripts. Each StringTie run in this step will take as input the sorted read alignments (BAM file) obtained in step 1 for the 
    #                                          corresponding sample and the -G option with the merged transcripts (GTF file) generated by stringtie --merge in step 3. Please note that this 
    #                                          is the only case where the -G option is not used with a reference annotation, but with the global, merged set of transcripts as observed 
    #                                          across all samples. (This step is the equivalent of the Tablemaker step described in the original Ballgown pipeline.)
    #                                       -> Ballgown can now be used to load the coverage tables generated in the previous step and perform various statistical analyses for differential 
    #                                          expression, generate plots etc.
    # Alternate, faster workflow;           -> if there is no interest in novel isoforms (i.e. assembled transcripts present in the samples but missing from the reference annotation), or 
    #                                          if only a well known set of transcripts of interest are targeted by the analysis. This simplified protocol has only 3 steps (depicted below) 
    #                                          as it bypasses the individual assembly of each RNA-Seq sample and the "transcript merge" step. DE workflow This simplified workflow attempts 
    #                                          to directly estimate and analyze the expression of a known set of transcripts as given in the reference annotation file.
    # Note: Starting with version 1.0.1 StringTie supports the -e option which is strongly recommended for this alternate workflow, as it limits the processing of read alignments to only 
    # those overlapping the targeted reference loci. 
    
    
    # eXpress
    # --------------------------------------------------------------------------------------------------
    # Usage:                                -> express [options]* <target_seqs.fasta> <aligned_reads.(sam/bam)>
    # Options:  <target_seqs.fasta>         -> A file of target sequences in multi-FASTA format. See Input Files for more details.
    #           <lib_1.sam,...,lib_N.sam>   -> A comma-separated list of filenames for reads aligned to the target sequences in SAM format. See Input Files for more details.
    #           -h/--help                   -> Prints the help message and exits
    #           -o/--output-dir <string>    -> Sets the name of the directory in which eXpress will write all of its output. The default is "./".
    #           -B/--additional-batch <int> -> Specifies the number of additional batch EM rounds to perform on the data using the initial results from the online EM as a seed. Can 
    #                                          improve accuracy at the cost of time.
    #           -O/--additional-online <int>-> Specifies the number of additional online EM rounds to perform on the data after the initial online round. Can improve accuracy at 
    #                                          the cost of time.
    #           -m/--frag-len-mean <int>    -> Specifies the mean fragment length. While the empirical distribution is estimated from paired-end reads on-the-fly, this value 
    #                                          paramaterizes the prior distribution. If only single-end reads are available, this prior distribution is also used to determine the 
    #                                          effective length. Default is 200.
    #           -s/--frag-len-stddev <int>  -> Specifies the fragment length standard deviation. While the empirical distribution is estimated from paired-end reads on-the-fly, this 
    #                                          value paramaterizes the prior distribution. If only single-end reads are available, this prior distribution is also used to determine 
    #                                          the effective length. Default is 60.
    #           -H/--haplotype-file <string>-> Specifies the location of a comma-separated file of sets of target IDs (one set per line) specifying which targets represent multiple 
    #                                          haplotypes of a single feature (ie, transcript). Useful for allele-specific expression.
    #           --output-align-prob         -> With this option, eXpress outputs an additional file called hits.prob.(sam/bam) containing identical copies of all input alignments 
    #                                          with an additional XP tag that contains the estimated probability that each alignment of the read (pair) is the "correct" one. The 
    #                                          XP values for all alignments of of the same read (pair) will sum to 1.
    #           --output-align-samp         -> With this option, eXpress outputs an additional file called hits.samp.(sam/bam) containing a single alignment for each fragment 
    #                                          sampled at random based on the alignment likelihoods calculated by eXpress.
    #           --fr-stranded               -> With this option, eXpress only accepts alignments (single-end or paired) where the first (or only) read is aligned to the forward 
    #                                          target sequence and the second read is aligned to the reverse-complemented target sequence. In directional sequencing, this is 
    #                                          equivalent to second-strand only. If all reads are single-end, --f-stranded should be used instead. Disabled by default.
    #           --rf-stranded               -> With this option, eXpress only accepts alignments (single-end or paired) where the first (or only) read is aligned to the 
    #                                          reverse-completemented target sequence and the second read is aligned to the forward target sequence. In directional sequencing, 
    #                                          this is equivalent to first-strand only. If all reads are single-end, --r-stranded should be used instead. Disabled by default.
    #           --f-stranded                -> With this option, eXpress only accepts single-end alignments to the forward target sequence. In directional sequencing, this is 
    #                                          equivalent to second-strand only. Disabled by default.
    #           --r-stranded                -> With this option, eXpress only accepts single-end alignments to the reverse target sequence. In directional sequencing, this is 
    #                                          equivalent to second-strand only. Disabled by default.
    #           --no-update-check           -> With this option, eXpress will not ping our server to see if a newer version is available.
    # Advanced: -f/--forget-param <float>   -> A parameter specifying the rate at which the prior is "forgotten" by increasing the mass of fragments during online processing. 
    #                                          Larger numbers (max of 1) mean a slower rate, which decreases convergence but improves stability. Smaller numbers (minumum of 0.5) 
    #                                          increase the rate, which may lead to faster convergence but can also lead to instability.
    #           --library-size <int>        -> Specifies the number of fragments in the library to be used in the FPKM calculation. If left unspecified, this number will be computed 
    #                                          from the input.
    #           --max-indel-size <int>      -> A parameter specifying the maximum allowed size of a single indel. Alignments with larger indels will be ignored. A geometric prior for 
    #                                          indel length is fit so that all but 10e-6 of the probability mass lies within the allowed region. The default is 10.
    #           --calc-covar                -> With this option, eXpress calculates the covariance between targets and outputs them for use in differential expression analysis. This 
    #                                          calculation requires slightly more time and memory.
    #           --expr-alpha <float>        -> A parameter specifying the weight of uniform the target abundance prior, in pseudo-counts per bp. The default is 0.01.
    #           --stop-at <int>             -> A parameter specifying the number of fragments to process before quitting.
    #           --burn-out <int>            -> A parameter specifying the number of fragments after which to stop learning the auxiliary parameters (fragment length, bias, error).
    #           --no-bias-correct           -> With this option, eXpress will not measure and account for sequence-specific biases. Will lead to a slight initial increase in speed 
    #                                         at the expense of accuracy.
    #           --no-error-model            -> With this option, eXpress will not measure and account for errors in alignments. Will lead to an increase in speed, but may greatly 
    #                                          decrease accuracy.
    #           --aux-param-file <string>   -> Specifies an auxiliary parameter file output by a different run of eXpress to be used as the auxiliary parameters for this round. 
    #                                          Greatly improves speed and should be used when a subset of the targets or fragments are being used in a second estimation.
    #
    # Input Files:
    # Target Sequences (FASTA):             -> eXpress requires a multi-FASTA file of target sequences for which the abundances will be measured. In the case of RNA-Seq, these are 
    #                                          the transcript sequences. If the transcriptome of your organism is not annotated, you can generate this file from your sequencing 
    #                                          reads using a de novo transcriptome assembler such as Trinity, Oases, or Trans-ABySS. If your organism has a reference genome, you 
    #                                          can assemble transcripts directly from mapped reads using Cufflinks. If your genome is already annotated (in GTF/GFF), you can generate 
    #                                          a multi-FASTA file using the UCSC Genome Browser by uploading your annotation as a track and downloading the sequences under the "Tables" tab.
    # Read Alignments (SAM/BAM):            -> eXpress also requires a file, multiple files, or a piped stream of SAM or binary SAM (BAM) alignments as input. The SAM alignments 
    #                                          should be generated by mapping your sequencing reads to the target sequences specified in the multi-FASTA input file described above. 
    #                                          For more details on the SAM format, see the specification. Many short read mappers including Bowtie, Bowtie2, BWA, and MAQ can produce 
    #                                          output in this format. It is important that you allow many multi-mappings (preferably unlimited) in order to allow eXpress to select 
    #                                          the correct alignment instead of the mapper. See Getting Started for an example using Bowtie in both streaming and file input modes.
    #                                       -> If using paired-end reads, the read names must match for each pair, excluding '/1' and '/2' suffix identifiers. Also, the SAM file 
    #                                          supplied to eXpress should be grouped by read id. If you aligned your reads with Bowtie, your alignments will be properly ordered 
    #                                          already. If you used another tool, you should ensure that they are properly sorted. You can sort your SAM using the following command:
    #                                       -> sort -k 1 hits.sam > hits.sam.sorted
    # Sort your BAM using this command:     -> samtools sort -n hits.bam hits.sorted
    #                                       -> If multiple libraries were prepared for the same sample or multiple read lengths were used in different sequencing runs, the alignments 
    #                                          for each should be grouped in separate SAM files so that auxiliary parameters can be estimated independently. The filenames can then be 
    #                                          input into eXpress as a comma-separated (with no spaces) list of SAM files. See above for an example. When this feature is used, separate 
    #                                          parameter estimates will be output for each library, but only a single abundance file will be produced.
    # Output Files:
    # Target Abundances (results.xprs):     -> This file is always output and contains the target abundances and other values calculated based on the input sequences and read alignments. 
    #                                          The file has 10 tab-delimited columns, sorted by the bundle_id (column 1).
    # The columns are defined as follows:      #   Column Name         Example         Description
    #                                          1   bundle_id           10              ID of bundle the target belongs to. A bundle is defined as the transitive closure of targets that 
    #                                                                                  share multi-mapping reads.
    #                                          2   target_id           NM_016467       The ID given to the target in the input multi-FASTA file.
    #                                          3   length              2182            The number of base pairs in the target sequence given in the input multi-FASTA file.
    #                                          4   eff_length          783.136288      The length of the target adjusted for fragment biases (length, sequence-specificity, and relative 
    #                                                                                  position). This number is what the fragment counts are normalized by to calculate FPKM, not the true length.
    #                                          5   tot_counts          99              The number of fragments mapping (uniquely or ambiguously) to this target.
    #                                          6   uniq_counts         7               The number of fragments uniquely mapping to this target.
    #                                          7   est_counts          26.702456       The estimated number of fragments generated from this target in the sequencing experiment.
    #                                          8   eff_counts          74.399258       The estimated number of fragments generated from this target in the sequencing experiment, adjusted 
    #                                                                                  for fragment and length biases. In other words, his is the expected number of reads from the experiment 
    #                                                                                  if these biases did not exist. This is the value recommended for input to count-biased differential 
    #                                                                                  expression tools.
    #                                          9   ambig_distr_alpha   3.154652        The alpha parameter for the posterior beta-binomial distribution fit to the ambiguous reads.
    #                                          10  ambig_distr_beta    2.293653        The beta parameter for the posterior beta-binomial distribution fit to the ambiguous reads.
    #                                          11  fpkm                3.514176        The estimated relative abundance of this target in the sample in units of fragments per kilobase per 
    #                                                                                  million mapped. This value is proportional to est_counts divided by eff_length.
    #                                          12  fpkm_conf_low       2.119151        The lower bound of the 95% confidence interval for the FPKM.
    #                                          13  fpkm_conf_high      4.909200        The upper bound of the 95% confidence interval for the FPKM.
    #                                          14  solvable            T               A binary (T/F) value indicating whether the likelihood function has a unique maximum. If false (F), 
    #                                                                                  the reported posterior distribution is uniform.
    #                                          15  tpm                 2.347222e+05    Transcripts per million. See description.
    #                                       -> See the Methods for more details on how these values are calculated.
    # Parameter Estimates (params.xprs):    -> This file contains the values of the other parameters (besides abundances and counts) estimated by eXpress. The file is separated 
    #                                          into sections for each parameter type, beginning with a '>' symbol. Following this symbol is the section header containing a name 
    #                                          for the parameter type followed by the values on subsequent lines. All values belong to this parameter field until the next '>' or 
    #                                          the end of the file. The following parameter types are output to this file:
    #                                       #   Parameter Type                  Description                                         Output Format
    #                                       1   Fragment Length Distribution    The empirical distribution on fragments lengths.    The fragment length range is listed next to the 
    #                                                                                                                               section header in paranthesis (0-800 by default). 
    #                                                                                                                               The next line contains a tab-delimited list of 
    #                                                                                                                               probabilities for these lengths in order.
    #                                       2   First Read Mismatch             The first-order Markov model for mismatches         Each line begins with the nucleotide position in the 
    #                                                                           between the reference and observed nucleotides      read followed by a colon (0-indexed). The column header 
    #                                                                           for the first read sequenced in a pair.             denotes the which "substitution" the probability is for. 
    #                                                                                                                               For example, a value in the column labeled "CG->*T" in 
    #                                                                                                                               the row labeled 10 is the conditional probability that 
    #                                                                                                                               a read has a 'T' at the 11th position given it is mapped 
    #                                                                                                                               to a reference having a 'C' in the 10th position and a 
    #                                                                                                                               'G' in the 11th. Note that since this is a conditional 
    #                                                                                                                               probability, CG->*A, CG->*C, CG->*G, CG->*T will sum to 1.
    #                                       3   Second Read Mismatch            The first-order Markov model for mismatches         Same as above.
    #                                                                           between the reference and observed nucleotides      
    #                                                                           for the second read sequenced in a pair.            
    #                                       4   5' Sequence-Specific Bias       Parameters relating to the likelihood of the        This section is divided into 3 subsections. First is a 
    #                                                                           sequence surrounding the 5' end of a fragment in    matrix of the empirical nucleotide distribution for observed 
    #                                                                           transcript coordinates. See Roberts, et al.         fragments ("Observed Marginal Distribution") at each position 
    #                                                                           (2010a) for more details.                           in a window surrounding the 5' end of the fragment. The 
    #                                                                                                                               column headers give the 0-indexed position number with 
    #                                                                                                                               negatives being upstream in the target sequence. Each row 
    #                                                                                                                               gives the probability for a different nucleotide, which is 
    #                                                                                                                               specified in the first column followed by a colon. Note 
    #                                                                                                                               that since this is a probability distribution, each column 
    #                                                                                                                               will sum to 1. The second subsection contains the "Observed 
    #                                                                                                                               Conditional Probabilities". These are the conditional 
    #                                                                                                                               probabilities for the 3rd order Markov model, the columns 
    #                                                                                                                               specifying the conditional event in the observed fragments 
    #                                                                                                                               and the row specifying the window position. The third matrix 
    #                                                                                                                               is the "Expected Conditional Probabilities". This matrix is 
    #                                                                                                                               similar to the previous, except the probabilities are 
    #                                                                                                                               calculated assuming target sampling based only on fragment 
    #                                                                                                                               length and relative abundance, and fragment sampling within 
    #                                                                                                                               a target dependent only on length (no sequence biases). 
    #                                                                                                                               Bias weights in eXpress are calculated by taking the ratio 
    #                                                                                                                               of obesrved to expected probability.
    #                                       5   3' Sequence-Specific Bias       Parameters relating to the likelihood of the        Same as above, except for the 3' fragment end.
    #                                                                           sequence surrounding the 3' end of a fragment in    
    #                                                                           transcript coordinates. See Roberts, et al.         
    #                                                                           (2010a) for more details.                           
    #                                       -> If multiple alignment files were provided, a separate parameter output will be output for each with a unique index identifying its position 
    #                                          in the command-line argument given by the user (ie, the second SAM file in the argument list will be named 'params.2.xprs').
    # Count Variance-Covar (varcov.xprs)    -> This file is produced only when the --calc-covars option flag is used as described above. The file contains the estimated variances and 
    #                                          covariances on the counts between pairs of targets that shared multi-mapped reads, primarily to be used in differential expression analysis. 
    #                                          Since the covariance between targets in different bundles is always 0, the full sparse matrix is broken up into smaller tab-delimited matrices 
    #                                          for each bundle. An example of this output for the sample dataset used in the Getting Starting tutorial is shown below:
    #                                       -> Each bundle's matrix is headed by an identifier line that begins with a greater than symbol (>) followed by the bundle id and a 
    #                                          comma-separated list of targets in the bundle. The ordering of this list provides the indices for the matrix that is to follow. For 
    #                                          example, in bundle 1 of the output above, the fifth value in the second row (-2.862072e+02) is the covariance between NM_153633 and 
    #                                           NM_014620. Notice that an identical value is also in the second column of the fifth row, as the variance-covariance matrix will always be symmetric.
    #                                          >5: NM_014620, NM_153693, NR_003084, NM_153633, NM_018953, NM_004503
    #                                          2.067753e+02	-0.000000e+00	-0.000000e+00	-0.000000e+00	-0.000000e+00	-0.000000e+00
    #                                       



