#!/usr/bin/Rscript
# Script: .R
# Author: 
# 
# Description: Compare different gene counts produced the Eqp pipeline
#     
# Usage:
#     
# Source:
#     
# Reference:
# 
library('DESeq')
library(reshape2)
library(plyr)



###############################################################################
## returns string w/o leading or trailing whitespace
###############################################################################

trim <- function (x) gsub('^\\s+|\\s+$', '', x)





###############################################################################
## get upper or lower triangle of a matrix
###############################################################################

get.triangle <- function (input.matrix, tri.corner='upper') {
    
    if (tri.corner == 'upper') {
        # get upper
        input.matrix[lower.tri(input.matrix)] <- NA
    } else {
        # get lower
        input.matrix[upper.tri(input.matrix)] <- NA
    }
    return(input.matrix)
}



###############################################################################
## Capitalize words
###############################################################################

capitalize <- function(s, strict = FALSE) {
    cap <- function(s) paste(toupper(substring(s, 1, 1)),
            {s <- substring(s, 2); if(strict) tolower(s) else s},
            sep = '', collapse = ' ')
    sapply(strsplit(s, split = ' '), cap, USE.NAMES = !is.null(names(s)))
}




###############################################################################
## Read project counts files
###############################################################################

readProjectCounts <- function (dataset.dir, projects.frame, arg.id.type='gene', samples.names, dataset.name, gene.trans.frame, count.type='COUNT', methods=c(), count.list=list()) {
    
    if (length(count.list) == 0 || ! 'Aligners' %in% names(count.list)) {
        count.list <- list (`Aligners`=c())
    }
    
    if (count.type == 'COUNT') {
        count.file.extension <- 'cnt'
    } else if (count.type == 'FPKM') {
        count.file.extension <- 'fpkm'
    } else {
        message('Unknown count type: ', count.type)
        stopifnot (count.type == 'COUNT' || count.type == 'FPKM')
    }
    
    for (i in 1:nrow(projects.frame)) {
        
        # get data from projects-gene.txt
        quant.dir       <- file.path(dataset.dir)
        aligner         <- projects.frame[['Aligner']][i]         # e.g.: HISAT
        quantifier      <- projects.frame[['Quantifier']][i]      # e.g.: EQPQM
        setting         <- projects.frame[['Setting name']][i]    # e.g.: ambig
        combi.string    <- projects.frame[['Combinations']][i]    # e.g.: 1,2,3;2;3;...
        id.type         <- projects.frame[['Count type']][i]      # e.g.: gene
        combi.list      <- strsplit(combi.string, ';')[[1]]       # e.g.: '1,2,3' '2' '3'
        method          <- paste(paste(quantifier, setting, sep='-'), aligner, sep='/')    # e.g.: EQPQM/HISAT
        
        # set colClasses, comment.char and header
        if (quantifier %in% c('EQPQM','CUFFL','BTSEQ')) {
            column.setting = c('character', 'numeric')
            header.setting = FALSE
            comment.setting = ''
        } else if (quantifier == 'FEATC') {
            column.setting = c('character', rep('NULL', 5), 'numeric')
            header.setting = TRUE
            comment.setting = '#'
        } else if (quantifier == 'HTSEQ') {
            column.setting = c('character', 'numeric')
            header.setting = FALSE
            comment.setting = '_'
        } else if (quantifier == 'FLXCP') {
            column.setting = c(rep('NULL', 9), 'character', rep('NULL', 8) , 'character', rep('NULL', 4))
            header.setting = FALSE
            comment.setting = ''
        } else if (quantifier == 'RSEM') {
            column.setting = c('NULL', 'character', rep('NULL', 2), 'numeric', rep('NULL', 2))
            header.setting = TRUE
            comment.setting = ''
        } else if (quantifier == 'KLLST') {
            column.setting = c('character', rep('NULL', 2), 'numeric', 'NULL')
            header.setting = TRUE
            comment.setting = ''
        } else if (quantifier == 'FLXSM') {
            column.setting = c('NULL', 'character', rep('NULL', 7), 'numeric', 'NULL', 'NULL', 'NULL')
            header.setting = FALSE
            comment.setting = ''
        }
        
        
        if ((length(methods) == 0 || method %in% methods)) {# && (id.type == arg.id.type)) {
            # make data frame containing all quantifications
            quant.frame <- lapply(1:length(combi.list), function(j) {
                combi.align <- strsplit(combi.list[j], ',')[[1]]
                sample.name = ifelse(length(combi.align) > 1, 
                    paste(dataset.name, paste(combi.align, collapse='_'), sep='_'),
                    samples.names[strtoi(combi.align)])
                # create count file name
                count.file.name <- paste(aligner, paste(combi.align, collapse='_'), paste(setting, count.file.extension, sep='.'), sep='_')
                # add count file path to count file name
                count.file  <- file.path(quant.dir, paste(quantifier, 'quant', sep='_'), count.file.name)
                # get talbe from count.file
                input.frame <- read.table(count.file, header=header.setting, colClasses=column.setting, comment.char=comment.setting)
                if (quantifier == 'FLXCP') {
                    input.frame <- data.frame(input.frame[,1], as.numeric(unlist(strsplit(input.frame[,2],';'))))
                }
                # change table header information
                names(input.frame) <- c('Id', sample.name)
                # return table
                return(input.frame)
            })
            name.frame <- unlist(lapply(1:length(combi.list), function(j) {
                combi.align <- strsplit(combi.list[j], ',')[[1]]
                sample.name = ifelse(length(combi.align) > 1, 
                    paste(dataset.name, paste(combi.align, collapse='_'), sep='_'),
                    samples.names[strtoi(combi.align)])
                return(sample.name)
            }))
            # combine tables from quant.frame into one single table
            quant.table <- Reduce(function(x,y) merge(x, y, by = c('Id'), all=TRUE), quant.frame)
            # replace NA entries in quant.table
            quant.table[is.na(quant.table)] <- 0
            message(nrow(quant.table), ' rows read for method ... ', method)
            # change 'Id' to 'Gene Id'
            names(quant.table) <- gsub('Id', paste(capitalize(id.type), 'Id'), names(quant.table))
            # handle transcript data
            if(id.type == 'trans') {
                quant.table <- merge(quant.table, gene.trans.frame, by=('Trans Id'), all.y=TRUE)
                quant.table[is.na(quant.table)] <- 0
                quant.table <- quant.table[,c(ncol(quant.table),2:(ncol(quant.table)-1))]
            }
            quant.table <- quant.table[order(quant.table[['Gene Id']]),]
            quant.table <- aggregate(quant.table[-1],by=quant.table['Gene Id'],sum)
            # add quant.table into count.list
            count.list[[method]] <- quant.table
            # update 'Aligners' list inside count.list
            count.list[['Aligners']] <- unique(c(count.list[['Aligners']], aligner))
            count.list[['Quantifiers']] <- unique(c(count.list[['Quantifiers']], paste(quantifier, setting, sep='-')))
            count.list[['Samples']] <- unique(c(count.list[['Samples']], name.frame))
            count.list[['Settings']] <- unique(c(count.list[['Settings']], paste(paste(quantifier, setting, sep='-'), aligner, sep='/')))
        }
    }
    return (count.list)
}




###############################################################################
## Make biorad taqman data list
###############################################################################

makeBioTaqCntList <- function (mean.cnt.list, qpcr.method, out.dir, gene.id.list) {
    
    # make ne wlist
    biotaq.norm.cnt.list <- list()
    biotaq.file <- file.path(out.dir, paste('Mean_Gene_Count_Summary_for_SEQC_', qpcr.method, '.txt', sep=''))
    qpcr.string <- paste(qpcr.method, '-default/QPCR', sep='')
    biotaq.norm.cnt.list[[qpcr.string]] <- read.table(biotaq.file, header=TRUE, colClasses=c('factor', rep('numeric', 2)), sep='\t')
    colnames(biotaq.norm.cnt.list[[qpcr.string]]) <- c('Gene Id', 'SEQC-A-samples', 'SEQC-B-samples')
    biotaq.norm.cnt.list[[qpcr.string]] <- biotaq.norm.cnt.list[[qpcr.string]][order(biotaq.norm.cnt.list[[qpcr.string]][['Gene Id']]),]
    biotaq.norm.cnt.list[[qpcr.string]] <- aggregate(biotaq.norm.cnt.list[[qpcr.string]][-1],by=biotaq.norm.cnt.list[[qpcr.string]]['Gene Id'],sum)
    # filter by gene ids (obtained from gtf file)
    gene.id.list.frame <- data.frame(sort(intersect(gene.id.list, biotaq.norm.cnt.list[[qpcr.string]][['Gene Id']])))
    colnames(gene.id.list.frame) <- 'Gene Id'
    biotaq.norm.cnt.list[[qpcr.string]] <- merge(biotaq.norm.cnt.list[[qpcr.string]], gene.id.list.frame, by=('Gene Id'), all.y=TRUE)
    biotaq.norm.cnt.list[[qpcr.string]][is.na(biotaq.norm.cnt.list[[qpcr.string]])] <- 0
    biotaq.sum.A <- 0.0
    biotaq.sum.B <- 0.0
    for (method in names(mean.cnt.list)) {
        message('Status: ', 'Adding ', method, ' to ', qpcr.method, ' list ... ')
        biotaq.norm.cnt.list[[method]] <- data.frame(biotaq.norm.cnt.list[[qpcr.string]][,1])
        colnames(biotaq.norm.cnt.list[[method]]) <- 'Gene Id'
        biotaq.norm.cnt.list[[method]] <- merge(biotaq.norm.cnt.list[[method]], mean.cnt.list[[method]], by=('Gene Id'), all.x=TRUE)
        biotaq.norm.cnt.list[[method]][is.na(biotaq.norm.cnt.list[[method]])] <- 0
        biotaq.sum.A <- biotaq.sum.A + sum(biotaq.norm.cnt.list[[method]][,2])
        biotaq.sum.B <- biotaq.sum.B + sum(biotaq.norm.cnt.list[[method]][,3])
    }
    biotaq.A <- sum(biotaq.norm.cnt.list[[qpcr.string]][,2])
    biotaq.B <- sum(biotaq.norm.cnt.list[[qpcr.string]][,3])
    biotaq.scale.A <- biotaq.sum.A / length(mean.cnt.list) / biotaq.A
    biotaq.scale.B <- biotaq.sum.B / length(mean.cnt.list) / biotaq.B
    biotaq.norm.cnt.list[[qpcr.string]][,2] <- biotaq.norm.cnt.list[[qpcr.string]][,2] * biotaq.scale.A
    biotaq.norm.cnt.list[[qpcr.string]][,3] <- biotaq.norm.cnt.list[[qpcr.string]][,3] * biotaq.scale.B
    return(biotaq.norm.cnt.list)
}



###############################################################################
## Compute summary
###############################################################################

computeSummary <- function (count.frame, sample.types, gene.subset=c()) {
    
    if (length(gene.subset) == 0) {
        genes.ind = 1:nrow(count.frame)
    } else {
        genes.ind <- match(gene.subset, count.frame[['Gene Id']], nomatch=0)
        genes.ind <- genes.ind [genes.ind != 0]
    }
    
    detected.genes.num     <- list()
    expression.genes.total <- list()
    expression.genes.mean  <- list()
    for (sample.type in sample.types) {
        sample.type.ind <- grep(sample.type, names(count.frame))
        if (length(sample.type.ind) > 0) {
            # get data as matrix -> e.g.: (107,12,4,117)
            sample.mat <- matrix(as.numeric(unlist(count.frame[genes.ind, sample.type.ind])), ncol=length(sample.type.ind))
            # calculate number of genes -> e.g.: (|x,x,x,x| = 4)
            detected.genes.num[[sample.type]]     <- apply(sample.mat, 2, function (x) { return(sum(! is.na(x))) } )
            # calculate total number of expressed genes -> e.g.: (107 + 12 + 4 + 117 = 240)
            expression.genes.total[[sample.type]] <- apply(sample.mat, 2, function (x) { return(sum(x, na.rm=TRUE)) } )
            # calculate mean of expressed genes -> e.g.: (240 / 4 = 60.00)
            expression.genes.mean[[sample.type]]  <- apply(sample.mat, 2, function (x) { return(sum(x, na.rm=TRUE) / sum(! is.na(x))) } )
        }
    }
    return(list(detected.genes.num=detected.genes.num, expression.genes.total=expression.genes.total, expression.genes.mean=expression.genes.mean))
}




###############################################################################
## Get real tool names
###############################################################################

getNames <- function (tool.list) {
    
    # transforms id names into real tool names
    tool.names <- lapply(tool.list, function(x) {
        prefix <- strsplit(x, '-')[[1]][1]
        if (prefix %in% c('EQPQM', 'FEATC', 'HTSEQ')) {suffix <- strsplit(x, '-')[[1]][2]}
        if (prefix == 'EQPQM') {return(paste('EQP-QM', suffix, sep='-'))}
        else if (prefix == 'FEATC') {return(paste('featureCounts', suffix, sep='-'))}
        else if (prefix == 'HTSEQ') {return(paste('HTSeq', suffix, sep='-'))}
        else if (prefix == 'STAR') {return('STAR')}
        else if (prefix == 'HISAT') {return('HISAT')}
        else if (prefix == 'TOPHAT') {return('TopHat2')}
        else if (prefix == 'CUFFL') {return('Cufflinks')}
        else if (prefix == 'FLXCP') {return('Flux Capacitor')}
        else if (prefix == 'BTSEQ') {return('BitSeq')}
        else if (prefix == 'RSEM') {return('RSEM')}
        else if (prefix == 'KLLST') {return('kallisto')}
        else if (prefix == 'FLXSM') {return('Flux Simulator')}
        else {return(x)}
        })
    return(tool.names)
}




###############################################################################
## Normalize counts
###############################################################################

normalize.counts <- function (count.frame, count.name, normalization.frame=NULL, gene.ids.frame=NULL, sample.types=c('SEQC-A', 'SEQC-B'),
      normalization.method='CPM', sample.condition=NULL, normalization.scale=1000) {
    
    message('Normalizing ', count.name, ' with method ', normalization.method, ' ...')
    # check if gene ids in the count.frame are the same as the gene.ids.frame
    if (! is.null(gene.ids.frame)) {
        input.gene.ids <- as.character(count.frame[[1]])
        gene.ids       <- as.character(gene.ids.frame[[1]])
        gene.ids.ind   <- match(input.gene.ids, gene.ids, nomatch=0)
    } else {
        gene.ids.ind <- 1:nrow(count.frame)
        gene.ids     <- as.character(count.frame[[1]])
    }
    
    # use count.frame as normalization.frame, if normalization.frame is empty
    if (is.null(normalization.frame)) {
        normalization.frame <- count.frame
        if (length(normalization.frame) != length(count.frame)) {
            message('Normalization frame has a different number of columns (', length(normalization.frame), ') than the count.frame (', length(count.frame), ').')
            message('Normalization frame:')
            print(head(normalization.frame))
            message('Count frame:')
            print(head(count.frame))
            stopifnot(length(normalization.frame) == length(count.frame))
        }
    }
    
    # make new counts matrix, get from normalization frame
    counts.mat <- array(0, dim=c(length(gene.ids), length(count.frame) - 1))
    counts.mat[gene.ids.ind, ] <- as.matrix(subset(count.frame, gene.ids.ind != 0)[,-1])
    normalization.mat <- as.matrix(normalization.frame[,-1])
    
    # apply DESeq mormalization method
    if (normalization.method == 'DESeq') {
        sample.names <- names(normalization.frame)[-1]
        # define sample type to sample name e.g.: SEQC-A-BC01-s_1 --(condition)-->   SEQC-A
        if (is.null(sample.condition)) {
            sample.condition <- sample.names
            for (sample.type in sample.types) {
                sample.condition[grep(sample.type, sample.names)] <- sample.type
            }
        }
        # add all conditions into one frame
        design.frame <- data.frame(row.names = sample.names, condition = sample.condition)
        # create count data set out of 
        normalization.cds <- newCountDataSet(round(normalization.mat), design.frame[['condition']])
        
        ## Normalize counts
        normalization.cds <- estimateSizeFactors(normalization.cds)
        # scaling factor for each sample
        scaling.factors   <- sizeFactors(normalization.cds)
        #print(scaling.factors)
        
        if (sum(is.na(scaling.factors)) > 0) {
            message('Size factors: ', paste(sizeFactors (normalization.cds), collapse=', '))
            normalization.method <- 'CPM'
        } else {
            # make cpm factor from sum of all genes from each sample divided by scale factor
            cpm.factor <- mean(colSums(normalization.mat) / normalization.scale)
            #print(cpm.factor)
            scaling.factor.mat <- matrix(rep(scaling.factors * cpm.factor, nrow(counts.mat)), ncol=ncol(counts.mat), byrow=TRUE)
            #print(scaling.factor.mat[1:2,1:12])
        }
    }
    
    # only use column sums for normalisation
    if (normalization.method == 'CPM') {
        sum.counts <- colSums(normalization.mat) / normalization.scale
        scaling.factor.mat <- matrix(rep(sum.counts, nrow(counts.mat)), ncol=ncol(counts.mat), byrow=TRUE)
    } else if (normalization.method != 'DESeq') {
        message('Unknown normalization type: ', normalization.method)
        stopifnot(normalization.method == 'DESeq' || normalization.method == 'CPM')
    }
    
    normalized.frame <- data.frame(gene.ids, counts.mat / scaling.factor.mat)
    #normalized.frame <- data.frame(gene.ids, counts(normalization.cds) / scaling.factors)
    #normalized.frame <- data.frame(gene.ids, counts.mat / scaling.factors)
    names(normalized.frame) <- names(count.frame)
    #message('##### count.frame #####')
    #print(count.frame)
    #message('##### normalized.frame #####')
    #print(normalized.frame)
    return (normalized.frame)
}




###############################################################################
## Denormalize by length (FPKM or TPM)
###############################################################################

length.denormalize <- function (input.frame, gene.length.frame) {

  common.gene.ids    <- intersect(input.frame[[1]], gene.length.frame[[1]])
  input.ind          <- match(common.gene.ids, input.frame[[1]], nomatch=0)
  input.filter.frame <- input.frame[input.ind, ]

  gen.len.ind              <- match(common.gene.ids, gene.length.frame[[1]], nomatch=0)
  gene.length.filter.frame <- gene.length.frame[gen.len.ind, ]

  result.frame <- data.frame(lapply(input.filter.frame,
                                    function (x, y) {
                                      if (! is.numeric(x[1])) return (x)
                                      return (1.0 * x * y)
                                    }, gene.length.filter.frame[["Length"]]/1000),
                             check.names=FALSE)
  
  return (result.frame)
}




###############################################################################
## Correct results
## -> compare stdev of gene against mean stdev of gene
##      - assume Poisson -> Var(x) = E(x))
## -> don't use upper and lower most samples for stdv and mean
##      - k = max(|samples| mod 8, 3)
##      - don't use: samples 1 -> k and n-k+1 -> |samples|
## -> for each gene and correct the genes for which there is a large deviation.
###############################################################################

correctResults <- function (input.frame, input.name, sample.types=c('SEQC-A', 'SEQC-B'), correct.pseudo_count=1) {
    
    message('Correcting ', input.name, ' ...')
    result.frame <- input.frame
    genes.corrected <- list(sample.types)
    
    # do for SEQC-A and SEQC-B
    for (sample.type in sample.types) {
        
        # get sample type name position in the result.frame
        sample.type.ind <- grep (sample.type, names(result.frame))
        
        ## check if there are more than 12 samples
        ## ??? n <- 13: [1]  -1  -2  -3 -11 -12 -13
        ## ??? n <- 12: [1]  1  2  3  4  5  6  7  8  9 10 11 12
        ## ??? n <- 32: [1]  -1  -2  -3  -4 -29 -30 -31 -32
        sample.cols.ind <- 1:length(sample.type.ind)
        if (length(sample.type.ind) > 12)  {
            n <- length(sample.type.ind)
            k <- max(n %/% 8, 3)
            sample.cols.ind <- c(-1 * 1:k, -1 * ((n-k+1):n))
        }
        
        ## Compute a trimmed mean on the middle 3/4 of the data (top and lower part of data)
        ## over rows(1) of matrix of all sample.type sampes in sample.type.ind
        ## compute mean (do for each gene (row)) e.g.: [1] 39.698644  7.507281  1.618229 48.081564
        gene_count.mean <- apply(as.matrix(result.frame[,sample.type.ind]), 1, function (x) return (mean(sort(x)[sample.cols.ind])))
        ## compute standard deviation (do for each gene (row)) e.g.: [1] 2.3950863 1.5765115 0.9594312 4.5822255
        gene_count.sd   <- apply(as.matrix(result.frame[,sample.type.ind]), 1, function (x) return (sd(sort(x)[sample.cols.ind])))
        
        ## check if genes (rows) need to be corrected
        ## We assume most Poisson distributed data so the variance equals the mean
        ## for all rows (genes)
        sd.factor <- 20
        gene_count.excl.bool <- apply(cbind(gene_count.mean, gene_count.sd, result.frame[,sample.type.ind]), 1, function (x) {
            ## row x :  mean       sd          seqA-1     seqA-2      seqA-3     seqA-4     seqA-5     seqA-6
            ## x e.g.:  39.698644   2.3950863  42.276122  38.5777587  38.416048  40.799690  41.986205  36.136043
            if (length(x) < 4) return (FALSE) # check if there is only one sample (e.g.:     mean    sd  seqA-1)
            mean.x <- as.numeric(x[1]) # get mean
            stdv.x <- max(as.numeric(x[2]), sqrt(mean.x + correct.pseudo_count)) # get stdev
            ## calculate difference threshold e.g.: max(39.698644 / 2, 20 * 2.3950863) -> max(19.849322, 47.901726)
            diff.thresh <- max(mean.x/2, sd.factor * stdv.x)
            x.sort <- as.numeric(sort(x[c(-1, -2)], decreasing=TRUE)) # sort gene entries -> descending 
            ## check if ! ( mean(genes) - diff.thresh < min(genes)  and  max(genes) < mean(genes) + diff.thresh )
            ## e.g.:    ! ( 39.698644   - 47.901726   < 36.136043   ...  42.276122  < 39.698644   + 47.901726   )
            ## e.g.:    ! (               -8.203082   < 36.136043   ...  42.276122  < 87.60037                  )
            ## -> check if genes are not in between mean(genes) +/- threshold
            ## -> return TRUE if one or both are outside
            #message(' ##### Name: ', input.name, ' ##### ')
            #message(' # ', mean.x, ' - ', diff.thresh, ' < ', x.sort[length(x.sort)])
            #message(' && ', x.sort[1], ' < ', mean.x, ' + ', diff.thresh)
            return (!(mean.x - diff.thresh < x.sort[length(x.sort)] && x.sort[1] < mean.x + diff.thresh))
        })
        
        # correcting rows (genes) for which genes are not in between mean(genes) +/- threshold
        if (sum(gene_count.excl.bool) > 0) {
            message('Correcting ', sum(gene_count.excl.bool), ' rows of gene counts for sample type ', sample.type, ' of ', input.name)
            if (sum(gene_count.excl.bool) <= 20) {
                message('Gene(s) with corrected counts: ', paste(result.frame[[1]][gene_count.excl.bool], collapse=', '))
            }
            # access all data entire rows(genes) which are marked as TRUE and correct results
            result.frame[gene_count.excl.bool, sample.type.ind] <-
            matrix(apply(cbind(gene_count.mean[gene_count.excl.bool], gene_count.sd[gene_count.excl.bool], result.frame[gene_count.excl.bool, sample.type.ind]), 1, function (x) {
                x.result    <- as.numeric(x[c(-1, -2)])                # data in row (gene) x
                x.leng      <- length(x.result)                        # length of row x
                x.order     <- order(x.result, decreasing=TRUE)        # order of data in x (descending) (e.g.: 100,50,70 -> 1,3,2)
                x.sorted    <- x.result[x.order]                       # sorted data of x according to x.order
                mean.x      <- as.numeric(x[1])                        # get mean from x
                stdv.x      <- max(as.numeric(x[2]), sqrt(mean.x + correct.pseudo_count)) # get max between stdev and sqrt(mean+1) from x
                diff.thresh <- max(mean.x/2, sd.factor * stdv.x)       # calculate diff.threshold max(man/2, stdv*factor)
                over.thresh <- sum(!(x.result < mean.x + diff.thresh)) # number of data entries which are not smaller than mean + diff.threshold
                undr.thresh <- sum(!(x.result > mean.x - diff.thresh)) # number of data entries which are not greater than mean - diff.threshold
                if (over.thresh > 0) {                                 # access entries found in over.thresh and shift the data one entrie in decreasing direction (removes highest entry)
                    x.result[x.order][1:over.thresh] <- x.sorted[over.thresh+1]
                }
                if (undr.thresh > 0) {                                 # access entries found in undr.thresh and shfit the data one entrie in increasing direction (removes lowest entry)
                    x.result[x.order][(x.leng-undr.thresh+1):x.leng] <- x.sorted[x.leng-undr.thresh]
                }
                return(x.result)
             }), ncol=length(sample.type.ind), byrow=TRUE)
        }
        genes.corrected[[sample.type]] <- result.frame[[1]][gene_count.excl.bool]
        # print out corrected genes
        if (sum(gene_count.excl.bool) > 0) {
            message('Method ', input.name, ' corrected genes ... ')
            print(genes.corrected[[sample.type]])
        }
    }
    return (result.frame)
}




###############################################################################
# compute CV, remove 'NA' values
###############################################################################
coefficientOfVariation <- function(x, na.rm=TRUE) {
    
    if (sum(! is.na(x)) > 1) {
        return (sd(x, na.rm=na.rm) / mean(x, na.rm=na.rm))
    } else {
        return (NA)
    }
}




###############################################################################
## Compute the coefficient of variation for each gene
## relative standard deviation -> stdev/mean
###############################################################################

computeCVs <- function (input.frame, input.name, sample.types) {
    
    message('Coefficient of variation for ', input.name, ' ...')
    
    # create CV matrix
    gene.CVs.mat <- array(NA, dim=c(nrow(input.frame), length(sample.types)))
    # calculete coefficient of variation for each sample type
    for (i in 1:length(sample.types)) {
        sample.type      <- sample.types[i]
        # get samples names from sample.type
        sample.type.ind  <- grep(sample.type, names(input.frame))
        # get samples data from samples names
        sample.type.mat  <- as.matrix(input.frame[,sample.type.ind])
        gene.expr.bool   <- rowSums(sample.type.mat) > 0
        # filter empty genes and only get non 0 rows
        gene.CVs.mat[gene.expr.bool,i] <- apply(sample.type.mat[gene.expr.bool, ], 1, coefficientOfVariation)
    }
    
    # calculate CV mean of all sample types for each gene
    CV.list.mean <- apply(gene.CVs.mat, 1, mean, na.rm=TRUE)
    CV.list.A <- gene.CVs.mat[,1]
    CV.list.B <- gene.CVs.mat[,2]
    names(CV.list.mean) <- input.frame[[1]]
    names(CV.list.A) <- input.frame[[1]]
    names(CV.list.B) <- input.frame[[1]]
    
    return(list(cv.mean.list=CV.list.mean, cv.A.list=CV.list.A, cv.B.list=CV.list.B))
}




###############################################################################
## Compute Fold Changes
###############################################################################

computeFCs <- function (input.frame, input.name, samples.matched.A, samples.matched.B, pseudo_count=0.5) {
    
    message('Fold change ', input.name, ' ... ')
    
    # create FC matrix
    fc.mat <- array(NA, dim=c(nrow(input.frame), length(samples.matched.A)))
    colnames(fc.mat) <- paste(samples.matched.A, samples.matched.B, sep='/')
    rownames(fc.mat) <- input.frame[['Gene Id']]
    # calculate fold change for each sample A - sample B combination e.g.: fcCV for SEQC-A-BC01-s_1 and SEQC-B-BC02-s_1 for each gene
    for (i in 1:min(length(samples.matched.A), length(samples.matched.B))) {
        sample.A <- samples.matched.A[i]
        sample.B <- samples.matched.B[i]
        # compute fold change vector for each gene and save into column i
        fold.change.frame <- computeFoldChangeVector(input.frame, sample.A, sample.B, pseudo_count)
        fc.mat[fold.change.frame[['Gene Id']],i] <- fold.change.frame[['fc']]
    }
    
    return(fc.mat)
}




###############################################################################
## Compute the coefficient of variation for FC (Fold Change) of each gene
###############################################################################

computeFcCVs <- function (fc.mat, input.name, log.fc.threshold) {
    
    message('Fold change coefficient of variation for ', input.name, ' ... ')
    
    # log2.thresh.func
    log2.thresh.func <- function(x) {
        if (abs(log2(x)) < log.fc.threshold) {
            if (log2(x) >= 0) {return(log.fc.threshold)}
            else {return(-log.fc.threshold)}
        } else {
            return(log2(x))
        }
    }
    
    # apply function
    log2.fc.mat <- matrix(mapply(log2.thresh.func, fc.mat), nrow=nrow(fc.mat), ncol=ncol(fc.mat))
    colnames(log2.fc.mat) <- colnames(fc.mat)
    rownames(log2.fc.mat) <- rownames(fc.mat)
    
    # calculate CV of the log2(FC matrix) for each gene
    fcCV.fc.norm <- apply(fc.mat, 1, coefficientOfVariation, na.rm=TRUE)
    fcCV.fc.log <- apply(log2.fc.mat, 1, coefficientOfVariation, na.rm=TRUE)
    names(fcCV.fc.norm) <- rownames(fc.mat)
    names(fcCV.fc.log) <- rownames(fc.mat)
    
    return(list(intra.gene.CV.FC=fcCV.fc.norm, intra.gene.CV.log.FC=fcCV.fc.log))
}





###############################################################################
## Compute fold change vector
###############################################################################

computeFoldChangeVector <- function (input.frame, sample.A, sample.B, pseudo_count=0.5) {
    
    # prepare A and B frame
    prepAB <- function (sample.X) {
        # get gene name and sample data
        sample.X.ind <- grep(sample.X, names(input.frame))
        result.frame.X <- input.frame[,c(1, sample.X.ind)]
        names(result.frame.X) <- c('Gene Id', 'Count')
        # set NA values to 0
        counts.X <- as.numeric(result.frame.X[['Count']])
        counts.X[is.na(counts.X)] <- 0
        # combine non-NA count vector with gene name vector
        result.frame.X <- data.frame(result.frame.X[['Gene Id']], counts.X)
        names(result.frame.X) <- c('Gene Id', 'Count')
        return(result.frame.X)
    }
    
    # prepare sample A
    result.frame.A <- prepAB(sample.A)
    # prepare sample B
    result.frame.B <- prepAB(sample.B)
    # get matching A and B gene positions (set to 0 if no match was found)
    result.A.B.ind <- match(result.frame.A[['Gene Id']], result.frame.B[['Gene Id']], nomatch=0)

    # check if there were any matching genes
    if (sum(result.A.B.ind != 0) == 0) {
        print(paste('Warning no matching gene ids found for', sample.A, 'and', sample.B))
        return (list(`Gene Id`=c(), fc=c()))
    }
    
    # get matching gene names
    result.id <- as.character(result.frame.A[['Gene Id']])[result.A.B.ind != 0]
    # get matching gene data
    result.A <- as.numeric(result.frame.A[['Count']])[result.A.B.ind != 0]
    result.B <- as.numeric(result.frame.B[['Count']])[result.A.B.ind != 0]
    # get A-B fold change vector
    result.fc <- (result.B + pseudo_count) / (result.A + pseudo_count)

    # make A-B fold change frame
    result.frame <- data.frame(result.id, result.fc)
    names(result.frame) <- c('Gene Id', 'fc')
    
    return(result.frame)
}




###############################################################################
## Compute correlation within methods (reproducibility)
###############################################################################

computeIntraMethodCorrelations <- function (input.frame, input.name, samples.list, sample.types, id.name='Gene Id', num.samples.thresh=-1) {
    
    message('Computing intra method count correlations for ', input.name, ' ... ')
    
    # initialize lists
    samples.seqc <- list()
    pearson.cor  <- list()
    CD.cor <- list()
    spearman.cor <- list()
    
    #separate A and B samples
    for (sample.type in sample.types) {
        samples.seqc[[sample.type]] <- samples.list[grep(sample.type, samples.list)]
        if (length(samples.seqc[[sample.type]]) > 0) {
            pearson.cor[[sample.type]]  <- list()
            CD.cor[[sample.type]] <- list()
            spearman.cor[[sample.type]] <- list()
        }
    }
    
    # set sample.thresh to length of samples.list if num.samples.thresh was specified as <= 0
    if (num.samples.thresh <= 0) {
        num.samples.thresh <- length(samples.list)
    }
    
    # for A and B samples
    for (sample.type in sample.types) {
        # get number of samples
        sample.type.num <- length(samples.seqc[[sample.type]])
        if (sample.type.num > 0) {
            # for each sample compute correlation to other samples
            
            # first sample
            for (s1 in 1:min(sample.type.num, num.samples.thresh)) {
                sample.1 <- samples.seqc[[sample.type]][s1] # get sample name from first sample
                sample.1.ind <- grep(sample.1, names(input.frame)) # get sample position
                result.1.frame <- input.frame[!is.na(input.frame[[sample.1.ind]]), c(1, sample.1.ind)] # get sample data
                names(result.1.frame) <- c(id.name, 'Count')
                
                # second sample
                for (s2 in 1:min(sample.type.num, num.samples.thresh)) {
                    sample.2 <- samples.seqc[[sample.type]][s2] # get sample name from second sample
                    sample.2.ind <- grep(sample.2, names(input.frame)) # get sample position
                    result.2.frame <- input.frame[!is.na(input.frame[[sample.2.ind]]), c(1, sample.2.ind)] # get sample data
                    names(result.2.frame) <- c(id.name, 'Count')
                    
                    matched.id.ind <- match(result.1.frame[[id.name]], result.2.frame[[id.name]], nomatch=0) # get matching positions e.g.: 1,2,3,4
                    
                    # calculate pearson over genes for samples
                    pearson.cor.sample <- cor(result.1.frame[['Count']][matched.id.ind != 0], result.2.frame[['Count']][matched.id.ind], use='pairwise.complete.obs')
                    pearson.cor[[sample.type]][[sample.1]] <- c(pearson.cor[[sample.type]][[sample.1]], pearson.cor.sample)
                    names(pearson.cor[[sample.type]][[sample.1]])[length(pearson.cor[[sample.type]][[sample.1]])] <- sample.2
                    
                    # calculate CD over genes for samples
                    CD.cor.sample <- cor(result.1.frame[['Count']][matched.id.ind != 0], result.2.frame[['Count']][matched.id.ind], use='pairwise.complete.obs')^2
                    CD.cor[[sample.type]][[sample.1]] <- c(CD.cor[[sample.type]][[sample.1]], CD.cor.sample)
                    names(CD.cor[[sample.type]][[sample.1]])[length(CD.cor[[sample.type]][[sample.1]])] <- sample.2
                    
                    # calculate spearman over genes for samples
                    spearman.cor.sample <- cor(result.1.frame[['Count']][matched.id.ind != 0], result.2.frame[['Count']][matched.id.ind], method='spearman', use='pairwise.complete.obs')
                    spearman.cor[[sample.type]][[sample.1]] <- c(spearman.cor[[sample.type]][[sample.1]], spearman.cor.sample)
                    names(spearman.cor[[sample.type]][[sample.1]])[length(spearman.cor[[sample.type]][[sample.1]])] <- sample.2
                }
            }
        }
    }
    return (list(intra.met.pearson=pearson.cor, intra.met.spearman=spearman.cor, intra.met.CD=CD.cor))
}




###############################################################################
## Compute fc-correlation within methods (reproducibility)
###############################################################################

computeIntraMethodFcCorrelations <- function (input.frame, input.name, samples.matched.A, samples.matched.B, pseudo_count=0.5, id.name='Gene Id') {
    
    message('Computing intra method (log) fold change correlations for ', input.name, ' ... ')
    
    # removes leading and trailing whitespace characters
    samples <- trim(gsub(paste(' ', '', sep=''), '', names(input.frame)[2:length(input.frame)]))
    
    # initialize lists
    pearson.fc.cor <- list()
    pearson.log.fc.cor <- list()
    spearman.fc.cor <- list()
    spearman.log.fc.cor <- list()
    CD.fc.cor <- list()
    CD.log.fc.cor <- list()
    # set number of correlations
    sample.num <- min(length(samples.matched.A), length(samples.matched.B))
    # get first samples
    for (i in 1:sample.num) {
        # get first A and B sample
        sample.A.1 <- samples.matched.A[i]
        sample.B.1 <- samples.matched.B[i]
        # set first correlation name
        sample.cor.name.1 <- paste(sample.A.1, sample.B.1, sep='/')
        # initialize vectors
        pearson.fc.cor[[sample.cor.name.1]] <- c()
        pearson.log.fc.cor[[sample.cor.name.1]] <- c()
        spearman.fc.cor[[sample.cor.name.1]] <- c()
        spearman.log.fc.cor[[sample.cor.name.1]] <- c()
        CD.fc.cor[[sample.cor.name.1]] <- c()
        CD.log.fc.cor[[sample.cor.name.1]] <- c()
        # initialize name vector
        combined.sample.names <- c()
        
        # compute fold change vector of first A and B sample
        result.1.fc.list <- computeFoldChangeVector(input.frame, sample.A.1, sample.B.1, pseudo_count)
        
        # get second samples
        for (j in 1:sample.num) {
            # get second A and B samples
            sample.A.2 <- samples.matched.A[j]
            sample.B.2 <- samples.matched.B[j]
            # set second correlation name and append to name vector
            sample.cor.name.2 <- paste(sample.A.2, sample.B.2, sep='/')
            combined.sample.names <- append(combined.sample.names, sample.cor.name.2)
            
            # compute fold change vector of second A and B sample
            result.2.fc.list <- computeFoldChangeVector(input.frame, sample.A.2, sample.B.2, pseudo_count)
            
            # get matching gene id positions
            result.id.ind <- match(as.character(result.2.fc.list[['Gene Id']]), as.character(result.1.fc.list[['Gene Id']]), nomatch=0)
            
            for (method in c('pearson', 'spearman')) {
                # compute log fold change correlation
                sample.log.cor <- cor( log2(result.2.fc.list[['fc']][result.id.ind != 0]), 
                                       log2(result.1.fc.list[['fc']][result.id.ind]),
                                       method=method, 
                                       use='pairwise.complete.obs')
                # compute fold change correlation
                sample.cor <- cor( result.2.fc.list[['fc']][result.id.ind != 0], 
                                   result.1.fc.list[['fc']][result.id.ind], 
                                   method=method, 
                                   use='pairwise.complete.obs')
                
                # save results
                if (method == 'pearson') {
                    pearson.fc.cor[[sample.cor.name.1]] <- c(pearson.fc.cor[[sample.cor.name.1]], sample.cor)
                    pearson.log.fc.cor[[sample.cor.name.1]] <- c(pearson.log.fc.cor[[sample.cor.name.1]], sample.log.cor)
                    CD.fc.cor[[sample.cor.name.1]] <- c(CD.fc.cor[[sample.cor.name.1]], sample.cor^2)
                    CD.log.fc.cor[[sample.cor.name.1]] <- c(CD.log.fc.cor[[sample.cor.name.1]], sample.log.cor^2)
                } else {
                    spearman.fc.cor[[sample.cor.name.1]] <- c(spearman.fc.cor[[sample.cor.name.1]], sample.cor)
                    spearman.log.fc.cor[[sample.cor.name.1]] <- c(spearman.log.fc.cor[[sample.cor.name.1]], sample.log.cor)
                }
            }
        }
        # set names
        names(pearson.fc.cor[[sample.cor.name.1]]) <- combined.sample.names
        names(pearson.log.fc.cor[[sample.cor.name.1]]) <- combined.sample.names
        names(spearman.fc.cor[[sample.cor.name.1]]) <- combined.sample.names
        names(spearman.log.fc.cor[[sample.cor.name.1]]) <- combined.sample.names
        names(CD.fc.cor[[sample.cor.name.1]]) <- combined.sample.names
        names(CD.log.fc.cor[[sample.cor.name.1]]) <- combined.sample.names
    }
    return (list(intra.met.pearson.fc=pearson.fc.cor,         intra.met.spearman.fc=spearman.fc.cor,         intra.met.CD.fc=CD.fc.cor,
                 intra.met.pearson.log.fc=pearson.log.fc.cor, intra.met.spearman.log.fc=spearman.log.fc.cor, intra.met.CD.log.fc=CD.log.fc.cor))
}




###############################################################################
## Compute correlation of fold-change between methods
###############################################################################

computeFoldChangeCorrelations <- function (result1.frame, result2.frame, method.name.1, method.name.2, 
      samples.matched.A, samples.matched.B, pseudo_count.1=0.5, pseudo_count.2=0.5, 
      permute.samples=FALSE, filter.genes=NULL) {
    
    message('Computing inter method (log) fold change correlations between ', method.name.2, ' and ', method.name.1, ' ... ')
    # <removed filterMatchedSamples> -> get intersection of result1.frame and result2.frame samples
    
    # initialize parameters
    pearson.fc.cor  <- c()
    pearson.log.fc.cor <- c()
    spearman.fc.cor <- c()
    spearman.log.fc.cor <- c()
    CD.fc.cor <- c()
    CD.log.fc.cor <- c()
    false.negative.rate <- c()
    false.positive.rate <- c()
    sensitivity <- c() # TPR
    specificity <- c() # TNR
    false.discovery.rate <- c()
    precision <- c() # PPV
    
    # get common genes
    common.genes <- intersect(result1.frame[[1]], result2.frame[[1]])
    # filter genes according to 'filter.genes' parameter
    if (! is.null(filter.genes)) {
        common.genes <- intersect(common.genes, filter.genes)
    }
    #stopifnot(length(common.genes) >= 100)
    # get gene positions and corresponding filtered result1.frame and result2.frame
    common.genes.ind.1 <- match(common.genes, result1.frame[[1]], nomatch=0)
    result1.filter.frame <- result1.frame[common.genes.ind.1, ]
    common.genes.ind.2 <- match(common.genes, result2.frame[[1]], nomatch=0)
    result2.filter.frame <- result2.frame[common.genes.ind.2, ]
    
    # make calculations
    leng.A <- length(samples.matched.A)
    leng.B <- length(samples.matched.B)
    for (i in 1:min(leng.A, leng.B)) {
        # get A.1, A.2, B.1, B.2 samples
        sample.A.1 <- samples.matched.A[i]
        sample.B.1 <- samples.matched.B[i]
        if (permute.samples) {
            sample.A.2 <- samples.matched.A[i %% leng.A + 1]
            sample.B.2 <- samples.matched.B[i %% leng.B + 1]
        } else {
            sample.A.2 <- sample.A.1
            sample.B.2 <- sample.B.1
        }
        
        # compute fold change for A.1 and B.1 in result1.frame (intra method 1 FC)
        result1.fc.list <- computeFoldChangeVector(result1.filter.frame, sample.A.1, sample.B.1, pseudo_count.1)
        names(result1.fc.list[['fc']]) <- result1.fc.list[['Gene Id']]
        # compute fold change for A.2 and B.2 in result2.frame (intra method 2 FC)
        result2.fc.list <- computeFoldChangeVector(result2.filter.frame, sample.A.2, sample.B.2, pseudo_count.2)
        names(result2.fc.list[['fc']]) <- result2.fc.list[['Gene Id']]
        # make correlations for (log) pearson, spearman and CD
        for (method in c('pearson', 'spearman')) {
            # compute log fold change correlation
            sample.log.cor <- cor(log2(result2.fc.list[['fc']]), log2(result1.fc.list[['fc']]), method=method, use='pairwise.complete.obs')
            # compute fold change correlation
            sample.cor <- cor(result2.fc.list[['fc']], result1.fc.list[['fc']], method=method, use='pairwise.complete.obs')
            # save results
            if (method == 'pearson') {
                pearson.fc.cor <- c(pearson.fc.cor, sample.cor)
                pearson.log.fc.cor <- c(pearson.log.fc.cor, sample.log.cor)
                CD.fc.cor <- c(CD.fc.cor, sample.cor^2)
                CD.log.fc.cor <- c(CD.log.fc.cor, sample.log.cor^2)
            } else {
                spearman.fc.cor <- c(spearman.fc.cor, sample.cor)
                spearman.log.fc.cor <- c(spearman.log.fc.cor, sample.log.cor)
            }
        }
        # set names
        fc.pair.name <- paste(sample.A.1, sample.B.1, sep='/')
        if (permute.samples) {
            fc.pair.name <- paste(paste(sample.A.1, sample.B.1, sep='-'), paste(sample.A.2, sample.B.2, sep='-'), sep='/')
        }
        names(pearson.fc.cor)[length(pearson.fc.cor)] <- fc.pair.name
        names(pearson.log.fc.cor)[length(pearson.log.fc.cor)] <- fc.pair.name
        names(spearman.fc.cor)[length(spearman.fc.cor)] <- fc.pair.name
        names(spearman.log.fc.cor)[length(spearman.log.fc.cor)] <- fc.pair.name
        names(CD.fc.cor)[length(CD.fc.cor)] <- fc.pair.name
        names(CD.log.fc.cor)[length(CD.log.fc.cor)] <- fc.pair.name
        
        
        ## calculate sensitivity and specificity
        # genes whose (FC is not NA) and (|log2(FC)| >= 1)
        # -> 2-fold expression increase -> log2(FC) = +1 (over  expressed gene)
        # -> 2-fold expression decrease -> log2(FC) = -1 (under expressed gene)
        # -> ABS gets over and under expressed genes
        # get names of over/under expressed genes in result1.fc.list (intra method 1 FC)
        result1.sig.id     <- result1.fc.list[['Gene Id']][! is.na (result1.fc.list[['fc']]) & abs(log2(result1.fc.list[['fc']])) >= 1]
        # get names of over/under expressed genes in result2.fc.list (intra method 2 FC)
        result2.sig.id     <- result2.fc.list[['Gene Id']][! is.na (result2.fc.list[['fc']]) & abs(log2(result2.fc.list[['fc']])) >= 1]
        
        ##                   | null hyp. true | null hyp. false |
        ## ------------------+----------------+-----------------+-----------
        ## decl. significant |     V (FP)     |      S (TP)     |    R
        ## ------------------+----------------+-----------------+-----------
        ## not decl. sign.   |     U (TN)     |      T (FN)     |   m-R
        ## ------------------+----------------+-----------------+-----------
        ##                   |     m_0 (P)    |    m-m_0 (N)    |  m (all)
        
        # m     -> total number hypotheses tested
        # m_0   -> number of true null hypotheses
        # m-m_0 -> number of true alternative hypotheses
        # V     -> number of false positives (Type I error) (false discoveries)
        # S     -> number of true positives (true discoveries)
        # T     -> number of false negatives (Type II error)
        # U     -> number of true negatives
        # R     -> number of rejected null hypotheses (also called "discoveries")
        # In m hypothesis tests of which m_0 are true null hypotheses, R is an observable random variable, and S, T, U, and V are unobservable random variables.
        
        # sensitivity or true positive rate        -> TPR  =  TP / (TP + FN)  =  TP / P  =  1 - FNR
        # specificity or true negative rate        -> SPC  =  TN / (TN + FP)  =  TN / N  =  1 - FPR
        # fall-out or false positive rate          -> FPR  =  FP / (FP + TN)  =  FP / N  =  1 - SPC 
        # false negative rate                      -> FNR  =  FN / (TP + FN)  =  FN / P  =  1 - TPR
        # precision or positive predictive value   -> PPV  =  TP / (TP + FP)  =             1 - FDR
        # negative predictive value                -> NPV  =  TN / (TN + FN)
        # false discovery rate                     -> FDR  =  FP / (TP + FP)  =             1 - PPV
        
        # H_0    -> gene is over/under expressed
        # H_a    -> gene is NOT over/under expressed
        # m      -> common.genes                                 -> all.id
        all.id            <- common.genes
        # m_0    -> result2.sig.id                               -> all.positives.id
        all.positives.id  <- result2.sig.id
        # m-m_0  -> setdiff(all.id, all.positives.id)            -> all.negatives.id
        all.negatives.id  <- setdiff(all.id, all.positives.id)
        # V (FP) -> setdiff(result1.sig.id, all.positives.id)    -> false.positives.id
        false.positives.id <- setdiff(result1.sig.id, all.positives.id)
        # U (TN) -> setdiff(all.negatives.id, result1.sig.id)    -> true.negatives.id
        true.negatives.id <- setdiff(all.negatives.id, result1.sig.id)
        # S (TP) -> intersect(result1.sig.id, all.positives.id)  -> true.positives.id
        true.positives.id <- intersect(result1.sig.id, all.positives.id)
        # T (FN) -> setdiff(all.positives.id, result1.sig.id)    -> false.negatives.id
        false.negatives.id <- setdiff(all.positives.id, result1.sig.id)
        # R      -> union(result1.sig.id, all.positives.id)      -> dec.significant.id
        dec.significant.id <- union(result1.sig.id, all.positives.id)
        
        # FPR = FP / N
        false.positive.rate <- c(false.positive.rate, length(false.positives.id) / length(all.negatives.id))
        names(false.positive.rate)[length(false.positive.rate)] <- fc.pair.name
        # FNR = FN / P
        false.negative.rate <- c(false.negative.rate, length(false.negatives.id) / length(all.positives.id))
        names(false.negative.rate)[length(false.negative.rate)] <- fc.pair.name
        # sensitivity -> TPR = TP / P
        sensitivity <- c(sensitivity, length(true.positives.id) / length(all.positives.id))
        names(sensitivity)[length(sensitivity)] <- fc.pair.name
        # specificity -> TNR = TN / N
        specificity <- c(specificity, length(true.negatives.id) / length(all.negatives.id))
        names(specificity)[length(specificity)] <- fc.pair.name
        # FDR = FP / (TP + FP)
        false.discovery.rate <- c(false.discovery.rate, length(false.positives.id) / length(dec.significant.id))
        names(false.discovery.rate)[length(false.discovery.rate)] <- fc.pair.name
        # precision -> PPV = TP / (TP + FP)
        precision <- c(precision, length(true.positives.id) / length(dec.significant.id))
        names(precision)[length(precision)] <- fc.pair.name
    }
    # return values
    return (list( inter.met.pearson.fc           = pearson.fc.cor, 
                  inter.met.pearson.log.fc       = pearson.log.fc.cor, 
                  inter.met.spearman.fc          = spearman.fc.cor, 
                  inter.met.spearman.log.fc      = spearman.log.fc.cor, 
                  inter.met.CD.fc                = CD.fc.cor, 
                  inter.met.CD.log.fc            = CD.log.fc.cor,                   
                  inter.met.false.negative.rate  = false.negative.rate, 
                  inter.met.sensitivity          = sensitivity,
                  inter.met.false.positive.rate  = false.positive.rate, 
                  inter.met.specificity          = specificity, 
                  inter.met.false.discovery.rate = false.discovery.rate, 
                  inter.met.precision            = precision             ))
}




###############################################################################
## Compute correlation between methods
###############################################################################

computeInterMethodCorrelations <- function (result1.frame, result2.frame, method.name.1, method.name.2, gene.length.frame=NA, 
                                            pseudo_count=0.5, normalize.counts=FALSE, permute.samples=FALSE, length.normalization.type="FPKM") {
    
    message('Computing inter method (log) gene expression correlations between ', method.name.2, ' and ', method.name.1, ' ... ')
    
    # get sample names
    samples1 <- names(result1.frame)[-1]
    samples2 <- names(result2.frame)[-1]
    samples <- intersect(samples1, samples2)
    
    # initiate parameters
    pearson.cor      <- array(0, dim=length(samples))
    spearman.cor     <- array(0, dim=length(samples))
    CD.cor           <- array(0, dim=length(samples))
    pearson.log.cor  <- array(0, dim=length(samples))
    spearman.log.cor <- array(0, dim=length(samples))
    CD.log.cor       <- array(0, dim=length(samples))
    
    # get common gene ids
    common.gene.ids <- intersect(result1.frame[[1]], result2.frame[[1]])
    # intersect with gene length frame
    if (length(gene.length.frame) > 1) {
        common.gene.ids <- intersect(common.gene.ids, gene.length.frame[[1]])
    }
    # test number of common genes
    if (length(common.gene.ids) < 0.8 * min(c(nrow(result1.frame), nrow(result2.frame)))) {
        message ("Warning: gene id mismatch - we can only use ", length(common.gene.ids), " of ", min(c(nrow(result1.frame), nrow(result2.frame))), " genes")
    }
    
    # get positions of common genes
    result1.ind   <- match(common.gene.ids, result1.frame[[1]], nomatch=0)
    result2.ind   <- match(common.gene.ids, result2.frame[[1]], nomatch=0)
    # make frame with only common gene ids
    result1.filter.frame <- result1.frame[result1.ind,]
    result2.filter.frame <- result2.frame[result2.ind,]
    names(result1.filter.frame) <- names(result1.frame)
    names(result2.filter.frame) <- names(result2.frame)
    
    # initialize normalized frame
    result1.normalized.frame <- result1.filter.frame
    result2.normalized.frame <- result2.filter.frame
    # normalize FPKM or TPM
    if (length(gene.length.frame) > 1) {
        if (length(grep("FPKM", method.name.1)) == 0) {
            ## message("Length normalizing data for ", method.name.1)
            result1.normalized.frame <- length.normalize (result1.filter.frame, gene.length.frame, length.normalization.type)
        }
        if (normalize.counts) {
            result2.normalized.frame <- length.normalize (result2.filter.frame, gene.length.frame, length.normalization.type)
        }
    }
    
    # make count correlations
    for (i in 1:length(samples)) {
        sample.1 <- samples[i]
        sample.2.ind <- i
        if (permute.samples) {
            sample.2.ind <- i %% length(samples) + 1
        }
        sample.2 <- samples[sample.2.ind]
        
        # make correlations for (log) pearson, spearman and CD
        for (method in c('pearson', 'spearman')) {
            # compute log correlation
            sample.log.cor <- cor(log2(result1.normalized.frame[[sample.1]]+pseudo_count), log2(result2.normalized.frame[[sample.2]] + pseudo_count), method=method, use='pairwise.complete.obs')
            # compute correlation
            sample.cor <- cor(result1.normalized.frame[[sample.1]], result2.normalized.frame[[sample.2]], method=method, use='pairwise.complete.obs')
            # save results
            if (method == 'pearson') {
                pearson.cor[i] <- sample.cor
                pearson.log.cor[i] <- sample.log.cor
                CD.cor[i] <- sample.cor^2
                CD.log.cor[i] <- sample.log.cor^2
            } else {
                spearman.cor[i] <- sample.cor
                spearman.log.cor[i] <- sample.log.cor
            }
        }
    }
    
    # get sample names
    if (permute.samples) {
        sample.pairs <- paste(samples, samples[1:length(samples) %% length(samples) + 1], sep="/")
    } else {
        sample.pairs <- samples
    }
    names(pearson.cor)  <- sample.pairs
    names(pearson.log.cor) <- sample.pairs
    names(CD.cor) <- sample.pairs
    names(CD.log.cor) <- sample.pairs
    names(spearman.cor) <- sample.pairs
    names(spearman.log.cor) <- sample.pairs
    
    # return correlations
    return (list( inter.met.pearson      = pearson.cor, 
                  inter.met.spearman     = pearson.log.cor, 
                  inter.met.CD           = CD.cor, 
                  inter.met.log.pearson  = CD.log.cor, 
                  inter.met.log.spearman = spearman.cor, 
                  inter.met.log.CD       = spearman.log.cor))
}




###############################################################################
## Normalize by length (FPKM or TPM)
##
## Originally TPM is defined as:
##     TPM_g = r_g × rl × 10^6 / (fl_g × T)
##     r_g   = total number of reads mapped to g
##     T     = sum_{g in G} r_g × rl / fl_g
##     fl_g  = feature length (i.e. transcript length)
##     rl    = read length
##
## We are working with CPMs which works in the same way.
##     T_g = r_g × rl × 10^6 / (fl_g × T)
##         = (r_g × 10^6 / R) × rl / (fl_g × T / R)
##         = cpm_g × rl / (fl_g × T_CPM / 10^6) = cpm_g × rl × 10^6 / (fl_g × T_CPM)
##         = cpm_g / fl_g × 10^6 / (T_CPM / rl)
##         = cpm_g / fl_g  / ((sum_{g in G} cpm_g / fl_g) / 10^6))
##
##   T_CPM = sum_{g in G} (r_g × 10^6 / R) × rl / fl_g = sum_{g in G} cpm_g × rl / fl_g
##   R     = sum_{g in G} r_g
##
## Note that TPMs just use a different scaling factor (T instead of R) and, hence,
## correlations are not affected by switching to TPMs.
###############################################################################

length.normalize <- function (input.frame, gene.length.frame, length.normalization.type="FPKM") {

  common.gene.ids    <- intersect(input.frame[[1]], gene.length.frame[[1]])
  input.ind          <- match(common.gene.ids, input.frame[[1]], nomatch=0)
  input.filter.frame <- input.frame[input.ind, ]

  gen.len.ind           <- match(common.gene.ids, gene.length.frame[[1]], nomatch=0)
  gene.length.filter.frame <- gene.length.frame[gen.len.ind, ]

  result.frame <- data.frame(lapply(input.filter.frame,
                                    function (x, y) {
                                      if (! is.numeric(x[1])) return (x)
                                      return (1.0 * x / y)
                                    }, gene.length.filter.frame[["Length"]]/1000),
                             check.names=FALSE)
  
  if (length.normalization.type == "TPM") {
    result.frame <- data.frame(lapply(input.filter.frame,
                                      function (x) {
                                        if (! is.numeric(x[1])) return (x)
                                        return (x / (sum(x) / 10^6))
                                      }), check.names=FALSE)
  } else if (length.normalization.type != "FPKM") {
    message("Unknown length normalization type: ", length.normalization.type)
    stopifnot (length.normalization.type %in% c("FPKM", "TPM"))
  }
  
  return (result.frame)
}
















###############################################################################
## Write correlation data to file
###############################################################################

write.matrix.data <- function (input.frame, method.name, table.name, column.names) {
    
    # make matrix
    input.matrix <- make.matrix(input.frame)
    mat.row.names <- rownames(input.matrix)
    
    if (length(column.names) > 1) {
        # make sub frames
        barcode.frame  <- sapply(strsplit(mat.row.names, '-'), function (x) return (x[[3]]))
        lane.frame     <- sapply(strsplit(mat.row.names, '-'), function (x) return (x[[4]]))
        method.frame   <- unlist(lapply(method.name, function(x) rep(x, length(mat.row.names))))
        type.frame     <- sapply(strsplit(mat.row.names, '-'), function (x) return (paste(x[[1]], x[[2]], sep='-')))
        
        # make main frame
        matrix.frame <- data.frame(method.frame, type.frame, barcode.frame, lane.frame, mat.row.names, input.matrix)
    } else {
        # make main frame
        matrix.frame <- data.frame(mat.row.names, input.matrix)
    }
    names(matrix.frame) <- c(column.names, colnames(input.matrix))
    
    # write text file
    if (make.txt.bool) {
        txt.file <- paste(gsub(' ', '_', table.name), 'txt', sep='.')
        matrix.file <- file.path(output.dir, txt.file)
        write.table(matrix.frame, matrix.file, row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
        message('File \'', txt.file, '\' with data values created ...')
    }
}




###############################################################################
## Make matrix out of nested n x n list
###############################################################################

make.matrix <- function (input.list, round.to=-1) {
    
    if (round.to < 0) {
        map.mat <- matrix(unlist(input.list), nrow=length(input.list[[1]]), ncol=length(input.list), dimnames=list(names(input.list[[1]]), names(input.list)))
    } else {
        map.mat <- matrix(round(unlist(input.list), round.to), nrow=length(input.list[[1]]), ncol=length(input.list), dimnames=list(names(input.list[[1]]), names(input.list)))
    }
    return(map.mat)
}




###############################################################################
## Plot heatmap
## input: symmetrical list of list with samples
###############################################################################

heatMap <- function (input.list, plot.title.name, legend.title, x.lab.name, y.lab.name, use.triangle.bool=TRUE, make.new.matrix=TRUE) {
    
    # create matrix for plotting
    if (make.new.matrix) {
        map.mat <- make.matrix(input.list, 2)
    } else {
        map.mat <- round(input.list,2)
    }
    if (use.triangle.bool) {
        map.mat.upper <- get.triangle(map.mat, 'upper')
        legend.setting <- list(c(0.3, 0.85),'horizontal',20,4,0.5,element_blank(),element_blank())
        l.char.space <- 0
        t.char.space <- 0
        map.mat.melt <- melt(map.mat.upper, na.rm=TRUE)
    } else {
        map.mat.upper <- map.mat
        legend.setting <- list('right', 'vertical', 2, 10, 0.0,element_text(size=16, margin=margin(20,00,00,00)),element_text(size=16, margin=margin(00,20,00,00)))
        l.char.space <- nchar(legend.title)/6
        t.char.space <- 1
        map.mat.melt <- melt(map.mat.upper, na.rm=FALSE)
    }
    
    p <- ggplot(data = map.mat.melt, aes(Var2, Var1, fill = value)) + 
            scale_fill_gradient2(low='blue', high='red', mid='white', midpoint=0.5, limit=c(0,1), space='Lab', name=legend.title) + 
            geom_tile(color='white') +                                                                      # make white background tiles
            xlab(x.lab.name) + ylab(y.lab.name) + 
            #ggtitle(plot.title.name) + 
            coord_fixed() + 
            geom_text(aes(Var2, Var1, label=value), color='black', size=4) +                                # draw correlation numbers
            theme(plot.title=element_text(size=20), 
                #axis.title.x=element_text(size=16, margin=margin(20,00,00,00)), 
                #axis.title.x = element_blank(),                                                             # remove x axis title
                #axis.title.y=element_text(size=16, margin=margin(00,20,00,00)),
                #axis.title.y = element_blank(),                                                             # remove y axis title
                axis.title.x = legend.setting[[6]],                                                             # remove y axis title
                axis.title.y = legend.setting[[7]],                                                             # remove y axis title
                plot.margin=margin(10,10,20,10),                                                            # top, right, bottom, left
                axis.text.x=element_text(size=16, angle=45, vjust=0.9, hjust=1.0), 
                axis.text.y=element_text(size=16, margin=margin(00,10,00,00)), 
                panel.grid.major=element_blank(),                                                           # remove grid lines
                panel.grid.minor=element_blank(), 
                panel.border=element_blank(),                                                               # remove border lines
                panel.background=element_blank(),                                                           # remove background
                axis.ticks=element_blank(),                                                                 # remove axis ticks
                legend.justification='center',
                legend.position=legend.setting[[1]],                                                               # horizontal, vertical
                legend.direction=legend.setting[[2]],
                legend.text=element_text(size=16),
                legend.title=element_text(size=16)) + 
            guides(fill=guide_colorbar(barwidth=legend.setting[[3]], barheight=legend.setting[[4]], title.position='top', title.hjust=legend.setting[[5]]))  # adjust legend position
    
    # print data
    if (make.print.bool) {
        print(p)
        message('Plot printed ...')
    }
    
    # print pdf
    if (make.pdf.bool) {
        pdf.file <- paste(gsub(' ', '_', plot.title.name), 'pdf', sep='.')
        w.char.space <- 2 #(nchar(max(names(input.list)))/2.5)
        h.char.space <- 2 #w.char.space/sqrt(2)
        w.bars <- 30 #2*length(unique(color.group))*length(unique(shape.group))
        h.bars <- 30 # 20
        l.char.space <- 2
        ggsave(file.path(output.dir, pdf.file), width=(1+w.char.space+l.char.space+w.bars), height=(1+h.char.space+h.bars), scale=1.0, units='cm')
        message('File \'', pdf.file, '\' created ...')
    }
}




###############################################################################
## Plot a boxplot for method
###############################################################################

boxPlot <- function (input.list, plot.title.name, x.lab.name='counting method', y.lab.name='gene count', reduce.jitter=TRUE, use.jitter=FALSE, use.axis.scale=FALSE, dataset.name='SEQC', use.scale=c(0.0,0.0), 
                     x.angle.setting=45, v.just.setting=0.9, h.just.setting=1.0, scale.factor.setting=1.0) {
    
    # create data frame for plotting
    #sum.tool <- unlist(lapply(names(input.list),function(x) return(rep(x, 6426))))#length(input.list[[x]])))))
    #sum.data <- unlist(lapply(input.list, function(x) return(unlist(sort(x)[seq(1,64251,10)], use.names=FALSE))), use.names=FALSE)
    i.len <- length(input.list[[1]])
    sum.tool <- unlist(lapply(names(input.list),function(x) return(rep(x, length(input.list[[x]])))))
    sum.data <- unlist(lapply(input.list, function(x) return(unlist(x, use.names=FALSE))), use.names=FALSE)
    max.data <- unlist(lapply(input.list, function(x) return(unlist(c(max(x),rep(NA, i.len-1)), use.names=FALSE))), use.names=FALSE)
    min.data <- unlist(lapply(input.list, function(x) return(unlist(c(min(x),rep(NA, i.len-1)), use.names=FALSE))), use.names=FALSE)
    sum.frame <- data.frame(sum.tool, sum.data, max.data, min.data)#)
    
    # set colors and shape groups
    color.group <- as.factor(sapply(sum.frame$sum.tool, function(x) strsplit(as.character(x), '/')[[1]][1]))
    shape.group <- as.factor(sapply(sum.frame$sum.tool, function(x) strsplit(as.character(x), '/')[[1]][2]))
    sum.frame <- cbind(sum.frame, color.group, shape.group)
    print(unique(color.group))
    print(unique(shape.group))
    
    # make frame
    #sum.point <- unlist(lapply(names(input.list),function(x) return(x)))
    #max.point <- unlist(lapply(input.list, function(x) return(unlist(max(x), use.names=FALSE))), use.names=FALSE)
    #min.point <- unlist(lapply(input.list, function(x) return(unlist(min(x), use.names=FALSE))), use.names=FALSE)
    #point.frame <- data.frame(sum.point, max.point, min.point)#)
    #color.point <- unique(as.factor(sapply(sum.frame$sum.tool, function(x) strsplit(as.character(x), '/')[[1]][1])))
    #shape.point <- unique(as.factor(sapply(sum.frame$sum.tool, function(x) strsplit(as.character(x), '/')[[1]][2])))
    #point.frame <- cbind(point.frame, color.point, shape.point)
    #geom_pointrange(data=point.frame, aes(x=color.point, ymin=min.point, ymax=max.point)) + 
    
    
    # set jitter alpha
    if (use.jitter) {
        if (reduce.jitter) {
            jitter.data <- data.frame()
            xth.data <- round(((length(sum.frame[['sum.data']]) / (length(unique(color.group))*length(unique(shape.group))))/200),0)
            for (i in 1:length(sum.frame[['sum.data']])) {
                if ((i %% xth.data) == 0) {
                    jitter.data <- rbind(jitter.data,sum.frame[i,])
                }
            }
            jitter.shape <- as.factor(sapply(jitter.data$sum.tool, function(x) strsplit(as.character(x), '/')[[1]][2]))
        } else {
            jitter.data <- sum.frame
            jitter.shape <- shape.group
        }
    }
    # print data as boxplots
    p <- ggplot(data=sum.frame, aes(x=color.group, y=sum.data))+#, ymin=min.data, ymax=max.data)) + 
            stat_boxplot(geom='errorbar') + 
            geom_pointrange(aes(x=color.group, ymin=min.data, ymax=max.data), size=0.5, linetype='dashed') + 
            geom_boxplot(aes(fill=color.group), alpha=1.0, outlier.size=NA) + 
            #geom_crossbar(width = 0.2) +
            #scale_shape_discrete(name='aligner')+#, labels=names(c.list.aligners)) + 
            #scale_fill_brewer(name='quantifier', palette='Paired')+#, labels=names(c.list.quantifiers)) + 
            #scale_colour_discrete(name='aligner', palette='Set3')+#, labels=names(c.list.aligners)) + 
            xlab(x.lab.name) + ylab(y.lab.name) + 
            #ggtitle(plot.title.name) + 
            theme(plot.title=element_text(size=20), 
                axis.title.x=element_text(size=16, margin=margin(20,00,00,00)), 
                axis.title.y=element_text(size=16, margin=margin(00,20,00,00)),
                axis.ticks.x = element_blank(), 
                plot.margin=margin(10,10,10,10), # top, right, bottom, left
                axis.text.x=element_text(size=12, angle=x.angle.setting, vjust=v.just.setting, hjust=h.just.setting), 
                axis.text.y=element_text(size=12, margin=margin(00,10,00,00)))
    if (dataset.name == 'SEQC') {
        # manual colors for BTSEQ−default, CUFFL−default, EQPQM−ambig, EQPQM−unambig, FEATC−default, FLXCP−default, HTSEQ−inonempt, HTSEQ−istrict, HTSEQ−union, KLLST−default, RSEM−default
        p <- p + scale_fill_manual(name='quantifier', values=c('tomato','wheat','skyblue','royalblue','tan1','darkorchid1','palegreen','seagreen1','seagreen4','gold1','tan4')) +
                 scale_color_manual(name='quantifier', values=c('tomato','wheat','skyblue','royalblue','tan1','darkorchid1','palegreen','seagreen1','seagreen4','gold1','tan4'))
    }
    if (dataset.name == 'SEQS') {
        # manual colors for BTSEQ−default, CUFFL−default, EQPQM−ambig, EQPQM−unambig, FEATC−default, FLXCP−default, FLXSM-default, HTSEQ−inonempt, HTSEQ−istrict, HTSEQ−union, KLLST−default, RSEM−default
        p <- p + scale_fill_manual(name='quantifier', values=c('tomato','wheat','skyblue','royalblue','tan1','darkorchid1','hotpink','palegreen','seagreen1','seagreen4','gold1','tan4')) +
                 scale_color_manual(name='quantifier', values=c('tomato','wheat','skyblue','royalblue','tan1','darkorchid1','hotpink','palegreen','seagreen1','seagreen4','gold1','tan4'))
    }
    if (dataset.name == 'TaqMan') {
        # manual colors for BTSEQ−default, CUFFL−default, EQPQM−ambig, EQPQM−unambig, FEATC−default, FLXCP−default, FLXSM-default, HTSEQ−inonempt, HTSEQ−istrict, HTSEQ−union, KLLST−default, RSEM−default
        p <- p + scale_fill_manual(name='quantifier', values=c('tomato','wheat','skyblue','royalblue','tan1','darkorchid1','palegreen','seagreen1','seagreen4','gold1','tan4','hotpink')) +
                 scale_color_manual(name='quantifier', values=c('tomato','wheat','skyblue','royalblue','tan1','darkorchid1','palegreen','seagreen1','seagreen4','gold1','tan4','hotpink'))
    }
    if (dataset.name == 'BioRad') {
        # manual colors for BTSEQ−default, CUFFL−default, EQPQM−ambig, EQPQM−unambig, FEATC−default, FLXCP−default, FLXSM-default, HTSEQ−inonempt, HTSEQ−istrict, HTSEQ−union, KLLST−default, RSEM−default
        p <- p + scale_fill_manual(name='quantifier', values=c('hotpink','tomato','wheat','skyblue','royalblue','tan1','darkorchid1','palegreen','seagreen1','seagreen4','gold1','tan4')) +
                 scale_color_manual(name='quantifier', values=c('hotpink','tomato','wheat','skyblue','royalblue','tan1','darkorchid1','palegreen','seagreen1','seagreen4','gold1','tan4'))
    }
    if (use.scale[[1]] != 0.0 | use.scale[[2]] != 0.0) {
        p <- p + coord_cartesian(ylim = use.scale)
    }
    if (use.jitter) {
        p <- p + geom_jitter(aes(color=color.group), alpha=0.8, size=2, width = 1.0, height = 1.0, na.rm=TRUE)
        #p <- p + geom_jitter(data=sum.frame, aes(shape=shape.group), alpha=0.8, size=2, width = 1.0, height = 1.0, na.rm=TRUE)
    }
    if (use.axis.scale) {
        #p <- p + coord_trans(y='log2',limy=5)
        p <- p + scale_size_manual(values=c(-10,-5,0,5,10))
    }
    # add facets
    p <- p + facet_grid(. ~ shape.group, scale='free_x', space='free_x')
    
    # print data
    if (make.print.bool) {
        print(p)
        message('Plot printed ...')
    }
    
    # print pdf
    if (make.pdf.bool) {
        pdf.file <- paste(gsub(' ', '_', plot.title.name), 'pdf', sep='.')
        w.char.space <- 0 #nchar(max(names(input.list)))/2
        h.char.space <- ((nchar(max(names(input.list)))/2.5)/sqrt(2))
        w.bars <- 5*length(c.list.settings)
        h.bars <- 20
        l.char.space <- 6
        ggsave(file.path(output.dir, pdf.file), width=(1+w.char.space+w.bars+l.char.space)*0.3, height=(1+h.char.space+h.bars), scale=scale.factor.setting, units='cm', limitsize=FALSE)
        message('File \'', pdf.file, '\' created ...')
    }
    
    
    # sum.method <- unlist(lapply(c('expression.genes.total', 'expression.genes.mean'), function(x) return(rep(x, (samples.counts.A + samples.counts.B) * c.list.length))))
    # sum.sample <- rep(c(rep('SEQC-A', samples.counts.A*c.list.length), rep('SEQC-B', samples.counts.B*c.list.length)), 2)
    # sum.tool.A <- unlist(lapply(names(c.list.counts),function(x) return(rep(x, samples.counts.A))))
    # sum.tool.B <- unlist(lapply(names(c.list.counts),function(x) return(rep(x, samples.counts.B))))
    # sum.tool <- rep(c(sum.tool.A, sum.tool.B), 2)
    # sum.data <- unlist(lapply(c(summary.list$expression.genes.total, summary.list$expression.genes.mean), function(x) return(unlist(x, use.names=FALSE))), use.names=FALSE)
    # sum.frame <- data.frame(sum.method, sum.sample, sum.tool, sum.data)
    # # set colors and shape groups
    # color.group <- as.factor(sapply(sum.frame$sum.tool, function(x) strsplit(as.character(x), '/')[[1]][1]))
    # shape.group <- as.factor(sapply(sum.frame$sum.tool, function(x) strsplit(as.character(x), '/')[[1]][2]))
    # # set facet naming
    # facet.label.list <- list('sum.method'=c('expression.genes.mean'='mean gene expression', 'expression.genes.total'='total gene expression'), 
          # 'sum.sample'=c('SEQC-A'='SEQC-A', 'SEQC-B'='SEQC-B'))
    # facet.label.func <- function(variable,value) {return(facet.label.list[[variable]][value])}
    # # print data as 2x2 boxplot
    # p <- ggplot(data=sum.frame, aes(x=sum.tool, y=sum.data)) +
            # stat_boxplot(geom='errorbar') +
            # xlab('counting method') + ylab('sample type') + 
            # geom_boxplot(aes(fill=color.group), alpha=1.0) +
            # geom_jitter(aes(shape=shape.group), alpha=0.7, size=3) +
            # facet_grid(sum.method ~ sum.sample, scale='free', labeller=facet.label.func) +
            # scale_fill_brewer(name='quantifier', palette='Set1', labels=names(c.list.quantifiers)) +
            # scale_shape_discrete(name='aligner', labels=names(c.list.aligners)) + 
            # theme(text=element_text(size=12), axis.text.x=element_text(angle=45, vjust=0.9, hjust=1.0))
            # #scale_fill_manual(values=rainbow(3))
    # ggsave(file.path(output.dir, plot.titles.list['summary.title'][[1]]), width=29.7, height=21, scale=1.4, units='cm')
    # print(p)
}



###############################################################################
## Plot inter method correlation heatmap
## input: list of inter method correlations
###############################################################################

multiHeatMap <- function (input.list, plot.title.name, legend.title, lab.name) {
    
    # create data frame for
    data.vector <- unlist(input.list)
    tool.2.frame <- unlist(lapply(names(data.vector), function(x) return(strsplit(x,'.',fixed=TRUE)[[1]][[1]])))
    tool.1.frame <- unlist(lapply(names(data.vector), function(x) return(strsplit(x,'.',fixed=TRUE)[[1]][[2]])))
    sample.frame <- unlist(lapply(names(data.vector), function(x) return(strsplit(x,'.',fixed=TRUE)[[1]][[3]])))
    item.frame <- round(unname(data.vector),4)
    
    plot.frame <- data.frame(sample.frame, tool.2.frame, tool.1.frame, item.frame)
    
    p <- ggplot(data=plot.frame, aes(x=tool.2.frame, y=tool.1.frame, fill=item.frame)) + 
            scale_fill_gradient2(low='blue', high='red', mid='white', midpoint=0.5, limit=c(0,1), space='Lab', name=legend.title) + 
            scale_x_discrete(expand=c(0,0)) + 
            scale_y_discrete(expand=c(0,0)) + 
            geom_tile(color='white') +                                                                      # make white background tiles
            xlab(lab.name) + ylab(lab.name) + 
            #ggtitle(plot.title.name) + 
            coord_fixed() + 
            geom_text(aes(x=tool.2.frame, y=tool.1.frame, label=item.frame), color='black', size=4) +                                # draw correlation numbers
            theme(plot.title=element_text(size=20), 
                #axis.title.x=element_text(size=16, margin=margin(20,00,00,00)), 
                axis.title.x = element_blank(),                                                             # remove x axis title
                #axis.title.y=element_text(size=16, margin=margin(00,20,00,00)),
                axis.title.y = element_blank(),                                                             # remove y axis title
                plot.margin=margin(10,10,20,10),                                                            # top, right, bottom, left
                axis.text.x=element_text(size=12, angle=45, vjust=0.9, hjust=1.0), 
                axis.text.y=element_text(size=12, margin=margin(00,10,00,00)), 
                panel.grid.major=element_blank(),                                                           # remove grid lines
                panel.grid.minor=element_blank(), 
                panel.border=element_blank(),                                                               # remove border lines
                panel.background=element_blank(),                                                           # remove background
                axis.ticks=element_blank(),                                                                 # remove axis ticks
                legend.justification='center',
                legend.position='right',                                                               # horizontal, vertical
                legend.direction='vertical',
                legend.text=element_text(size=10),
                legend.title=element_text(size=10)) + 
            facet_grid(. ~ sample.frame, scale='free') +
            guides(fill=guide_colorbar(barwidth=2, barheight=10, title.position='top', title.hjust=0.0))  # adjust legend position
    
    # print data
    if (make.print.bool) {
        print(p)
        message('Plot printed ...')
    }
    
    # print pdf
    if (make.pdf.bool) {
        pdf.file <- paste(gsub(' ', '_', plot.title.name), 'pdf', sep='.')
        w.bars <- 2*length(unique(tool.2.frame))*length(unique(sample.frame))
        h.bars <- 2*length(unique(tool.1.frame))
        w.char.space <- max(nchar(unique(tool.2.frame)))/2.5
        h.char.space <- w.char.space/sqrt(2)
        l.char.space <- nchar(legend.title)/6
        ggsave(file.path(output.dir, pdf.file), width=(1+w.char.space+l.char.space+w.bars), height=(1+h.char.space+h.bars), scale=0.9, units='cm', limitsize=FALSE)
        message('File \'', pdf.file, '\' created ...')
    }
}




###############################################################################
## Make Method vs Method Scatterplot
###############################################################################

scatterPlot <- function(input.list.A, input.list.B, plot.title.name, method.1, method.2, x.lab.name='counting method', y.lab.name='gene count', alpha.val=1.0, dataset.name='SEQC', single.data=FALSE) {
    
    if (!single.data) {
        d.set.A <- paste(dataset.name,'-A',sep='')
        d.set.B <- paste(dataset.name,'-B',sep='')
        input.list <- list(input.list.A,input.list.B)
        names(input.list) <- list(d.set.A,d.set.B)
        item.frame <- data.frame()
        # create data frame (sample.type, xvar, yvar)
        #for (method.1 in names(input.list.A)) {
        #for (method.2 in names(input.list.A)) {
        #method.1 <- 'EQPQM-ambig/HISAT'
        #method.2 <- 'FEATC-default/HISAT'
        for (sample.type in c(d.set.A,d.set.B)) {
            data.m1 <- unname(unlist(input.list[[sample.type]][[method.1]]))
            data.m2 <- unname(unlist(input.list[[sample.type]][[method.2]]))
            method.m1 <- rep(method.1, length(data.m1))
            method.m2 <- rep(method.2, length(data.m2))
            sample.type <- rep(as.factor(sample.type),length(data.m1))
            d.frame <- data.frame(method.m1, method.m2, sample.type, data.m1, data.m2)
            d.frame <- na.omit(d.frame)
            d.frame$data.m1 <- round((d.frame$data.m1), 1)
            d.frame$data.m2 <- round((d.frame$data.m2), 1)
            d.frame <- unique(d.frame)
            item.frame <- rbind(item.frame, d.frame)
        }
    } else {
        sample.type <- paste(dataset.name,'-A-B',sep='')
        item.frame <- data.frame()
        data.m1 <- unname(unlist(input.list.A[[method.1]]))
        data.m2 <- unname(unlist(input.list.A[[method.2]]))
        method.m1 <- rep(method.1, length(data.m1))
        method.m2 <- rep(method.2, length(data.m2))
        sample.type <- rep(as.factor(sample.type),length(data.m1))
        d.frame <- data.frame(method.m1, method.m2, sample.type, data.m1, data.m2)
        d.frame <- na.omit(d.frame)
        d.frame$data.m1 <- round(d.frame$data.m1, 2)
        d.frame$data.m2 <- round(d.frame$data.m2, 2)
        d.frame <- unique(d.frame)
        item.frame <- rbind(item.frame, d.frame)
    }
    #}
    #}
    
    
    # Extend the regression lines beyond the domain of the data
    p <- ggplot(data=item.frame, aes(x=data.m1, y=data.m2, color=sample.type)) + 
                geom_point(shape=1, alpha=alpha.val) +
                scale_colour_hue(l=50) + # Use a slightly darker palette than normal
                geom_smooth(method=lm,   # Add linear regression lines
                            se=FALSE,    # Don't add shaded confidence region
                            fullrange=TRUE) + # Extend regression lines
                xlab(method.1) + ylab(method.2) +
                #scale_x_continuous(limits=c(-15,15), breaks=c(-10,-5,0,5,10)) +
                #scale_y_continuous(limits=c(-15,15), breaks=c(-10,-5,0,5,10)) +
                coord_cartesian(xlim=c(-15,15), ylim = c(-15,15))
                theme(legend.title=element_blank())
    if (x.lab.name != 'counting method' & y.lab.name != 'gene count') {
        p <- p + xlab(x.lab.name) + ylab(y.lab.name)
    }
    #facet_grid(method.m1 ~ method.m2, scale='free')
    # print data
    if (make.print.bool) {
        print(p)
        message('Plot printed ...')
    }
    
    # print pdf
    if (make.pdf.bool) {
        pdf.file <- paste(gsub(' ', '_', plot.title.name), 'pdf', sep='.')
        w.char.space <- 2 #nchar(max(names(input.list)))/2
        h.char.space <- 2 #(nchar(method.1)/2.5)/sqrt(2)
        w.bars <- 18 #2*length(unique(color.group))*length(unique(shape.group))
        h.bars <- 18 # 20
        l.char.space <- 2
        ggsave(file.path(output.dir, pdf.file), width=(1+w.char.space+w.bars+l.char.space), height=(1+h.char.space+h.bars), scale=1.0, units='cm')
        message('File \'', pdf.file, '\' created ...')
    }
}




###############################################################################
# Scatter Plot of FC(A,B) vs mean(A,B)
###############################################################################

meanFCScatterPlot <- function(m.FC.list, m.norm.count.list, m.sample.type.A, m.sample.type.B, 
            method.occurence, color.vector, plot.title.name, round.to.xth=100, method.names, sample.thresh, plot.variance=FALSE) {
    
    # make sample combination
    sample.AB <- paste(m.sample.type.A, m.sample.type.B, sep='/')
    # make data vector
    item.frame <- data.frame()
    for (method.n in method.names) {
        tmp.fc.list <- m.FC.list[[method.n]][,sample.AB]
        tmp.mean.list <- (m.norm.count.list[[method.n]][,m.sample.type.A] + m.norm.count.list[[method.n]][,m.sample.type.B]) / 2
        data.m1 <- tmp.fc.list
        data.m2 <- tmp.mean.list
        #data.m2 <- tmp.mean.list
        method.type <- rep(as.factor(method.n),length(data.m1))
        gene.names <- m.norm.count.list[[1]][['Gene Id']]
        d.frame <- data.frame(gene.names, method.type, data.m1, data.m2)
        d.frame <- na.omit(d.frame)
        # d.frame$data.m2 <- unlist(lapply(d.frame$data.m2, function(x) {
                               # if (x < 1) return(round(x, 4))
                               # else if ((x >= 1) && (x < 10)) return(round(x, 3))
                               # else if ((x >= 10) && (x < 100)) return(round(x, 2))
                               # else if ((x >= 100) && (x < 1000)) return(round(x, 1))
                               # else return(round(x, 0))}))
        item.frame <- rbind(item.frame, d.frame)
    }
    rownames(item.frame) <- 1:length(item.frame[,1])
    # check if log2 FC and log2 mean is the same (rounded) for genes
    g.num <- length(m.norm.count.list[[1]][['Gene Id']])
    p.vec <- c()
    for (i in 0:(method.occurence-1)) {
        p.vec <- c(p.vec, i*g.num)
    }
    multi.occ.list <- list()
    multi.occ.list <- lapply(c(1:g.num), function(i, p.vec, sample.thresh, item.frame) {
                           gene.positions <- p.vec+i
                           FC.i <- item.frame[gene.positions,3]
                           mn.i <- item.frame[gene.positions,4]
                           s.t <- sample.thresh/100
                           FC.i
                           # check if the FC and mn values are in between 5% of the mean of each other
                           #if ((c.neg.FC <= FC.i) && (FC.i <= c.pos.FC) && (c.neg.mn <= mn.i) && (mn.i <= c.pos.mn)) {
                           if ((mean(FC.i)*(1+s.t) >= FC.i && mean(FC.i)*(1-s.t) <= FC.i) && (mean(mn.i)*(1+s.t) >= mn.i && mean(mn.i)*(1-s.t) <= mn.i)) {
                               gene.names <- as.factor(as.character(item.frame[i,1]))
                               if (plot.variance) {
                                   method.type <- as.factor(paste('log2(FC(A,B)) and log2(var(log2(FC(A,B)))) within ', as.character(sample.thresh), '%'))
                               } else {
                                   method.type <- as.factor(paste('log2(FC(A,B)) and log2(mean(A,B)) within ', as.character(sample.thresh), '%'))
                               }
                               data.m1 <- mean(FC.i)
                               data.m2 <- mean(mn.i)
                               o.frame <- data.frame(gene.names, method.type, data.m1, data.m2)
                               #o.list <- list(gene.names=gene.names, method.type=method.type, data.m1=data.m1, data.m2=data.m2)
                               return(o.frame)
                           }
                       }, p.vec, sample.thresh, item.frame)
    
    
    multi.occ.list <- multi.occ.list[!sapply(multi.occ.list, is.null)]
    gene.names <- sapply(multi.occ.list, function(x) return(x$gene.names))
    method.type <- sapply(multi.occ.list, function(x) return(x$method.type))
    data.m1 <- c()
    data.m2 <- c()
    data.m1 <- sapply(multi.occ.list, function(x) return(x$data.m1))
    data.m2 <- sapply(multi.occ.list, function(x) return(x$data.m2))
    multi.occ.frame <- data.frame(gene.names, method.type, data.m1, data.m2)
    
    # make plot frame - m.item.frame
    set.alpha <- rep(1.0,length(multi.occ.frame[,1]))
    m.item.frame <- cbind(multi.occ.frame, set.alpha)
    m.item.frame$data.m1 <- log2(m.item.frame$data.m1)
    m.item.frame$data.m2 <- log2(m.item.frame$data.m2)
    m.item.frame$data.m1 <- round(m.item.frame$data.m1*round.to.xth)/round.to.xth
    m.item.frame$data.m2 <- round(m.item.frame$data.m2*round.to.xth)/round.to.xth
    
    # make plot frame - a.item.frame
    a.item.frame <- item.frame
    a.item.frame$data.m1 <- log2(a.item.frame$data.m1)
    a.item.frame$data.m2 <- log2(a.item.frame$data.m2)
    alpha.list <- list()
    for (i in c(100,90,80,70,60,50,40,30,20,10,9,8,7,6,5,4,3,2,1)) {
        a.item.frame$data.m1 <- round(a.item.frame$data.m1*i)/i
        a.item.frame$data.m2 <- round(a.item.frame$data.m2*i)/i
        alpha.list[[as.character(i)]] <- !duplicated(a.item.frame[c('method.type', 'data.m1', 'data.m2')], fromLast=TRUE)
    }
    
    # make plot.frame - p.item.frame
    set.alpha <- rep(0.01,length(item.frame[,1]))
    p.item.frame <- cbind(item.frame, set.alpha)
    p.item.frame$data.m1 <- log2(p.item.frame$data.m1)
    p.item.frame$data.m2 <- log2(p.item.frame$data.m2)
    p.item.frame$data.m1 <- round(p.item.frame$data.m1*round.to.xth)/round.to.xth
    p.item.frame$data.m2 <- round(p.item.frame$data.m2*round.to.xth)/round.to.xth
    for (i in names(alpha.list)) {
        p.item.frame[alpha.list[[i]],'set.alpha'] <- 1/as.integer(i)
    }
    p.item.frame <- rbind(p.item.frame, m.item.frame)
    p.item.frame <- p.item.frame[!duplicated(p.item.frame[c('method.type', 'data.m1', 'data.m2')], fromLast=TRUE),]
    set.alpha <- p.item.frame[,'set.alpha']
    
    # sort according to alpha values
    p.item.frame <- p.item.frame[order(p.item.frame$set.alpha),]
    p.item.frame <- p.item.frame[!duplicated(p.item.frame[c('data.m1', 'data.m2')], fromLast=TRUE),]
    
    # make plot
    p <- ggplot(data=p.item.frame, aes(x=data.m2, y=data.m1, color=method.type, alpha=set.alpha)) + 
                scale_alpha(guide = 'none', range=c(0.0,1.0)) +
                geom_point() +
                scale_color_manual(name='method', values=color.vector) + 
                #scale_x_continuous(trans='log2')
                #scale_x_continuous(limits=c(-0.05,2.0), breaks=c(0.0,0.5,1.0,1.5)) +
                #scale_y_continuous(limits=c(-3.5,3.5), breaks=c(-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0)) +
                #scale_colour_hue(l=50) +
                theme(plot.title=element_text(size=20), 
                      axis.title.x=element_text(size=16, margin=margin(20,00,00,00)), 
                      axis.title.y=element_text(size=16, margin=margin(00,20,00,00)),
                      axis.ticks.x = element_blank(), 
                      plot.margin=margin(10,10,10,10), # top, right, bottom, left
                      axis.text.x=element_text(size=12, angle=0, vjust=0.9, hjust=1.0), 
                      axis.text.y=element_text(size=12, margin=margin(00,10,00,00)), 
                      legend.text=element_text(size=10),
                      legend.title=element_text(size=10))
    if (plot.variance) {
        p <- p + scale_x_continuous(limits = c(-18.0, 8.0), breaks=c(-15.0,-10.0,-5.0,-0.0,5.0)) + 
                 scale_y_continuous(limits = c(-13.0, 13.0), breaks=c(-10.0,-5.0,-0.0,5.0,10.0)) + 
                 xlab('log2(var(log2(FC(A,B))))') + ylab('log2(FC(A,B))')
    } else {
        p <- p + scale_x_continuous(limits = c(-18.0, 18.0), breaks=c(-15.0,-10.0,-5.0,-0.0,5.0,10.0,15.0)) + 
                 scale_y_continuous(limits = c(-13.0, 13.0), breaks=c(-10.0,-5.0,-0.0,5.0,10.0)) + 
                 xlab('log2(mean(A,B))') + ylab('log2(FC(A,B))')
    }
    
    # create pdf
    pdf.file <- paste(gsub(' ', '_', plot.title.name), 'pdf', sep='.')
    w.char.space <- 2 #nchar(max(names(input.list)))/2
    h.char.space <- 2 #(nchar(method.1)/2.5)/sqrt(2)
    w.bars <- 24 #2*length(unique(color.group))*length(unique(shape.group))
    h.bars <- 24 # 20
    l.char.space <- 16
    ggsave(file.path(output.dir, pdf.file), width=(1+w.char.space+w.bars+l.char.space), height=(1+h.char.space+h.bars), scale=1.0, units='cm')
    message('File \'', pdf.file, '\' created ...')
    return(multi.occ.frame)
}





###############################################################################
# Confidence interval
###############################################################################

confInv <- function(data.v, alpha) {
    
    # Konfidenzintervall:
    # Erwartungswert eines normalverteilten Merkmals mit bekannter Varianz \sigma^2:
    # z_(1-α/2) -> α = 0.05
    # z_(0.975) -> 97,5-%-Quantil
    # [mean - z_(1-α/2) * sd / sqrt(n)  ;  mean + z_(1-α/2) * sd / sqrt(n)]
    
    m <- mean(data.v)
    z <- qnorm(1-alpha/2)
    s <- sd(data.v)
    n <- sqrt(length(data.v))
    c.neg <- m - z * s/n
    c.pos <- m + z * s/n
    
    return(c(c.neg,c.pos))
}




###############################################################################
# ROC calculations
###############################################################################

N.c <- function(gold.vec,g.thresh) {
    N <- sum(abs(gold.vec) < g.thresh)
    return(N)
}
TN.c <- function(FC.vec,thresh,gold.vec,g.thresh) {
    #TN <- sum((FC.vec < thresh) & (gold.vec < g.thresh))
    TN <- sum((abs(FC.vec) < thresh) & (abs(gold.vec) < g.thresh) & ((gold.vec < 0.0) & (FC.vec < 0.0) | (gold.vec >= 0.0) & (FC.vec >= 0.0)))
    return(TN)
}
P.c <- function(gold.vec,g.thresh) {
    P <- sum(abs(gold.vec) >= g.thresh)
    return(P)
}
TP.c <- function(FC.vec,thresh,gold.vec,g.thresh) {
    TP <- sum((abs(FC.vec) >= thresh) & (abs(gold.vec) >= g.thresh) & ((gold.vec < 0.0) & (FC.vec < 0.0) | (gold.vec >= 0.0) & (FC.vec >= 0.0)))
    return(TP)
}
TPR.c <- function(FC.vec,gold.vec,thresh,g.thresh,pseudo) {
    P <- P.c(gold.vec,g.thresh)
    TP <- TP.c(FC.vec,thresh,gold.vec,g.thresh)
    TPR <- (TP + pseudo)/(P + pseudo)
    return(TPR)
}
FPR.c <- function(FC.vec,gold.vec,thresh,g.thresh,pseudo) {
    N <- N.c(gold.vec,g.thresh)
    TN <- TN.c(FC.vec,thresh,gold.vec,g.thresh)
    FPR <- 1-((TN + pseudo)/(N + pseudo))
    return(FPR)
}


plot.ROC <- function(input.list, gold.name, density.setting, pseudo, g.thresh, plot.title.name) {
    
    # make data
    data.list <- list()
    data.list[['TPR']] <- list()
    data.list[['FPR']] <- list()
    for (method in names(input.list)) {
        if (method != gold.name) {
            data.vec <- input.list[[method]]
            gold.vec <- input.list[[gold.name]]
            min.val <- min(data.vec)
            max.val <- max(data.vec)
            thresh.vec <- seq(min.val,max.val,length.out=density.setting)
            TPR.vec <- c()
            FPR.vec <- c()
            for (thresh in thresh.vec) {
                TPR.vec <- c(TPR.vec,TPR.c(data.vec,gold.vec,thresh,g.thresh,pseudo))
                FPR.vec <- c(FPR.vec,FPR.c(data.vec,gold.vec,thresh,g.thresh,pseudo))
            }
            data.list[['TPR']][[method]] <- TPR.vec
            data.list[['FPR']][[method]] <- FPR.vec
        }
    }
    
    # make frame
    data.TPR <- unname(unlist(data.list[['TPR']]))
    data.FPR <- unname(unlist(data.list[['FPR']]))
    method.m <- rep(names(data.list[['TPR']]), density.setting)
    item.frame <- data.frame(method.m, data.TPR, data.FPR)
    
    # make plot
    p <- ggplot(data=item.frame, aes(x=data.FPR, y=data.TPR, color=method.m)) + 
                geom_line(size=1.0, alpha=1.0) +
                scale_colour_hue(l=50) + # Use a slightly darker palette than normal
                xlab('False Positive Rate') + ylab('True Positive Rate') +
                #scale_x_continuous(limits=c(-15,15), breaks=c(-10,-5,0,5,10)) +
                #scale_y_continuous(limits=c(-15,15), breaks=c(-10,-5,0,5,10)) +
                coord_cartesian(xlim=c(0.0,1.0), ylim = c(0.0,1.0), expand=FALSE) +
                scale_color_manual(name='quantifier', values=c('tomato','wheat','skyblue','royalblue','tan1','darkorchid1','palegreen','seagreen1','seagreen4','gold1','tan4'))
    
    # print pdf
    if (make.pdf.bool) {
        pdf.file <- paste(gsub(' ', '_', plot.title.name), 'pdf', sep='.')
        w.char.space <- 2 #nchar(max(names(input.list)))/2
        h.char.space <- 2 #(nchar(method.1)/2.5)/sqrt(2)
        w.bars <- 18 #2*length(unique(color.group))*length(unique(shape.group))
        h.bars <- 18 # 20
        l.char.space <- 2
        ggsave(file.path(output.dir, pdf.file), width=(1+w.char.space+w.bars+l.char.space), height=(1+h.char.space+h.bars), scale=1.0, units='cm')
        message('File \'', pdf.file, '\' created ...')
    }
}



make.TPR.FPR <- function(input.list, gold.name, gold.list, pseudo, g.thresh) {
    print(paste(gold.name,'method','TPR','FPR',sep=' '))
    for (method in names(input.list)) {
        if (method != gold.name) {
            data.vec <- c(log2(input.list[[method]]+pseudo))
            gold.vec <- c(log2(gold.list[[gold.name]]+pseudo))
            thresh <- g.thresh
            TPR.vec <- TPR.c(data.vec,gold.vec,thresh,g.thresh,pseudo)
            FPR.vec <- FPR.c(data.vec,gold.vec,thresh,g.thresh,pseudo)
            print(paste(method,round(TPR.vec,2),' & ',round(FPR.vec,2),sep=' '))
        }
    }
}
