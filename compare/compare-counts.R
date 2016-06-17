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

#library('parallel')
library(ggplot2)
library(RColorBrewer)
#library(grid)








###############################################################################
###############################################################################
###############################################################################
######
######  Initialize parameters
######
###############################################################################
###############################################################################
###############################################################################




###############################################################################
## set directories and get arguments
###############################################################################

## print program starting message
message('')
message('###############################################################################')
message('Status: ', 'Script: ', 'compare-counts.R ... ')
message(rep('-', 79))

## get script path and set new working directory
old.dir <- getwd() # at the end set old one
if (length(commandArgs()) < 6) {
    dataset.name <- 'SEQC'
    script.dir <- 'D:/uni/Masterarbeit/statistics_scripts'
    quant.dir <- 'D:/uni/Masterarbeit/quantifications'
    data.dir <- 'D:/uni/Masterarbeit/statistics_scripts'
    output.dir <- 'D:/uni/Masterarbeit/statistics_out'
    gene.lengths.file <- 'ng-human-GRCh38-Ensembl-glen.txt'
    gene.2.trans.file <- 'ng-human-GRCh38-Ensembl-trans2gene.txt'
    projects.gene.file <- paste(paste('make_statistics', dataset.name, sep='_'), 'conf', sep='.')
} else {
    #script.dir <- dirname(normalizePath(strsplit(commandArgs()[4], '=')[[1]][2]))
    #input.dir <- commandArgs()[6]
    #output.dir <- commandArgs()[7]
    #dataset.name <- commandArgs()[8]
}

message('Script directory ... ')
print(script.dir)
message('Quantification directory ... ')
print(quant.dir)
message('Data directory ... ')
print(data.dir)
message('Output directory ... ')
print(output.dir)
message('Dataset name ... ')
print(dataset.name)

# get local libraries
#source(file.path(script.dir, 'exon-util-lib.R'))
#source(file.path(script.dir, 'compare-lib.R'))
source(file.path(script.dir, 'util-lib.R'))

# set parameters
num.threads <- 1
logFC.thresh <- 0.2
VClogFC.print.thresh <- 50000.0
normalized.pseudo_count  <- 0.0001
normalization.scale <- 1000000
normalization.method <- 'DESeq'
#normalization.method <- 'CPM'
save.SEQC.SEQS <- FALSE
load.SEQC.SEQS <- TRUE
SEQCxSEQS <- FALSE

# loading
write.count.matrix.bool <- FALSE
write.normalized.matrix.bool <- FALSE
read.results.bool <- FALSE
read.cnt.matrix <- FALSE
read.norm.cnt.matrix <- FALSE
normalized.results.bool <- FALSE
correct.results.bool <- FALSE
make.biotaq.norm.cnt.list <- FALSE

# modifying
create.summary.bool <- FALSE
create.norm.summary.bool <- FALSE
calc.mean.norm.cnt.list <- FALSE
calc.sd.norm.cnt.list <- FALSE
calc.mean.FC.bool <- FALSE

#used
fold.change.bool <- FALSE
calc.sd.mean.fc.list <- FALSE
stat.correlation.list.bool <- FALSE

# statistics
coefficient.variation.bool <- FALSE
fold.change.cv.bool <- FALSE
intra.method.cnt.cor.bool <- FALSE
intra.method.cnt.cor.thresh <- 6 # threshold for number of samples to compare
intra.method.fc.cor.bool <- FALSE
inter.method.fc.cor.bool <- FALSE
inter.method.cnt.cor.bool <- FALSE

# printing
sample.naming.vector <- c('SEQC-A'='(A)', 'SEQC-B'='(B)')

make.txt.bool <- FALSE
write.stat.bool <- FALSE
intra.method.fc.cor.write.bool <- FALSE
intra.method.cnt.cor.write.bool <- FALSE
inter.method.fc.cor.write.bool <- FALSE
inter.method.cnt.cor.write.bool <- FALSE

make.pdf.bool <- TRUE
make.print.bool <- FALSE
mean.sd.scatter.plot.bool <- FALSE
biotaq.scatter.plot.bool <- FALSE
expression.genes.mean.plot.bool <- FALSE
mean.fc.scatter.plot.bool <- FALSE
cnt.cnt.scatter.plot.bool <- FALSE
fc.fc.scatter.plot.bool <- FALSE
biotaq.fc.scatter.plot.bool <- FALSE
cor.heatmap.plot.bool <- FALSE
FC.gold.plot.bool <- FALSE
boxplot.FC.plot.bool <- FALSE

rocplot.FC.bool <- TRUE
cnt.fc.scatter.plot.bool <- FALSE
FC.genes.mean.plot.bool <- FALSE
intra.fc.box.plot.bool <- FALSE
CV.list.plot.bool <- FALSE
fcCV.list.plot.bool <- FALSE
intra.method.cnt.cor.list.plot.bool <- FALSE
intra.method.fc.cor.list.plot.bool <- FALSE
inter.method.fc.cor.list.plot.bool <- FALSE
inter.method.cnt.cor.list.plot.bool <- FALSE
CV.list.scatter.plot.bool <- FALSE
FC.mean.scatter.plot.bool <- FALSE
FC.mean.scatter.all.plot.bool <- FALSE
mean.sd.fc.scatter.plot.bool <- FALSE


# list for pdf naming
# CV -> Coefficient of Variation (relative standard deviation (RSD))
# FC -> Fold Change
# CC -> Correlation Coefficient
# CD -> Coefficient of Determination ((Pearson CC)^2)
# SG -> Significant Gene Expression Change (Over and Under Expressed Genes) (Log2(FC) >= +/- 1)
plot.titles.list <- c(
    'gene.sample.count.matrix'       =  'Gene Count Matrix', 
    'gene.sample.norm.matrix'        =  'Normalized Gene Count Matrix', 
    'mean.cnt.sum'                   =  'Mean Gene Count Summary', 
    'mean.norm.cnt.sum'              =  'Mean Normalized Gene Count Summary', 
    'mean.taq.norm.cnt.sum'          =  'Mean TaqMan Normalized Gene Count', 
    'mean.bio.norm.cnt.sum'          =  'Mean Bio-rad Normalized Gene Count', 
    'taq.fc.sum'                     =  'TaqMan FC', 
    'bio.fc.sum'                     =  'Bio-rad FC', 
    'mean.all.norm.cnt'              =  'Mean Normalized Gene Count', 
    'cnt.sum'                        =  'Gene Count Summary',
    'norm.cnt.sum'                   =  'Normalized Gene Count Summary',
    'norm.cnt.sd'                    =  'Normalized Gene Count Standard Deviation',
    'norm.cnt.cnt.plot'              =  'Count vs Count',
    'fc.fc.plot'                     =  'FC vs FC',
    'norm.cnt.sim.plot'              =  'Count vs Simulated',
    'mean.sd.norm.plot'              =  'SD vs Mean',
    'mean.fc.plot'                   =  'FC vs Mean',
    'cnt.fc.plot'                    =  'FC vs Count',
    'mean.sd.fc.plot'                =  'SD FC vs Mean FC',
    'mean.log.FC'                    =  'Log2 FC of Mean Gene Counts',
    'box.FC.taq'                     =  'Log2 FC TaqMan',
    'box.FC.bio'                     =  'Log2 FC Bio-Rad',
    'box.FC.SEQC'                    =  'Log2 FC SEQC',
    'box.FC.SEQS'                    =  'Log2 FC SEQS',
    'box.FC.taq.abs'                 =  'ABS Log2 FC TaqMan',
    'box.FC.bio.abs'                 =  'ABS Log2 FC Bio-Rad',
    'box.FC.SEQC.abs'                =  'ABS Log2 FC SEQC',
    'box.FC.SEQS.abs'                =  'ABS Log2 FC SEQS',
    'intra.gene.CV'                  =  'Intra Gene CV', 
    'intra.gene.CV.FC'               =  'Intra Gene CV of FC', 
    'intra.gene.CV.log.FC'           =  'Intra Gene CV of Log2 FC', 
    'intra.gene.FC'                  =  'Intra Gene FC', 
    'intra.gene.log.FC'              =  'Intra Gene Log2 FC', 
    'intra.met.pearson'              =  'Intra Method Pearson CC', 
    'intra.met.spearman'             =  'Intra Method Spearman CC', 
    'intra.met.CD'                   =  'Intra Method CD', 
    'intra.met.pearson.fc'           =  'Intra Method FC Pearson CC', 
    'intra.met.spearman.fc'          =  'Intra Method FC Spearman CC', 
    'intra.met.CD.fc'                =  'Intra Method FC CD', 
    'intra.met.pearson.log.fc'       =  'Intra Method Log2 FC Pearson CC', 
    'intra.met.spearman.log.fc'      =  'Intra Method Log2 FC Spearman CC', 
    'intra.met.CD.log.fc'            =  'Intra Method Log2 FC CD', 
    'inter.met.pearson'              =  'Inter Method Pearson CC', 
    'inter.met.spearman'             =  'Inter Method Spearman CC', 
    'inter.met.CD'                   =  'Inter Method CD', 
    'inter.met.log.pearson'          =  'Inter Method Log2 Pearson CC', 
    'inter.met.log.spearman'         =  'Inter Method Log2 Spearman CC', 
    'inter.met.log.CD'               =  'Inter Method Log2 CD', 
    'inter.met.pearson.fc'           =  'Inter Method FC Pearson CC', 
    'inter.met.spearman.fc'          =  'Inter Method FC Spearman CC', 
    'inter.met.CD.fc'                =  'Inter Method FC CD', 
    'inter.met.pearson.log.fc'       =  'Inter Method Log2 FC Pearson CC', 
    'inter.met.spearman.log.fc'      =  'Inter Method Log2 FC Spearman CC', 
    'inter.met.CD.log.fc'            =  'Inter Method Log2 FC CD', 
    'inter.met.false.negative.rate'  =  'Inter Method SG FNR', 
    'inter.met.sensitivity'          =  'Inter Method SG Sensitivity',
    'inter.met.false.positive.rate'  =  'Inter Method SG FPR', 
    'inter.met.specificity'          =  'Inter Method SG Specificity',  
    'inter.met.false.discovery.rate' =  'Inter Method SG FDR', 
    'inter.met.precision'            =  'Inter Method SG Precision'
)




###############################################################################
## Set project directory and pseudo counts
###############################################################################

#project     <- 'SEQC-NG00008.1-EQP2.0-Ensembl'
#project.dir <- file.path('/dlab', 'ldrive', 'PHCHBS-I21605', 'schuisv1', 'ngs', 'RNA-seq-test', project)
#result.dir  <- file.path(project.dir, 'benchmark-files', 'result-files')
#image.dir   <- file.path(project.dir, 'benchmark-files', 'image-files')

# set main directories -> only use output.dir
image.dir <- output.dir
dataset.dir <- file.path(quant.dir, dataset.name)

# get ref file
samples.file <- file.path(dataset.dir, paste0(dataset.name, '_fastq.ref'))
samples.frame <- read.table(samples.file, colClasses = c('integer', 'character', 'integer'), header = FALSE)
samples.id <- samples.frame[[1]]
message('Sample ids ... ')
print(samples.id)
samples.names <- samples.frame[[2]]
message('Sample names ... ')
print(samples.names)
samples.count <- samples.frame[[3]]
message('Sample count ... ')
print(samples.count)

id.type      <- Sys.getenv('COUNT_TYPE', 'gene')
message('Id type ... ')
print(id.type)
id.type.name <- paste(capitalize(id.type), 'Id')
message('Id name ... ')
print(id.type.name)










###############################################################################
###############################################################################
###############################################################################
######
######  Get data
######
###############################################################################
###############################################################################
###############################################################################




###############################################################################
## Read project descriptions and data
###############################################################################

# get data from projects-gene file
projects.file  <- file.path(data.dir, projects.gene.file)
# assign reading frame
projects.frame <- read.delim(projects.file, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, comment.char='#')

# get gene ids and gene length
gene.data.file  <- file.path(data.dir, gene.lengths.file)
gene.trans.file  <- file.path(data.dir, gene.2.trans.file)
gene.data.frame <- read.delim(gene.data.file, stringsAsFactors=FALSE, fill=FALSE, check.names=FALSE)
gene.trans.frame <- read.delim(gene.trans.file, stringsAsFactors=FALSE, fill=FALSE, check.names=FALSE, header=FALSE)
message('Gene data frame (first 10) ... ')
print(gene.data.frame[1:10,])
names(gene.trans.frame) <- c('Trans Id', 'Gene Id')
message('Gene trans frame(first 10) ... ')
print(gene.trans.frame[1:10,])

# call readProjectCounts to read project quantification data
if (read.results.bool) {
    count.list <- readProjectCounts(dataset.dir, projects.frame, id.type, samples.names, dataset.name, gene.trans.frame)
    # separate count list, aligner list and quantifier list
    c.list.aligners    <- count.list[['Aligners']]
    names(c.list.aligners) <- getNames(c.list.aligners)
    c.list.quantifiers <- count.list[['Quantifiers']]
    names(c.list.quantifiers) <- getNames(c.list.quantifiers)
    c.list.samples    <- count.list[['Samples']]
    c.list.settings    <- count.list[['Settings']]
    c.list.counts  <- count.list[!(names(count.list) %in% c('Aligners', 'Quantifiers', 'Samples', 'Settings'))]
    c.list.length <- length(c.list.counts)
    message('Aligner list ... ')
    print(c.list.aligners)
    message('Quantifier list ... ')
    print(c.list.quantifiers)
    message('Sample list ... ')
    print(c.list.samples)
    message('Setting list ... ')
    print(c.list.settings)
    #message('Status: ', 'Count list ... ')
    #print(c.list.counts)
    message('Length of count list ... ')
    print(c.list.length)

    # round count lists
    message('Round count list ... ')
    for (m in names(c.list.counts)) {
        c.list.counts[[m]][,2:ncol(c.list.counts[[m]])] <- round(c.list.counts[[m]][,2:ncol(c.list.counts[[m]])])
    }
}


###############################################################################
# read count matrix
###############################################################################

if (read.cnt.matrix) {
    # reading count matrices from output dir
    message('Status: ', 'Reading count matrices ... ')
    message(rep('-', 79))
    c.list.counts <- vector('list', nrow(projects.frame))
    n.name.vector <- c()
    for (i in 1:nrow(projects.frame)) {
        # get data from projects-gene.txt
        aligner         <- projects.frame[['Aligner']][i]         # e.g.: HISAT
        quantifier      <- projects.frame[['Quantifier']][i]      # e.g.: EQPQM
        setting         <- projects.frame[['Setting name']][i]    # e.g.: ambig
        combi.string    <- projects.frame[['Combinations']][i]    # e.g.: 1,2,3;2;3;...
        combi.list      <- strsplit(combi.string, ';')[[1]]       # e.g.: '1,2,3' '2' '3'
        matrix.file <- file.path(output.dir, paste(paste(dataset.name, 'Gene_Count_Matrix', paste(quantifier, setting, aligner, sep='-'), sep='_'), 'txt', sep='.'))
        matrix.name <- paste(quantifier, '-', setting, '/', aligner, sep='')
        n.name.vector[[i]] <- matrix.name
        message('Status: ', 'Reading ', matrix.name, ' ... ')
        c.list.counts[[i]] <- read.table(matrix.file, header=TRUE, colClasses=c('character', rep('numeric', length(combi.list))), sep='\t')
        colnames(c.list.counts[[i]])[-1] <- gsub('.', '-', colnames(c.list.counts[[i]])[-1], fixed=TRUE)
        colnames(c.list.counts[[i]])[1] <- gsub('.', ' ', colnames(c.list.counts[[i]])[1], fixed=TRUE)
    }
    names(c.list.counts) <- n.name.vector
    message(rep('#', 79))
}


###############################################################################
# read normalized count matrix
###############################################################################

if (read.norm.cnt.matrix) {
    # reading count matrices from output dir
    message('Status: ', 'Reading normalized count matrices ... ')
    message(rep('-', 79))
    normalized.count.list <- vector('list', nrow(projects.frame))
    n.name.vector <- c()
    c.list.aligners <- c()
    c.list.quantifiers <- c()
    c.list.settings <- c()
    for (i in 1:nrow(projects.frame)) {
        # get data from projects-gene.txt
        aligner         <- projects.frame[['Aligner']][i]         # e.g.: HISAT
        quantifier      <- projects.frame[['Quantifier']][i]      # e.g.: EQPQM
        setting         <- projects.frame[['Setting name']][i]    # e.g.: ambig
        combi.string    <- projects.frame[['Combinations']][i]    # e.g.: 1,2,3;2;3;...
        combi.list      <- strsplit(combi.string, ';')[[1]]       # e.g.: '1,2,3' '2' '3'
        matrix.file <- file.path(output.dir, paste(paste(dataset.name, 'Normalized_Gene_Count_Matrix', paste(quantifier, setting, aligner, sep='-'), sep='_'), 'txt', sep='.'))
        matrix.name <- paste(quantifier, '-', setting, '/', aligner, sep='')
        n.name.vector[[i]] <- matrix.name
        c.list.aligners[[i]] <- aligner
        c.list.quantifiers[[i]] <- paste(quantifier, setting, sep='-')
        c.list.settings[[i]] <- matrix.name
        message('Status: ', 'Reading ', matrix.name, ' ... ')
        normalized.count.list[[i]] <- read.table(matrix.file, header=TRUE, colClasses=c('character', rep('numeric', length(combi.list))), sep='\t')
        colnames(normalized.count.list[[i]])[-1] <- gsub('.', '-', colnames(normalized.count.list[[i]])[-1], fixed=TRUE)
        colnames(normalized.count.list[[i]])[1] <- gsub('.', ' ', colnames(normalized.count.list[[i]])[1], fixed=TRUE)
    }
    names(normalized.count.list) <- n.name.vector
    c.list.length <- length(c.list.settings)
    c.list.samples <- gsub('.', '-', colnames(normalized.count.list[[1]])[-1], fixed=TRUE)
    c.list.aligners <- unique(c.list.aligners)
    c.list.quantifiers <- unique(c.list.quantifiers)
    names(c.list.aligners) <- getNames(c.list.aligners)
    names(c.list.quantifiers) <- getNames(c.list.quantifiers)
    message('Aligner list ... ')
    print(c.list.aligners)
    message('Quantifier list ... ')
    print(c.list.quantifiers)
    message('Sample list ... ')
    print(c.list.samples)
    message('Setting list ... ')
    print(c.list.settings)
    #message('Status: ', 'Count list ... ')
    #print(c.list.counts)
    message('Length of count list ... ')
    print(c.list.length)
    message(rep('#', 79))
}




###############################################################################
## Define samples
###############################################################################

# -> get fastq sample names
#samples <- list.files(file.path(project.dir, 'samples-scratch'), 'SEQC-.*')
if (!load.SEQC.SEQS) {
    samples <- c.list.samples

    sample.types <- c('SEQC-A', 'SEQC-B')
    names(sample.types) <- c('A', 'B')
    # split um A and B samples
    samples.matched.A <- samples[grep(sample.types['A'], samples)]
    samples.matched.B <- samples[grep(sample.types['B'], samples)]
    message('Samples A ... ')
    print(samples.matched.A)
    message('Samples B ... ')
    print(samples.matched.B)
    samples.matched <- list()
    samples.matched[['SEQC-A']] <- samples.matched.A
    samples.matched[['SEQC-B']] <- samples.matched.B

    # get total numbe rof samples
    samples.counts.A <- length(samples.matched.A)
    samples.counts.B <- length(samples.matched.B)
    message('Count samples A ... ')
    print(samples.counts.A)
    message('Count samples B ... ')
    print(samples.counts.B)

    message(rep('#', 79))
}



###############################################################################
## write count matrix
###############################################################################

if (write.count.matrix.bool) {
    message('Status: ', 'Creating count matrix ... ')
    message(rep('-', 79))
    for (method.name in names(c.list.counts)) {
        # prepare data
        print.list <- list()
        for (sample.name in names(c.list.counts[[method.name]])) {
            if (sample.name == 'Gene Id') {
                gene.name.frame <- c.list.counts[[method.name]][[sample.name]]
            } else {
                print.list[[sample.name]] <- c.list.counts[[method.name]][[sample.name]]
                names(print.list[[sample.name]]) <- gene.name.frame
            }
        }
        # write data
        write.matrix.data( input.frame  = print.list, 
                           method.name  = method.name,  
                           table.name   = paste(dataset.name, plot.titles.list['gene.sample.count.matrix'][[1]], gsub('/', '-', method.name), sep=' '), 
                           column.names = 'Gene Id')
    }
    message(rep('#', 79))
}




###############################################################################
## Create summary
###############################################################################

if (create.summary.bool) {
    message('Status: ', 'Creating summary ... ')
    message(rep('-', 79))
    # calculate number of detected genes, totaly expressed genes detected.genes.num
    result.summary.list <- lapply (c.list.counts, computeSummary, sample.types)
    
    # get names of result.summary.list
    summaryFields <- names(result.summary.list[[1]])
    # rearrange result.summary.list into -> field $ sample.type -> e.g.: detected.genes.num $ SEQC-A
    summary.list <- list ()
    for (field in summaryFields) {
        summary.list[[field]] <- list ()
        for (sample.type in sample.types) {
            summary.list[[field]][[sample.type]] <- lapply(result.summary.list, function (x) return (x[[field]][[sample.type]]))
        }
    }
    
    # rearrange result.summary.list into -> method.name $ sample.type $ sample.name (named with summary.type)
    print.list <- list()
    for (method.name in names(result.summary.list)) {
        for (summary.type in names(result.summary.list[[method.name]])) {
            print.list[[method.name]][[summary.type]] <- c(result.summary.list[[method.name]][[summary.type]][[sample.types[[1]]]], 
                                                           result.summary.list[[method.name]][[summary.type]][[sample.types[[2]]]])
            names(print.list[[method.name]][[summary.type]]) <- unname(unlist(samples.matched))
        }
        names(print.list[[method.name]]) <- c('Gene Count', 'Total Expression', 'Mean Expression')
        # write data
        message('Making summary matrix ... ')
        write.matrix.data( input.frame  = print.list[[method.name]], 
                           method.name  = method.name,  
                           table.name   = paste(dataset.name, plot.titles.list['cnt.sum'][[1]], 'for', gsub('/', '-', method.name)), 
                           column.names = 'Sample Name')
    }
    message(rep('#', 79))
}




###############################################################################
## normalize data
###############################################################################

if (normalized.results.bool) {
    # extract cufflinks
    #cuffl.cnt.list <- list()
    #cuffl.pos <- grep('CUFFL', names(c.list.counts))
    #cuffl.cnt.list[[names(c.list.counts)[cuffl.pos]]] <- c.list.counts[[cuffl.pos]]
    #c.list.counts <- c.list.counts[[!cuffl.pos]]
    message('Status: ', 'Normalizing with ', normalization.method, ' ... ')
    message(rep('-', 79))
    #single core lapply call
    normalized.count.list <- vector('list', length(c.list.counts))
    for (i in 1:length(c.list.counts)) {
        ##normalized.count.list[[i]] <- data.frame(c.list.counts[[i]][['Gene Id']])
        ##names(normalized.count.list[[i]]) <- 'Gene Id'
        # for each quantifier - alignment combination
        ##for (seqc.name in c('SEQC-A', 'SEQC-B')) {
            ##seqc.pos <- c(1,grep(seqc.name, names(c.list.counts[[i]])))
            ##tmp.norm.cnt.frame <- normalize.counts( count.frame          = c.list.counts[[i]][,seqc.pos],
        normalized.count.list[[i]] <- normalize.counts( count.frame          = c.list.counts[[i]], 
                                                        count.name           = names(c.list.counts)[i], 
                                                        gene.ids.frame       = gene.data.frame['Gene Id'], 
                                                        normalization.method = normalization.method, 
                                                        normalization.scale  = normalization.scale )
        ##     normalized.count.list[[i]] <- merge(normalized.count.list[[i]], tmp.norm.cnt.frame, by = c('Gene Id'), all=TRUE)
        ##}
        ##normalized.count.list[[i]] <- normalized.count.list[[i]][,order(order(names(c.list.counts[[i]])))]
    }
    names(normalized.count.list) <- names(c.list.counts)
    message(rep('#', 79))
}




###############################################################################
## correct data
###############################################################################

if (correct.results.bool) {
    message('Status: ', 'Correcting ', 'results ... ')
    message(rep('-', 79))
    corrected.count.list <- vector('list', length(normalized.count.list))
    for (i in 1:length(normalized.count.list)) {
        corrected.count.list[[i]] <- correctResults( input.frame = normalized.count.list[[i]], 
                                                     input.name  = names(normalized.count.list)[i] )
    }

    normalized.count.list <- corrected.count.list
    names(normalized.count.list) <- names(c.list.counts)
    message(rep('#', 79))
}



###############################################################################
## Compare to Taqman and Biorad data
###############################################################################

# TODO:
#source(file.path(Sys.getenv('HOME'), 'ngs', 'pipelines', 'exon-pipeline', 'R', 'benchmark-files', 'compare-counts-taqman.R'))
#source(file.path(Sys.getenv('HOME'), 'ngs', 'pipelines', 'exon-pipeline', 'R', 'benchmark-files', 'compare-counts-biorad.R'))




###############################################################################
## write normalized matrix
###############################################################################

if (write.normalized.matrix.bool) {
    message('Status: ', 'Creating normalized count matrix ... ')
    message(rep('-', 79))
    for (method.name in names(normalized.count.list)) {
        # prepare data
        print.list <- list()
        for (sample.name in names(normalized.count.list[[method.name]])) {
            if (sample.name == 'Gene Id') {
                gene.name.frame <- normalized.count.list[[method.name]][[sample.name]]
            } else {
                print.list[[sample.name]] <- normalized.count.list[[method.name]][[sample.name]]
                names(print.list[[sample.name]]) <- gene.name.frame
            }
        }
        # write data
        write.matrix.data( input.frame  = print.list, 
                           method.name  = method.name,  
                           table.name   = paste(dataset.name, plot.titles.list['gene.sample.norm.matrix'][[1]], gsub('/', '-', method.name), sep=' '), 
                           column.names = 'Gene Id')
    }
    message(rep('#', 79))
}




###############################################################################
## Create Normalized Gene Count Summary
###############################################################################

if (create.norm.summary.bool) {
    message('Status: ', 'Creating normalized summary ... ')
    message(rep('-', 79))
    # calculate number of detected genes, totaly expressed genes detected.genes.num
    norm.result.summary.list <- lapply (normalized.count.list, computeSummary, sample.types)
    
    # get names of norm.result.summary.list
    normSummaryFields <- names(norm.result.summary.list[[1]])
    # rearrange norm.result.summary.list into -> field $ sample.type -> e.g.: detected.genes.num $ SEQC-A
    norm.summary.list <- list ()
    for (field in normSummaryFields) {
        norm.summary.list[[field]] <- list ()
        for (sample.type in sample.types) {
            norm.summary.list[[field]][[sample.type]] <- lapply(norm.result.summary.list, function (x) return (x[[field]][[sample.type]]))
        }
    }
    
    # rearrange norm.result.summary.list into -> method.name $ sample.type $ sample.name (named with summary.type)
    print.list <- list()
    for (method.name in names(norm.result.summary.list)) {
        for (summary.type in names(norm.result.summary.list[[method.name]])) {
            print.list[[method.name]][[summary.type]] <- c(norm.result.summary.list[[method.name]][[summary.type]][[sample.types[[1]]]], 
                                                           norm.result.summary.list[[method.name]][[summary.type]][[sample.types[[2]]]])
            names(print.list[[method.name]][[summary.type]]) <- unname(unlist(samples.matched))
        }
        names(print.list[[method.name]]) <- c('Gene Count', 'Total Expression', 'Mean Expression')
        # write data
        message('Making normalized summary matrix ... ')
        write.matrix.data( input.frame  = print.list[[method.name]], 
                           method.name  = method.name,  
                           table.name   = paste(dataset.name, plot.titles.list['norm.cnt.sum'][[1]], 'for', gsub('/', '-', method.name)), 
                           column.names = 'Sample Name')
    }
    message(rep('#', 79))
}




###############################################################################
# calculate mean normalized count list over all samples A and B
###############################################################################

if (calc.mean.norm.cnt.list) {
    # calculate mean for each gene
    message('Status: ', 'Creating mean normalized count list ... ')
    message(rep('-', 79))
    mean.all.norm.cnt.list <- list()
    mean.inner.norm.cnt.list <- list()
    for (method.n in names(normalized.count.list)) {
        message('Computing mean count for method ', method.n,  ' ... ')
        # calculate overall mean
        SEQC.A.samples <- apply(normalized.count.list[[method.n]][,samples.matched.A], 1, function(x) return(mean(x)))
        SEQC.B.samples <- apply(normalized.count.list[[method.n]][,samples.matched.B], 1, function(x) return(mean(x)))
        gene.id.frame <- normalized.count.list[[method.n]][,'Gene Id']
        mean.all.norm.cnt.list[[method.n]] <- data.frame(gene.id.frame, SEQC.A.samples, SEQC.B.samples)
        names(mean.all.norm.cnt.list[[method.n]]) <- as.factor(c('Gene Id', 'SEQC-A-samples', 'SEQC-B-samples'))
        # calculate mean of inner 50% (remove 25% top and bottom values)
        lenA <- length(samples.matched.A)
        lenB <- length(samples.matched.B)
        posAs <- (lenA/4)+1
        posAe <- (lenA/4)*3
        posBs <- (lenB/4)+1
        posBe <- (lenB/4)*3
        SEQC.A.samples <- apply(normalized.count.list[[method.n]][,samples.matched.A], 1, function(x) return(mean(sort(x)[posAs:posAe])))
        SEQC.B.samples <- apply(normalized.count.list[[method.n]][,samples.matched.B], 1, function(x) return(mean(sort(x)[posBs:posBe])))
        gene.id.frame <- normalized.count.list[[method.n]][,'Gene Id']
        mean.inner.norm.cnt.list[[method.n]] <- data.frame(gene.id.frame, SEQC.A.samples, SEQC.B.samples)
        names(mean.inner.norm.cnt.list[[method.n]]) <- as.factor(c('Gene Id', 'SEQC-A-samples', 'SEQC-B-samples'))
    }
    # write data
    # message('Making inner mean summary matrix ... ')
    # write.matrix.data( input.frame  = mean.inner.norm.cnt.list[[method.name]], 
                       # method.name  = method.name,  
                       # table.name   = paste(plot.titles.list['norm.cnt.sum'][[1]], 'for', gsub('/', '-', method.name)), 
                       # column.names = 'Sample Name')
    message(rep('#', 79))
}



###############################################################################
# calculate standard deviation A and B
###############################################################################

if (calc.sd.norm.cnt.list) {
    # calculate mean for each gene
    message('Status: ', 'Creating standard deviation normalized list ... ')
    message(rep('-', 79))
    sd.norm.cnt.list <- list()
    for (method.n in names(normalized.count.list)) {
        message('Computing standard deviation for method ', method.n,  ' ... ')
        # calculate overall mean
        SEQC.A.samples <- apply(normalized.count.list[[method.n]][,samples.matched.A], 1, function(x) return(sd(x)))
        SEQC.B.samples <- apply(normalized.count.list[[method.n]][,samples.matched.B], 1, function(x) return(sd(x)))
        gene.id.frame <- normalized.count.list[[method.n]][,'Gene Id']
        sd.norm.cnt.list[[method.n]] <- data.frame(gene.id.frame, SEQC.A.samples, SEQC.B.samples)
        names(sd.norm.cnt.list[[method.n]]) <- as.factor(c('Gene Id', 'SEQC-A-samples', 'SEQC-B-samples'))
    }
    # write data
    # message('Making inner mean summary matrix ... ')
    # write.matrix.data( input.frame  = sd.norm.cnt.list[[method.name]], 
                       # method.name  = method.name,  
                       # table.name   = paste(plot.titles.list['norm.cnt.sum'][[1]], 'for', gsub('/', '-', method.name)), 
                       # column.names = 'Sample Name')
    message(rep('#', 79))
}




###############################################################################
# make taqman and biorad data
###############################################################################

if (make.biotaq.norm.cnt.list) {
    message('Status: ', 'Reading Bio-rad and TaqMan data ... ')
    message(rep('-', 79))
    biorad.norm.cnt.list <- makeBioTaqCntList(mean.cnt.list = mean.all.norm.cnt.list, qpcr.method  = 'BIORD', out.dir = output.dir, gene.id.list = gene.id.frame)
    taqman.norm.cnt.list <- makeBioTaqCntList(mean.cnt.list = mean.all.norm.cnt.list, qpcr.method  = 'TQMAN', out.dir = output.dir, gene.id.list = gene.id.frame)
    # FC
    sample.type.A <- 'SEQC-A-samples'
    sample.type.B <- 'SEQC-B-samples'
    biorad.FC.list <- list()
    for (method.n in names(biorad.norm.cnt.list)) {
        FC.matrix <- computeFCs( input.frame       = biorad.norm.cnt.list[[method.n]], 
                                 input.name        = names(biorad.norm.cnt.list)[method.n], 
                                 samples.matched.A = c(sample.type.A), 
                                 samples.matched.B = c(sample.type.B), 
                                 pseudo_count      = normalized.pseudo_count )
        biorad.FC.list[[method.n]] <- FC.matrix
    }
    taqman.FC.list <- list()
    for (method.n in names(taqman.norm.cnt.list)) {
        FC.matrix <- computeFCs( input.frame       = taqman.norm.cnt.list[[method.n]], 
                                 input.name        = names(taqman.norm.cnt.list)[method.n], 
                                 samples.matched.A = c(sample.type.A), 
                                 samples.matched.B = c(sample.type.B), 
                                 pseudo_count      = normalized.pseudo_count )
        taqman.FC.list[[method.n]] <- FC.matrix
    }
    message(rep('#', 79))
}




###############################################################################
## Compute FC of mean over A and B samples
###############################################################################

if (calc.mean.FC.bool) {
    mean.normalized.list <- list()
    for (method.n in names(normalized.count.list)) {
        SEQC.A.samples <- apply(normalized.count.list[[method.n]][,samples.matched.A], 1, function(x) return(sum(x)/length(x)))
        SEQC.B.samples <- apply(normalized.count.list[[method.n]][,samples.matched.B], 1, function(x) return(sum(x)/length(x)))
        gene.id.frame <- normalized.count.list[[method.n]][,'Gene Id']
        mean.normalized.list[[method.n]] <- data.frame(gene.id.frame, SEQC.A.samples, SEQC.B.samples)
        names(mean.normalized.list[[method.n]]) <- as.factor(c('Gene Id', 'SEQC-A-samples', 'SEQC-B-samples'))
    }
    sample.type.A <- 'SEQC-A-samples'
    sample.type.B <- 'SEQC-B-samples'
    mean.FC.list <- list()
    for (method.n in names(mean.normalized.list)) {
        FC.matrix <- computeFCs( input.frame       = mean.normalized.list[[method.n]], 
                                 input.name        = names(mean.normalized.list)[method.n], 
                                 samples.matched.A = c(sample.type.A), 
                                 samples.matched.B = c(sample.type.B), 
                                 pseudo_count      = normalized.pseudo_count )
        mean.FC.list[[method.n]] <- FC.matrix
    }
}




###############################################################################
# manage SEQC x SEQS
###############################################################################

if (save.SEQC.SEQS) {
    message('Status: ', 'Managing SEQC and SEQS ... ')
    message(rep('-', 79))
    if (dataset.name == 'SEQC') {
        # info data
        SEQC.c.list.length <- c.list.length
        SEQC.c.list.samples <- c.list.samples
        SEQC.c.list.aligners <- c.list.aligners
        SEQC.c.list.quantifiers <- c.list.quantifiers
        # count lists
        SEQC.c.list.counts <- c.list.counts
        SEQC.normalized.count.list <- normalized.count.list
        # summary lists
        SEQC.summary.list <- summary.list
        SEQC.norm.summary.list <- norm.summary.list
        # mean and sd lists
        SEQC.mean.all.norm.cnt.list <- mean.all.norm.cnt.list
        SEQC.sd.norm.cnt.list <- sd.norm.cnt.list
        SEQC.mean.FC.list <- mean.FC.list
    }
    if (dataset.name == 'SEQS') {
        # info data
        SEQS.c.list.length <- c.list.length
        SEQS.c.list.samples <- c.list.samples
        SEQS.c.list.aligners <- c.list.aligners
        SEQS.c.list.quantifiers <- c.list.quantifiers
        # count lists
        SEQS.c.list.counts <- c.list.counts
        SEQS.normalized.count.list <- normalized.count.list
        # summary lists
        SEQS.summary.list <- summary.list
        SEQS.norm.summary.list <- norm.summary.list
        # mean and sd lists
        SEQS.mean.all.norm.cnt.list <- mean.all.norm.cnt.list
        SEQS.sd.norm.cnt.list <- sd.norm.cnt.list
        SEQS.mean.FC.list <- mean.FC.list
    } 
    message(rep('#', 79))
}




###############################################################################
# LOAD SEQC x SEQS
###############################################################################

if (load.SEQC.SEQS) {
    message('Status: ', 'Loading SEQC and SEQS ... ')
    message(rep('-', 79))
    if (dataset.name == 'SEQC') {
        # info data
        c.list.length <- SEQC.c.list.length
        c.list.samples <- SEQC.c.list.samples
        c.list.aligners <- SEQC.c.list.aligners
        c.list.quantifiers <- SEQC.c.list.quantifiers
        # count lists
        c.list.counts <- SEQC.c.list.counts
        normalized.count.list <- SEQC.normalized.count.list
        # summary lists
        summary.list <- SEQC.summary.list
        norm.summary.list <- SEQC.norm.summary.list
        # mean and sd lists
        mean.all.norm.cnt.list <- SEQC.mean.all.norm.cnt.list
        sd.norm.cnt.list <- SEQC.sd.norm.cnt.list
        mean.FC.list <- SEQC.mean.FC.list
    }
    if (dataset.name == 'SEQS') {
        # info data
        c.list.length <- SEQS.c.list.length
        c.list.samples <- SEQS.c.list.samples
        c.list.aligners <- SEQS.c.list.aligners
        c.list.quantifiers <- SEQS.c.list.quantifiers
        # count lists
        c.list.counts <- SEQS.c.list.counts
        normalized.count.list <- SEQS.normalized.count.list
        # summary lists
        summary.list <- SEQS.summary.list
        norm.summary.list <- SEQS.norm.summary.list
        # mean and sd lists
        mean.all.norm.cnt.list <- SEQS.mean.all.norm.cnt.list
        sd.norm.cnt.list <- SEQS.sd.norm.cnt.list
        mean.FC.list <- SEQS.mean.FC.list
    } 
    samples <- c.list.samples

    sample.types <- c('SEQC-A', 'SEQC-B')
    names(sample.types) <- c('A', 'B')
    # split um A and B samples
    samples.matched.A <- samples[grep(sample.types['A'], samples)]
    samples.matched.B <- samples[grep(sample.types['B'], samples)]
    message('Samples A ... ')
    print(samples.matched.A)
    message('Samples B ... ')
    print(samples.matched.B)
    samples.matched <- list()
    samples.matched[['SEQC-A']] <- samples.matched.A
    samples.matched[['SEQC-B']] <- samples.matched.B

    # get total numbe rof samples
    samples.counts.A <- length(samples.matched.A)
    samples.counts.B <- length(samples.matched.B)
    message('Count samples A ... ')
    print(samples.counts.A)
    message('Count samples B ... ')
    print(samples.counts.B)
    
    message(rep('#', 79))
}


## write statistics
if (write.stat.bool) {
    # lists
    stat.name.list <- list()
    stat.dataset.list <- list()
    # stats
    stat.name.list[['SEQC']]   <- c('Mean Count A', 'Mean Count B', 'Mean Log2 FC', 'SD Mean A', 'SD Mean B', 'SD Log2 FC', 'Mean SD A', 'Mean SD B', 'Mean CV A', 'Mean CV B')
    stat.name.list[['SEQS']]   <- c('Mean Count A', 'Mean Count B', 'Mean Log2 FC', 'SD Mean A', 'SD Mean B', 'SD Log2 FC', 'Mean SD A', 'Mean SD B', 'Mean CV A', 'Mean CV B')
    stat.name.list[['taqman']] <- c('Mean Count A', 'Mean Count B', 'Mean Log2 FC', 'SD Mean A', 'SD Mean B', 'SD Log2 FC')
    stat.name.list[['biorad']] <- c('Mean Count A', 'Mean Count B', 'Mean Log2 FC', 'SD Mean A', 'SD Mean B', 'SD Log2 FC')
    # comparisons
    stat.name.list[['SEQCxTQMAN']] <- c('Pearson A', 'Pearson B', 'Spearman A', 'Spearman B', 'ABS Log2 SD/SD A', 'ABS Log2 SD/SD B', 'Pearson Log2 FC', 'Spearman Log2 FC', 'ABS Log2 SD/SD FC')
    stat.name.list[['SEQCxBIORD']] <- c('Pearson A', 'Pearson B', 'Spearman A', 'Spearman B', 'ABS Log2 SD/SD A', 'ABS Log2 SD/SD B', 'Pearson Log2 FC', 'Spearman Log2 FC', 'ABS Log2 SD/SD FC')
    stat.name.list[['SEQSxFLXSM']] <- c('Pearson A', 'Pearson B', 'Spearman A', 'Spearman B', 'ABS Log2 SD/SD A', 'ABS Log2 SD/SD B', 'Pearson Log2 FC', 'Spearman Log2 FC', 'ABS Log2 SD/SD FC')
    stat.name.list[['SEQCxSEQS']]  <- c('Pearson A', 'Pearson B', 'Spearman A', 'Spearman B', 'ABS Log2 SD/SD A', 'ABS Log2 SD/SD B', 'Pearson Log2 FC', 'Spearman Log2 FC', 'ABS Log2 SD/SD FC')
    stat.name.list[['SEQCxFLXSM']] <- c('Pearson A', 'Pearson B', 'Spearman A', 'Spearman B', 'ABS Log2 SD/SD A', 'ABS Log2 SD/SD B', 'Pearson Log2 FC', 'Spearman Log2 FC', 'ABS Log2 SD/SD FC')
    # data SEQC.mean.all.norm.cnt.list
    stat.data.list <- list()
    for (stat.name in names(stat.name.list)) {
        stat.data.list[[stat.name]] <- list()
        for (x.name in stat.name.list[[stat.name]]) {
            stat.data.list[[stat.name]][[x.name]] <- list()
        }
    }
    for (method in names(SEQC.mean.all.norm.cnt.list)) {
        stat.data.list[['SEQC']][['Mean Count A']][[method]] <- mean(SEQC.mean.all.norm.cnt.list[[method]][['SEQC-A-samples']])
        stat.data.list[['SEQC']][['Mean Count B']][[method]] <- mean(SEQC.mean.all.norm.cnt.list[[method]][['SEQC-B-samples']])
        stat.data.list[['SEQC']][['Mean Log2 FC']][[method]] <- mean(SEQC.mean.FC.list[[method]])
        stat.data.list[['SEQC']][['SD Mean A']][[method]] <- sd(SEQC.mean.all.norm.cnt.list[[method]][['SEQC-A-samples']])
        stat.data.list[['SEQC']][['SD Mean B']][[method]] <- sd(SEQC.mean.all.norm.cnt.list[[method]][['SEQC-B-samples']])
        stat.data.list[['SEQC']][['SD Log2 FC']][[method]] <- sd(SEQC.mean.FC.list[[method]])
        stat.data.list[['SEQC']][['Mean SD A']][[method]] <- mean(SEQC.sd.norm.cnt.list[[method]][['SEQC-A-samples']])
        stat.data.list[['SEQC']][['Mean SD B']][[method]] <- mean(SEQC.sd.norm.cnt.list[[method]][['SEQC-B-samples']])
        stat.data.list[['SEQC']][['Mean CV A']][[method]] <- mean((SEQC.sd.norm.cnt.list[[method]][['SEQC-A-samples']]+normalized.pseudo_count)/(SEQC.mean.all.norm.cnt.list[[method]][['SEQC-A-samples']]+normalized.pseudo_count))
        stat.data.list[['SEQC']][['Mean CV B']][[method]] <- mean((SEQC.sd.norm.cnt.list[[method]][['SEQC-B-samples']]+normalized.pseudo_count)/(SEQC.mean.all.norm.cnt.list[[method]][['SEQC-B-samples']]+normalized.pseudo_count))
    }
    for (method in names(SEQS.mean.all.norm.cnt.list)) {
        stat.data.list[['SEQS']][['Mean Count A']][[method]] <- mean(SEQS.mean.all.norm.cnt.list[[method]][['SEQC-A-samples']])
        stat.data.list[['SEQS']][['Mean Count B']][[method]] <- mean(SEQS.mean.all.norm.cnt.list[[method]][['SEQC-B-samples']])
        stat.data.list[['SEQS']][['Mean Log2 FC']][[method]] <- mean(SEQS.mean.FC.list[[method]])
        stat.data.list[['SEQS']][['SD Mean A']][[method]] <- sd(SEQS.mean.all.norm.cnt.list[[method]][['SEQC-A-samples']])
        stat.data.list[['SEQS']][['SD Mean B']][[method]] <- sd(SEQS.mean.all.norm.cnt.list[[method]][['SEQC-B-samples']])
        stat.data.list[['SEQS']][['SD Log2 FC']][[method]] <- sd(SEQS.mean.FC.list[[method]])
        stat.data.list[['SEQS']][['Mean SD A']][[method]] <- mean(SEQS.sd.norm.cnt.list[[method]][['SEQC-A-samples']])
        stat.data.list[['SEQS']][['Mean SD B']][[method]] <- mean(SEQS.sd.norm.cnt.list[[method]][['SEQC-B-samples']])
        stat.data.list[['SEQS']][['Mean CV A']][[method]] <- mean((SEQS.sd.norm.cnt.list[[method]][['SEQC-A-samples']]+normalized.pseudo_count)/(SEQS.mean.all.norm.cnt.list[[method]][['SEQC-A-samples']]+normalized.pseudo_count))
        stat.data.list[['SEQS']][['Mean CV B']][[method]] <- mean((SEQS.sd.norm.cnt.list[[method]][['SEQC-B-samples']]+normalized.pseudo_count)/(SEQS.mean.all.norm.cnt.list[[method]][['SEQC-B-samples']]+normalized.pseudo_count))
    }
    for (method in names(taqman.norm.cnt.list)) {
        stat.data.list[['taqman']][['Mean Count A']][[method]] <- mean(taqman.norm.cnt.list[[method]][['SEQC-A-samples']])
        stat.data.list[['taqman']][['Mean Count B']][[method]] <- mean(taqman.norm.cnt.list[[method]][['SEQC-B-samples']])
        stat.data.list[['taqman']][['Mean Log2 FC']][[method]] <- mean(taqman.FC.list[[method]])
        stat.data.list[['taqman']][['SD Mean A']][[method]] <- sd(taqman.norm.cnt.list[[method]][['SEQC-A-samples']])
        stat.data.list[['taqman']][['SD Mean B']][[method]] <- sd(taqman.norm.cnt.list[[method]][['SEQC-B-samples']])
        stat.data.list[['taqman']][['SD Log2 FC']][[method]] <- sd(taqman.FC.list[[method]])
    }
    for (method in names(biorad.norm.cnt.list)) {
        stat.data.list[['biorad']][['Mean Count A']][[method]] <- mean(biorad.norm.cnt.list[[method]][['SEQC-A-samples']])
        stat.data.list[['biorad']][['Mean Count B']][[method]] <- mean(biorad.norm.cnt.list[[method]][['SEQC-B-samples']])
        stat.data.list[['biorad']][['Mean Log2 FC']][[method]] <- mean(biorad.FC.list[[method]])
        stat.data.list[['biorad']][['SD Mean A']][[method]] <- sd(biorad.norm.cnt.list[[method]][['SEQC-A-samples']])
        stat.data.list[['biorad']][['SD Mean B']][[method]] <- sd(biorad.norm.cnt.list[[method]][['SEQC-B-samples']])
        stat.data.list[['biorad']][['SD Log2 FC']][[method]] <- sd(biorad.FC.list[[method]])
    }
    for (method in names(SEQC.mean.all.norm.cnt.list)) {
        stat.data.list[['SEQCxTQMAN']][['Pearson A']][[method]] <- stat.correlation.list[['taq.cnt.A.pearson']]['TQMAN-default/QPCR',method]
        stat.data.list[['SEQCxTQMAN']][['Pearson B']][[method]] <- stat.correlation.list[['taq.cnt.B.pearson']]['TQMAN-default/QPCR',method]
        stat.data.list[['SEQCxTQMAN']][['Spearman A']][[method]] <- stat.correlation.list[['taq.cnt.A.spearman']]['TQMAN-default/QPCR',method]
        stat.data.list[['SEQCxTQMAN']][['Spearman B']][[method]] <- stat.correlation.list[['taq.cnt.B.spearman']]['TQMAN-default/QPCR',method]
        stat.data.list[['SEQCxTQMAN']][['ABS Log2 SD/SD A']][[method]] <- stat.correlation.list[['taq.cnt.A.sd']]['TQMAN-default/QPCR',method]
        stat.data.list[['SEQCxTQMAN']][['ABS Log2 SD/SD B']][[method]] <- stat.correlation.list[['taq.cnt.B.sd']]['TQMAN-default/QPCR',method]
        stat.data.list[['SEQCxTQMAN']][['Pearson Log2 FC']][[method]] <- stat.correlation.list[['taq.log.fc.pearson']]['TQMAN-default/QPCR',method]
        stat.data.list[['SEQCxTQMAN']][['Spearman Log2 FC']][[method]] <- stat.correlation.list[['taq.log.fc.spearman']]['TQMAN-default/QPCR',method]
        stat.data.list[['SEQCxTQMAN']][['ABS Log2 SD/SD FC']][[method]] <- stat.correlation.list[['taq.log.fc.sd']]['TQMAN-default/QPCR',method]
    }
    for (method in names(SEQC.mean.all.norm.cnt.list)) {
        stat.data.list[['SEQCxBIORD']][['Pearson A']][[method]] <- stat.correlation.list[['bio.cnt.A.pearson']]['BIORD-default/QPCR',method]
        stat.data.list[['SEQCxBIORD']][['Pearson B']][[method]] <- stat.correlation.list[['bio.cnt.B.pearson']]['BIORD-default/QPCR',method]
        stat.data.list[['SEQCxBIORD']][['Spearman A']][[method]] <- stat.correlation.list[['bio.cnt.A.spearman']]['BIORD-default/QPCR',method]
        stat.data.list[['SEQCxBIORD']][['Spearman B']][[method]] <- stat.correlation.list[['bio.cnt.B.spearman']]['BIORD-default/QPCR',method]
        stat.data.list[['SEQCxBIORD']][['ABS Log2 SD/SD A']][[method]] <- stat.correlation.list[['bio.cnt.A.sd']]['BIORD-default/QPCR',method]
        stat.data.list[['SEQCxBIORD']][['ABS Log2 SD/SD B']][[method]] <- stat.correlation.list[['bio.cnt.B.sd']]['BIORD-default/QPCR',method]
        stat.data.list[['SEQCxBIORD']][['Pearson Log2 FC']][[method]] <- stat.correlation.list[['bio.log.fc.pearson']]['BIORD-default/QPCR',method]
        stat.data.list[['SEQCxBIORD']][['Spearman Log2 FC']][[method]] <- stat.correlation.list[['bio.log.fc.spearman']]['BIORD-default/QPCR',method]
        stat.data.list[['SEQCxBIORD']][['ABS Log2 SD/SD FC']][[method]] <- stat.correlation.list[['bio.log.fc.sd']]['BIORD-default/QPCR',method]
    }
    for (method in names(SEQC.mean.all.norm.cnt.list)) {
        stat.data.list[['SEQSxFLXSM']][['Pearson A']][[method]] <- stat.correlation.list[['SEQS.cnt.A.pearson']]['FLXSM-default/SIMED',method]
        stat.data.list[['SEQSxFLXSM']][['Pearson B']][[method]] <- stat.correlation.list[['SEQS.cnt.B.pearson']]['FLXSM-default/SIMED',method]
        stat.data.list[['SEQSxFLXSM']][['Spearman A']][[method]] <- stat.correlation.list[['SEQS.cnt.A.spearman']]['FLXSM-default/SIMED',method]
        stat.data.list[['SEQSxFLXSM']][['Spearman B']][[method]] <- stat.correlation.list[['SEQS.cnt.B.spearman']]['FLXSM-default/SIMED',method]
        stat.data.list[['SEQSxFLXSM']][['ABS Log2 SD/SD A']][[method]] <- stat.correlation.list[['SEQS.cnt.A.sd']]['FLXSM-default/SIMED',method]
        stat.data.list[['SEQSxFLXSM']][['ABS Log2 SD/SD B']][[method]] <- stat.correlation.list[['SEQS.cnt.B.sd']]['FLXSM-default/SIMED',method]
        stat.data.list[['SEQSxFLXSM']][['Pearson Log2 FC']][[method]] <- stat.correlation.list[['SEQS.log.fc.pearson']]['FLXSM-default/SIMED',method]
        stat.data.list[['SEQSxFLXSM']][['Spearman Log2 FC']][[method]] <- stat.correlation.list[['SEQS.log.fc.spearman']]['FLXSM-default/SIMED',method]
        stat.data.list[['SEQSxFLXSM']][['ABS Log2 SD/SD FC']][[method]] <- stat.correlation.list[['SEQS.log.fc.sd']]['FLXSM-default/SIMED',method]
    }
    for (method in names(SEQC.mean.all.norm.cnt.list)) {
        stat.data.list[['SEQCxSEQS']][['Pearson A']][[method]] <- stat.correlation.list[['SEQC.SEQS.A.cnt.pearson']][method,method]
        stat.data.list[['SEQCxSEQS']][['Pearson B']][[method]] <- stat.correlation.list[['SEQC.SEQS.B.cnt.pearson']][method,method]
        stat.data.list[['SEQCxSEQS']][['Spearman A']][[method]] <- stat.correlation.list[['SEQC.SEQS.A.cnt.spearman']][method,method]
        stat.data.list[['SEQCxSEQS']][['Spearman B']][[method]] <- stat.correlation.list[['SEQC.SEQS.B.cnt.spearman']][method,method]
        stat.data.list[['SEQCxSEQS']][['ABS Log2 SD/SD A']][[method]] <- stat.correlation.list[['SEQC.SEQS.A.cnt.sd']][method,method]
        stat.data.list[['SEQCxSEQS']][['ABS Log2 SD/SD B']][[method]] <- stat.correlation.list[['SEQC.SEQS.B.cnt.sd']][method,method]
        stat.data.list[['SEQCxSEQS']][['Pearson Log2 FC']][[method]] <- stat.correlation.list[['SEQC.SEQS.log.fc.pearson']][method,method]
        stat.data.list[['SEQCxSEQS']][['Spearman Log2 FC']][[method]] <- stat.correlation.list[['SEQC.SEQS.log.fc.spearman']][method,method]
        stat.data.list[['SEQCxSEQS']][['ABS Log2 SD/SD FC']][[method]] <- stat.correlation.list[['SEQC.SEQS.log.fc.sd']][method,method]
    }
    for (method in names(SEQC.mean.all.norm.cnt.list)) {
        stat.data.list[['SEQCxFLXSM']][['Pearson A']][[method]] <- stat.correlation.list[['SEQC.SEQS.A.cnt.pearson']][method,'FLXSM-default/SIMED']
        stat.data.list[['SEQCxFLXSM']][['Pearson B']][[method]] <- stat.correlation.list[['SEQC.SEQS.B.cnt.pearson']][method,'FLXSM-default/SIMED']
        stat.data.list[['SEQCxFLXSM']][['Spearman A']][[method]] <- stat.correlation.list[['SEQC.SEQS.A.cnt.spearman']][method,'FLXSM-default/SIMED']
        stat.data.list[['SEQCxFLXSM']][['Spearman B']][[method]] <- stat.correlation.list[['SEQC.SEQS.B.cnt.spearman']][method,'FLXSM-default/SIMED']
        stat.data.list[['SEQCxFLXSM']][['ABS Log2 SD/SD A']][[method]] <- stat.correlation.list[['SEQC.SEQS.A.cnt.sd']][method,'FLXSM-default/SIMED']
        stat.data.list[['SEQCxFLXSM']][['ABS Log2 SD/SD B']][[method]] <- stat.correlation.list[['SEQC.SEQS.B.cnt.sd']][method,'FLXSM-default/SIMED']
        stat.data.list[['SEQCxFLXSM']][['Pearson Log2 FC']][[method]] <- stat.correlation.list[['SEQC.SEQS.log.fc.pearson']][method,'FLXSM-default/SIMED']
        stat.data.list[['SEQCxFLXSM']][['Spearman Log2 FC']][[method]] <- stat.correlation.list[['SEQC.SEQS.log.fc.spearman']][method,'FLXSM-default/SIMED']
        stat.data.list[['SEQCxFLXSM']][['ABS Log2 SD/SD FC']][[method]] <- stat.correlation.list[['SEQC.SEQS.log.fc.sd']][method,'FLXSM-default/SIMED']
    }
    
    for (stat.name in names(stat.name.list)) {
        matrix.frame <- data.frame(names(stat.data.list[[stat.name]][[1]]))
        for (x.name in names(stat.data.list[[stat.name]])) {
            m.frame <- data.frame(round(unname(unlist(stat.data.list[[stat.name]][[x.name]])),2))
            matrix.frame <- cbind(matrix.frame,m.frame)
        }
        colnames(matrix.frame) <- c('Method', names(stat.data.list[[stat.name]]))
        txt.file <- paste(paste('table', stat.name, 'round_2', sep='_'), 'txt', sep='.')
        matrix.file <- file.path(output.dir, txt.file)
        write.table(matrix.frame, matrix.file, row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
        message('File \'', txt.file, '\' with data values created ...')
    }
}

###############################################################################
###############################################################################
###############################################################################
######
######  Make statistics
######
###############################################################################
###############################################################################
###############################################################################

###############################################################################
## Compute the coefficient of variation for genes
## relative standard deviation
## 'intra.gene.CV'             =  'Intra Gene CV'
###############################################################################

if (coefficient.variation.bool) {
    message('Status: ', 'Computing coefficients of variation ...')
    message(rep('-', 79))
    #single core lapply call
    CV.list <- vector('list', length(normalized.count.list))
    for (i in 1:length(normalized.count.list)) {
        CV.list[[i]] <- computeCVs( input.frame  = normalized.count.list[[i]], 
                                    input.name   = names(normalized.count.list)[i], 
                                    sample.types = sample.types )
    }

    names(CV.list) <- names(normalized.count.list)
    
    # write data
    for (method.name in names(CV.list)) {
        # change names for matrix printing
        names(CV.list[[method.name]]) <- c('Mean CV A B', 'CV Samples A', 'CV Samples B')
        write.matrix.data( input.frame  = CV.list[[method.name]], 
                           method.name  = method.name,  
                           table.name   = paste(plot.titles.list['intra.gene.CV'][[1]], 'for', gsub('/', '-', method.name)), 
                           column.names = 'Gene Id')
        # change to original names
        names(CV.list[[method.name]]) <- c('cv.mean.list', 'cv.A.list', 'cv.B.list')
    }
    message(rep('#', 79))
}




###############################################################################
## Compute Fold Change
## 'intra.gene.FC'          =  'Intra Gene FC'
## 'intra.gene.log.FC'      =  'Intra Gene Log2 FC'
###############################################################################

if (fold.change.bool) {
    message('Status: ', 'Computing fold change ...')
    message(rep('-', 79))
    #single core lapply call
    FC.list <- vector('list', length(normalized.count.list))
    for (i in 1:length(normalized.count.list)) {
        FC.list[[i]] <- computeFCs( input.frame       = normalized.count.list[[i]], 
                                    input.name        = names(normalized.count.list)[i], 
                                    samples.matched.A = samples.matched.A, 
                                    samples.matched.B = samples.matched.B, 
                                    pseudo_count      = normalized.pseudo_count )
    }
    names(FC.list) <- names(normalized.count.list)
    message(rep('#', 79))
}




###############################################################################
# calculate standard deviation A and B
###############################################################################

if (calc.sd.mean.fc.list) {
    # calculate mean for each gene
    message('Status: ', 'Creating standard deviation and mean of FC list ... ')
    message(rep('-', 79))
    sd.mean.fc.list <- list()
    for (method.n in names(normalized.count.list)) {
        message('Computing standard deviation and mean FC for method ', method.n,  ' ... ')
        # calculate overall mean
        sd.fc.sample <- apply(FC.list[[method.n]][,1:ncol(FC.list[[method.n]])], 1, function(x) return(sd(x)))
        mean.fc.sample <- apply(FC.list[[method.n]][,1:ncol(FC.list[[method.n]])], 1, function(x) return(mean(x)))
        gene.id.frame <- normalized.count.list[[method.n]][,'Gene Id']
        sd.mean.fc.list[[method.n]] <- data.frame(gene.id.frame, sd.fc.sample, mean.fc.sample)
        colnames(sd.mean.fc.list[[method.n]]) <- as.factor(c('Gene Id', 'SD FC', 'mean FC'))
        rownames(sd.mean.fc.list[[method.n]]) <- as.factor(1:length(sd.mean.fc.list[[method.n]][,1]))
    }
    # write data
    # message('Making inner mean summary matrix ... ')
    # write.matrix.data( input.frame  = sd.norm.cnt.list[[method.name]], 
                       # method.name  = method.name,  
                       # table.name   = paste(plot.titles.list['norm.cnt.sum'][[1]], 'for', gsub('/', '-', method.name)), 
                       # column.names = 'Sample Name')
    message(rep('#', 79))
}



###############################################################################
# calculate correlation
###############################################################################

if (stat.correlation.list.bool) {
    # calculate mean for each gene
    message('Status: ', 'Creating correlation list ... ')
    message(rep('-', 79))
    data.correlation.list <- list()
    data.correlation.list[['SEQC.cnt']] <- SEQC.mean.all.norm.cnt.list
    data.correlation.list[['SEQC.fc']] <- SEQC.mean.FC.list
    data.correlation.list[['SEQS.cnt']] <- SEQS.mean.all.norm.cnt.list
    data.correlation.list[['SEQS.fc']] <- SEQS.mean.FC.list
    data.correlation.list[['taq.cnt']] <- taqman.norm.cnt.list
    data.correlation.list[['taq.fc']] <- taqman.FC.list
    data.correlation.list[['bio.cnt']] <- biorad.norm.cnt.list
    data.correlation.list[['bio.fc']] <- biorad.FC.list
    data.correlation.list[['SEQC.log.fc']] <- lapply(names(SEQC.mean.FC.list), function(x) return(log2(SEQC.mean.FC.list[[x]]+normalized.pseudo_count)))
    data.correlation.list[['SEQS.log.fc']] <- lapply(names(SEQS.mean.FC.list), function(x) return(log2(SEQS.mean.FC.list[[x]]+normalized.pseudo_count)))
    data.correlation.list[['taq.log.fc']] <- lapply(names(taqman.FC.list), function(x) return(log2(taqman.FC.list[[x]]+normalized.pseudo_count)))
    data.correlation.list[['bio.log.fc']] <- lapply(names(biorad.FC.list), function(x) return(log2(biorad.FC.list[[x]]+normalized.pseudo_count)))
    names(data.correlation.list[['SEQC.log.fc']]) <- names(SEQC.mean.FC.list)
    names(data.correlation.list[['SEQS.log.fc']]) <- names(SEQS.mean.FC.list)
    names(data.correlation.list[['taq.log.fc']]) <- names(taqman.FC.list)
    names(data.correlation.list[['bio.log.fc']]) <- names(biorad.FC.list)
    combi.type.list <- list()
    # count
    combi.type.list[['SEQC.cnt.A.spearman']] <- c('SEQC.cnt', 'SEQC.cnt', 'spearman', 'SEQC-A-samples')
    combi.type.list[['SEQC.cnt.A.pearson']] <- c('SEQC.cnt', 'SEQC.cnt', 'pearson', 'SEQC-A-samples')
    combi.type.list[['SEQC.cnt.B.spearman']] <- c('SEQC.cnt', 'SEQC.cnt', 'spearman', 'SEQC-B-samples')
    combi.type.list[['SEQC.cnt.B.pearson']] <- c('SEQC.cnt', 'SEQC.cnt', 'pearson', 'SEQC-B-samples')
    combi.type.list[['SEQS.cnt.A.spearman']] <- c('SEQS.cnt', 'SEQS.cnt', 'spearman', 'SEQC-A-samples')
    combi.type.list[['SEQS.cnt.A.pearson']] <- c('SEQS.cnt', 'SEQS.cnt', 'pearson', 'SEQC-A-samples')
    combi.type.list[['SEQS.cnt.B.spearman']] <- c('SEQS.cnt', 'SEQS.cnt', 'spearman', 'SEQC-B-samples')
    combi.type.list[['SEQS.cnt.B.pearson']] <- c('SEQS.cnt', 'SEQS.cnt', 'pearson', 'SEQC-B-samples')
    combi.type.list[['SEQC.SEQS.A.cnt.spearman']] <- c('SEQC.cnt', 'SEQS.cnt', 'spearman', 'SEQC-A-samples')
    combi.type.list[['SEQC.SEQS.A.cnt.pearson']] <- c('SEQC.cnt', 'SEQS.cnt', 'pearson', 'SEQC-A-samples')
    combi.type.list[['SEQC.SEQS.B.cnt.spearman']] <- c('SEQC.cnt', 'SEQS.cnt', 'spearman', 'SEQC-B-samples')
    combi.type.list[['SEQC.SEQS.B.cnt.pearson']] <- c('SEQC.cnt', 'SEQS.cnt', 'pearson', 'SEQC-B-samples')
    combi.type.list[['taq.cnt.A.spearman']] <- c('taq.cnt', 'taq.cnt', 'spearman', 'SEQC-A-samples')
    combi.type.list[['taq.cnt.A.pearson']] <- c('taq.cnt', 'taq.cnt', 'pearson', 'SEQC-A-samples')
    combi.type.list[['taq.cnt.B.spearman']] <- c('taq.cnt', 'taq.cnt', 'spearman', 'SEQC-B-samples')
    combi.type.list[['taq.cnt.B.pearson']] <- c('taq.cnt', 'taq.cnt', 'pearson', 'SEQC-B-samples')
    combi.type.list[['bio.cnt.A.spearman']] <- c('bio.cnt', 'bio.cnt', 'spearman', 'SEQC-A-samples')
    combi.type.list[['bio.cnt.A.pearson']] <- c('bio.cnt', 'bio.cnt', 'pearson', 'SEQC-A-samples')
    combi.type.list[['bio.cnt.B.spearman']] <- c('bio.cnt', 'bio.cnt', 'spearman', 'SEQC-B-samples')
    combi.type.list[['bio.cnt.B.pearson']] <- c('bio.cnt', 'bio.cnt', 'pearson', 'SEQC-B-samples')
    # FC
    combi.type.list[['SEQC.fc.spearman']] <- c('SEQC.fc', 'SEQC.fc', 'spearman', 'FC')
    combi.type.list[['SEQC.fc.pearson']] <- c('SEQC.fc', 'SEQC.fc', 'pearson', 'FC')
    combi.type.list[['SEQS.fc.spearman']] <- c('SEQS.fc', 'SEQS.fc', 'spearman', 'FC')
    combi.type.list[['SEQS.fc.pearson']] <- c('SEQS.fc', 'SEQS.fc', 'pearson', 'FC')
    combi.type.list[['SEQC.SEQS.fc.spearman']] <- c('SEQC.fc', 'SEQS.fc', 'spearman', 'FC')
    combi.type.list[['SEQC.SEQS.fc.pearson']] <- c('SEQC.fc', 'SEQS.fc', 'pearson', 'FC')
    combi.type.list[['taq.fc.spearman']] <- c('taq.fc', 'taq.fc', 'spearman', 'FC')
    combi.type.list[['taq.fc.pearson']] <- c('taq.fc', 'taq.fc', 'pearson', 'FC')
    combi.type.list[['bio.fc.spearman']] <- c('bio.fc', 'bio.fc', 'spearman', 'FC')
    combi.type.list[['bio.fc.pearson']] <- c('bio.fc', 'bio.fc', 'pearson', 'FC')
    # log2 FC
    combi.type.list[['SEQC.log.fc.spearman']] <- c('SEQC.log.fc', 'SEQC.log.fc', 'spearman', 'FC')
    combi.type.list[['SEQC.log.fc.pearson']] <- c('SEQC.log.fc', 'SEQC.log.fc', 'pearson', 'FC')
    combi.type.list[['SEQS.log.fc.spearman']] <- c('SEQS.log.fc', 'SEQS.log.fc', 'spearman', 'FC')
    combi.type.list[['SEQS.log.fc.pearson']] <- c('SEQS.log.fc', 'SEQS.log.fc', 'pearson', 'FC')
    combi.type.list[['SEQC.SEQS.log.fc.spearman']] <- c('SEQC.log.fc', 'SEQS.log.fc', 'spearman', 'FC')
    combi.type.list[['SEQC.SEQS.log.fc.pearson']] <- c('SEQC.log.fc', 'SEQS.log.fc', 'pearson', 'FC')
    combi.type.list[['taq.log.fc.spearman']] <- c('taq.log.fc', 'taq.log.fc', 'spearman', 'FC')
    combi.type.list[['taq.log.fc.pearson']] <- c('taq.log.fc', 'taq.log.fc', 'pearson', 'FC')
    combi.type.list[['bio.log.fc.spearman']] <- c('bio.log.fc', 'bio.log.fc', 'spearman', 'FC')
    combi.type.list[['bio.log.fc.pearson']] <- c('bio.log.fc', 'bio.log.fc', 'pearson', 'FC')
    # log2 sd
    combi.type.list[['SEQC.cnt.A.sd']] <- c('SEQC.cnt', 'SEQC.cnt', 'SD', 'SEQC-A-samples')
    combi.type.list[['SEQC.cnt.B.sd']] <- c('SEQC.cnt', 'SEQC.cnt', 'SD', 'SEQC-B-samples')
    combi.type.list[['SEQS.cnt.A.sd']] <- c('SEQS.cnt', 'SEQS.cnt', 'SD', 'SEQC-A-samples')
    combi.type.list[['SEQS.cnt.B.sd']] <- c('SEQS.cnt', 'SEQS.cnt', 'SD', 'SEQC-B-samples')
    combi.type.list[['SEQC.SEQS.A.cnt.sd']] <- c('SEQC.cnt', 'SEQS.cnt', 'SD', 'SEQC-A-samples')
    combi.type.list[['SEQC.SEQS.B.cnt.sd']] <- c('SEQC.cnt', 'SEQS.cnt', 'SD', 'SEQC-B-samples')
    combi.type.list[['taq.cnt.A.sd']] <- c('taq.cnt', 'taq.cnt', 'SD', 'SEQC-A-samples')
    combi.type.list[['taq.cnt.B.sd']] <- c('taq.cnt', 'taq.cnt', 'SD', 'SEQC-B-samples')
    combi.type.list[['bio.cnt.A.sd']] <- c('bio.cnt', 'bio.cnt', 'SD', 'SEQC-A-samples')
    combi.type.list[['bio.cnt.B.sd']] <- c('bio.cnt', 'bio.cnt', 'SD', 'SEQC-B-samples')
    combi.type.list[['SEQC.fc.sd']] <- c('SEQC.fc', 'SEQC.fc', 'SD', 'FC')
    combi.type.list[['SEQS.fc.sd']] <- c('SEQS.fc', 'SEQS.fc', 'SD', 'FC')
    combi.type.list[['SEQC.SEQS.fc.sd']] <- c('SEQC.fc', 'SEQS.fc', 'SD', 'FC')
    combi.type.list[['taq.fc.sd']] <- c('taq.fc', 'taq.fc', 'SD', 'FC')
    combi.type.list[['bio.fc.sd']] <- c('bio.fc', 'bio.fc', 'SD', 'FC')
    combi.type.list[['SEQC.log.fc.sd']] <- c('SEQC.log.fc', 'SEQC.log.fc', 'SD', 'FC')
    combi.type.list[['SEQS.log.fc.sd']] <- c('SEQS.log.fc', 'SEQS.log.fc', 'SD', 'FC')
    combi.type.list[['SEQC.SEQS.log.fc.sd']] <- c('SEQC.log.fc', 'SEQS.log.fc', 'SD', 'FC')
    combi.type.list[['taq.log.fc.sd']] <- c('taq.log.fc', 'taq.log.fc', 'SD', 'FC')
    combi.type.list[['bio.log.fc.sd']] <- c('bio.log.fc', 'bio.log.fc', 'SD', 'FC')
    # make calculations
    stat.correlation.list <- list()
    for (correlation.type in names(combi.type.list)) {
        message('Status: ', 'Creating correlation for ', correlation.type, ' ... ')
        combi.list <- combi.type.list[[correlation.type]]
        data.cor.list.1 <- data.correlation.list[[combi.list[[1]]]]
        data.cor.list.2 <- data.correlation.list[[combi.list[[2]]]]
        stat.correlation.list[[correlation.type]] <- matrix(NA, nrow=length(names(data.cor.list.1)), ncol=length(names(data.cor.list.2)))
        rownames(stat.correlation.list[[correlation.type]]) <- names(data.cor.list.1)
        colnames(stat.correlation.list[[correlation.type]]) <- names(data.cor.list.2)
        for (method.1 in names(data.cor.list.1)) {
            for (method.2 in names(data.cor.list.2)) {
                if (combi.list[[3]] == 'SD') {
                    if (combi.list[[4]] == 'FC') {
                        stat.correlation.list[[correlation.type]][method.1,method.2] <- abs(log2((sd(data.cor.list.1[[method.1]])+normalized.pseudo_count)/(sd(data.cor.list.2[[method.2]])+normalized.pseudo_count)))
                    } else {
                        stat.correlation.list[[correlation.type]][method.1,method.2] <- abs(log2((sd(data.cor.list.1[[method.1]][[combi.list[[4]]]]+normalized.pseudo_count)/(sd(data.cor.list.2[[method.2]][[combi.list[[4]]]]+normalized.pseudo_count)))))
                    }
                } else {
                    if (combi.list[[4]] == 'FC') {
                        stat.correlation.list[[correlation.type]][method.1,method.2] <- cor(data.cor.list.1[[method.1]],data.cor.list.2[[method.2]],method=combi.list[[3]])
                    } else {
                        stat.correlation.list[[correlation.type]][method.1,method.2] <- cor(data.cor.list.1[[method.1]][[combi.list[[4]]]],data.cor.list.2[[method.2]][[combi.list[[4]]]],method=combi.list[[3]])
                    }
                }
            }
        }
    }
    message(rep('#', 79))
}


###############################################################################
## Compute the coefficient of variation of FC and log(FC) for genes
## 'intra.gene.CV.FC'          =  'Intra Gene CV of FC'
## 'intra.gene.CV.log.FC'      =  'Intra Gene CV of Log2 FC'
###############################################################################

if (fold.change.cv.bool) {
    message('Status: ', 'Computing log FC coefficients of variation ...')
    message(rep('-', 79))
    #single core lapply call
    fcCV.list <- vector('list', length(normalized.count.list))
    for (i in 1:length(normalized.count.list)) {
        fcCV.list[[i]] <- computeFcCVs( fc.mat            = FC.list[[i]], 
                                        input.name        = names(normalized.count.list)[i],
                                        log.fc.threshold  = logFC.thresh)
    }

    names(fcCV.list) <- names(normalized.count.list)
    
    # write data
    for (method.name in names(fcCV.list)) {
        # change names for matrix printing
        names(fcCV.list[[method.name]]) <- c('CV(FC)', 'CV(log2(FC))')
        write.matrix.data( input.frame  = fcCV.list[[method.name]], 
                           method.name  = method.name,  
                           table.name   = paste(plot.titles.list['intra.gene.CV.FC'][[1]], 'for', gsub('/', '-', method.name), 'with min(abs(log2(FC))) of', as.character(logFC.thresh)), 
                           column.names = 'Gene Id')
        # change to original names
        names(fcCV.list[[method.name]]) <- c('intra.gene.CV.FC', 'intra.gene.CV.log.FC')
    }
    message(rep('#', 79))
}




###############################################################################
## Reproducibility correlations (precision)
## sample vs sample
## 'intra.met.pearson'         =  'Intra Method Pearson CC'
## 'intra.met.spearman'        =  'Intra Method Spearman CC'
## 'intra.met.CD'              =  'Intra Method CD'
###############################################################################

if (intra.method.cnt.cor.bool) {
    message('Status: ', 'Computing reproducibility - intra method count correlations ...')
    message(rep('-', 79))
    #single core lapply call
    intra.method.cnt.cor.list <- vector('list', length(normalized.count.list))
    for (i in 1:length(normalized.count.list)) {
        intra.method.cnt.cor.list[[i]] <- computeIntraMethodCorrelations( input.frame        = normalized.count.list[[i]], 
                                                                          input.name         = names(normalized.count.list)[i], 
                                                                          samples.list       = c.list.samples, 
                                                                          sample.types       = sample.types, 
                                                                          id.name            = id.type.name, 
                                                                          num.samples.thresh = intra.method.cnt.cor.thresh)
    }
    names(intra.method.cnt.cor.list) <- names(normalized.count.list)
    
    # write results
    if (intra.method.cnt.cor.write.bool) {
        for (method.name in names(intra.method.cnt.cor.list)) {
            for (precision.method in names(intra.method.cnt.cor.list[[method.name]])) {
                for (sample.type in names(intra.method.cnt.cor.list[[method.name]][[precision.method]])) {
                    # write pearson correlation matrix to file
                    write.matrix.data( input.frame  = intra.method.cnt.cor.list[[method.name]][[precision.method]][[sample.type]], 
                                       method.name  = method.name, 
                                       table.name   = paste(plot.titles.list[precision.method][[1]], 'for', gsub('/', '-', method.name), sample.naming.vector[[sample.type]], sep=' '), 
                                       column.names = 'Sample Name')
                }
            }
        }
    }
    message(rep('#', 79))
}




###############################################################################
## Reproducibility: log fold change correlations
## 'intra.met.pearson.fc'      =  'Intra Method FC Pearson CC'
## 'intra.met.spearman.fc'     =  'Intra Method FC Spearman CC'
## 'intra.met.CD.fc'           =  'Intra Method FC CD'
## 'intra.met.pearson.log.fc'  =  'Intra Method Log2 FC Pearson CC'
## 'intra.met.spearman.log.fc' =  'Intra Method Log2 FC Spearman CC'
## 'intra.met.CD.log.fc'       =  'Intra Method Log2 FC CD'
###############################################################################

## log fold change correlations
if (intra.method.fc.cor.bool) {
    message('Status: ', 'Computing (log) fold change correlations (reproducibility) ...')
    message(rep('-', 79))
    #single core lapply call
    intra.method.fc.cor.list <- vector('list', length(normalized.count.list))
    for (i in 1:length(normalized.count.list)) {
        intra.method.fc.cor.list[[i]] <- computeIntraMethodFcCorrelations( input.frame       = normalized.count.list[[i]], 
                                                                           input.name        = names(normalized.count.list)[i], 
                                                                           samples.matched.A = samples.matched.A, 
                                                                           samples.matched.B = samples.matched.B, 
                                                                           pseudo_count      = normalized.pseudo_count)
    }
    names(intra.method.fc.cor.list) <- names(normalized.count.list)
    
    # write results
    if (intra.method.fc.cor.write.bool) {
        for (method.name in names(intra.method.fc.cor.list)) {
            for (precision.method in names(intra.method.fc.cor.list[[method.name]])) {
                # write pearson correlation matrix to file
                write.matrix.data( input.frame  = intra.method.fc.cor.list[[method.name]][[precision.method]], 
                                   method.name  = method.name, 
                                   table.name   = paste(plot.titles.list[precision.method][[1]], 'for', gsub('/', '-', method.name), sep=' '), 
                                   column.names = 'Sample Name')
            }
        }
    }
    message(rep('#', 79))
}


###############################################################################
##  Comparisons between methods (agreement)
###############################################################################


###############################################################################
## Fold change correlations among methods
# 'inter.met.pearson.fc'           =  'Inter Method FC Pearson CC'
# 'inter.met.spearman.fc'          =  'Inter Method FC Spearman CC'
# 'inter.met.CD.fc'                =  'Inter Method FC CD'
# 'inter.met.pearson.log.fc'       =  'Inter Method Log2 FC Pearson CC'
# 'inter.met.spearman.log.fc'      =  'Inter Method Log2 FC Spearman CC'
# 'inter.met.CD.log.fc'            =  'Inter Method Log2 FC CD'
# 'inter.met.false.negative.rate'  =  'Inter Method SG FNR'
# 'inter.met.sensitivity'          =  'Inter Method SG Sensitivity
# 'inter.met.false.positive.rate'  =  'Inter Method SG FPR'
# 'inter.met.specificity'          =  'Inter Method SG Specificity',
# 'inter.met.false.discovery.rate' =  'Inter Method SG FDR'
# 'inter.met.precision'            =  'Inter Method SG Precision'
###############################################################################

## (log) fold change correlations among methods
if (inter.method.fc.cor.bool) {
    message('Status: ', 'Computing (log) fold change correlations among methods ...')
    message(rep('-', 79))
    inter.method.fc.cor.list <- list()
    for (i in 1:length(normalized.count.list)) {
        method.name.2 <- names(normalized.count.list)[i]
        #single core lapply call
        inter.method.fc.cor.list[[method.name.2]] <- vector('list', length(normalized.count.list))
        for (j in 1:length(normalized.count.list)) {
            inter.method.fc.cor.list[[method.name.2]][[j]] <- computeFoldChangeCorrelations( 
                                                                result1.frame     = normalized.count.list[[j]], 
                                                                result2.frame     = normalized.count.list[[i]], 
                                                                method.name.1     = names(normalized.count.list)[j], 
                                                                method.name.2     = method.name.2, 
                                                                samples.matched.A = samples.matched.A, 
                                                                samples.matched.B = samples.matched.B, 
                                                                pseudo_count.1    = normalized.pseudo_count, 
                                                                pseudo_count.2    = normalized.pseudo_count)
        }
        names(inter.method.fc.cor.list[[method.name.2]]) <- names(normalized.count.list)
    }
    
    # write results
    inter.method.type.naming <- c( 'inter.met.pearson.fc'          = 'for',    'inter.met.spearman.fc'          = 'for',    'inter.met.CD.fc'               = 'for',
                                   'inter.met.pearson.log.fc'      = 'for',    'inter.met.spearman.log.fc'      = 'for',    'inter.met.CD.log.fc'           = 'for',
                                   'inter.met.false.negative.rate' = 'for H0', 'inter.met.sensitivity'          = 'for H0', 'inter.met.false.positive.rate' = 'for H0',
                                   'inter.met.specificity'         = 'for H0', 'inter.met.false.discovery.rate' = 'for H0', 'inter.met.precision'           = 'for H0' )
    if (inter.method.fc.cor.write.bool) {
        for (method.name.2 in names(inter.method.fc.cor.list)) {
            for (inter.method.type in names(inter.method.fc.cor.list[[method.name.2]][[1]])) {
                # write pearson correlation matrix to file
                write.matrix.data( input.frame  = lapply(inter.method.fc.cor.list[[method.name.2]], function(method.name.1) return(method.name.1[[inter.method.type]])), 
                                   method.name  = '', 
                                   table.name   = paste(plot.titles.list[inter.method.type][[1]], inter.method.type.naming[[inter.method.type]], gsub('/', '-', method.name.2), sep=' '), 
                                   column.names = 'Sample FC Name')
            }
        }
    }
    message(rep('#', 79))
}




###############################################################################
## Expression correlations among methods
# 'inter.met.pearson'              =  'Inter Method Pearson CC'
# 'inter.met.spearman'             =  'Inter Method Spearman CC'
# 'inter.met.CD'                   =  'Inter Method CD'
# 'inter.met.log.pearson'          =  'Inter Method Log2 Pearson CC'
# 'inter.met.log.spearman'         =  'Inter Method Log2 Spearman CC'
# 'inter.met.log.CD'               =  'Inter Method Log2 CD'
###############################################################################

## expression correlations among methods
if (inter.method.cnt.cor.bool) {
    message('Status: ', 'Computing gene expression correlations among methods ...')
    message(rep('-', 79))
    inter.method.cnt.cor.list <- list()
    for (i in 1:length(normalized.count.list)) {
        method.name.2 <- names(normalized.count.list)[i]
        #single core lapply call
        inter.method.cnt.cor.list[[method.name.2]] <- vector('list', length(normalized.count.list))
        for (j in 1:length(normalized.count.list)) {
            inter.method.cnt.cor.list[[method.name.2]][[j]] <- computeInterMethodCorrelations( 
                                                                 result1.frame     = normalized.count.list[[j]], 
                                                                 result2.frame     = normalized.count.list[[i]], 
                                                                 method.name.1     = names(normalized.count.list)[j],
                                                                 method.name.2     = method.name.2, 
                                                                 pseudo_count      = normalized.pseudo_count)
        }
        names(inter.method.cnt.cor.list[[method.name.2]]) <- names(normalized.count.list)
    }
    
    # write results
    if (inter.method.cnt.cor.write.bool) {
        for (method.name.2 in names(inter.method.cnt.cor.list)) {
            for (inter.method.type in names(inter.method.cnt.cor.list[[method.name.2]][[1]])) {
                # write pearson correlation matrix to file
                write.matrix.data( input.frame  = lapply(inter.method.cnt.cor.list[[method.name.2]], function(method.name.1) return(method.name.1[[inter.method.type]])), 
                                   method.name  = '', 
                                   table.name   = paste(plot.titles.list[inter.method.type][[1]], 'for', gsub('/', '-', method.name.2), sep=' '), 
                                   column.names = 'Sample Name')
            }
        }
    }
    message(rep('#', 79))
}





###############################################################################
###############################################################################
###############################################################################
######
######  Print Results
######
###############################################################################
###############################################################################
###############################################################################

###############################################################################
# plot ROC
###############################################################################
if (rocplot.FC.bool) {
    # plot mean gene expression
    expression.plot.list <- list()
    d.set.setting <- list()
    
    #"EQPQM-ambig/HISAT"     "EQPQM-ambig/STAR"      "EQPQM-ambig/TOPHAT"    "EQPQM-unambig/HISAT"  
    #"EQPQM-unambig/STAR"    "EQPQM-unambig/TOPHAT"  "FEATC-default/HISAT"   "FEATC-default/STAR"   
    #"FEATC-default/TOPHAT"  "HTSEQ-union/HISAT"     "HTSEQ-union/STAR"      "HTSEQ-union/TOPHAT"   
    #"HTSEQ-istrict/HISAT"   "HTSEQ-istrict/STAR"    "HTSEQ-istrict/TOPHAT"  "HTSEQ-inonempt/HISAT" 
    #"HTSEQ-inonempt/STAR"   "HTSEQ-inonempt/TOPHAT" "CUFFL-default/HISAT"   "CUFFL-default/STAR"   
    #"CUFFL-default/TOPHAT"  "FLXCP-default/HISAT"   "FLXCP-default/STAR"    "FLXCP-default/TOPHAT" 
    #"BTSEQ-default/BWTRS"   "RSEM-default/BRSEM"    "KLLST-default/AFREE"   "KLLST-default/AFREE"  
    expression.plot.list[['roc.FC.taq.2.0']] <- list(log2(taqman.FC.list[['BTSEQ-default/BWTRS']]+normalized.pseudo_count),
                                                     log2(taqman.FC.list[['CUFFL-default/HISAT']]+normalized.pseudo_count),
                                                     log2(taqman.FC.list[['EQPQM-ambig/HISAT']]+normalized.pseudo_count),
                                                     log2(taqman.FC.list[['EQPQM-unambig/HISAT']]+normalized.pseudo_count),
                                                     log2(taqman.FC.list[['FEATC-default/HISAT']]+normalized.pseudo_count),
                                                     log2(taqman.FC.list[['FLXCP-default/HISAT']]+normalized.pseudo_count),
                                                     log2(taqman.FC.list[['HTSEQ-inonempt/HISAT']]+normalized.pseudo_count),
                                                     log2(taqman.FC.list[['HTSEQ-istrict/HISAT']]+normalized.pseudo_count),
                                                     log2(taqman.FC.list[['HTSEQ-union/HISAT']]+normalized.pseudo_count),
                                                     log2(taqman.FC.list[['KLLST-default/AFREE']]+normalized.pseudo_count),
                                                     log2(taqman.FC.list[['RSEM-default/BRSEM']]+normalized.pseudo_count),
                                                     log2(taqman.FC.list[['TQMAN-default/QPCR']]+normalized.pseudo_count))
    names(expression.plot.list[['roc.FC.taq.2.0']]) <- c('EQPQM-ambig/HISAT',   'CUFFL-default/HISAT', 'EQPQM-ambig/HISAT',    'EQPQM-unambig/HISAT',
                                                         'FEATC-default/HISAT', 'FLXCP-default/HISAT', 'HTSEQ-inonempt/HISAT', 'HTSEQ-istrict/HISAT',
                                                         'HTSEQ-union/HISAT',   'KLLST-default/AFREE', 'RSEM-default/BRSEM',   'TQMAN-default/QPCR')
    expression.plot.list[['roc.FC.SEQS.2.0']] <- list(log2(SEQS.mean.FC.list[['BTSEQ-default/BWTRS']]+normalized.pseudo_count),
                                                     log2(SEQS.mean.FC.list[['CUFFL-default/HISAT']]+normalized.pseudo_count),
                                                     log2(SEQS.mean.FC.list[['EQPQM-ambig/HISAT']]+normalized.pseudo_count),
                                                     log2(SEQS.mean.FC.list[['EQPQM-unambig/HISAT']]+normalized.pseudo_count),
                                                     log2(SEQS.mean.FC.list[['FEATC-default/HISAT']]+normalized.pseudo_count),
                                                     log2(SEQS.mean.FC.list[['FLXCP-default/HISAT']]+normalized.pseudo_count),
                                                     log2(SEQS.mean.FC.list[['HTSEQ-inonempt/HISAT']]+normalized.pseudo_count),
                                                     log2(SEQS.mean.FC.list[['HTSEQ-istrict/HISAT']]+normalized.pseudo_count),
                                                     log2(SEQS.mean.FC.list[['HTSEQ-union/HISAT']]+normalized.pseudo_count),
                                                     log2(SEQS.mean.FC.list[['KLLST-default/AFREE']]+normalized.pseudo_count),
                                                     log2(SEQS.mean.FC.list[['RSEM-default/BRSEM']]+normalized.pseudo_count),
                                                     log2(SEQS.mean.FC.list[['FLXSM-default/SIMED']]+normalized.pseudo_count))
    names(expression.plot.list[['roc.FC.SEQS.2.0']]) <- c('EQPQM-ambig/HISAT',   'CUFFL-default/HISAT', 'EQPQM-ambig/HISAT',    'EQPQM-unambig/HISAT',
                                                         'FEATC-default/HISAT', 'FLXCP-default/HISAT', 'HTSEQ-inonempt/HISAT', 'HTSEQ-istrict/HISAT',
                                                         'HTSEQ-union/HISAT',   'KLLST-default/AFREE', 'RSEM-default/BRSEM',   'FLXSM-default/SIMED')
    #expression.plot.list[['roc.FC.taq.1.0']] <- lapply(taqman.FC.list, function(x) return(abs(log2(x+normalized.pseudo_count))))
    #names(expression.plot.list[['roc.FC.taq.1.0']]) <- names(taqman.FC.list)
    d.set.setting[['roc.FC.taq.2.0']] <- list('TQMAN-default/QPCR',100,normalized.pseudo_count,2.0,'TaqMan ROC DEG 2.0')
    d.set.setting[['roc.FC.SEQS.2.0']] <- list('FLXSM-default/SIMED',100,normalized.pseudo_count,2.0,'SEQS ROC DEG 2.0')
    d.set.setting[['roc.FC.SEQS.2.0']] <- list('FLXSM-default/SIMED',100,normalized.pseudo_count,2.0,'SEQS ROC DEG 2.0')
    d.set.setting[['roc.FC.SEQS.2.0']] <- list('FLXSM-default/SIMED',100,normalized.pseudo_count,2.0,'SEQS ROC DEG 2.0')
    #BIORD-default/QPCR
    #FLXSM-default/SIMED
    for (plot.name in names(expression.plot.list)) {
        plot.ROC( input.list      = expression.plot.list[[plot.name]], 
                  gold.name       = d.set.setting[[plot.name]][[1]],
                  density.setting = d.set.setting[[plot.name]][[2]], 
                  pseudo          = d.set.setting[[plot.name]][[3]], 
                  g.thresh        = d.set.setting[[plot.name]][[4]], 
                  plot.title.name = d.set.setting[[plot.name]][[5]])
    }
}

###############################################################################
# Mean Gene Expression
###############################################################################
if (expression.genes.mean.plot.bool) {
    # plot mean gene expression
    plot.scales <- list()
    scales.max.A <- max(unlist(summary.list[['expression.genes.total']][['SEQC-A']], function(x) return(max(x))))
    scales.min.A <- min(unlist(summary.list[['expression.genes.total']][['SEQC-A']], function(x) return(min(x))))
    scales.max.B <- max(unlist(summary.list[['expression.genes.total']][['SEQC-B']], function(x) return(max(x))))
    scales.min.B <- min(unlist(summary.list[['expression.genes.total']][['SEQC-B']], function(x) return(min(x))))
    plot.scales[['cnt.sum']] <- c(min(scales.min.A, scales.min.B), max(scales.max.A, scales.max.B))
    scales.max.A <- max(unlist(norm.summary.list[['expression.genes.total']][['SEQC-A']], function(x) return(max(x))))
    scales.min.A <- min(unlist(norm.summary.list[['expression.genes.total']][['SEQC-A']], function(x) return(min(x))))
    scales.max.B <- max(unlist(norm.summary.list[['expression.genes.total']][['SEQC-B']], function(x) return(max(x))))
    scales.min.B <- min(unlist(norm.summary.list[['expression.genes.total']][['SEQC-B']], function(x) return(min(x))))
    plot.scales[['norm.cnt.sum']] <- c(min(scales.min.A, scales.min.B), max(scales.max.A, scales.max.B))
    scales.max.A <- max(unlist(summary.list[['expression.genes.mean']][['SEQC-A']], function(x) return(max(x))))
    scales.min.A <- min(unlist(summary.list[['expression.genes.mean']][['SEQC-A']], function(x) return(min(x))))
    scales.max.B <- max(unlist(summary.list[['expression.genes.mean']][['SEQC-B']], function(x) return(max(x))))
    scales.min.B <- min(unlist(summary.list[['expression.genes.mean']][['SEQC-B']], function(x) return(min(x))))
    plot.scales[['mean.cnt.sum']] <- c(min(scales.min.A, scales.min.B), max(scales.max.A, scales.max.B))
    scales.max.A <- max(unlist(norm.summary.list[['expression.genes.mean']][['SEQC-A']], function(x) return(max(x))))
    scales.min.A <- min(unlist(norm.summary.list[['expression.genes.mean']][['SEQC-A']], function(x) return(min(x))))
    scales.max.B <- max(unlist(norm.summary.list[['expression.genes.mean']][['SEQC-B']], function(x) return(max(x))))
    scales.min.B <- min(unlist(norm.summary.list[['expression.genes.mean']][['SEQC-B']], function(x) return(min(x))))
    plot.scales[['mean.norm.cnt.sum']] <- c(min(scales.min.A, scales.min.B), max(scales.max.A, scales.max.B))
    for (sample.type in sample.types) {
        expression.plot.list <- list()
        expression.plot.list[['cnt.sum']] <- summary.list[['expression.genes.total']][[sample.type]]
        expression.plot.list[['norm.cnt.sum']] <- norm.summary.list[['expression.genes.total']][[sample.type]]
        expression.plot.list[['mean.cnt.sum']] <- summary.list[['expression.genes.mean']][[sample.type]]
        expression.plot.list[['mean.norm.cnt.sum']] <- norm.summary.list[['expression.genes.mean']][[sample.type]]
        for (plot.name in names(expression.plot.list)) {
            boxPlot( input.list      = expression.plot.list[[plot.name]], 
                     plot.title.name = paste(dataset.name, plot.titles.list[plot.name][[1]], sample.naming.vector[[sample.type]]), 
                     x.lab.name      = 'Counting Method', 
                     y.lab.name      = paste(dataset.name, plot.titles.list[plot.name][[1]], sample.naming.vector[[sample.type]]),
                     reduce.jitter   = FALSE,
                     use.jitter      = TRUE,
                     dataset.name    = dataset.name,
                     use.scale       = plot.scales[[plot.name]])
         }
    }
}


###############################################################################
# box plot FC
###############################################################################
if (boxplot.FC.plot.bool) {
    # plot mean gene expression
    expression.plot.list <- list()
    expression.plot.list[['box.FC.taq']] <- lapply(taqman.FC.list, function(x) return(log2(x+normalized.pseudo_count)))
    expression.plot.list[['box.FC.bio']] <- lapply(biorad.FC.list, function(x) return(log2(x+normalized.pseudo_count)))
    expression.plot.list[['box.FC.SEQC']] <- lapply(SEQC.mean.FC.list, function(x) return(log2(x+normalized.pseudo_count)))
    expression.plot.list[['box.FC.SEQS']] <- lapply(SEQS.mean.FC.list, function(x) return(log2(x+normalized.pseudo_count)))
    names(expression.plot.list[['box.FC.taq']]) <- names(taqman.FC.list)
    names(expression.plot.list[['box.FC.bio']]) <- names(biorad.FC.list)
    names(expression.plot.list[['box.FC.SEQC']]) <- names(SEQC.mean.FC.list)
    names(expression.plot.list[['box.FC.SEQS']]) <- names(SEQS.mean.FC.list)
    expression.plot.list[['box.FC.taq.abs']] <- lapply(taqman.FC.list, function(x) return(abs(log2(x+normalized.pseudo_count))))
    expression.plot.list[['box.FC.bio.abs']] <- lapply(biorad.FC.list, function(x) return(abs(log2(x+normalized.pseudo_count))))
    expression.plot.list[['box.FC.SEQC.abs']] <- lapply(SEQC.mean.FC.list, function(x) return(abs(log2(x+normalized.pseudo_count))))
    expression.plot.list[['box.FC.SEQS.abs']] <- lapply(SEQS.mean.FC.list, function(x) return(abs(log2(x+normalized.pseudo_count))))
    names(expression.plot.list[['box.FC.taq.abs']]) <- names(taqman.FC.list)
    names(expression.plot.list[['box.FC.bio.abs']]) <- names(biorad.FC.list)
    names(expression.plot.list[['box.FC.SEQC.abs']]) <- names(SEQC.mean.FC.list)
    names(expression.plot.list[['box.FC.SEQS.abs']]) <- names(SEQS.mean.FC.list)
    d.set.naming <- list()
    d.set.naming[['box.FC.taq']] <- 'TaqMan'
    d.set.naming[['box.FC.bio']] <- 'BioRad'
    d.set.naming[['box.FC.SEQC']] <- 'SEQC'
    d.set.naming[['box.FC.SEQS']] <- 'SEQS'
    d.set.naming[['box.FC.taq.abs']] <- 'TaqMan'
    d.set.naming[['box.FC.bio.abs']] <- 'BioRad'
    d.set.naming[['box.FC.SEQC.abs']] <- 'SEQC'
    d.set.naming[['box.FC.SEQS.abs']] <- 'SEQS'
    for (plot.name in names(expression.plot.list)) {
        boxPlot( input.list      = expression.plot.list[[plot.name]], 
                 plot.title.name = paste('BoxPlot', plot.titles.list[plot.name][[1]], sep=' '), 
                 x.lab.name      = 'Counting Method', 
                 y.lab.name      = plot.titles.list[plot.name][[1]],
                 reduce.jitter   = FALSE,
                 use.jitter      = FALSE,
                 dataset.name    = d.set.naming[[plot.name]])
    }
}

###############################################################################
# Scatter Plot of counts biorad and taqman data  mean.taq.norm.cnt.sum
###############################################################################
if (biotaq.scatter.plot.bool) {
    # plot results
    tmp.taq.A.list <- lapply(taqman.norm.cnt.list, function(x) return(log2(x[['SEQC-A-samples']]+normalized.pseudo_count)))
    tmp.taq.B.list <- lapply(taqman.norm.cnt.list, function(x) return(log2(x[['SEQC-B-samples']]+normalized.pseudo_count)))
    scatter.plot.list <- list()
    for (taq.num in 1:length(taqman.norm.cnt.list)) {
        scatter.plot.list[[taq.num]] <- c('TQMAN-default/QPCR', names(taqman.norm.cnt.list)[[taq.num]])
    }
    for (taq.num in 1:length(scatter.plot.list)) {
        message('Status: ', 'Plot TaqMan Scatter Plot for ... ', scatter.plot.list[[taq.num]][[2]])
        scatterPlot( input.list.A    = tmp.taq.A.list, 
                     input.list.B    = tmp.taq.B.list, 
                     plot.title.name = paste(dataset.name, 'Scatter Plot for', plot.titles.list['mean.taq.norm.cnt.sum'][[1]], 'for', gsub('/', '-', scatter.plot.list[[taq.num]][[2]]),sep=' '), 
                     method.1        = scatter.plot.list[[taq.num]][[1]],
                     method.2        = scatter.plot.list[[taq.num]][[2]],
                     x.lab.name      = 'TaqMan', 
                     y.lab.name      = scatter.plot.list[[taq.num]][[2]],
                     alpha.val       = 1.0,
                     dataset.name    = dataset.name)
    }
    tmp.bio.A.list <- lapply(biorad.norm.cnt.list, function(x) return(log2(x[['SEQC-A-samples']]+normalized.pseudo_count)))
    tmp.bio.B.list <- lapply(biorad.norm.cnt.list, function(x) return(log2(x[['SEQC-B-samples']]+normalized.pseudo_count)))
    scatter.plot.list <- list()
    for (bio.num in 1:length(biorad.norm.cnt.list)) {
        scatter.plot.list[[bio.num]] <- c('BIORD-default/QPCR', names(biorad.norm.cnt.list)[[bio.num]])
    }
    for (bio.num in 1:length(scatter.plot.list)) {
        message('Status: ', 'Plot Bio-rad Scatter Plot for ... ', scatter.plot.list[[bio.num]][[2]])
        scatterPlot( input.list.A    = tmp.bio.A.list, 
                     input.list.B    = tmp.bio.B.list, 
                     plot.title.name = paste(dataset.name, 'Scatter Plot for', plot.titles.list['mean.bio.norm.cnt.sum'][[1]], 'for', gsub('/', '-', scatter.plot.list[[bio.num]][[2]]),sep=' '), 
                     method.1        = scatter.plot.list[[bio.num]][[1]],
                     method.2        = scatter.plot.list[[bio.num]][[2]],
                     x.lab.name      = 'Bio-Rad', 
                     y.lab.name      = scatter.plot.list[[bio.num]][[2]],
                     alpha.val       = 1.0,
                     dataset.name    = dataset.name)
    }
}


###############################################################################
# Scatter Plot of FC biorad and taqman data  mean.taq.norm.cnt.sum
###############################################################################
if (biotaq.fc.scatter.plot.bool) {
    # plot results taqman.FC.list
    tmp.taq.fc.list <- lapply(taqman.FC.list, function(x) return(log2(x+normalized.pseudo_count)))
    scatter.plot.list <- list()
    for (taq.num in 1:length(taqman.FC.list)) {
        scatter.plot.list[[taq.num]] <- c('TQMAN-default/QPCR', names(taqman.FC.list)[[taq.num]])
    }
    for (taq.num in 1:length(scatter.plot.list)) {
        message('Status: ', 'Plot TaqMan Scatter FC for ... ', scatter.plot.list[[taq.num]][[2]])
        scatterPlot( input.list.A    = tmp.taq.fc.list, 
                     input.list.B    = FALSE, 
                     plot.title.name = paste(dataset.name, 'Scatter Plot for', plot.titles.list['taq.fc.sum'][[1]], 'for', gsub('/', '-', scatter.plot.list[[taq.num]][[2]]),sep=' '), 
                     method.1        = scatter.plot.list[[taq.num]][[1]],
                     method.2        = scatter.plot.list[[taq.num]][[2]],
                     x.lab.name      = 'TaqMan', 
                     y.lab.name      = scatter.plot.list[[taq.num]][[2]],
                     alpha.val       = 1.0,
                     dataset.name    = dataset.name,
                     single.data     = TRUE)
    }
    tmp.bio.fc.list <- lapply(biorad.FC.list, function(x) return(log2(x+normalized.pseudo_count)))
    scatter.plot.list <- list()
    for (bio.num in 1:length(biorad.FC.list)) {
        scatter.plot.list[[bio.num]] <- c('BIORD-default/QPCR', names(biorad.FC.list)[[bio.num]])
    }
    for (bio.num in 1:length(scatter.plot.list)) {
        message('Status: ', 'Plot Bio-rad Scatter FC for ... ', scatter.plot.list[[bio.num]][[2]])
        scatterPlot( input.list.A    = tmp.bio.fc.list, 
                     input.list.B    = FALSE, 
                     plot.title.name = paste(dataset.name, 'Scatter Plot for', plot.titles.list['bio.fc.sum'][[1]], 'for', gsub('/', '-', scatter.plot.list[[bio.num]][[2]]),sep=' '), 
                     method.1        = scatter.plot.list[[bio.num]][[1]],
                     method.2        = scatter.plot.list[[bio.num]][[2]],
                     x.lab.name      = 'Bio-Rad', 
                     y.lab.name      = scatter.plot.list[[bio.num]][[2]],
                     alpha.val       = 1.0,
                     dataset.name    = dataset.name,
                     single.data     = TRUE)
    }
}


###############################################################################
# Scatter Plot of sd x mean counts
###############################################################################
if (mean.sd.scatter.plot.bool) {
    for (method in names(mean.all.norm.cnt.list)) {
        tmp.mean.A.list <- lapply(mean.all.norm.cnt.list, function(x) return(log2(x[['SEQC-A-samples']]+normalized.pseudo_count)))
        tmp.mean.B.list <- lapply(mean.all.norm.cnt.list, function(x) return(log2(x[['SEQC-B-samples']]+normalized.pseudo_count)))
        tmp.sd.A.list <- lapply(sd.norm.cnt.list, function(x) return(log2(x[['SEQC-A-samples']]+normalized.pseudo_count)))
        tmp.sd.B.list <- lapply(sd.norm.cnt.list, function(x) return(log2(x[['SEQC-B-samples']]+normalized.pseudo_count)))
        tmp.comb.A.list <- list(tmp.mean.A.list[[method]],tmp.sd.A.list[[method]])
        tmp.comb.B.list <- list(tmp.mean.B.list[[method]],tmp.sd.B.list[[method]])
        names(tmp.comb.A.list) <- c('mean','standard deviation')
        names(tmp.comb.B.list) <- c('mean','standard deviation')
        scatterPlot( input.list.A    = tmp.comb.A.list, 
                     input.list.B    = tmp.comb.B.list, 
                     plot.title.name = paste(dataset.name, 'Scatter Plot', plot.titles.list['mean.sd.norm.plot'][[1]], 'for', gsub('/','-',method), sep=' '), 
                     method.1        = 'mean',
                     method.2        = 'standard deviation',
                     x.lab.name      = 'log2 mean', 
                     y.lab.name      = 'log2 standard deviation',
                     alpha.val       = 1.0,
                     dataset.name    = dataset.name)
    }
}


###############################################################################
# Scatter Plot of count x count
###############################################################################
if (cnt.cnt.scatter.plot.bool) {
    if (SEQCxSEQS) {
        for (pos in 1:(length(names(SEQS.mean.all.norm.cnt.list))-1)) {
            method <- names(SEQS.mean.all.norm.cnt.list)[[pos]]
            tmp.SEQC.A.list <- lapply(SEQC.mean.all.norm.cnt.list, function(x) return(log2(x[['SEQC-A-samples']]+normalized.pseudo_count)))
            tmp.SEQC.B.list <- lapply(SEQC.mean.all.norm.cnt.list, function(x) return(log2(x[['SEQC-B-samples']]+normalized.pseudo_count)))
            tmp.SEQS.A.list <- lapply(SEQS.mean.all.norm.cnt.list, function(x) return(log2(x[['SEQC-A-samples']]+normalized.pseudo_count)))
            tmp.SEQS.B.list <- lapply(SEQS.mean.all.norm.cnt.list, function(x) return(log2(x[['SEQC-B-samples']]+normalized.pseudo_count)))
            tmp.comb.A.list <- list(tmp.SEQC.A.list[[method]],tmp.SEQS.A.list[[method]])
            tmp.comb.B.list <- list(tmp.SEQC.B.list[[method]],tmp.SEQS.B.list[[method]])
            names(tmp.comb.A.list) <- c('SEQC','SEQS')
            names(tmp.comb.B.list) <- c('SEQC','SEQS')
            scatterPlot( input.list.A    = tmp.comb.A.list, 
                         input.list.B    = tmp.comb.B.list, 
                         plot.title.name = paste('SEQCxSEQS', 'Scatter Plot', plot.titles.list['norm.cnt.cnt.plot'][[1]], 'for', gsub('/','-',method), sep=' '), 
                         method.1        = 'SEQC',
                         method.2        = 'SEQS',
                         x.lab.name      = paste('SEQC', method, sep=' '), 
                         y.lab.name      = paste('SEQS', method, sep=' '),
                         alpha.val       = 1.0,
                         dataset.name    = dataset.name)
        }
    } else {
        for (pos.1 in 1:(length(names(mean.all.norm.cnt.list))-1)) {
            for (pos.2 in (pos.1+1):length(names(mean.all.norm.cnt.list))) {
                method.1 <- names(mean.all.norm.cnt.list)[[pos.1]]
                method.2 <- names(mean.all.norm.cnt.list)[[pos.2]]
                tmp.cnt.A.list <- lapply(mean.all.norm.cnt.list, function(x) return(log2(x[['SEQC-A-samples']]+normalized.pseudo_count)))
                tmp.cnt.B.list <- lapply(mean.all.norm.cnt.list, function(x) return(log2(x[['SEQC-B-samples']]+normalized.pseudo_count)))
                tmp.comb.A.list <- list(tmp.cnt.A.list[[method.1]],tmp.cnt.A.list[[method.2]])
                tmp.comb.B.list <- list(tmp.cnt.B.list[[method.1]],tmp.cnt.B.list[[method.2]])
                names(tmp.comb.A.list) <- c(method.1,method.2)
                names(tmp.comb.B.list) <- c(method.1,method.2)
                scatterPlot( input.list.A    = tmp.comb.A.list, 
                             input.list.B    = tmp.comb.B.list, 
                             plot.title.name = paste(dataset.name, 'Scatter Plot', plot.titles.list['norm.cnt.cnt.plot'][[1]], 'for', gsub('/','-',method.1), 'vs', gsub('/','-',method.2), sep=' '), 
                             method.1        = method.1,
                             method.2        = method.2,
                             x.lab.name      = method.1, 
                             y.lab.name      = method.2,
                             alpha.val       = 1.0,
                             dataset.name    = dataset.name)
            }
        }
    }
}




###############################################################################
# Scatter Plot of FC x FC
###############################################################################
if (fc.fc.scatter.plot.bool) {
    if (SEQCxSEQS) {
        for (pos in 1:(length(names(SEQS.mean.all.norm.cnt.list))-1)) {
            method <- names(SEQS.mean.all.norm.cnt.list)[[pos]]
            tmp.SEQC.fc.list <- lapply(SEQC.mean.FC.list, function(x) return(log2(x+normalized.pseudo_count)))
            tmp.SEQS.fc.list <- lapply(SEQS.mean.FC.list, function(x) return(log2(x+normalized.pseudo_count)))
            tmp.comb.fc.list <- list(tmp.SEQC.fc.list[[method]],tmp.SEQS.fc.list[[method]])
            names(tmp.comb.fc.list) <- c('SEQC','SEQS')
            scatterPlot( input.list.A    = tmp.comb.fc.list, 
                         input.list.B    = FALSE, 
                         plot.title.name = paste('SEQCxSEQS', 'Scatter Plot', plot.titles.list['fc.fc.plot'][[1]], 'for', gsub('/','-',method), sep=' '), 
                         method.1        = 'SEQC',
                         method.2        = 'SEQS',
                         x.lab.name      = paste('SEQC', method, 'log2 FC', sep=' '), 
                         y.lab.name      = paste('SEQS', method, 'log2 FC', sep=' '),
                         alpha.val       = 1.0,
                         dataset.name    = dataset.name,
                         single.data     = TRUE)
        }
    } else {
        for (pos.1 in 1:(length(names(mean.FC.list))-1)) {
            for (pos.2 in (pos.1+1):length(names(mean.FC.list))) {
        #pos.1 <- length(names(mean.FC.list))
        #    for (pos.2 in 1:(length(names(mean.FC.list))-1)) {
                method.1 <- names(mean.FC.list)[[pos.1]]
                method.2 <- names(mean.FC.list)[[pos.2]]
                tmp.fc.list <- lapply(mean.FC.list, function(x) return(log2(x+normalized.pseudo_count)))
                tmp.comb.fc.list <- list(tmp.fc.list[[method.1]],tmp.fc.list[[method.2]])
                names(tmp.comb.fc.list) <- c(method.1,method.2)
                scatterPlot( input.list.A    = tmp.comb.fc.list, 
                             input.list.B    = FALSE, 
                             plot.title.name = paste(dataset.name, 'Scatter Plot', plot.titles.list['fc.fc.plot'][[1]], 'for', gsub('/','-',method.1), 'vs', gsub('/','-',method.2), sep=' '), 
                             method.1        = method.1,
                             method.2        = method.2,
                             x.lab.name      = paste(method.1, 'log2 FC', sep=' '), 
                             y.lab.name      = paste(method.2, 'log2 FC', sep=' '),
                             alpha.val       = 1.0,
                             dataset.name    = dataset.name,
                             single.data     = TRUE)
            }
        }
    }
}



###############################################################################
# Scatter Plot of fc x mean counts
###############################################################################
if (mean.fc.scatter.plot.bool) {
    for (method in names(mean.all.norm.cnt.list)) {
        tmp.mean.list <- lapply(mean.all.norm.cnt.list, function(x) return(0.5*(log2(x[['SEQC-A-samples']]+normalized.pseudo_count) + log2(x[['SEQC-B-samples']]+normalized.pseudo_count))))
        tmp.fc.list <- lapply(mean.FC.list, function(x) return(log2(x+normalized.pseudo_count)))
        tmp.comb.list <- list(tmp.mean.list[[method]],tmp.fc.list[[method]])
        names(tmp.comb.list) <- c('mean','fold change')
        scatterPlot( input.list.A    = tmp.comb.list, 
                     input.list.B    = FALSE, 
                     plot.title.name = paste(dataset.name, 'Scatter Plot', plot.titles.list['mean.fc.plot'][[1]], 'for', gsub('/','-',method), sep=' '), 
                     method.1        = 'mean',
                     method.2        = 'fold change',
                     x.lab.name      = '1/2 log2 A*B', 
                     y.lab.name      = 'log2 fold change',
                     alpha.val       = 1.0,
                     dataset.name    = dataset.name,
                     single.data     = TRUE)
    }
}



###############################################################################
# Scatter Plot of fc x mean counts
###############################################################################
if (cnt.fc.scatter.plot.bool) {
    for (method in names(mean.all.norm.cnt.list)) {
        tmp.cnt.A.list <- lapply(mean.all.norm.cnt.list, function(x) return(log2(x[['SEQC-A-samples']]+normalized.pseudo_count)))
        tmp.cnt.B.list <- lapply(mean.all.norm.cnt.list, function(x) return(log2(x[['SEQC-B-samples']]+normalized.pseudo_count)))
        tmp.fc.A.list <- lapply(mean.FC.list, function(x) return(log2(x+normalized.pseudo_count)))
        tmp.fc.B.list <- lapply(mean.FC.list, function(x) return(log2(x+normalized.pseudo_count)))
        tmp.comb.A.list <- list(tmp.cnt.A.list[[method]],tmp.fc.A.list[[method]])
        tmp.comb.B.list <- list(tmp.cnt.B.list[[method]],tmp.fc.B.list[[method]])
        names(tmp.comb.A.list) <- c('count','fold change')
        names(tmp.comb.B.list) <- c('count','fold change')
        scatterPlot( input.list.A    = tmp.comb.A.list, 
                     input.list.B    = tmp.comb.B.list, 
                     plot.title.name = paste(dataset.name, 'Scatter Plot', plot.titles.list['cnt.fc.plot'][[1]], 'for', gsub('/','-',method), sep=' '), 
                     method.1        = 'count',
                     method.2        = 'fold change',
                     x.lab.name      = 'log2 count', 
                     y.lab.name      = 'log2 fold change',
                     alpha.val       = 1.0,
                     dataset.name    = dataset.name,
                     single.data     = FALSE)
    }
}



###############################################################################
# Heatmaps
###############################################################################
if (cor.heatmap.plot.bool) {
    #
    legend.title.list <- list()
    legend.title.list[['SEQC.cnt.A.spearman']] <- c('HeatMap SEQC Count A','Spearman Correlation')
    legend.title.list[['SEQC.cnt.A.pearson']] <- c('HeatMap SEQC Count A','Pearson Correlation')
    legend.title.list[['SEQC.cnt.B.spearman']] <- c('HeatMap SEQC Count B','Spearman Correlation')
    legend.title.list[['SEQC.cnt.B.pearson']] <- c('HeatMap SEQC Count B','Pearson Correlation')
    legend.title.list[['SEQS.cnt.A.spearman']] <- c('HeatMap SEQS Count A','Spearman Correlation')
    legend.title.list[['SEQS.cnt.A.pearson']] <- c('HeatMap SEQS Count A','Pearson Correlation')
    legend.title.list[['SEQS.cnt.B.spearman']] <- c('HeatMap SEQS Count B','Spearman Correlation')
    legend.title.list[['SEQS.cnt.B.pearson']] <- c('HeatMap SEQS Count B','Pearson Correlation')
    legend.title.list[['SEQC.SEQS.A.cnt.spearman']] <- c('HeatMap SEQCxSEQS Count A','Spearman Correlation')
    legend.title.list[['SEQC.SEQS.A.cnt.pearson']] <- c('HeatMap SEQCxSEQS Count A','Pearson Correlation')
    legend.title.list[['SEQC.SEQS.B.cnt.spearman']] <- c('HeatMap SEQCxSEQS Count B','Spearman Correlation')
    legend.title.list[['SEQC.SEQS.B.cnt.pearson']] <- c('HeatMap SEQCxSEQS Count B','Pearson Correlation')
    legend.title.list[['taq.cnt.A.spearman']] <- c('HeatMap TaqMan Count A','Spearman Correlation')
    legend.title.list[['taq.cnt.A.pearson']] <- c('HeatMap TaqMan Count A','Pearson Correlation')
    legend.title.list[['taq.cnt.B.spearman']] <- c('HeatMap TaqMan Count B','Spearman Correlation')
    legend.title.list[['taq.cnt.B.pearson']] <- c('HeatMap TaqMan Count B','Pearson Correlation')
    legend.title.list[['bio.cnt.A.spearman']] <- c('HeatMap Bio-Rad Count A','Spearman Correlation')
    legend.title.list[['bio.cnt.A.pearson']] <- c('HeatMap Bio-Rad Count A','Pearson Correlation')
    legend.title.list[['bio.cnt.B.spearman']] <- c('HeatMap Bio-Rad Count B','Spearman Correlation')
    legend.title.list[['bio.cnt.B.pearson']] <- c('HeatMap Bio-Rad Count B','Pearson Correlation')
    legend.title.list[['SEQC.fc.spearman']] <- c('HeatMap SEQC FC','Spearman Correlation')
    legend.title.list[['SEQC.fc.pearson']] <- c('HeatMap SEQC FC','Pearson Correlation')
    legend.title.list[['SEQS.fc.spearman']] <- c('HeatMap SEQS FC','Spearman Correlation')
    legend.title.list[['SEQS.fc.pearson']] <- c('HeatMap SEQS FC','Pearson Correlation')
    legend.title.list[['SEQC.SEQS.fc.spearman']] <- c('HeatMap SEQCxSEQS FC','Spearman Correlation')
    legend.title.list[['SEQC.SEQS.fc.pearson']] <- c('HeatMap SEQCxSEQS FC','Pearson Correlation')
    legend.title.list[['taq.fc.spearman']] <- c('HeatMap TaqMan FC','Spearman Correlation')
    legend.title.list[['taq.fc.pearson']] <- c('HeatMap TaqMan FC','Pearson Correlation')
    legend.title.list[['bio.fc.spearman']] <- c('HeatMap Bio-Rad FC','Spearman Correlation')
    legend.title.list[['bio.fc.pearson']] <- c('HeatMap Bio-Rad FC','Pearson Correlation')
    legend.title.list[['SEQC.log.fc.spearman']] <- c('HeatMap SEQC log2 FC','Spearman Correlation')
    legend.title.list[['SEQC.log.fc.pearson']] <- c('HeatMap SEQC log2 FC','Pearson Correlation')
    legend.title.list[['SEQS.log.fc.spearman']] <- c('HeatMap SEQS log2 FC','Spearman Correlation')
    legend.title.list[['SEQS.log.fc.pearson']] <- c('HeatMap SEQS log2 FC','Pearson Correlation')
    legend.title.list[['SEQC.SEQS.log.fc.spearman']] <- c('HeatMap SEQCxSEQS log2 FC','Spearman Correlation')
    legend.title.list[['SEQC.SEQS.log.fc.pearson']] <- c('HeatMap SEQCxSEQS log2 FC','Pearson Correlation')
    legend.title.list[['taq.log.fc.spearman']] <- c('HeatMap TaqMan log2 FC','Spearman Correlation')
    legend.title.list[['taq.log.fc.pearson']] <- c('HeatMap TaqMan log2 FC','Pearson Correlation')
    legend.title.list[['bio.log.fc.spearman']] <- c('HeatMap Bio-Rad log2 FC','Spearman Correlation')
    legend.title.list[['bio.log.fc.pearson']] <- c('HeatMap Bio-Rad log2 FC','Pearson Correlation')
    legend.title.list[['SEQC.cnt.A.sd']] <- c('HeatMap SEQC Count A','ABS log2 SD - log2 SD')
    legend.title.list[['SEQC.cnt.B.sd']] <- c('HeatMap SEQC Count B','ABS log2 SD - log2 SD')
    legend.title.list[['SEQS.cnt.A.sd']] <- c('HeatMap SEQS Count A','ABS log2 SD - log2 SD')
    legend.title.list[['SEQS.cnt.B.sd']] <- c('HeatMap SEQS Count B','ABS log2 SD - log2 SD')
    legend.title.list[['SEQC.SEQS.A.cnt.sd']] <- c('HeatMap SEQCxSEQS Count A','ABS log2 SD - log2 SD')
    legend.title.list[['SEQC.SEQS.B.cnt.sd']] <- c('HeatMap SEQCxSEQS Count B','ABS log2 SD - log2 SD')
    legend.title.list[['taq.cnt.A.sd']] <- c('HeatMap TaqMan Count A','ABS log2 SD - log2 SD')
    legend.title.list[['taq.cnt.B.sd']] <- c('HeatMap TaqMan Count B','ABS log2 SD - log2 SD')
    legend.title.list[['bio.cnt.A.sd']] <- c('HeatMap Bio-Rad Count A','ABS log2 SD - log2 SD')
    legend.title.list[['bio.cnt.B.sd']] <- c('HeatMap Bio-Rad Count B','ABS log2 SD - log2 SD')
    legend.title.list[['SEQC.fc.sd']] <- c('HeatMap SEQC FC','ABS log2 SD - log2 SD')
    legend.title.list[['SEQS.fc.sd']] <- c('HeatMap SEQS FC','ABS log2 SD - log2 SD')
    legend.title.list[['SEQC.SEQS.fc.sd']] <- c('HeatMap SEQCxSEQS FC','ABS log2 SD - log2 SD')
    legend.title.list[['taq.fc.sd']] <- c('HeatMap TaqMan FC','ABS log2 SD - log2 SD')
    legend.title.list[['bio.fc.sd']] <- c('HeatMap Bio-Rad FC','ABS log2 SD - log2 SD')
    legend.title.list[['SEQC.log.fc.sd']] <- c('HeatMap SEQC log2 FC','ABS log2 SD - log2 SD')
    legend.title.list[['SEQS.log.fc.sd']] <- c('HeatMap SEQS log2 FC','ABS log2 SD - log2 SD')
    legend.title.list[['SEQC.SEQS.log.fc.sd']] <- c('HeatMap SEQCxSEQS log2 FC','ABS log2 SD - log2 SD')
    legend.title.list[['taq.log.fc.sd']] <- c('HeatMap TaqMan log2 FC','ABS log2 SD - log2 SD')
    legend.title.list[['bio.log.fc.sd']] <- c('HeatMap Bio-Rad log2 FC','ABS log2 SD - log2 SD')
    
    for (method.type in names(stat.correlation.list)) {
        heatMap( input.list        = stat.correlation.list[[method.type]],
                 plot.title.name   = paste(legend.title.list[[method.type]][[1]],legend.title.list[[method.type]][[2]],sep=' '), 
                 legend.title      = paste(legend.title.list[[method.type]][[1]],legend.title.list[[method.type]][[2]],sep='\n'), 
                 x.lab.name        = '',
                 y.lab.name        = '',
                 make.new.matrix   = FALSE)
    }
    
    
}



###############################################################################
# Scatter Plot of sd x mean counts
###############################################################################
if (mean.sd.fc.scatter.plot.bool) {
    for (method in names(mean.all.norm.cnt.list)) {
        tmp.mean.list <- lapply(sd.mean.fc.list, function(x) return(log2(x[['mean FC']]+normalized.pseudo_count)))
        tmp.sd.list <- lapply(sd.mean.fc.list, function(x) return(log2(x[['SD FC']]+normalized.pseudo_count)))
        tmp.comb.list <- list(tmp.mean.list[[method]],tmp.sd.list[[method]])
        names(tmp.comb.list) <- c('mean','standard deviation')
        scatterPlot( input.list.A    = tmp.comb.list, 
                     input.list.B    = FALSE, 
                     plot.title.name = paste(dataset.name, 'Scatter Plot', plot.titles.list['mean.sd.fc.plot'][[1]], 'for', gsub('/','-',method), sep=' '), 
                     method.1        = 'mean',
                     method.2        = 'standard deviation',
                     x.lab.name      = 'log2 mean fold change', 
                     y.lab.name      = 'log2 standard deviation fold change',
                     alpha.val       = 1.0,
                     dataset.name    = dataset.name)
    }
}


###############################################################################
# Scatter Plot of mean counts
###############################################################################
if (biotaq.scatter.plot.bool) {
    tmp.mean.A.list <- lapply(mean.all.norm.cnt.list, function(x) return(log2(x[['SEQC-A-samples']]+normalized.pseudo_count)))
    tmp.mean.B.list <- lapply(mean.all.norm.cnt.list, function(x) return(log2(x[['SEQC-B-samples']]+normalized.pseudo_count)))
    make.print.bool <- TRUE
    make.pdf.bool <- FALSE
    scatterPlot( input.list.A    = tmp.mean.A.list, 
                 input.list.B    = tmp.mean.B.list, 
                 plot.title.name = paste('Scatter Plot for', plot.titles.list['mean.all.norm.cnt'][[1]],sep=' '), 
                 method.1        = names(tmp.mean.A.list)[[1]],
                 method.2        = names(tmp.mean.B.list)[[1]],
                 x.lab.name      = '', 
                 y.lab.name      = '')
}

###############################################################################
# FC Gene Expression
###############################################################################
if (FC.genes.mean.plot.bool) {
    # plot mean gene expression
    boxPlot( input.list      = mean.FC.list, 
             plot.title.name = plot.titles.list['mean.log.FC'][[1]], 
             x.lab.name      = 'Counting Method', 
             y.lab.name      = plot.titles.list['mean.log.FC'][[1]],
             reduce.jitter   = FALSE,
             use.jitter      = FALSE)
}


###############################################################################
# FC gold
###############################################################################
if (FC.gold.plot.bool) {
    b.plot.list <- list()
    b.plot.list[['TQMAN-default/QPCR']] <- abs(log2(taqman.FC.list[['TQMAN-default/QPCR']]+normalized.pseudo_count))
    b.plot.list[['BIORD-default/QPCR']] <- abs(log2(biorad.FC.list[['BIORD-default/QPCR']]+normalized.pseudo_count))
    b.plot.list[['FLXSM-default/SIMED']] <- abs(log2(SEQS.mean.FC.list[['FLXSM-default/SIMED']]+normalized.pseudo_count))
    # plot mean gene expression
    boxPlot( input.list      = b.plot.list, 
             plot.title.name = 'Gold Standard log2 FC', 
             x.lab.name      = 'Gold Standard', 
             y.lab.name      = 'log2 FC',
             reduce.jitter   = FALSE,
             use.jitter      = FALSE,
             x.angle.setting = 0,
             v.just.setting  = 0,
             h.just.setting  = 0.5,
             scale.factor.setting = 0.6)
}


###############################################################################
## intra FC box plot
## 'intra.met.pearson.fc'      =  'Intra Method FC Pearson CC'
## 'intra.met.spearman.fc'     =  'Intra Method FC Spearman CC'
## 'intra.met.CD.fc'           =  'Intra Method FC CD'
## 'intra.met.pearson.log.fc'  =  'Intra Method Log2 FC Pearson CC'
## 'intra.met.spearman.log.fc' =  'Intra Method Log2 FC Spearman CC'
## 'intra.met.CD.log.fc'       =  'Intra Method Log2 FC CD'
###############################################################################
if (intra.fc.box.plot.bool) {
    # plot mean gene expression
    for (sample.type in sample.types) {
        expression.plot.list <- list()
        expression.plot.list[['intra.met.pearson.fc']] <- lapply(names(intra.method.fc.cor.list), function(x) return(intra.method.fc.cor.list[[x]][['intra.met.pearson.fc']]))
        print(names(expression.plot.list[['intra.met.pearson.fc']]))
        for (plot.name in names(expression.plot.list)) {
            boxPlot( input.list      = expression.plot.list[[plot.name]], 
                     plot.title.name = paste(plot.titles.list[plot.name][[1]], sample.naming.vector[[sample.type]]), 
                     x.lab.name      = 'Counting Method', 
                     y.lab.name      = paste(plot.titles.list[plot.name][[1]], sample.naming.vector[[sample.type]]),
                     reduce.jitter   = FALSE,
                     use.jitter      = TRUE)
         }
    }
}

###############################################################################
# Coefficient of Variation among Genes
###############################################################################
if (CV.list.plot.bool) {
    # plot results
    cv.naming.vector <- c(cv.mean.list='(mean(A,B))', cv.A.list='(A)', cv.B.list='(B)')
    for (cv.list.type in c('cv.mean.list', 'cv.A.list', 'cv.B.list')) {
        tmp.cv.list <- lapply(CV.list, function(x) return(x[[cv.list.type]]))
        boxPlot( input.list      = tmp.cv.list, 
                 plot.title.name = paste(plot.titles.list['intra.gene.CV'][[1]], cv.naming.vector[cv.list.type], sep=' '), 
                 x.lab.name      = 'Counting Method', 
                 y.lab.name      = paste(plot.titles.list['intra.gene.CV'][[1]], cv.naming.vector[cv.list.type], sep=' '))
    }
}

###############################################################################
# Coefficient of Variation among Genes of Fold Change of A and B samples
###############################################################################
if (fcCV.list.plot.bool) {
    # plot results
    fcCV.naming.vector <- c(intra.gene.CV.FC='(CV(FC))', intra.gene.CV.log.FC='(CV(log2(FC)))')
    #for (fcCV.list.type in c('intra.gene.CV.FC', 'intra.gene.CV.log.FC')) {
    for (fcCV.list.type in c('intra.gene.CV.log.FC')) {
        tmp.fcCV.list <- lapply(fcCV.list, function(x) return(x[[fcCV.list.type]]))
        boxPlot( input.list      = tmp.fcCV.list, 
                 plot.title.name = paste(plot.titles.list[fcCV.list.type][[1]]), 
                 x.lab.name      = 'counting method', 
                 y.lab.name      = paste(plot.titles.list[fcCV.list.type][[1]]),
                 use.axis.scale  = TRUE)
    }
     
    # remove big values
    tmp.fcCV.list <- lapply(fcCV.list, function(x) return(x[['intra.gene.CV.log.FC']]))
    norm.fcCV.list <- lapply(tmp.fcCV.list, function(tmp.vector) {x
        tmp.vector <- na.omit(tmp.vector)
        tmp.vector <- unlist(lapply(tmp.vector, function(tmp.entry) {
            if (abs(tmp.entry) < VClogFC.print.thresh) return(tmp.entry)}))
            #if (tmp.entry >= 0) return(tmp.entry)}))
        return(tmp.vector)})
    big.fcCV.list <- lapply(tmp.fcCV.list, function(tmp.vector) {
        tmp.vector <- na.omit(tmp.vector)
        tmp.vector <- unlist(lapply(tmp.vector, function(tmp.entry) {
            if (abs(tmp.entry) >= VClogFC.print.thresh) return(tmp.entry)}))
        return(tmp.vector)})
    boxPlot( input.list      = norm.fcCV.list,
             plot.title.name = paste(plot.titles.list['intra.gene.CV.log.FC'][[1]], 'reduced', sep='-'), 
             x.lab.name      = 'counting method', 
             y.lab.name      = paste(plot.titles.list['intra.gene.CV.log.FC'][[1]]),
             use.axis.scale  = TRUE)
    
}
 
###############################################################################
# Gene Count Correlation Coefficient for each method and sample type (sample x sample)
###############################################################################
 if (intra.method.cnt.cor.list.plot.bool) {
    # plot combined results
    legend.title.vector.short <- c(intra.met.pearson               = 'Pearson CC', 
                                   intra.met.spearman              = 'Spearman CC', 
                                   intra.met.CD                    = 'CD')
    for (precision.method in names(intra.method.cnt.cor.list[[1]])) {
        for (sample.type in names(intra.method.cnt.cor.list[[1]][[precision.method]])) {
            tmp.vector <- unlist(lapply(intra.method.cnt.cor.list, function(x) return(x[[precision.method]][[sample.type]])))
            names(tmp.vector) <- unlist(lapply(names(tmp.vector), function(x) return(paste(strsplit(x,'.',fixed=TRUE)[[1]][[2]],
                                                                                           strsplit(x,'.',fixed=TRUE)[[1]][[3]],
                                                                                           strsplit(x,'.',fixed=TRUE)[[1]][[1]], sep='.'))))
            multiHeatMap( input.list      = tmp.vector,
                          plot.title.name = paste('Combined', plot.titles.list[precision.method][[1]], sample.naming.vector[[sample.type]], sep=' '), 
                          legend.title    = legend.title.vector.short[[precision.method]], 
                          lab.name        = 'quantifier-option/aligner')
        }
    }
    # plot results
    legend.title.vector <- c(intra.met.pearson='Pearson Correlation\nCoefficient', intra.met.spearman='Spearman Correlation\n Coefficient', intra.met.CD='Coefficient of\nDetermination')
    for (method.name in names(intra.method.cnt.cor.list)) {
        for (precision.method in names(intra.method.cnt.cor.list[[method.name]])) {
            for (sample.type in names(intra.method.cnt.cor.list[[method.name]][[precision.method]])) {
                heatMap( input.list      = intra.method.cnt.cor.list[[method.name]][[precision.method]][[sample.type]], 
                         plot.title.name = paste(plot.titles.list[precision.method][[1]], 'for', gsub('/', '-', method.name), sample.naming.vector[[sample.type]]), 
                         legend.title    = legend.title.vector[[precision.method]], 
                         x.lab.name      = 'sample name',
                         y.lab.name      = 'sample name')
            }
        }
    }
}

###############################################################################
# Correlation Coefficient of Fold Change of A and B samples for each method (FC(AB) x FC(AB))
###############################################################################
if (intra.method.fc.cor.list.plot.bool) {
    # plot combined results
    legend.title.vector.short <- c(intra.met.pearson.fc         = 'FC Pearson CC', 
                             intra.met.spearman.fc              = 'FC Spearman CC', 
                             intra.met.CD.fc                    = 'FC CD', 
                             intra.met.pearson.log.fc           = 'Log2 FC Pearson CC', 
                             intra.met.spearman.log.fc          = 'Log2 FC Spearman CC', 
                             intra.met.CD.log.fc                = 'Log2 FC CD')
    for (precision.method in names(intra.method.fc.cor.list[[1]])) {
        tmp.vector <- unlist(lapply(intra.method.fc.cor.list, function(x) return(x[[precision.method]])))
        names(tmp.vector) <- unlist(lapply(names(tmp.vector), function(x) return(paste(strsplit(x,'.',fixed=TRUE)[[1]][[2]],
                                                                                       strsplit(x,'.',fixed=TRUE)[[1]][[3]],
                                                                                       strsplit(x,'.',fixed=TRUE)[[1]][[1]], sep='.'))))
        multiHeatMap( input.list      = tmp.vector,
                      plot.title.name = paste('Combined', plot.titles.list[precision.method][[1]], sep=' '), 
                      legend.title    = legend.title.vector.short[[precision.method]], 
                      lab.name        = 'quantifier-option/aligner')
    }
    # plot results
    legend.title.vector <- c(intra.met.pearson.fc      ='Fold Change Pearson\nCorrelation Coefficient', 
                             intra.met.spearman.fc     ='Fold Change Spearman\nCorrelation Coefficient', 
                             intra.met.CD.fc           ='Fold Change Coefficient\nof Determination', 
                             intra.met.pearson.log.fc  ='Log2 Fold Change Pearson\nCorrelation Coefficient', 
                             intra.met.spearman.log.fc ='Log2 Fold Change Spearman\nCorrelation Coefficient', 
                             intra.met.CD.log.fc       ='Log2 Fold Change Coefficient\nof Determination')
    for (method.name in names(intra.method.fc.cor.list)) {
        for (precision.method in names(intra.method.fc.cor.list[[method.name]])) {
            heatMap( input.list      = intra.method.fc.cor.list[[method.name]][[precision.method]], 
                     plot.title.name = paste(plot.titles.list[precision.method][[1]], 'for', gsub('/', '-', method.name)), 
                     legend.title    = legend.title.vector[[precision.method]], 
                     x.lab.name      = 'sample name',
                     y.lab.name      = 'sample.name')
        }
    }
}

###############################################################################
# Correlation Coefficient of Fold Change of A and B samples between methods (method(FC(AB)) x method(FC(AB)))
###############################################################################
if (inter.method.fc.cor.list.plot.bool) {
    # plot combined results
    legend.title.vector.short <- c(inter.met.pearson.fc      = 'FC Pearson CC', 
                             inter.met.spearman.fc           = 'FC Spearman CC', 
                             inter.met.CD.fc                 = 'FC CD', 
                             inter.met.pearson.log.fc        = 'Log2 FC Pearson CC', 
                             inter.met.spearman.log.fc       = 'Log2 FC Spearman CC', 
                             inter.met.CD.log.fc             = 'Log2 FC CD', 
                             inter.met.false.negative.rate   = 'SG FNR', 
                             inter.met.sensitivity           = 'SG Sensitivity', 
                             inter.met.false.positive.rate   = 'SG FPR', 
                             inter.met.specificity           = 'SG Specificity', 
                             inter.met.false.discovery.rate  = 'SG FDR', 
                             inter.met.precision             = 'SG Precision')
    for (inter.method.type in names(inter.method.fc.cor.list[[1]][[1]])) {
        multiHeatMap( input.list      = lapply(inter.method.fc.cor.list, function(x) return(lapply(x, function(y) return(y[[inter.method.type]])))), 
                      plot.title.name = paste('Combined', plot.titles.list[inter.method.type][[1]], sep=' '), 
                      legend.title    = legend.title.vector.short[[inter.method.type]], 
                      lab.name        = 'quantifier-option/aligner')
    }
    # plot single results
    legend.title.vector <- c(inter.met.pearson.fc      = 'Fold Change Pearson\nCorrelation Coefficient', 
                             inter.met.spearman.fc     = 'Fold Change Spearman\nCorrelation Coefficient', 
                             inter.met.CD.fc           = 'Fold Change Coefficient\nof Determination', 
                             inter.met.pearson.log.fc  = 'Log2 Fold Change Pearson\nCorrelation Coefficient', 
                             inter.met.spearman.log.fc = 'Log2 Fold Change Spearman\nCorrelation Coefficient', 
                             inter.met.CD.log.fc       = 'Log2 Fold Change Coefficient\nof Determination',
                             inter.met.false.negative.rate  = 'SG False Negative Rate', 
                             inter.met.sensitivity          = 'SG Sensitivity', 
                             inter.met.false.positive.rate  = 'SG False Positive Rate', 
                             inter.met.specificity          = 'SG Specificity', 
                             inter.met.false.discovery.rate = 'SG False Discovery Rate', 
                             inter.met.precision            = 'SG Precision')
    triangle.vector <- c( inter.met.pearson.fc      = TRUE, 
                          inter.met.spearman.fc     = TRUE, 
                          inter.met.CD.fc           = TRUE,
                          inter.met.pearson.log.fc  = TRUE,
                          inter.met.spearman.log.fc = TRUE,
                          inter.met.CD.log.fc       = TRUE,
                          inter.met.false.negative.rate  = FALSE,
                          inter.met.sensitivity          = FALSE,
                          inter.met.false.positive.rate  = FALSE,
                          inter.met.specificity          = FALSE, 
                          inter.met.false.discovery.rate = FALSE,
                          inter.met.precision            = FALSE)
    for (inter.method.type in names(inter.method.fc.cor.list[[1]][[1]])) {
        for (sample.name in names(inter.method.fc.cor.list[[1]][[1]][[1]])) {
            heatMap( input.list        = lapply(inter.method.fc.cor.list, function(x) 
                                           return(unlist(lapply(x, function(y) 
                                             return(y[[inter.method.type]][[sample.name]]))))), 
                     plot.title.name   = paste(plot.titles.list[[inter.method.type]], 'for', gsub('/', '#', sample.name), sep=' '), 
                     legend.title      = legend.title.vector[[inter.method.type]], 
                     x.lab.name        = 'Null Hypotheses',
                     y.lab.name        = 'Alternative Hypotheses',
                     use.triangle.bool = triangle.vector[[inter.method.type]])
        }
    }
}

###############################################################################
# Correlation Coefficient of Gene Expression between methods for each sample (method(sample)) x method(sample))
###############################################################################
if (inter.method.cnt.cor.list.plot.bool) {
    # plot combined results
    legend.title.vector.short <- c(inter.met.pearson         = 'Pearson CC', 
                             inter.met.spearman              = 'Spearman CC', 
                             inter.met.CD                    = 'CD', 
                             inter.met.log.pearson           = 'Log2 Pearson CC', 
                             inter.met.log.spearman          = 'Log2 Spearman CC', 
                             inter.met.log.CD                = 'Log2 CD')
    for (inter.method.type in names(inter.method.cnt.cor.list[[1]][[1]])) {
        multiHeatMap( input.list      = lapply(inter.method.cnt.cor.list, function(x) return(lapply(x, function(y) return(y[[inter.method.type]])))), 
                      plot.title.name = paste('Combined', plot.titles.list[inter.method.type][[1]], sep=' '), 
                      legend.title    = legend.title.vector.short[[inter.method.type]], 
                      lab.name        = 'quantifier-option/aligner')
    }
    # plot single results
    legend.title.vector <- c(inter.met.pearson      = 'Gene Expression Count Pearson\nCorrelation Coefficient', 
                             inter.met.spearman     = 'Gene Expression Count Spearman\nCorrelation Coefficient', 
                             inter.met.CD           = 'Gene Expression Count\nCoefficient of Determination', 
                             inter.met.log.pearson  = 'Log2 Gene Expression Count\nPearson Correlation Coefficient', 
                             inter.met.log.spearman = 'Log2 Gene Expression Count\nSpearman Correlation Coefficient', 
                             inter.met.log.CD       = 'Log2 Gene Expression Count\nCoefficient of Determination')
    for (inter.method.type in names(inter.method.cnt.cor.list[[1]][[1]])) {
        for (sample.name in names(inter.method.cnt.cor.list[[1]][[1]][[1]])) {
            heatMap( input.list        = lapply(inter.method.cnt.cor.list, function(x) 
                                           return(unlist(lapply(x, function(y) 
                                             return(y[[inter.method.type]][[sample.name]]))))), 
                     plot.title.name   = paste(plot.titles.list[[inter.method.type]], 'for', gsub('/', '#', sample.name), sep=' '), 
                     legend.title      = legend.title.vector[[inter.method.type]], 
                     x.lab.name        = '',
                     y.lab.name        = '')
        }
    }
}

###############################################################################
# Scatter Plot of CV.list
###############################################################################
if (CV.list.scatter.plot.bool) {
    # plot results
    tmp.cv.A.list <- lapply(CV.list, function(x) return(x[['cv.A.list']]))
    tmp.cv.B.list <- lapply(CV.list, function(x) return(x[['cv.B.list']]))
    scatter.plot.list <- list(scat.pl.1=c('HTSEQ-union/TOPHAT', 'FEATC-default/TOPHAT'),
                              scat.pl.2=c('EQPQM-ambig/TOPHAT', 'HTSEQ-union/TOPHAT'),
                              scat.pl.3=c('FEATC-default/TOPHAT', 'EQPQM-ambig/TOPHAT'),
                              scat.pl.4=c('HTSEQ-union/HISAT', 'FEATC-default/HISAT'),
                              scat.pl.5=c('EQPQM-ambig/HISAT', 'HTSEQ-union/HISAT'),
                              scat.pl.6=c('FEATC-default/HISAT', 'EQPQM-ambig/HISAT'),
                              scat.pl.7=c('HTSEQ-union/STAR', 'FEATC-default/STAR'),
                              scat.pl.8=c('EQPQM-ambig/STAR', 'HTSEQ-union/STAR'),
                              scat.pl.9=c('FEATC-default/STAR', 'EQPQM-ambig/STAR'))
    for (sc.pl in names(scatter.plot.list)) {
        scatterPlot( input.list.A    = tmp.cv.A.list, 
                     input.list.B    = tmp.cv.B.list, 
                     plot.title.name = paste('Scatter Plot for', plot.titles.list['intra.gene.CV'][[1]], 'for', gsub('/', '-', scatter.plot.list[[sc.pl]][[1]]) , 'and' , gsub('/', '-', scatter.plot.list[[sc.pl]][[2]]),sep=' '), 
                     method.1        = scatter.plot.list[[sc.pl]][[1]],
                     method.2        = scatter.plot.list[[sc.pl]][[2]],
                     x.lab.name      = '', 
                     y.lab.name      = '')
    }
}

###############################################################################
# Scatter Plot of FC(A,B) vs mean(A,B)
###############################################################################

if (FC.mean.scatter.plot.bool) {
    # plot results
    for (i in 1:length(samples.matched.A)) {
        sample.type.A <- samples.matched.A[[i]]
        sample.type.B <- samples.matched.B[[i]]
        sample.AB <- paste(sample.type.A, sample.type.B, sep='/')
        color.vector <- c(brewer.pal(9,'Set1'), "#000000")
        plot.title.name <- paste('Scatter Plot for', gsub('/', '#', sample.AB), sep=' ')
        # make plot
        meanFCScatterPlot( FC.list          = FC.list, 
                           count.list       = normalized.count.list, 
                           sample.type.A    = sample.type.A, 
                           sample.type.B    = sample.type.B, 
                           method.occurence = 9, 
                           color.vector     = color.vector, 
                           plot.title.name  = plot.title.name)
    }
}

###############################################################################
# Scatter Plot of mean over all samples of FC(A,B) vs mean(A,B)
###############################################################################

if (FC.mean.scatter.all.plot.bool) {
    # calculate mean for each gene
    mean.normalized.list <- list()
    for (method.n in names(normalized.count.list)) {
        SEQC.A.samples <- apply(normalized.count.list[[method.n]][,samples.matched.A], 1, function(x) return(sum(x)/length(x)))
        SEQC.B.samples <- apply(normalized.count.list[[method.n]][,samples.matched.B], 1, function(x) return(sum(x)/length(x)))
        gene.id.frame <- normalized.count.list[[method.n]][,'Gene Id']
        mean.normalized.list[[method.n]] <- data.frame(gene.id.frame, SEQC.A.samples, SEQC.B.samples)
        names(mean.normalized.list[[method.n]]) <- as.factor(c('Gene Id', 'SEQC-A-samples', 'SEQC-B-samples'))
    }
    sample.type.A <- 'SEQC-A-samples'
    sample.type.B <- 'SEQC-B-samples'
    mean.FC.list <- list()
    for (method.n in names(mean.normalized.list)) {
        FC.matrix <- computeFCs( input.frame       = mean.normalized.list[[method.n]], 
                                 input.name        = names(mean.normalized.list)[method.n], 
                                 samples.matched.A = c(sample.type.A), 
                                 samples.matched.B = c(sample.type.B), 
                                 pseudo_count      = normalized.pseudo_count )
        mean.FC.list[[method.n]] <- FC.matrix
    }
    # make plot list
    scatter.plot.list <- list()
    for (i in c.list.aligners) {
        scatter.names <- c()
        for (j in c.list.quantifiers) {
            scatter.names <- c(scatter.names, paste(j,i,sep='/'))
        }
        scatter.plot.list[[i]] <- scatter.names
    }
    for (i in c.list.quantifiers) {
        scatter.names <- c()
        for (j in c.list.aligners) {
            scatter.names <- c(scatter.names, paste(i,j,sep='/'))
        }
        scatter.plot.list[[i]] <- scatter.names
    }
    # plot results
    sample.AB <- paste(sample.type.A, sample.type.B, sep='/')
    multi.occ.mean.frame.list <- list()
    for (scatter.name in names(scatter.plot.list)) {
        plot.title.name <- paste('Scatter Plot FC x mean for', gsub('/', '#', sample.AB), 'with', scatter.name, sep=' ')
        color.vector <- c(brewer.pal(length(scatter.plot.list[[scatter.name]]),'Set1'), "#000000")
        # make plots
        multi.occ.mean.frame.list[[scatter.name]] <- meanFCScatterPlot( m.FC.list         = mean.FC.list, 
                                                                        m.norm.count.list = mean.normalized.list, 
                                                                        m.sample.type.A   = sample.type.A, 
                                                                        m.sample.type.B   = sample.type.B, 
                                                                        method.occurence  = length(scatter.plot.list[[scatter.name]]), 
                                                                        color.vector      = color.vector, 
                                                                        plot.title.name   = plot.title.name,
                                                                        round.to.xth      = 100,
                                                                        method.names      = scatter.plot.list[[scatter.name]],
                                                                        sample.thresh     = 5)
    }
    
    # variance
    var.FC.list <- list()
    for (method.n in names(mean.normalized.list)) {
        log.FC.list <- log2(FC.list[[method.n]])
        SEQC.A.samples <- apply(log.FC.list, 1, function(x) return(var(x)))
        SEQC.B.samples <- SEQC.A.samples
        gene.id.frame <- normalized.count.list[[method.n]][,'Gene Id']
        var.FC.list[[method.n]] <- data.frame(gene.id.frame, SEQC.A.samples, SEQC.B.samples)
        names(var.FC.list[[method.n]]) <- as.factor(c('Gene Id', 'SEQC-A-samples', 'SEQC-B-samples'))
    }
    multi.occ.var.frame.list <- list()
    for (scatter.name in names(scatter.plot.list)) {
        plot.title.name <- paste('Scatter Plot FC x variance for', gsub('/', '#', sample.AB), 'with', scatter.name, sep=' ')
        color.vector <- c(brewer.pal(length(scatter.plot.list[[scatter.name]]),'Set1'), "#000000")
        # make plots
        multi.occ.var.frame.list[[scatter.name]] <- meanFCScatterPlot( m.FC.list         = mean.FC.list, 
                                                                       m.norm.count.list = var.FC.list, 
                                                                       m.sample.type.A   = sample.type.A, 
                                                                       m.sample.type.B   = sample.type.B, 
                                                                       method.occurence  = length(scatter.plot.list[[scatter.name]]), 
                                                                       color.vector      = color.vector, 
                                                                       plot.title.name   = plot.title.name,
                                                                       round.to.xth      = 100,
                                                                       method.names      = scatter.plot.list[[scatter.name]],
                                                                       sample.thresh     = 5,
                                                                       plot.variance     = TRUE)
    }
}









###############################################################################
## Accuracy
###############################################################################
## Accuracy fold change
#boxPlot(accuracy.fc.pearson.list)
#boxPlot(accuracy.fc.spearman.list)
#boxPlot(sensitivity.fc.pearson.list)
#boxPlot(specificity.fc.pearson.list)
## Accuracy expression
#boxPlot(accuracy.pearson.list)
#boxPlot(accuracy.spearman.list)

###############################################################################
## Precision
###############################################################################
## Precision orig vs corrected Cufflinks and other (with transcripts)
#boxPlot(precision.orig.pearson.transcript.list)
#boxPlot(precision.pearson.transcript.list)
#boxPlot(precision.orig.spearman.transcript.list)
#boxPlot(precision.spearman.transcript.list)
## Precision orig vs corrected Cufflinks and other (without transcripts)
#boxPlot(precision.orig.pearson.list)
#boxPlot(precision.pearson.lis)
#boxPlot(precision.orig.spearman.list)
#boxPlot(precision.spearman.list)
## Fold change precision
#boxPlot(precision.fc.pearson.list)
#boxPlot(precision.fc.spearman.list)

###############################################################################
## Agreement
###############################################################################
## Agreement fold change
#boxPlot(agreement.fc.pearson.list)
#boxPlot(agreement.fc.spearman.list)



# check big fcCV list items -> only for 13 items, not dynamic
check.CVFC <- function(ps.count) {
    tmp.fcCV.list <- lapply(fcCV.list, function(x) return(x[['intra.gene.CV.log.FC']]))
    big.fcCV.list <- lapply(tmp.fcCV.list, function(tmp.vector) {
        tmp.vector <- na.omit(tmp.vector)
        tmp.vector <- unlist(lapply(tmp.vector, function(tmp.entry) {
            if (abs(tmp.entry) >= 5000) return(tmp.entry)}))
        return(tmp.vector)})
    big.method.names <- unname(unlist(lapply(names(big.fcCV.list), function(x) return(rep(x, length(big.fcCV.list[[x]]))))))
    big.gene.names <- unname(unlist(lapply(big.fcCV.list, function(x) return(names(x)))))
    id.list <- c()
    for (i in 1:length(big.method.names)) {
        id.list <- c(id.list,grep(big.gene.names[[i]],normalized.count.list[[big.method.names[[i]]]][[1]]))
    }
    problem.mat <- matrix(ncol=18,nrow=length(id.list))
    for (i in 1:length(big.method.names)) {
        problem.mat[i,] <- c(big.method.names[[i]],big.gene.names[[i]],unlist(normalized.count.list[[big.method.names[[i]]]][id.list[[i]],-1]))
    }
    colnames(problem.mat) <- c('Method', 'Gene Id', colnames(normalized.count.list[[1]])[-1])
    problem.frame <- data.frame(problem.mat)
    colnames(problem.frame)[1:length(colnames(problem.mat))] <- as.character(colnames(problem.mat))
    for (i in 1:min(length(samples.matched.A), length(samples.matched.B))) {
        # compute fold change vector for each gene and save into column i
        frame.A <- as.numeric(levels(problem.frame[[samples.matched.A[i]]]))[problem.frame[[samples.matched.A[i]]]]
        frame.B <- as.numeric(levels(problem.frame[[samples.matched.B[i]]]))[problem.frame[[samples.matched.B[i]]]]
        problem.frame[,18+i] <- (frame.A+ps.count) / (frame.B+ps.count)
        colnames(problem.frame)[length(colnames(problem.frame))] <- paste(samples.matched.A[i],samples.matched.B[i],sep='/')
    }
    colnames(problem.frame)[1:length(colnames(problem.mat))] <- as.character(colnames(problem.mat))
    problem.frame[,27] <- unlist(apply(problem.frame[,19:26], 1,function(x) return(sd(x))))
    colnames(problem.frame)[27] <- 'sd fc'
    problem.frame[,28] <- unlist(apply(log2(problem.frame[,19:26]), 1,function(x) return(sd(x))))
    colnames(problem.frame)[28] <- 'sd log fc'
    problem.frame[,29] <- unlist(apply(problem.frame[,19:26], 1,function(x) return(mean(x))))
    colnames(problem.frame)[29] <- 'mean fc'
    problem.frame[,30] <- unlist(apply(log2(problem.frame[,19:26]), 1,function(x) return(mean(x))))
    colnames(problem.frame)[30] <- 'mean log fc'
    problem.frame[,31] <- unlist(apply(problem.frame[,19:26], 1,function(x) return(sd(x)/mean(x))))
    colnames(problem.frame)[31] <- 'cv fc'
    problem.frame[,32] <- unlist(apply(log2(problem.frame[,19:26]), 1,function(x) return(sd(x)/mean(x))))
    colnames(problem.frame)[32] <- 'cv log fc'
    
    txt.file <- paste(paste('problem_file',as.character(ps.count),sep='_'), 'txt', sep='.')
    matrix.file <- file.path(output.dir, txt.file)
    write.table(problem.frame[,1:18], matrix.file, row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
    write.table(problem.frame[,19:26], matrix.file, row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE, append=TRUE)
    write.table(problem.frame[,27:32], matrix.file, row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE, append=TRUE)
    message('File \'', txt.file, '\' with data values created ...')
    
    return(problem.frame[,32])
}

#
if (FALSE) {
    pseudo.c <- c(1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10)
    method.gene <- paste(big.method.names,big.gene.names,sep='/')
    yvar <- c()
    for (i in pseudo.c) {
        yvar <- c(yvar, check.CVFC(i))
    }
    cond <- rep(method.gene,length(pseudo.c))
    xvar <- unlist(lapply(pseudo.c, function(x) return(rep(-log10(x),length(method.gene)))))
    scatter.frame <- data.frame(cond,xvar,yvar)
    p <- ggplot(data=scatter.frame, aes(x=xvar, y=yvar, color=cond)) + 
                    geom_point(shape=1) +
                    scale_colour_hue(l=50) + # Use a slightly darker palette than normal
                    geom_smooth(se=FALSE,    # Don't add shaded confidence region
                                fullrange=TRUE) + # Extend regression lines
                    xlab('-log10(pseudo count)') + ylab('CV(log(FC))')
    pdf.file <- paste('pseudo_count_vs_CV(log2(FC))', 'pdf', sep='.')
    ggsave(file.path(output.dir, pdf.file), width=4*length(pseudo.c)+4, height=length(method.gene)+1, scale=1.0, units='cm')
    message('File \'', pdf.file, '\' created ...')
}



# problem cuff tophat
#
#24634
if (FALSE) {
    for (i in 1:length(normalized.count.list[["EQPQM-ambig/HISAT"]][["Gene Id"]])) {
        test.x <- unlist(normalized.count.list[["CUFFL-default/TOPHAT"]][i ,-1])
        if ((max(test.x) > (sum(test.x)/10)) && (sum(test.x) > 1000)) {
            print(test.x)
        }
    }
    
    test.x.list <- list()
    for (i in 1:length(normalized.count.list[["EQPQM-ambig/HISAT"]][["Gene Id"]])) {
        test.x <- unlist(normalized.count.list[["CUFFL-default/TOPHAT"]][i ,-1])
        if ((max(test.x) > (sum(test.x)/10)) && (sum(test.x) > 1000)) {
            test.x.list[[length(test.x.list)+1]] <- test.x
        }
    }
}

