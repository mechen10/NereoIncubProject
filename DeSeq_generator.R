#!/bin/bash

# Deseq
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq")
# source("http://www.bioinformatics.nl/courses/RNAseq/DEseqExercise.R")
# internode_data = read.table(
#   "http://www.bioinformatics.nl/courses/RNAseq/maize_e3.table", row.names = 1,
#   header = T, sep = "\t" )

# Load my data
# This script is for plotting enrichment
library("optparse")
########################### OPT PARSE #################################
option_list = list(
  make_option(c("-i", "--OTUTableFP"), type="character",
              help="OTU Table input"),
  make_option(c("-m", "--metadata"), type="character",
              help="Metadata fp"),
  make_option(c("-c", "--category"), type="character",
              help="Category in metadata to collapse by (groups)"),
  make_option(c("-g", "--groups"), type="character",
              help="Groups that you want to compare in anaylsis. Comma separated. Note that there is NO DEFAULT"),
  make_option(c("-t", "--minthreshold"), type="character",
              help="minthreshold for filtering OTUTable"),
  make_option(c("-a", "--annotationsIncluded"), type="logical", action = "store_true",
              help="Include this flag if there are annotations in last column of data table"),
  make_option(c("-r", "--relativeAbund"), type="logical", action = "store_true",
              help="Include this flag if table is relative abundance, and what rarefaction depth it was"),
  make_option(c("-d", "--depth"), type="character", default = NULL,
              help="Include rarefaction depth IF the relativeAbund flag is TRUE")
);

make_option()
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

OTUTableFP = opt$OTUTableFP
metadataFP = opt$metadata
category = opt$category
groups = opt$groups
groups <- strsplit(groups, split = ",", fixed = TRUE)
minthreshold = opt$minthreshold
annotations = opt$annotationsIncluded
relativeAbund = opt$relativeAbund
depth = opt$depth

########################### LOAD DATA #################################
# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library("DESeq2")
setwd("/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis")
OTUTableFP <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/ANALYSIS_ALPHABETATAXA/summarize_taxa/rarefied_OTU_Table_sorted_L6.txt"
# OTUTableFP <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/OTU_MP_filt/OTU_Table_text.txt"
MFFP <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/ANALYSIS_ALPHABETATAXA/OTU_Tables_and_MP/MF_withalpha.txt"
minthreshold <- 10
category <- "ColRep"
annotations <- FALSE
relativeAbund <- TRUE
depth <- 1000
groups2 <- "NereotestH2OWater,NereotestExNWater,NereotestNereoWater,NereotestMastWater,NereotestNereoMastWater"
groups2 <- unlist(strsplit(groups2, split = ",", fixed = TRUE))
groups <- "NereotestExNExN,NereotestMastExN,NereotestNereoExN,NereotestNereoMastExN"
groups <- unlist(strsplit(groups, split = ",", fixed = TRUE))

MF <- read.delim(paste0(MFFP)
                 , header = TRUE
                 , row.names = 1
                 , stringsAsFactors = FALSE)
rownames(MF) <- gsub("-",".", rownames(MF))

OTUTable <- read.delim(paste0(OTUTableFP)
                      , header = TRUE
                      , skip = 1
                      , row.names = 1
                      , stringsAsFactors = FALSE)

############################ ******DESEQ***** ##########################
system("mkdir DESEQ")
system("echo 'LOG FOR DESEQ' > ./DESEQ/LOG.txt")

# # Make relative abundances
# colSumOTUTable <- colSums(OTUTable[1:(ncol(OTUTable)-1)])
# OTUTable.RelAbund <- OTUTable
# for (i in 1:nrow(OTUTable.RelAbund)) {
#   for (j in 1:(ncol(OTUTable.RelAbund)-1)) {
#     OTUTable.RelAbund[i,j] <- OTUTable.RelAbund[i,j]/colSumOTUTable[j]
#   }
# }

############ SETTING UP DATA ############
# DESEQ, first, remove OTUs that have less than 0.005 abundances in max sample
if (relativeAbund) {
  OTUTable.deseq <- OTUTable*depth
  OTUTable.deseq <- apply(OTUTable.deseq, c(1,2), as.integer)
} else {
  OTUTable.deseq <- OTUTable
}

if (annotations) {
  lastcolumn <- (ncol(OTUTable.deseq)-1)
} else {
  lastcolumn <- ncol(OTUTable.deseq)
}

# GET ONLY TREATMENTS I WANT
# Gets max of each row
mx = apply(OTUTable.deseq[,1:lastcolumn],1,FUN = sum)
# Gets rid of things that are less than 10
OTUTable.deseq = OTUTable.deseq[ mx > minthreshold,]

#### LOG info; record OTUs lost----------------------------------------------
system(paste0("echo 'Original dimensions of OTU table: "
              ,dim(OTUTable)[1], "x", dim(OTUTable)[2]
              ,"' >> ./DESEQ/LOG.txt"))
system(paste0("echo 'Filtering threshold: "
              , minthreshold
              ,"' >> ./DESEQ/LOG.txt"))
OTUsDELETED <- dim(OTUTable)-dim(OTUTable.deseq)[1]
system(paste0("echo 'Number of OTUs filtered out: "
              ,OTUsDELETED
              ,"' >> ./DESEQ/LOG.txt"))
system(paste0("echo 'New filtered OTU table dimensions: "
              ,dim(OTUTable.deseq)[1], "x", dim(OTUTable.deseq)[2]
              ,"' >> ./DESEQ/LOG.txt"))
####----------------------------------------------

# Now, create a new DESeq count dataset with replicates labelled

# Add 1 to every OTU
OTUTable.deseq <- OTUTable.deseq + 1

# For our data, we need to collapse by ColRep
categoryList = MF[unlist(lapply(colnames(OTUTable.deseq), function(x) grep(x, rownames(MF)))), paste0(category)]

system("mkdir ./DESEQ/RAW")
system("mkdir ./DESEQ/FILT")
for (groupOne in 1:(length(groups)-1)) {
  for (groupTwo in (groupOne+1):length(groups)) {
    # What am I comparing?
    currentCompare <- paste0(groups[groupOne],"vs",groups[groupTwo])
    # make colData for comparisonand filter OTU table by it
    colData <- cbind("condition" = categoryList[grep(paste0(groups[groupOne],"|",groups[groupTwo]), categoryList)])
    OTUTable.filt <- OTUTable.deseq[,grep(paste0(groups[groupOne],"|",groups[groupTwo]), categoryList)]
    # filter out small counts
    # Gets max of each row
    mx <- apply(OTUTable.filt,1,FUN = sum)
    # Gets rid of things that are less than 10
    OTUTable.filt <- OTUTable.filt[ mx > minthreshold,]
    
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(OTUTable.filt)
                                  , colData = colData
                                  , design = ~condition
    )
    
    dds <- DESeq(dds)
    res <- results(dds)
    res$taxonomy <- rownames(OTUTable.filt)
    res <- res[order(res$padj, na.last = TRUE),]
    
    write.table(res, col.names = NA, row.names = T
                , file = paste0("./DESEQ/RAW/",currentCompare, ".txt"), sep = "\t")

        # FILTER OUT MEANINGFUL ONES
    filt.res <- res[which(res$padj <= 0.05),]
    write.table(filt.res, col.names = NA, row.names = T
                , file = paste0("./DESEQ/FILT/",currentCompare, ".txt"), sep = "\t")
    
    
  }
}

for (groupOne in 1:(length(groups2)-1)) {
  for (groupTwo in (groupOne+1):length(groups2)) {
    # What am I comparing?
    currentCompare <- paste0(groups2[groupOne],"vs",groups2[groupTwo])
    # make colData for comparisonand filter OTU table by it
    colData <- cbind("condition" = categoryList[grep(paste0(groups2[groupOne],"|",groups2[groupTwo]), categoryList)])
    OTUTable.filt <- OTUTable.deseq[,grep(paste0(groups2[groupOne],"|",groups2[groupTwo]), categoryList)]
    # filter out small counts
    # Gets max of each row
    mx <- apply(OTUTable.filt,1,FUN = sum)
    # Gets rid of things that are less than 10
    OTUTable.filt <- OTUTable.filt[ mx > minthreshold,]
    
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(OTUTable.filt)
                                  , colData = colData
                                  , design = ~condition
    )
    
    dds <- DESeq(dds)
    res <- results(dds)
    res$taxonomy <- rownames(OTUTable.filt)
    res <- res[order(res$padj, na.last = TRUE),]
    
    write.table(res, col.names = NA, row.names = T
                , file = paste0("./DESEQ/RAW/",currentCompare, ".txt"), sep = "\t")
    
    # FILTER OUT MEANINGFUL ONES
    filt.res <- res[which(res$padj <= 0.05),]
    write.table(filt.res, col.names = NA, row.names = T
                , file = paste0("./DESEQ/FILT/",currentCompare, ".txt"), sep = "\t")
    
    
  }
}

