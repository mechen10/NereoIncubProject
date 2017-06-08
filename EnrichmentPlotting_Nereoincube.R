#!bin/bash

# This script is for plotting enrichment
library("optparse")
########################### OPT PARSE #################################
option_list = list(
  make_option(c("-i", "--inputFolder"), type="character",
              help="Folder with reduced OTU tables that are only core OTUs for each treatment type"),
  make_option(c("-t", "--otuTable"), type="character",
              help="OTU Table input"),
  make_option(c("-m", "--metadata"), type="character",
              help="Metadata fp")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

inputFolder = opt$inputFolder
otuTableFP = opt$otuTableFP
metadataFP = opt$metadata

########################### FOR TESTING #################################
setwd('/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/')
inputFolder <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/TEMP_frombotclust/WATER_BatchCoreAnalysis'
otuTableFP <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/TEMP_frombotclust/OTU_Table_textformat.txt'
metadataFP <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/TEMP_frombotclust/MF_nochlpmito_m1000.txt'

########################### LOAD DATA #################################
system("mkdir OTUHEATMAP")

OTUTable <- read.delim(paste0(otuTableFP)
                       , header = TRUE
                       , skip = 1
                       , row.names = 1)
taxonomyNames <- as.data.frame(OTUTable[,ncol(OTUTable)])
OTUTable <- OTUTable[,-ncol(OTUTable)]

# Make OTU Table relative abundance
OTUTable.RelAbund <- OTUTable
colSumsOTUTable <- colSums(OTUTable)
for (i in 1:ncol(OTUTable)) {
  for (j in 1:nrow(OTUTable)) {
    OTUTable.RelAbund[j,i] <- OTUTable[j,i]/colSumsOTUTable[i]
  }
}


MF <- read.delim(paste0(metadataFP)
                 , header = TRUE
                 , row.names = 1
                 , stringsAsFactors = FALSE)
rownames(MF) <- gsub("-",".", rownames(MF))

ListOfFiles <- system(paste0('ls ',inputFolder), intern = TRUE)

Exp.Files <- ListOfFiles[-grep("taxaIDLegend", ListOfFiles)]

for (i in Exp.Files) {
  assign(i, read.delim(paste0(inputFolder,"/",i)
                       , stringsAsFactors = FALSE
                       , header = TRUE
                       , row.names = 1
                       , strip.white = TRUE))
}
get(i)
############### START ##############
# Reorder ExpFiles into order I want
Exp.Files <- Exp.Files["NereotestExNWater.txt.txt","NereotestNereoWater.txt.txt","Nereotest"]


# Get all unique OTUs
allOTUs <- c()
for (i in Exp.Files) {
  allOTUs <- c(allOTUs, rownames(get(i)))
}

uniqueOTUs <- unique(allOTUs)

# Get sample names in order

allSamples <- c()
for (i in Exp.Files) {
  allSamples <- c(allSamples, rownames(MF)[grep(gsub(".txt*", "", i), MF$ColRep)])
}



