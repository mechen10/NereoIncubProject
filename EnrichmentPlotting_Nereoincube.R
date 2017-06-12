#!bin/bash

# This script is for plotting enrichment
library("optparse")
########################### OPT PARSE #################################
option_list = list(
  make_option(c("-f", "--first"), type="character",
              help="Folder with reduced OTU tables that are only core OTUs for each treatment type (water)"),
  make_option(c("-s", "--second"), type="character",
              help="Folder with reduced OTU tables that are only core OTUs for each treatment type (ExN)"),
  make_option(c("-t", "--otuTable"), type="character",
              help="OTU Table input"),
  make_option(c("-m", "--metadata"), type="character",
              help="Metadata fp")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

inputFolder = opt$first
inputFolder2 = opt$second
otuTableFP = opt$otuTableFP
metadataFP = opt$metadata

########################### FOR TESTING #################################
setwd('/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/')
inputFolder <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/CORE_ANALYSIS/WATER_coreanalysis/'
inputFolder2 <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/CORE_ANALYSIS/EXN_coreanalysis/'
otuTableFP <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/OTU_MP_filt/OTU_Table_text.txt'
metadataFP <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/OTU_MP_filt/MF_nochlpmito_m1000.txt'

########################### LOAD DATA #################################
library(gplots)
system("mkdir OTUHEATMAP")

OTUTable <- read.delim(paste0(otuTableFP)
                       , header = TRUE
                       , skip = 1
                       , row.names = 1)
taxonomyNames <- as.data.frame(OTUTable[,ncol(OTUTable)], row.names = rownames(OTUTable))
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
ListOfFiles2 <- system(paste0('ls ',inputFolder2), intern = TRUE)

Exp.Files <- ListOfFiles[-grep("taxaIDLegend", ListOfFiles)]
Exp.Files2 <- ListOfFiles2[-grep("taxaIDLegend", ListOfFiles2)]


for (i in Exp.Files) {
  assign(i, read.delim(paste0(inputFolder,"/",i)
                       , stringsAsFactors = FALSE
                       , header = TRUE
                       , row.names = 1
                       , strip.white = TRUE))
}
for (i in Exp.Files2) {
  assign(i, read.delim(paste0(inputFolder2,"/",i)
                       , stringsAsFactors = FALSE
                       , header = TRUE
                       , row.names = 1
                       , strip.white = TRUE))
}

############### START ##############
# Reorder ExpFiles into order I want
# Exp.Files <- Exp.Files["NereotestExNWater.txt","NereotestNereoWater.txt","NereotestMastWater.txt", "NereotestNereoMastWater.txt"]
Exp.Files <- Exp.Files[c(1,4,2,3)]

# Exp.Files2 <- Exp.Files2["NereotestNereoExN.txt", "NereotestMastExN.txt", "NereotestNereoMastExN.txt"]
Exp.Files2 <- Exp.Files2[c(2, 1, 3)]


# Get all unique OTUs
allOTUs <- c()
for (i in Exp.Files) {
  allOTUs <- c(allOTUs, rownames(get(i)))
}

uniqueOTUs <- unique(allOTUs)


allOTUs2 <- c()
for (i in Exp.Files2) {
  allOTUs2 <- c(allOTUs2, rownames(get(i)))
}

uniqueOTUs2 <- unique(allOTUs2)

# Get sample names in order

allSamples <- c()
for (i in Exp.Files) {
  allSamples <- c(allSamples, rownames(MF)[grep(gsub(".txt*", "", i), MF$ColRep)])
}

allSamples2 <- c()
for (i in Exp.Files2) {
  allSamples2 <- c(allSamples2, rownames(MF)[grep(gsub(".txt*", "", i), MF$ColRep)])
}

# Filter OTU table based on all core OTUs
allOTUsCombo <- unique(c(allOTUs, allOTUs2))
OTUTable.filt <- OTUTable.RelAbund[allOTUsCombo,c(allSamples,allSamples2)]

# Collapse OTU Table by sample 
MF.filt <- MF[sapply(colnames(OTUTable.filt), function(x) grep(x,rownames(MF))),]
collapsedTemp <- aggregate(t(OTUTable.filt), by = list(MF.filt$ColRep), FUN = mean )
OTUTable.filt.col <- t(data.frame(collapsedTemp, row.names = 1))
rownames(OTUTable.filt.col) <- gsub("X","",rownames(OTUTable.filt.col))

# Calculate collapse of Control (ExN and Water)
ExNWater.filt.col <- OTUTable.filt.col[,grep("NereotestExNWater",colnames(OTUTable.filt.col))]

ExNExN.filt <- OTUTable.RelAbund[allOTUsCombo,grep("ExN.Nereotest.ExN", colnames(OTUTable.RelAbund))]
ExNExN.filt.col <- rowMeans(ExNExN.filt)

# Calculate Fold Change
# First, get rid of the control in the OTU Table
OTUTable.filt.col.treatonly <- OTUTable.filt.col[,-grep("NereotestExNWater", colnames(OTUTable.filt.col))]

FoldChangeTable1 <- log(OTUTable.filt.col.treatonly[,grep("Water", colnames(OTUTable.filt.col.treatonly))]/ExNWater.filt.col,2)
FoldChangeTable2 <- log(OTUTable.filt.col.treatonly[,-grep("Water", colnames(OTUTable.filt.col.treatonly))]/ExNExN.filt.col,2)
FoldChangeTable <- cbind(FoldChangeTable1, FoldChangeTable2)

# Get rid of NA's and Inf's
for (i in 1:nrow(FoldChangeTable)) {
  for (j in 1:ncol(FoldChangeTable)) {
    if (is.na(FoldChangeTable[i,j])) {
      FoldChangeTable[i,j] <- -10
    } else if (FoldChangeTable[i,j] == Inf) {
      FoldChangeTable[i,j] <- 10
    } else if (FoldChangeTable[i,j] == -Inf){
      FoldChangeTable[i,j] <- -10
    } #else if (FoldChangeTable[i,j] >= 10) {
      # FoldChangeTable[i,j] <- 10
    #}
  }
}

# Reorder Samples
orderedSamples <- c("NereotestNereoWater"
                    , "NereotestNereoExN"
                    
                    , "NereotestMastWater"
                    , "NereotestMastExN"
                    
                    , "NereotestNereoMastWater"
                    , "NereotestNereoMastExN")


FoldChangeTable <- FoldChangeTable[,sapply(orderedSamples, function(x) grep(x,colnames(FoldChangeTable)))]

# Now, filter out fold changes that are less than 5
todelete <- c()
for (i in 1:nrow(FoldChangeTable)) {
  if (max(abs(FoldChangeTable[i,])) < 2) {
    todelete <- c(todelete, i)
  }
}

FoldChangeTable.filt <- FoldChangeTable[-todelete,]

# THIS IS FOLD CHANGE MAPPED
hc <- hclust(as.dist(1-cor(t(FoldChangeTable.filt))))
hc.v <- hclust(as.dist(1-cor(FoldChangeTable.filt)))
redgreenRamp <- colorRampPalette(c("orange","white","blue"))

# MAKE LABEL for ROWS
taxonomyLabels <- matrix(nrow = nrow(taxonomyNames), ncol = 1)
rownames(taxonomyLabels) <- rownames(taxonomyNames)
for (i in 1:nrow(taxonomyNames)){
  splitName <- strsplit(as.character(taxonomyNames[i,1]), split = "; __", fixed = TRUE)
  newName <- paste0(splitName[[1]][c(3)], ": ",splitName[[1]][c(5)], "_",splitName[[1]][c(6)])
  taxonomyLabels[i,1] <- newName
}

taxonomyLabels.filt <- taxonomyLabels[sapply(rownames(FoldChangeTable.filt), function(x) grep(paste0("^",x,"$"), rownames(taxonomyLabels))),]
taxonomyLabels.filt <- gsub("NA: NA_NA","Unidentified", taxonomyLabels.filt)

######### CORE OR NOT ############

OTUTable.Core <- OTUTable.RelAbund[allOTUsCombo,c(allSamples,allSamples2)]
for (i in 1:nrow(OTUTable.Core)) {
  for (j in 1:ncol(OTUTable.Core)) {
    for (k in Exp.Files) {
      if (colnames(OTUTable.Core)[j] %in% colnames(get(k))) {
        if (rownames(OTUTable.Core)[i] %in% rownames(get(k))) {
          OTUTable.Core[i,j] <- 1
        } else {
          OTUTable.Core[i,j] <- 0
        }
      }
      
    }
    for (l in Exp.Files2) {
      if (colnames(OTUTable.Core)[j] %in% colnames(get(l))) {
        if (rownames(OTUTable.Core)[i] %in% rownames(get(l))) {
          OTUTable.Core[i,j] <- 1
        } else {
          OTUTable.Core[i,j] <- 0
        }
      }
    }
  }
}

# Collapse OTU Table by sample
MF.core.filt <- MF[sapply(colnames(OTUTable.Core), function(x) grep(x,rownames(MF))),]
collapsedTemp <- aggregate(t(OTUTable.Core), by = list(MF.core.filt$ColRep), FUN = mean )
OTUTable.core.filt.col <- t(data.frame(collapsedTemp, row.names = 1))
rownames(OTUTable.core.filt.col) <- gsub("X","",rownames(OTUTable.core.filt.col))

# Reorder Samples and filter to match above


OTUTable.core.filt.col <- OTUTable.core.filt.col[sapply(paste0("^",rownames(FoldChangeTable.filt),"$"), function(x) grep(x, rownames(OTUTable.core.filt.col)))
                                                 ,sapply(paste0("^",colnames(FoldChangeTable.filt),"$"), function(x) grep(x, colnames(OTUTable.core.filt.col)))]

for (i in 1:nrow(OTUTable.core.filt.col)) {
  for (j in 1:ncol(OTUTable.core.filt.col)) {
    if (OTUTable.core.filt.col[i,j] == 1) {
      OTUTable.core.filt.col[i,j] <- "*"
    } else {
      OTUTable.core.filt.col[i,j] <- ""
      
    }
  }
}



###### PLOT #######
pdf("OTUHEATMAP/Heatmap.pdf", pointsize = 14, height = 14, width = 7)
heatmap.2(as.matrix(FoldChangeTable.filt)
        , Rowv = hc
        , labCol= c("WATER", "NMF","WATER","NMF","WATER","NMF")
        , labRow = taxonomyLabels.filt
        , Colv = NA
        , cexCol = 1
        , ColSideColors = c("green","green","red","red","brown","brown")
        , cexRow = 0.5
        , scale = "none"
        
        # , margins = c(7,10)
        
        , col = redgreenRamp(17)
        , breaks = c(-10,-9,-8,-7,-6,-5,-4,-3,-2,2,3,4,5,6,7,8,9,10)
        
        , cellnote = as.matrix(OTUTable.core.filt.col[hc$order,])
        , notecol = "black"
        , colsep = c(2,4)
        , dendrogram = "none"
        , trace = "none"
        , density.info = "none"
        
        , lmat = rbind(c(8, 5, 6), c(7 , 1, 3), c( 4, 2, 9))
        , lhei = c(2.5,1.5,10)
        , lwid = c(1,5,3)
        , key.xlab = "Fold-change"
        )
dev.off()

# quartz()
# heatmap.2(as.matrix(OTUTable.core.filt.col)
#         , Rowv = as.dendrogram(hc)
#         , labRow = NA
#         , Colv = NA
#         , labCol = c("WATER","NMF","WATER","NMF","WATER","NMF")
#         , scale = "none"
#         , col = c("white","darkgreen")
#         , trace = "none"
#         , density.info = "none"
#         , dendrogram = "none"
#         , ColSideColors = c("green","green","red","red","brown","brown")
#         )

