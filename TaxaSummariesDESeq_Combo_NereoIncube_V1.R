#!/bin/bash


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
  make_option(c("-o","--groups"), type="character",
              help="Groups that you want to compare in anaylsis. Comma separated. Note that there is NO DEFAULT"),
  make_option(c("-T","--groupstwo"), type="character",
              help="Groups that you want to compare in anaylsis. Comma separated. Note that there is NO DEFAULT"),
  make_option(c("-t", "--minthreshold"), type="character",
              help="minthreshold for filtering OTUTable"),
  make_option(c("-a", "--annotationsIncluded"), type="logical", action = "store_true", default = FALSE,
              help="Include this flag if there are annotations in last column of data table"),
  make_option(c("-r", "--relativeAbund"), type="logical", action = "store_true", default = FALSE,
              help="Include this flag if table is relative abundance, and what rarefaction depth it was"),
  make_option(c("-d", "--depth"), type="character", default = NULL,
              help="Include rarefaction depth IF the relativeAbund flag is TRUE")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

OTUTableFP = opt$OTUTableFP
MFFP = opt$metadata
category = opt$category
groups = opt$groups
groups <- unlist(strsplit(groups, split = ",", fixed = TRUE))
groups2 <- opt$groupstwo
groups2 <- unlist(strsplit(groups2, split = ",", fixed = TRUE))
minthreshold = opt$minthreshold
annotations = opt$annotationsIncluded
relativeAbund = opt$relativeAbund
depth = opt$depth


#################### set temp info ####################

#### For genus, used this:
# summarize_taxa.py -i /Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/OTU_MP_filt/OTU_Table_nochlpmito_m1000.biom -d ';' -L 6 --suppress_biom_table_output -o taxa_sum_lvl6 -a
### For OTU level, just used OTU table in text format
# biom convert -i OTUTable --to-tsv --header-key taxonomy -o OTU_Table_text.txt

# OTU
setwd("/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/OTULEVEL")
OTUTableFP <- '/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/OTU_MP_filt/OTU_Table_text.txt'
MFFP <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/ANALYSIS_ALPHABETATAXA/OTU_Tables_and_MP/MF_withalpha.txt"
minthreshold <- 10
category <- "ColRep"
annotations <- TRUE
relativeAbund <- FALSE
depth <- NULL
groups2 <- "NereotestH2OWater,NereotestExNWater,NereotestNereoWater,NereotestMastWater,NereotestNereoMastWater"
groups2 <- unlist(strsplit(groups2, split = ",", fixed = TRUE))
groups <- "NereotestExNExN,NereotestMastExN,NereotestNereoExN,NereotestNereoMastExN"
groups <- unlist(strsplit(groups, split = ",", fixed = TRUE))
delimiter <- "..__" # For OTU

# GENUS
setwd("/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/GENUSLEVEL")
OTUTableFP <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/ANALYSIS_ALPHABETATAXA/summarize_taxa/OTU_Table_nochlpmito_m1000_sorted_L6.txt"
MFFP <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/ANALYSIS_ALPHABETATAXA/OTU_Tables_and_MP/MF_withalpha.txt"
minthreshold <- 10
category <- "ColRep"
annotations <- FALSE
relativeAbund <- FALSE
depth <- NULL
groups2 <- "NereotestH2OWater,NereotestExNWater,NereotestNereoWater,NereotestMastWater,NereotestNereoMastWater"
groups2 <- unlist(strsplit(groups2, split = ",", fixed = TRUE))
groups <- "NereotestExNExN,NereotestMastExN,NereotestNereoExN,NereotestNereoMastExN"
groups <- unlist(strsplit(groups, split = ",", fixed = TRUE))
delimiter <- ";__" # For genus

########################### LOAD DATA #################################
library("DESeq2")
# library("RColorBrewer")
library("MASS")
library("gplots")
### THIS IS FOR THIS PROJECT:
ExN.toLoadPostSeq <- c("NereotestExNExNvsNereotestNereoExN"
                       , "NereotestExNExNvsNereotestMastExN"
                       , "NereotestExNExNvsNereotestNereoMastExN"
                       )
ExNWater.toLoadPostSeq <- c("NereotestExNWatervsNereotestNereoWater"
                            , "NereotestExNWatervsNereotestMastWater"
                            , "NereotestExNWatervsNereotestNereoMastWater")
######

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


########## START ##############

system("mkdir TAXASUMMARIES")
system("mkdir DESEQ")
system("echo 'LOG FOR DESEQ' > ./DESEQ/LOG.txt")


######### ********** DESeq First ************ ########

############ SETTING UP DATA ############
# Check if relabund or counts; if relAbund, then change to counts.
# If there's taxonomy at the end, then remove this.
if (annotations) {
  lastcolumn <- (ncol(OTUTable)-1)
  rownames(OTUTable)<- make.names(OTUTable[,ncol(OTUTable)], unique = TRUE)
  OTUTable <- OTUTable[,1:lastcolumn]
  
} else {
  lastcolumn <- ncol(OTUTable)
}

if (relativeAbund) {
  OTUTable.deseq <- OTUTable*depth
  OTUTable.deseq <- apply(OTUTable.deseq, c(1,2), as.integer)
} else {
  OTUTable.deseq <- OTUTable
}

## Sort from most to least abundant
OTUTable.deseq <- OTUTable.deseq[order(rowSums(OTUTable.deseq)),]

## Filter OTUs that are <100 counts
mxOTU <- unlist(lapply(1:nrow(OTUTable.deseq), function(x) max(OTUTable.deseq[x,])))
OTUTable.deseq <- OTUTable.deseq[-which(mxOTU < 100), ]


## Change taxasum into relative abundance
taxasum <- data.frame(sapply(1:ncol(OTUTable.deseq), function(x) {
  OTUTable.deseq[,x]/sum(OTUTable.deseq[,x])
}), row.names = rownames(OTUTable.deseq))

colnames(taxasum) <- colnames(OTUTable.deseq)

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
    
    dds <- DESeq(dds, test = "Wald")
    res <- results(dds, cooksCutoff = FALSE)
    res$taxonomy <- rownames(OTUTable.filt)
    
    res <- res[order(res$padj, na.last = TRUE),]
    
    write.table(res, col.names = NA, row.names = T
                , file = paste0("./DESEQ/RAW/",currentCompare, ".txt"), sep = "\t", quote = FALSE)
    
    # FILTER OUT MEANINGFUL ONES
    filt.res <- res[which(res$padj <= 0.05),]
    write.table(filt.res, col.names = NA, row.names = T
                , file = paste0("./DESEQ/FILT/",currentCompare, ".txt"), sep = "\t", quote = FALSE)
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
    
    dds <- DESeq(dds, test = "Wald")
    res <- results(dds, cooksCutoff = FALSE)
    res$taxonomy <- rownames(OTUTable.filt)
    res <- res[order(res$padj, na.last = TRUE),]
    
    write.table(res, col.names = NA, row.names = T
                , file = paste0("./DESEQ/RAW/",currentCompare, ".txt"), sep = "\t", quote = FALSE)
    
    # FILTER OUT MEANINGFUL ONES
    filt.res <- res[which(res$padj <= 0.05),]
    write.table(filt.res, col.names = NA, row.names = T
                , file = paste0("./DESEQ/FILT/",currentCompare, ".txt"), sep = "\t", quote = FALSE)
    
    
  }
}

########### Take output; combine ############
# Load outputfiles from DEseq
ExN.deseq <- c()
for (file in ExN.toLoadPostSeq) {
  newName <- gsub("Nereotest","",file)
  newName <- gsub("ExNExNvs","",newName)
  assign(paste0(newName), read.delim(paste0("./DESEQ/RAW/",file,".txt")
                                     , header = TRUE
                                     , row.names = 1
                                     , strip.white = TRUE
                                     , stringsAsFactors = FALSE))
  
  ExN.deseq <- c(ExN.deseq, paste0(newName))
}

Water.deseq <- c()
for (file in ExNWater.toLoadPostSeq) {
  newName <- gsub("Nereotest","",file)
  newName <- gsub("ExNWatervs","",newName)
  assign(paste0(newName), read.delim(paste0("./DESEQ/RAW/",file,".txt")
                                     , header = TRUE
                                     , row.names = 1
                                     , strip.white = TRUE
                                     , stringsAsFactors = FALSE))
  
  Water.deseq <- c(Water.deseq, paste0(newName))
}
# Now, combine them and filter P
temp <- merge(NereoExN[,c("log2FoldChange","stat","padj")], MastExN[,c("log2FoldChange","stat","padj")], by = 0)
colnames(temp) <- gsub(".x",".Nereo", colnames(temp))
colnames(temp) <- gsub(".y", ".Mast", colnames(temp))
temp <- data.frame(temp, row.names = 1)
temp <- merge(temp, NereoMastExN[,c("log2FoldChange","stat","padj")], by = 0)
temp <- data.frame(temp, row.names = 1)
colnames(temp)[7:9] <- c(paste0(colnames(temp)[7:9], ".NereoMast"))


temp2 <- merge(NereoWater[,c("log2FoldChange","stat","padj")], MastWater[,c("log2FoldChange","stat","padj")], by = 0)
colnames(temp2) <- gsub(".x",".Nereo", colnames(temp2))
colnames(temp2) <- gsub(".y", ".Mast", colnames(temp2))
temp2 <- data.frame(temp2, row.names = 1)
temp2 <- merge(temp2, NereoMastWater[,c("log2FoldChange","stat","padj")], by = 0)
temp2 <- data.frame(temp2, row.names = 1)
colnames(temp2)[7:9] <- c(paste0(colnames(temp2)[7:9], ".NereoMast"))

deseq.combo <- merge(temp,temp2, by = 0)
colnames(deseq.combo) <- gsub(".x",".ExN", colnames(deseq.combo))
colnames(deseq.combo) <- gsub(".y",".water", colnames(deseq.combo))
deseq.combo <- data.frame(deseq.combo, row.names = 1)

# Filter by padj
colpadj <- grep("padj", colnames(deseq.combo))
colstat <- grep("stat", colnames(deseq.combo))
colfc <- grep("log2FoldChange", colnames(deseq.combo))
deseq.sig <- c()
deseq.stat <- c()
deseq.fc <- c()
for (r in 1:nrow(deseq.combo)) {
  tempLine <- deseq.combo[r,colpadj]
  tempLine2 <- deseq.combo[r,colfc]
  if (any(tempLine < 0.05, na.rm = TRUE)) {
    if (any(abs(tempLine2) > 2, na.rm = TRUE)) {
      deseq.sig <- rbind(deseq.sig, deseq.combo[r,colpadj])
      deseq.stat <- rbind(deseq.stat, deseq.combo[r,colstat])
      deseq.fc <- rbind(deseq.fc, deseq.combo[r,colfc])
    }
  }
}

write.table(deseq.fc, row.names = TRUE, col.names = NA, sep = "\t", file = "DESEQ/deseq.fc.txt")
write.table(deseq.sig, row.names = TRUE, col.names = NA, sep = "\t", file = "DESEQ/deseq.sig.txt")

######## *******Taxa summaries****** ############
# Make sure MF and taxa summaries are in each other
MF <- MF[(rownames(MF) %in% colnames(taxasum)),]
taxasum <- taxasum[,(colnames(taxasum) %in% rownames(MF))]

# Reorder
MF <- MF[order(match(rownames(MF),colnames(taxasum))),]

########## SPLIT INTO EXPERIMENTS #########

taxasum.Loneincube <- taxasum[,grep("Loneincube", MF$Type)]
MF.Loneincube <-  MF[rownames(MF) %in% colnames(taxasum.Loneincube),]
MF.Loneincube <- MF.Loneincube[order(match(rownames(MF.Loneincube),colnames(taxasum.Loneincube))),]

taxasum.Nereotest <- taxasum[,grep("Nereotest", MF$Type)]
MF.Nereotest <-  MF[rownames(MF) %in% colnames(taxasum.Nereotest),]
MF.Nereotest <-  MF.Nereotest[order(match(rownames(MF.Nereotest),colnames(taxasum.Nereotest))),]

taxasum.Water <- taxasum.Nereotest[,grep("(W|w)ater", MF.Nereotest$Substrate)]
MF.Water <-  MF[rownames(MF) %in% colnames(taxasum.Water),]
MF.Water <-  MF.Water[order(match(rownames(MF.Water),colnames(taxasum.Water))),]

taxasum.ExN <- taxasum.Nereotest[,grep("ExN", MF.Nereotest$Substrate)]
MF.ExN <-  MF[rownames(MF) %in% colnames(taxasum.ExN),]
MF.ExN <-  MF.ExN[order(match(rownames(MF.ExN),colnames(taxasum.ExN))),]

taxasum.Environ <- taxasum[,grep("Brockton", MF$ColRep)]
MF.Environ <- MF[rownames(MF) %in% colnames(taxasum.Environ),]
MF.Environ <- MF.Environ[order(match(rownames(MF.Environ), colnames(taxasum.Environ)))]

taxasum.Star <- taxasum[,grep("starfish", MF$Project)]
MF.Star <- MF[rownames(MF) %in% colnames(taxasum.Star),]
MF.Star <- MF.Star[order(match(rownames(MF.Star),colnames(taxasum.Star)))]
########## COLOURS ##########

# Filter deseq based on whether it is abundant at 5% or more
toDelete <- c()
for (r in 1:nrow(deseq.fc)) {
  rowAbund <- taxasum[match(rownames(deseq.fc)[r], rownames(taxasum)),]
  if (max(rowAbund) < 0.03) {
    toDelete <- c(toDelete, r)
  }
}
if (length(toDelete) > 0) {
  deseq.fc <- deseq.fc[-toDelete,]
  deseq.sig <- deseq.sig[-toDelete,]
  deseq.stat <- deseq.stat[-toDelete,]
}



set.seed(5)

# Alternative color generator
colAll <- c()
for (i in seq(0.2, 0.7,length.out = 4)) {
  for (j in seq(0.3, 0.8, length.out = 4)) {
    for (k in seq(0.3, 1, length.out = 3)) {
      if (!((i == j) & (i == k))) {
        colAll <- c(colAll, rgb(i,j,k))
      }
    }
  }
}
for (i in seq(0.2, 0.7,length.out = 4)) {
  for (j in seq(0.3, 0.8, length.out = 4)) {
    for (k in seq(0.3, 1, length.out = 3)) {
      if (!((i == j) & (i == k))) {
        colAll <- c(colAll, rgb(j,i,k))
      }
    }
  }
}
for (i in seq(0.2, 0.7,length.out = 4)) {
  for (j in seq(0.3, 0.8, length.out = 4)) {
    for (k in seq(0.3, 1, length.out = 3)) {
      if (!((i == j) & (i == k))) {
        colAll <- c(colAll, rgb(j,k,i))
      }
    }
  }
}
for (i in seq(0.2, 0.7,length.out = 4)) {
  for (j in seq(0.3, 0.8, length.out = 4)) {
    for (k in seq(0.3, 1, length.out = 3)) {
      if (!((i == j) & (i == k))) {
        colAll <- c(colAll, rgb(k,j,i))
      }
    }
  }
}
for (i in c(0.25,0.5,0.75,1)) {
  colAll <- c(colAll, rgb(i,0,0))
  colAll <- c(colAll, rgb(0,i,0))
  colAll <- c(colAll, rgb(0,0,i))
  colAll <- c(colAll, rgb(i,i,0))
  colAll <- c(colAll, rgb(i,0,i))
  colAll <- c(colAll, rgb(0,i,i))
}


colAll <- unique(colAll)


randomColors <- c(sample(c(colAll),nrow(deseq.fc)-1),"black")

colorLegend <- cbind(rownames(deseq.fc), c(randomColors))



############ LONE INCUBE #############

# Sort in order by colrep
MF.Loneincube$ColRep <- factor(MF.Loneincube$ColRep, levels = c("LoneincubeNereoNereo"
                                        ,"LoneincubeMastMast"
                                        , "LoneincubeNereowater"
                                        ,"LoneincubeMastwater")
)

# Sort taxasum by ColRep
order <- c()
for (i in levels(MF.Loneincube$ColRep)) {
  allOfType <- rownames(MF.Loneincube)[grep(paste0(i), MF.Loneincube$ColRep)]
  for (j in allOfType) {
    order <- c(order, grep(paste0(j), colnames(taxasum.Loneincube)))
  }
}
taxasum.Loneincube <- taxasum.Loneincube[,order]

# Align colors with LoneIncube
colorLegend.Loneincube <- data.frame(rep("#FFFFFF", nrow(taxasum.Loneincube)), row.names = rownames(taxasum.Loneincube)
                                     ,stringsAsFactors = FALSE)
for (r in 1:nrow(colorLegend.Loneincube)) {
  nameTemp <- rownames(colorLegend.Loneincube)[r]
  if ((nameTemp %in% colorLegend[,1]) ) {#& (tempRow >= 2)) {
    colorLegend.Loneincube[r,1] <- as.character(colorLegend[which(nameTemp == colorLegend[,1]),2])
  } 
}

# For legend
toDelete <- c()
for (r in 1:nrow(colorLegend.Loneincube)) {
  if (colorLegend.Loneincube[r,1] == "#FFFFFF") {
    toDelete <- c(toDelete, r)
  } 
}
if (length(toDelete) > 0){
  colorLegend.Loneincube.LEGEND <- data.frame(colorLegend.Loneincube)[-toDelete,]
  colorLegend.Loneincube.LEGEND <- data.frame(colorLegend.Loneincube.LEGEND)
  rownames(colorLegend.Loneincube.LEGEND) <- rownames(colorLegend.Loneincube)[-toDelete]
}

listNamesTemp <- c()
for (r in 1:nrow(colorLegend.Loneincube.LEGEND)) {
  newName <- rownames(colorLegend.Loneincube.LEGEND)[r]
  if (length(grep("Unassigned",newName)) > 0) {
    newName <- "Unidentified"
  } else {
    newName <- strsplit(newName, split = paste0(delimiter), fixed = TRUE)
    newName <- paste0(newName[[1]][3],": ",newName[[1]][5],"_",newName[[1]][6])
  }
  listNamesTemp <- c(listNamesTemp, newName)
}
rownames(colorLegend.Loneincube.LEGEND) <- make.names(listNamesTemp , unique = TRUE)
colorsToPlot.Loneincube <- as.vector(colorLegend.Loneincube.LEGEND[,1])

pdf( "./TAXASUMMARIES/LoneIncube.pdf", pointsize = 14)
par(fig = c(0,0.65,0,1), mar = c(4,4,4,4))
barplot(as.matrix(taxasum.Loneincube)
        , col = colorLegend.Loneincube[,1]
        , space = c(0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0)
        , border = FALSE
        , xaxt = "n"
        , ylab = "Relative Abundance"
        )
axis(side = 1
     # , las = 2
     , at = c(2,7,12.5,18.5)
     , labels = c("Nereo","Mast","Water(N)", "Water(M)")
     , tick = FALSE
    , line = -1
    , cex.axis = 0.6
    , las = 2
    )
par(fig = c(0.5,1,0,1), mar = c(2.1,0,4.1,0), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n")
legend("top"
       , legend = rev(rownames(colorLegend.Loneincube.LEGEND))[1:10]
      , pch = 22
      , pt.bg = rev(colorsToPlot.Loneincube)[1:10]
      , col = rev(colorsToPlot.Loneincube)[1:10]
      , cex = 0.5
      , pt.cex = 1.2
      , y.intersp = 1.5)

dev.off()

############ WATER #############

# Sort in order by colrep
MF.Water$ColRep <- factor(MF.Water$ColRep, levels = c("NereotestH2OWater"
                                                      ,"NereotestExNWater"
                                                      , "NereotestNereoWater"
                                                      , "NereotestMastWater"
                                                      ,"NereotestNereoMastWater")
)

# Sort taxasum by ColRep
order <- c()
for (i in levels(MF.Water$ColRep)) {
  allOfType <- rownames(MF.Water)[grep(paste0(i), MF.Water$ColRep)]
  for (j in allOfType) {
    order <- c(order, grep(paste0(j), colnames(taxasum.Water)))
  }
}
taxasum.Water <- taxasum.Water[,order]

# Align colors with LoneIncube
colorLegend.Water <- data.frame(rep("#FFFFFF", nrow(taxasum.Water)), row.names = rownames(taxasum.Water)
                                     ,stringsAsFactors = FALSE)
for (r in 1:nrow(colorLegend.Water)) {
  nameTemp <- rownames(colorLegend.Water)[r]
  if ((nameTemp %in% colorLegend[,1]) ) { 
    colorLegend.Water[r,1] <- as.character(colorLegend[which(nameTemp == colorLegend[,1]),2])
  } 
}

# For legend
toDelete <- c()
for (r in 1:nrow(colorLegend.Water)) {
  if (colorLegend.Water[r,1] == "#FFFFFF") {
    toDelete <- c(toDelete, r)
  } 
}
if (length(toDelete) > 0){
  colorLegend.Water.LEGEND <- data.frame(colorLegend.Water)[-toDelete,]
  colorLegend.Water.LEGEND <- data.frame(colorLegend.Water.LEGEND)
  rownames(colorLegend.Water.LEGEND) <- rownames(colorLegend.Water)[-toDelete]
}

listNamesTemp <- c()
for (r in 1:nrow(colorLegend.Water.LEGEND)) {
  newName <- rownames(colorLegend.Water.LEGEND)[r]
  if (length(grep("Unassigned", newName)) > 0) {
    newName <- "Unidentified"
  } else {
    newName <- strsplit(newName, split = paste0(delimiter), fixed = TRUE)
    newName <- paste0(newName[[1]][3],": ",newName[[1]][5],"_",newName[[1]][6])
  }
  listNamesTemp <- c(listNamesTemp, newName)
}
rownames(colorLegend.Water.LEGEND)<- make.names(listNamesTemp, unique = TRUE)

colorsToPlot.Water <- as.vector(colorLegend.Water.LEGEND[,1])

pdf( "./TAXASUMMARIES/Water.pdf", pointsize = 14)
par(fig = c(0,0.65,0,1), mar = c(6,4,4,4))
barplot(as.matrix(taxasum.Water)
        , col = colorLegend.Water[,1]
        , las = 2
        , space = c(0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0)
        , border = FALSE
        , xaxt = "n"
        , ylab = "Relative Abundance"
)
axis(side = 1
     , las = 2
     , at = c(2,7,8,12.5,13.5,18.5,19.5, 24,25)
     , labels = c("Water Alone","Water","with NMF only","Water", "with Nereo", "Water", "with Mast", "Water","with Nereo + Mast")
     , tick = FALSE
     , line = -1
     , cex.axis = 0.6
)
par(fig = c(0.5,1,0,1), mar = c(2.1,0,4.1,0), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n")
legend("top"
       , legend = rev(rownames(colorLegend.Water.LEGEND))[1:10]
       , pch = 22
       , pt.bg = rev(colorsToPlot.Water)[1:10]
       , col = rev(colorsToPlot.Water)[1:10]
       , cex = 0.5
       , pt.cex = 1.2)
dev.off()


############ EXN #############

# Sort in order by colrep
MF.ExN$ColRep <- factor(MF.ExN$ColRep, levels = c("NereotestExNExN"
                                                      ,"NereotestNereoExN"
                                                      , "NereotestMastExN"
                                                      ,"NereotestNereoMastExN")
)

# Sort taxasum by ColRep
order <- c()
for (i in levels(MF.ExN$ColRep)) {
  allOfType <- rownames(MF.ExN)[grep(paste0(i), MF.ExN$ColRep)]
  for (j in allOfType) {
    order <- c(order, grep(paste0(j), colnames(taxasum.ExN)))
  }
}
taxasum.ExN <- taxasum.ExN[,order]

# Align colors with ExN
colorLegend.ExN <- data.frame(rep("#FFFFFF", nrow(taxasum.ExN)), row.names = rownames(taxasum.ExN)
                                     ,stringsAsFactors = FALSE)
for (r in 1:nrow(colorLegend.ExN)) {
  nameTemp <- rownames(colorLegend.ExN)[r]
  if ((nameTemp %in% colorLegend[,1])){ 
    colorLegend.ExN[r,1] <- as.character(colorLegend[which(nameTemp == colorLegend[,1]),2])
  } 
}

# For legend
toDelete <- c()
for (r in 1:nrow(colorLegend.ExN)) {
  if (colorLegend.ExN[r,1] == "#FFFFFF") {
    toDelete <- c(toDelete, r)
  }
}
if (length(toDelete) > 0){
  colorLegend.ExN.LEGEND <- data.frame(colorLegend.ExN)[-toDelete,]
  colorLegend.ExN.LEGEND <- data.frame(colorLegend.ExN.LEGEND)
  rownames(colorLegend.ExN.LEGEND) <- rownames(colorLegend.ExN)[-toDelete]
}

listNamesTemp <- c()
for (r in 1:nrow(colorLegend.ExN.LEGEND)) {
  newName <- rownames(colorLegend.ExN.LEGEND)[r]
  if (length(grep("Unassigned", newName)) > 0) {
    newName <- "Unidentified"
  } else {
    newName <- strsplit(newName, split = paste0(delimiter), fixed = TRUE)
    newName <- paste0(newName[[1]][3],": ",newName[[1]][5],"_",newName[[1]][6])
  }
  listNamesTemp <- c(listNamesTemp, newName)
}
rownames(colorLegend.ExN.LEGEND) <- make.names(listNamesTemp, unique = TRUE)

colorsToPlot.ExN <- as.vector(colorLegend.ExN.LEGEND[,1])

pdf( "./TAXASUMMARIES/ExN.pdf", pointsize = 14)
par(fig = c(0,0.65,0,1), mar = c(6,4,4,4))
barplot(as.matrix(taxasum.ExN)
        , col = colorLegend.ExN[,1]
        , las = 2
        , space = c(0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0)
        , border = FALSE
        , xaxt = "n"
        , ylab = "Relative Abundance"
)
axis(side = 1
     , las = 2
     , at = c(2,8,9,13,14,18,19)
     , labels = c("NMF Alone","NMF", "with Nereo", "NMF", "with Mast", "NMF","with Nereo + Mast")
     , tick = FALSE
     , line = -1
     , cex.axis = 0.6
)
par(fig = c(0.5,1,0,1), mar = c(2.1,0,4.1,0), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n")
legend("top"
       , legend = rev(rownames(colorLegend.ExN.LEGEND))[1:10]
       , pch = 22
       , pt.bg = rev(colorsToPlot.ExN)[1:10]
       , col = rev(colorsToPlot.ExN)[1:10]
       , cex = 0.5
       , pt.cex = 1.2)

dev.off()


############ ENVIRON #############

# Sort in order by colrep
MF.Environ$ColRep <- factor(MF.Environ$ColRep, levels = c("EnvironmentalBrocktonYoungNereo"
                                                  ,"EnvironmentalBrocktonOldNereo"
                                                  , "EnvironmentalBrocktonMast"
                                                  )
)

# Sort taxasum by ColRep
order <- c()
for (i in levels(MF.Environ$ColRep)) {
  allOfType <- rownames(MF.Environ)[grep(paste0(i), MF.Environ$ColRep)]
  for (j in allOfType) {
    order <- c(order, grep(paste0(j), colnames(taxasum.Environ)))
  }
}
taxasum.Environ <- taxasum.Environ[,order]

# Align colors with LoneIncube
colorLegend.Environ <- data.frame(rep("#FFFFFF", nrow(taxasum.Environ)), row.names = rownames(taxasum.Environ)
                              ,stringsAsFactors = FALSE)
for (r in 1:nrow(colorLegend.Environ)) {
  nameTemp <- rownames(colorLegend.Environ)[r]
  tempRow <- sum(taxasum.Environ[r,]>0.03)
  if ((nameTemp %in% colorLegend[,1])){
    colorLegend.Environ[r,1] <- as.character(colorLegend[which(nameTemp == colorLegend[,1]),2])
  } 
}

# For legend
toDelete <- c()
for (r in 1:nrow(colorLegend.Environ)) {
  if (colorLegend.Environ[r,1] == "#FFFFFF") {
    toDelete <- c(toDelete, r)
  } 
}
if (length(toDelete) > 0){
  colorLegend.Environ.LEGEND <- data.frame(colorLegend.Environ)[-toDelete,]
  colorLegend.Environ.LEGEND <- data.frame(colorLegend.Environ.LEGEND)
  rownames(colorLegend.Environ.LEGEND) <- rownames(colorLegend.Environ)[-toDelete]
}

listNamesTemp <- c()
for (r in 1:nrow(colorLegend.Environ.LEGEND)) {
  newName <- rownames(colorLegend.Environ.LEGEND)[r]
  if (length(grep("Unassigned", newName)) > 0) {
    newName <- "Unidentified"
  } else {
    newName <- strsplit(newName, split = paste0(delimiter), fixed = TRUE)
    newName <- paste0(newName[[1]][3],": ",newName[[1]][5],"_",newName[[1]][6])
  }
  listNamesTemp <- c(listNamesTemp, newName)
}
rownames(colorLegend.Environ.LEGEND) <- make.names(listNamesTemp, unique = TRUE)

colorsToPlot.Environ <- as.vector(colorLegend.Environ.LEGEND[,1])


pdf( "./TAXASUMMARIES/Environ.pdf", pointsize = 14)
par(fig = c(0,0.65,0,1), mar = c(6,4,4,4))
barplot(as.matrix(taxasum.Environ)
        , col = colorLegend.Environ[,1]
        , las = 2
        , space = c(0,0,0,0,0,1,0,0,0,0,1,0,0,0,0)
        , border = FALSE
        , xaxt = "n"
        , ylab = "Relative Abundance"
)
axis(side = 1
     , las = 2
     , at = c(2.5,8.5,14)
     , labels = c("10cm Nereocystis","50cm Nereocystis", "Mastocarpus")
     , tick = FALSE
     , line = -1
     , cex.axis = 0.6
)
par(fig = c(0.5,1,0,1), mar = c(2.1,0,4.1,0), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n")
legend("top"
       , legend = rev(rownames(colorLegend.Environ.LEGEND))[1:10]
       , pch = 22
       , pt.bg = rev(colorsToPlot.Environ)[1:10]
       , col = rev(colorsToPlot.Environ)[1:10]
       , cex = 0.5
       , pt.cex = 1.2)
dev.off()

############ Star #############

# Sort in order by colrep
MF.Star$ColRep <- factor(MF.Star$ColRep, levels = c("EnvironmentalStarfishInnerNereocystis"
                                                          ,"EnvironmentalStarfishOuterNereocystis"
                                                          , "EnvironmentalStarfishInnerWater"
                                                    ,"EnvironmentalStarfishOuterWater"
)
)

# Sort taxasum by ColRep
order <- c()
for (i in levels(MF.Star$ColRep)) {
    allOfType <- rownames(MF.Star)[grep(paste0(i), MF.Star$ColRep)]
    for (j in allOfType) {
        order <- c(order, grep(paste0(j), colnames(taxasum.Star)))
    }
}
taxasum.Star <- taxasum.Star[,order]

# Align colors with LoneIncube
colorLegend.Star <- data.frame(rep("#FFFFFF", nrow(taxasum.Star)), row.names = rownames(taxasum.Star)
                                  ,stringsAsFactors = FALSE)
for (r in 1:nrow(colorLegend.Star)) {
    nameTemp <- rownames(colorLegend.Star)[r]
    tempRow <- sum(taxasum.Star[r,]>0.03)
    if ((nameTemp %in% colorLegend[,1])){
        colorLegend.Star[r,1] <- as.character(colorLegend[which(nameTemp == colorLegend[,1]),2])
    } 
}

# For legend
toDelete <- c()
for (r in 1:nrow(colorLegend.Star)) {
    if (colorLegend.Star[r,1] == "#FFFFFF") {
        toDelete <- c(toDelete, r)
    } 
}
if (length(toDelete) > 0){
    colorLegend.Star.LEGEND <- data.frame(colorLegend.Star)[-toDelete,]
    colorLegend.Star.LEGEND <- data.frame(colorLegend.Star.LEGEND)
    rownames(colorLegend.Star.LEGEND) <- rownames(colorLegend.Star)[-toDelete]
}

listNamesTemp <- c()
for (r in 1:nrow(colorLegend.Star.LEGEND)) {
    newName <- rownames(colorLegend.Star.LEGEND)[r]
    if (length(grep("Unassigned", newName)) > 0) {
        newName <- "Unidentified"
    } else {
        newName <- strsplit(newName, split = paste0(delimiter), fixed = TRUE)
        newName <- paste0(newName[[1]][3],": ",newName[[1]][5],"_",newName[[1]][6])
    }
    listNamesTemp <- c(listNamesTemp, newName)
}
rownames(colorLegend.Star.LEGEND) <- make.names(listNamesTemp, unique = TRUE)

colorsToPlot.Star <- as.vector(colorLegend.Star.LEGEND[,1])

pdf( "./TAXASUMMARIES/Star.pdf", pointsize = 14)
par(fig = c(0,0.65,0,1), mar = c(6,4,4,4))
barplot(as.matrix(taxasum.Star)
        , col = colorLegend.Star[,1]
        , las = 2
        , space = c(0,0,0,0,1,0,0,0,1,0,1,0)
        , border = FALSE
        , xaxt = "n"
        , ylab = "Relative Abundance"
)
axis(side = 1
     , las = 2
     , at = c(2.5,8.5,14)
     , labels = c("Inner Nereocystis","Outer Nereocystis", "Inner Water", "Outer Water")
     , tick = FALSE
     , line = -1
     , cex.axis = 0.6
)
par(fig = c(0.5,1,0,1), mar = c(2.1,0,4.1,0), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n")
legend("top"
       , legend = rev(rownames(colorLegend.Star.LEGEND))[1:10]
       , pch = 22
       , pt.bg = rev(colorsToPlot.Star)[1:10]
       , col = rev(colorsToPlot.Star)[1:10]
       , cex = 0.5
       , pt.cex = 1.2)
dev.off()

################ ALL PLOTS TOGETHER #################
# Want to put NMF-incubation experiment together and seaweed surfaces together
# Make new legend that combines everything
colorLegend.ExN.new <- cbind(rev(rownames(colorLegend.ExN.LEGEND))[1:20], rev(as.character(colorLegend.ExN.LEGEND[,1]))[1:20])
colorLegend.Water.new <- cbind(rev(rownames(colorLegend.Water.LEGEND))[1:20], rev(as.character(colorLegend.Water.LEGEND[,1]))[1:20])

duplicatedRows1 <- which(colorLegend.ExN.new[,1] %in% colorLegend.Water.new[,1])
comboLegend <- rbind(colorLegend.Water.new, colorLegend.ExN.new[-duplicatedRows1,])

pdf("TAXASUMMARIES/TaxaSummaries_combo.pdf", pointsize = 14, width = 10, height = 7)
par(oma = c(0,0,0,0))
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , bty = "n"
     , xlab = ""
     , ylab = "")
title(ylab = "Relative Abundance"
      , line = 3)
par(oma = c(0,2,0,0))
par(fig = c(0,0.6,0.4,1), mar = c(5.1,4.1,1,2.1), new = TRUE)
barplot(as.matrix(taxasum.Water)
        , col = colorLegend.Water[,1]
        , las = 2
        , space = c(0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0)
        , border = FALSE
        , xaxt = "n"
        , ylab = "WATER SAMPLES"
)
axis(side = 1
     , las = 1
     , at = c(2,7.5,13.5,19.5,25.5)
     , labels = c("WATER ONLY","NMF ALONE","WITH NEREO","WITH MAST","WITH N + M")
     , tick = FALSE
     , line = -1
     , cex.axis = 0.6
)
par(fig = c(0,0.6,0,0.6), mar = c(1,4.1,5.1,2.1), new = TRUE)
#Need to add empty plot to make spacing right
colnames(taxasum.ExN)
barplot(as.matrix(cbind(rep(0, nrow(taxasum.ExN)),rep(0, nrow(taxasum.ExN)),rep(0, nrow(taxasum.ExN)),rep(0, nrow(taxasum.ExN)),taxasum.ExN))
        , col = colorLegend.ExN[,1]
        , las = 2
        , space = c(0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,3,0,0,0,0)
        , border = FALSE
        , xaxt = "n"
        , ylab = "NMF SURFACE"
        
)

par(fig = c(0.55,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n")
legend(x = -1.1, y = 0
       , legend = (comboLegend[,1])
       , pch = 22
       , pt.bg = (comboLegend[,2])
       , col = (comboLegend[,2])
       , cex = 0.65
       , pt.cex = 1.2
       , bty = "n"
       , xjust = 0
       , yjust = 0.5
       )
par(fig = c(0,0.6,0.47,0.52), mar = c(0,4.1,0,2.1), new = TRUE)
barplot(cbind(rep(1, 24))
        , beside = TRUE
        , xlab = ""
        , xaxt = "n"
        , yaxt = "n"
        , col = c(rep("blue",4),rep("grey",5),rep("green",5),rep("purple",5),rep("brown",5))
        , space = c(0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0)
        , border = c(rep("blue",4),rep("grey",5),rep("green",5),rep("purple",5),rep("brown",5))
)
dev.off()

write.matrix(comboLegend, sep = "\t", file = "TAXASUMMARIES/comboLegend.txt")

######### ALL NEREO PLOT ##############
# Want to put NMF-incubation experiment together and seaweed surfaces together
# Make new legend that combines everything
taxasum.Environ.Nereoonly <- taxasum.Environ[,grep("Nereo", colnames(taxasum.Environ))]
taxasum.ExN.Nereoonly <- taxasum.ExN[,grep("Nereotest.ExN",colnames(taxasum.ExN))]
taxasum.Loneincube.Nereoonly <- taxasum.Loneincube[,grep("Nereo.Loneincube.Nereo", colnames(taxasum.Loneincube))]
taxasum.Star.Nereoonly <- taxasum.Star[,grep("Nereocystis",colnames(taxasum.Star))]

colorLegend.ExN.new <- cbind(rev(rownames(colorLegend.ExN.LEGEND))[1:20], rev(as.character(colorLegend.ExN.LEGEND[,1]))[1:20])
colorLegend.Loneincube.new <- cbind(rev(rownames(colorLegend.Loneincube.LEGEND))[1:20], rev(as.character(colorLegend.Loneincube.LEGEND[,1]))[1:20])
colorLegend.Star.new <- cbind(rev(rownames(colorLegend.Star.LEGEND))[1:20], rev(as.character(colorLegend.Star.LEGEND[,1]))[1:20])
colorLegend.Environ.new <- cbind(rev(rownames(colorLegend.Environ.LEGEND))[1:20], rev(as.character(colorLegend.Environ.LEGEND[,1]))[1:20])
comboLegend.Nereo <- rbind(colorLegend.ExN.new,colorLegend.Loneincube.new,colorLegend.Star.new,colorLegend.Environ.new)
duplicatedRows.nereo <- which(duplicated(comboLegend.Nereo[,1]))
comboLegend.Nereo <- comboLegend.Nereo[-duplicatedRows.nereo,]

pdf("TAXASUMMARIES/Nereo_combo.pdf", pointsize = 14, width = 10, height = 7)
par(oma = c(0,0,0,0))
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , bty = "n"
     , xlab = ""
     , ylab = "")
title(ylab = "Relative Abundance"
      , line = 3)
par(oma = c(0,2,0,0))
par(fig = c(0,0.6,0.4,1), mar = c(5.1,4.1,1,2.1), new = TRUE)
barplot(as.matrix(taxasum.Environ.Nereoonly)
        , col = colorLegend.Environ[,1]
        , las = 2
        , space = c(0,0,0,0,0,1,0,0,0,0)
        , border = FALSE
        , xaxt = "n"
        , ylab = "ENVIRONMENTAL SAMPLES"
)
axis(side = 1
     , las = 1
     , at = c(2.5,8.5)
     , labels = c("MERISTEM","MATURE BLADES")
     , tick = FALSE
     , line = -1
     , cex.axis = 0.6
)
par(fig = c(0,0.6,0,0.6), mar = c(1,4.1,5.1,2.1), new = TRUE)
#Need to add empty plot to make spacing right
barplot(as.matrix(cbind(taxasum.ExN.Nereoonly,taxasum.Loneincube.Nereoonly,rep(0,nrow(taxasum.Loneincube.Nereoonly))))
        , col = colorLegend.ExN[,1]
        , las = 2
        , space = c(0,0,0,0,0,1,0,0,0,0)
        , border = FALSE
        , xaxt = "n"
        , ylab = "LAB SAMPLES"
        
)

par(fig = c(0.55,1,0,0.5), mar = c(0,0,0,0), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n")
legend(x = -1.1, y = 0
       , legend = (comboLegend.Nereo[,1])
       , pch = 22
       , pt.bg = (comboLegend.Nereo[,2])
       , col = (comboLegend.Nereo[,2])
       , cex = 0.65
       , pt.cex = 1.2
       , bty = "n"
       , xjust = 0
       , yjust = 0.5
)
par(fig = c(0,0.6,0.47,0.52), mar = c(0,4.1,0,2.1), new = TRUE)
barplot(cbind(rep(1, 9))
        , beside = TRUE
        , xlab = ""
        , xaxt = "n"
        , yaxt = "n"
        , col = c(rep("green",5),rep("darkgreen",5))
        , space = c(0,0,0,0,0,1,0,0,0,0)
        , border = c(rep("green",5),rep("darkgreen",5))
)
colnames(taxasum.Star.Nereoonly)
par(fig = c(0.55,1,0.4,1), mar = c(5.1,4.1,1,2.1), new =TRUE)
barplot(as.matrix(taxasum.Star.Nereoonly)
        , col = colorLegend.Star[,1]
        , las = 2
        , space = c(0,0,0,0,1,0,0,0)
        , border = FALSE
        , xaxt = "n"
        , ylab = "REMOTE SAMPLES"
)
axis(side = 1
     , las = 1
     , at = c(2,8)
     , labels = c("INNER FOREST","OUTER FOREST")
     , tick = FALSE
     , line = -1
     , cex.axis = 0.6
)
dev.off()




############## *********PLOT DESEQ********* ################
# are fc, sig, and stat still the same?
all(rownames(deseq.fc) == rownames(deseq.sig))
all(rownames(deseq.fc) == rownames(deseq.stat))

# Change names of both sig and fc
deseq.fc.names <- deseq.fc
deseq.sig.names <- deseq.sig
deseq.stat.names <- deseq.stat
listNamesTemp <- c()
for (r in 1:nrow(deseq.fc.names)) {
  tempName <- rownames(deseq.fc.names)[r]
  tempName <- strsplit(tempName, split = paste0(delimiter), fixed = TRUE)
  newName <- paste0(tempName[[1]][3], ": ", tempName[[1]][5],"_",tempName[[1]][6])
  if (newName == "NA: NA_NA"){
    newName <- "Unidentified"
  }
  listNamesTemp <- c(listNamesTemp, newName)
}
rownames(deseq.fc.names) <- make.names(listNamesTemp, unique = TRUE)
rownames(deseq.sig.names) <- make.names(listNamesTemp, unique = TRUE)
rownames(deseq.stat.names) <- make.names(listNamesTemp, unique = TRUE)

# rownames(deseq.fc.names)

# Now filter deseq by names of comboExNWater

deseq.fc.filter <- deseq.fc.names[rownames(deseq.fc.names) %in% comboLegend[,1],]
deseq.sig.filter <- deseq.sig.names[rownames(deseq.sig.names) %in% comboLegend[,1],]
deseq.stat.filter <- deseq.stat.names[rownames(deseq.stat.names) %in% comboLegend[,1],]

# GET RID OF UNIDENTIFIED FIRST
if ("Unidentified" %in% rownames(deseq.fc.filter)) {
  deseq.fc.filter <- deseq.fc.filter[-grep("Unidentified", rownames(deseq.fc.filter)),]
  deseq.sig.filter <- deseq.sig.filter[-grep("Unidentified", rownames(deseq.sig.filter)),]
  deseq.stat.filter <- deseq.stat.filter[-grep("Unidentified", rownames(deseq.stat.filter)),]
  
}

# Now, get color
absMax <- ceiling(max(abs(deseq.fc.filter)))
ncolors <- absMax*2+1

colorRange <- colorRampPalette(c("blue","grey","red"))

# Order samples properly, just in case-- should already by done though
colnames(deseq.fc.filter) <- gsub("log2FoldChange.","",colnames(deseq.fc.filter))
colnames(deseq.sig.filter) <- gsub("padj.","",colnames(deseq.sig.filter))
colnames(deseq.stat.filter) <- gsub("stat.","",colnames(deseq.stat.filter))

deseq.fc.filter.ExN <- deseq.fc.filter[,unlist(lapply(c("Nereo.ExN","Mast.ExN","NereoMast.ExN"), function(x) grep(paste0("^",x,"$"), colnames(deseq.fc.filter))))]
deseq.fc.filter.Water <- deseq.fc.filter[,unlist(lapply(c("Nereo.water","Mast.water","NereoMast.water"), function(x) grep(paste0("^",x,"$"), colnames(deseq.fc.filter))))]

deseq.sig.filter.ExN <- deseq.sig.filter[,unlist(lapply(c("Nereo.ExN","Mast.ExN","NereoMast.ExN"), function(x) grep(paste0("^",x,"$"), colnames(deseq.sig.filter))))]
deseq.sig.filter.Water <- deseq.sig.filter[,unlist(lapply(c("Nereo.water","Mast.water","NereoMast.water"), function(x) grep(paste0("^",x,"$"), colnames(deseq.sig.filter))))]

deseq.stat.filter.ExN <- deseq.stat.filter[,unlist(lapply(c("Nereo.ExN","Mast.ExN","NereoMast.ExN"), function(x) grep(paste0("^",x,"$"), colnames(deseq.stat.filter))))]
deseq.stat.filter.Water <- deseq.stat.filter[,unlist(lapply(c("Nereo.water","Mast.water","NereoMast.water"), function(x) grep(paste0("^",x,"$"), colnames(deseq.stat.filter))))]

deseq.fc.Water.color <- deseq.fc.filter.Water
deseq.sig.Water.color <- deseq.sig.filter.Water
deseq.stat.Water.color <- deseq.stat.filter.Water


# Sort OTUs manually!
allenriched.Water <- c()
allreduced.Water <- c()
rest.Water <- c()
for (r in 1:nrow(deseq.fc.Water.color)) {
  if (all(deseq.fc.Water.color[r,] > 0)) {
    allenriched.Water <- c(allenriched.Water, r)
  } else if (all(deseq.fc.Water.color[r,] < 0)) {
    allreduced.Water <- c(allreduced.Water,r)
  } else {
    rest.Water <- rbind(rest.Water, deseq.fc.Water.color[r,])
  }
}
# Combine
rest.Water <- rest.Water[with(as.data.frame(rest.Water), order(-rest.Water[,"Nereo.water"]
                                                               ,-rest.Water[,"Mast.water"]
                                                               ,-rest.Water[,"NereoMast.water"])), ]
enrichedTemp <- deseq.fc.Water.color[allenriched.Water,]
enrichSums <- order(rowSums(enrichedTemp), decreasing = TRUE)
enrichedTemp <- enrichedTemp[enrichSums,]

reducedTemp <- deseq.fc.Water.color[allreduced.Water,]
reduceSums <- order(rowSums(reducedTemp), decreasing = TRUE)
reducedTemp <- reducedTemp[reduceSums,]

restSums <- order(rowSums(rest.Water), decreasing = TRUE)
rest.Water <- rest.Water[restSums,]

deseq.fc.Water.color2 <- rbind(enrichedTemp, rest.Water, reducedTemp)

# Now, reorder deseq.fc.ExN
deseq.fc.filter.ExN <- deseq.fc.filter.ExN[match(rownames(deseq.fc.Water.color2), rownames(deseq.fc.filter.ExN)),]

# Reorder other deseq files in accordance with deseq.fc
deseq.sig.filter.Water <- deseq.sig.filter.Water[match(rownames(deseq.fc.Water.color2), rownames(deseq.sig.filter.Water)),]
deseq.sig.filter.ExN <- deseq.sig.filter.ExN[match(rownames(deseq.fc.Water.color2), rownames(deseq.sig.filter.ExN)),]

deseq.stat.filter.Water <- deseq.stat.filter.Water[match(rownames(deseq.fc.Water.color2), rownames(deseq.stat.filter.Water)),]
deseq.stat.filter.ExN <- deseq.stat.filter.ExN[match(rownames(deseq.fc.Water.color2), rownames(deseq.stat.filter.ExN)),]


# Change ones that are not part of the group into NaNs
ExN.tokeep <- rownames(colorLegend.ExN.LEGEND)
Water.tokeep <- rownames(colorLegend.Water.LEGEND)

deseq.fc.ExN.NAs <- deseq.fc.filter.ExN
toReplace.ExN <- which(!(rownames(deseq.fc.ExN.NAs) %in% ExN.tokeep))

for (pos in toReplace.ExN) {
  # deseq.fc.ExN.NAs[pos,] <- c(NaN, NaN, NaN)
}

deseq.fc.Water.NAs <- deseq.fc.Water.color2
toReplace.Water <- which(!(rownames(deseq.fc.Water.NAs) %in% Water.tokeep))

for(pos in toReplace.Water) {
  # deseq.fc.Water.NAs[pos,] <- c(NaN,NaN,NaN)
}

# Combine into one

deseq.fc.ALL.NAs <- cbind(deseq.fc.Water.NAs,deseq.fc.ExN.NAs)

deseq.sig.ALL <- cbind(deseq.sig.filter.Water, deseq.sig.filter.ExN)
deseq.stat.ALL <- cbind(deseq.stat.filter.Water, deseq.stat.filter.ExN)

# Change deseq.sig into a matrix of stars for significance
deseq.sig.ALL.star <- deseq.sig.ALL
for (r in 1:nrow(deseq.sig.ALL.star)) {
  for (c in 1:ncol(deseq.sig.ALL.star)) {
    if (is.na(deseq.sig.ALL.star[r,c])) {
      deseq.sig.ALL.star[r,c] <- "-"
    } else if (as.numeric(deseq.sig.ALL.star[r,c]) > 0.05) {
      deseq.sig.ALL.star[r,c] <- ""
    } else if (as.numeric(deseq.sig.ALL.star[r,c]) > 0.01) {
      deseq.sig.ALL.star[r,c] <- "*"
    } else if (as.numeric(deseq.sig.ALL.star[r,c]) > 0.001) {
      deseq.sig.ALL.star[r,c] <- "**"
    } else if (as.numeric(deseq.sig.ALL.star[r,c]) <= 0.001) {
      deseq.sig.ALL.star[r,c] <- "***"
    }
  }
}

######### PLOT ############
pdf("./DESEQ/Heatmap_enrichedOTUs.pdf", pointsize = 14)
heatmap.2(as.matrix(deseq.fc.ALL.NAs)
        , Rowv = NA
        , Colv = NA
        , labCol = c("+Nereo","+Mast","+Nereo+Mast","+Nereo","+Mast","+Nereo+Mast")
        , col = colorRange(ncolors)
        , scale = "none"
        , dendrogram = "none"
        , offsetCol = 3
        
        , trace = "none"
        , density.info = "none"
        , key.xlab = "Fold-change"
        
        , na.color = "white"
        
        , margins = c(0,0)
        , lmat = rbind(c(0,4,5),c(2,1,6),c(0,3,7))
        , lhei = c(2,5,3)
        , lwid = c(0.5,3,5.5)
        , colsep = c(3)
        , sepwidth = c(0.1,0)
        
        , cellnote = as.matrix(deseq.sig.ALL.star)
        , notecol = "black"
        , notecex = 2.0
        )
par(fig = c(0,1,0,1), mar = c(1,1,1,1), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , bty = "n"
     , xlab = ""
     , ylab = "")
text(x = 0.05, y = c(1,0.95,0.90,0.85)
     , labels = c(" -","  *"," **","***")
     , pos = 2
     , cex = 0.8)
text(x = 0.1, y = c(1,0.95,0.90,0.85)
     , labels = c( "Absent","p <= 0.05", "p <= 0.01","p <= 0.001")
     , pos = 4
     , cex = 0.8)
text(x = -0.85, y = c(-0.44,-0.5)
     , labels = c("_____________","WATER")
     , pos = 1
     , cex = 0.8
     , xpd = "n")
text(x = -0.45, y = c(-0.44,-0.5)
     , labels = c("_____________","NMF SURFACE")
     , pos = 1
     , cex = 0.8
     , xpd = "n")

dev.off()

