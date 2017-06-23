#!/bin/bash


# This script is for plotting enrichment
library("optparse")
########################### OPT PARSE #################################
option_list = list(
  make_option(c("-t", "--taxasumFP"), type="character",
              help="Taxa summaries input"),
  make_option(c("-m", "--metadata"), type="character",
              help="Metadata fp")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

taxasumFP = opt$taxasumFP
metadataFP = opt$metadata

########################### LOAD DATA #################################

setwd("/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis")
taxasumFP <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/ANALYSIS_ALPHABETATAXA/summarize_taxa/rarefied_OTU_Table_sorted_L6.txt"
MFFP <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/ANALYSIS_ALPHABETATAXA/OTU_Tables_and_MP/MF_withalpha.txt"

MF <- read.delim(paste0(MFFP)
                 , header = TRUE
                 , row.names = 1
                 , stringsAsFactors = FALSE)
rownames(MF) <- gsub("-",".", rownames(MF))

taxasum <- read.delim(paste0(taxasumFP)
                      , header = TRUE
                      , skip = 1
                      , row.names = 1
                      , stringsAsFactors = FALSE)


##########START ##############

set.seed(4)

system("mkdir TAXASUMMARIES")
######### FILT LOW ABUND ###########

todelete <- c()
for (i in 1:nrow(taxasum)) {
  nLowAbund <- sum(taxasum[i,] >= 0.001)
  if (nLowAbund <= 5) {
    todelete <- c(todelete, i)
  }
}
taxasum <- taxasum[-todelete,]

# Make sure MF and taxa summaries are in each other
MF <- MF[(rownames(MF) %in% colnames(taxasum)),]
taxasum <- taxasum[,(colnames(taxasum) %in% rownames(MF))]

# Reorder
MF <- MF[order(match(rownames(MF),colnames(taxasum))),]

# nrow(MF)
# ncol(taxasum)
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

########## COLOURS ##########
listSplit <- strsplit(rownames(taxasum), split = ";__", fixed = TRUE)
correctedNames <- c()
for (i in listSplit) {
  correctedNames <- c(correctedNames, paste0(i[3],": ", i[5],"_",i[6]))
}
correctedNames <- gsub("NA: NA_NA", "Unidentified", correctedNames)


randomColors <- sample(colors()[-grep("white|gr(a|ey)", colors())], nrow(taxasum))
colorLegend <- cbind(correctedNames, randomColors)


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

# Get IDs for legend that are greater than 5%
deleteColor <- c()
deleteTaxa <- c()
colorPlot.Loneincube <- colorLegend
for (i in 1:nrow(taxasum.Loneincube)) {
  if (max(taxasum.Loneincube[i,]) < 0.05) {
    deleteColor <- c(deleteColor, i)
    # colorPlot.Loneincube[i,2] <- "white"
    deleteTaxa <- c(deleteTaxa, i)
  }
}
colorLegend.Loneincube <- colorLegend[-deleteColor,]
taxasum.Loneincube <- taxasum.Loneincube[-deleteTaxa,]

pdf( "./TAXASUMMARIES/LoneIncube.pdf", pointsize = 14)
par(fig = c(0,0.65,0,1), mar = c(4,4,4,4))
barplot(as.matrix(taxasum.Loneincube)
        , col = colorLegend.Loneincube[,2]
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
par(fig = c(0.5,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n")
legend("center"
       , legend = rev(colorLegend.Loneincube[,1])
      , pch = 22
      , pt.bg = rev(colorLegend.Loneincube[,2])
      , col = rev(colorLegend.Loneincube[,2])
      , cex = 0.5
      , pt.cex = 1.2)
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

# Get IDs for legend that are greater than 5%
deleteColor <- c()
deleteTaxa <- c()
# colorPlot.Water <- colorLegend
for (i in 1:nrow(taxasum.Water)) {
  if (max(taxasum.Water[i,]) < 0.05) {
    deleteColor <- c(deleteColor, i)
    # colorPlot.Water[i,2] <- "white"
    deleteTaxa <- c(deleteTaxa, i)
  }
}
colorLegend.Water <- colorLegend[-deleteColor,]
taxasum.Water <- taxasum.Water[-deleteTaxa,]

pdf( "./TAXASUMMARIES/Water.pdf", pointsize = 14)
par(fig = c(0,0.65,0,1), mar = c(6,4,4,4))
barplot(as.matrix(taxasum.Water)
        , col = colorLegend.Water[,2]
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
par(fig = c(0.5,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n")
legend("center"
       , legend = rev(colorLegend.Water[,1])
       , pch = 22
       , pt.bg = rev(colorLegend.Water[,2])
       , col = rev(colorLegend.Water[,2])
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

# Get IDs for legend that are greater than 5%
deleteColor <- c()
deleteTaxa <- c()
# colorPlot.ExN <- colorLegend
for (i in 1:nrow(taxasum.ExN)) {
  if (max(taxasum.ExN[i,]) < 0.05) {
    deleteColor <- c(deleteColor, i)
    # colorPlot.ExN[i,2] <- "white"
    deleteTaxa <- c(deleteTaxa, i)
  }
}
colorLegend.ExN <- colorLegend[-deleteColor,]
taxasum.ExN <- taxasum.ExN[-deleteColor,]

pdf( "./TAXASUMMARIES/ExN.pdf", pointsize = 14)
par(fig = c(0,0.65,0,1), mar = c(6,4,4,4))
barplot(as.matrix(taxasum.ExN)
        , col = colorLegend.ExN[,2]
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
par(fig = c(0.5,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n")
legend("center"
       , legend = rev(colorLegend.ExN[,1])
       , pch = 22
       , pt.bg = rev(colorLegend.ExN[,2])
       , col = rev(colorLegend.ExN[,2])
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

# Get IDs for legend that are greater than 5%
deleteColor <- c()
deleteTaxa <- c()
# colorPlot.ExN <- colorLegend
for (i in 1:nrow(taxasum.Environ)) {
  if (max(taxasum.Environ[i,]) < 0.05) {
    deleteColor <- c(deleteColor, i)
    # colorPlot.ExN[i,2] <- "white"
    deleteTaxa <- c(deleteTaxa, i)
  }
}
colorLegend.Environ <- colorLegend[-deleteColor,]
taxasum.Environ <- taxasum.Environ[-deleteColor,]


pdf( "./TAXASUMMARIES/Environ.pdf", pointsize = 14)
par(fig = c(0,0.65,0,1), mar = c(6,4,4,4))
barplot(as.matrix(taxasum.Environ)
        , col = colorLegend.Environ[,2]
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
par(fig = c(0.5,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n")
legend("center"
       , legend = rev(colorLegend.Environ[,1])
       , pch = 22
       , pt.bg = rev(colorLegend.Environ[,2])
       , col = rev(colorLegend.Environ[,2])
       , cex = 0.5
       , pt.cex = 1.2)
dev.off()

################ ALL PLOTS TOGETHER #################
### WORKING
# Want to put NMF-incubation experiment together and seaweed surfaces together
# Make new legend that combines everything

duplicatedRows1 <- which(colorLegend.ExN[,1] %in% colorLegend.Water[,1])
comboLegend <- rbind(colorLegend.ExN[-duplicatedRows1,], colorLegend.Water)
comboLegend <- comboLegend[order(comboLegend[,1]),]
# Put unidentified at the top
postemp <- grep("Unidentified", comboLegend[,1])
comboLegend <- comboLegend[c(1:(postemp-1),(postemp+1):(nrow(comboLegend)-1),postemp),]

pdf("TAXASUMMARIES/TaxaSummaries_combo.pdf", pointsize = 14, width = 9)
par(fig = c(0,0.65,0.4,1))
barplot(as.matrix(taxasum.Water)
        , col = colorLegend.Water[,2]
        , las = 2
        , space = c(0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0)
        , border = FALSE
        , xaxt = "n"
        , ylab = "Relative Abundance"
)
axis(side = 1
     , las = 2
     , at = c(2,7,8,12.5,13.5,18.5,19.5,25,26)
     , labels = c("Water Alone","Water","with NMF only","Water", "with Nereo", "Water", "with Mast", "Water","with Nereo + Mast")
     , tick = FALSE
     , line = -1
     , cex.axis = 0.6
)

par(fig = c(0,0.65,0,0.55), new = TRUE)
#Need to add empty plot to make spacing right
barplot(as.matrix(cbind(rep(0, nrow(taxasum.ExN)),taxasum.ExN))
        , col = colorLegend.ExN[,2]
        , las = 2
        , space = c(0,4,0,0,0,0,1,0,0,0,0,2,0,0,3,0,0,0,0)
        , border = FALSE
        , xaxt = "n"
        , ylab = "Relative Abundance"
)
axis(side = 1
     , las = 2
     , at = c(7,13,14,19,20,26,27)
     , labels = c("NMF Alone","NMF", "with Nereo", "NMF", "with Mast", "NMF","with Nereo + Mast")
     , tick = FALSE
     , line = -1
     , cex.axis = 0.6
)
par(fig = c(0.6,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n")
legend("center"
       , legend = rev(comboLegend[,1])
       , pch = 22
       , pt.bg = rev(comboLegend[,2])
       , col = rev(comboLegend[,2])
       , cex = 0.5
       , pt.cex = 1.2)
dev.off()

write.matrix(comboLegend, sep = "\t", file = "TAXASUMMARIES/comboLegend.txt")
