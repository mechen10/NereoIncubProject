#!/bin/bash

library("gplots")
## Plotting DEseq results (post combining, using OTUTraits.py)
setwd("/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/DESEQ")
fcNMFFP <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/DESEQ/ExNcompare/log2FoldChange.txt"
fcwaterFP <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/DESEQ/Watercompare/log2FoldChange.txt"
over5.FP <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/TAXASUMMARIES/comboLegend.txt"

################### LOAD DATA ####################
fcNMF <- read.delim(paste0(fcNMFFP)
                 , row.names = 1
                 , header = TRUE
                 , stringsAsFactors = FALSE
                 , na.strings = c("NA")
                 , strip.white = TRUE)
# Get rid of redundant second taxonomy
fcNMF <- fcNMF[,-grep('taxonomy', colnames(fcNMF))]

fcwater <- read.delim(paste0(fcwaterFP)
                    , row.names = 1
                    , header = TRUE
                    , stringsAsFactors = FALSE
                    , na.strings = c("NA")
                    , strip.white = TRUE)
# Get rid of redundant second taxonomy
fcwater <- fcwater[,-grep('taxonomy', colnames(fcwater))]


over5 <- read.delim(paste0(over5.FP)
                    , header = TRUE
                    , row.names = 1
                    , strip.white = TRUE)

################### BEGIN ######################
 # Processing #
# Switch NAs to zeros
for (r in 1:nrow(fcNMF)) {
  for (c in 1:ncol(fcNMF)) {
    if (is.na(fcNMF[r,c])) {
      fcNMF[r,c] <- 0
    }
  }
}

for (r in 1:nrow(fcwater)) {
  for (c in 1:ncol(fcwater)) {
    if (is.na(fcwater[r,c])) {
      fcwater[r,c] <- 0
    }
  }
}

#### Color ####

# get color gradient
orangegreen <- colorRampPalette(c("orange","white","green"))

#### Sort samples ####
# Sort samples
sampleOrder.NMF <- c(
                 "NereotestExNExNvsNereotestNereoExN"
                 ,"NereotestExNExNvsNereotestMastExN"
                 ,"NereotestExNExNvsNereotestNereoMastExN")
fcNMF.sorted <- fcNMF[,sapply(sampleOrder.NMF, function(x) grep(x, colnames(fcNMF)))]

sampleOrder.water <- c("NereotestExNWatervsNereotestNereoWater"
                     ,"NereotestExNWatervsNereotestMastWater"
                     ,"NereotestExNWatervsNereotestNereoMastWater"
                     )
fcwater.sorted <- fcwater[,sapply(sampleOrder.water, function(x) grep(x, colnames(fcwater)))]

#For labelling in plt
forLabCol.NMF <- c("NMF: with Nereo"
                  ,"NMF: with Mast"
                  ,"NMF: with Nereo + Mast")
forLabCol.water <- c("WATER: with Nereo"
                   ,"WATER: with Mast"
                   ,"WATER: with Nereo + Mast"
                  )
colnames(fcNMF.sorted) <- forLabCol.NMF
colnames(fcwater.sorted) <- forLabCol.water
# sort OTUs by dendrogram
# hc.NMF <- hclust(dist(1-cor(t(fcNMF.sorted))))
# hc.water <- hclust(dist(1-cor(t(fcwater.sorted))))

# Sort OTUs manually!
allenriched.NMF <- c()
allreduced.NMF <- c()
rest.NMF <- c()
for (r in 1:nrow(fcNMF.sorted)) {
  if (all(fcNMF.sorted[r,] > 0)) {
    allenriched.NMF <- c(allenriched.NMF, r)
  } else if (all(fcNMF.sorted[r,] < 0)) {
    allenriched.NMF <- c(allenriched.NMF,r)
  } else {
    rest.NMF <- rbind(rest.NMF, fcNMF.sorted[r,])
  }
}
# Combine
rest.NMF <- rest.NMF[with(as.data.frame(rest.NMF), order(-rest.NMF[,"NMF: with Nereo"]
                                                     ,-rest.NMF[,"NMF: with Mast"]
                                                     ,-rest.NMF[,"NMF: with Nereo + Mast"])), ]
fcNMF.sorted2 <- rbind(fcNMF.sorted[allenriched.NMF,], rest.NMF, fcNMF.sorted[allreduced.NMF,])

allenriched.water <- c()
allreduced.water <- c()
rest.water <- c()
for (r in 1:nrow(fcwater.sorted)) {
  if (all(fcwater.sorted[r,] > 0)) {
    allenriched.water <- c(allenriched.water, r)
  } else if (all(fcwater.sorted[r,] < 0)) {
    allenriched.water <- c(allenriched.water,r)
  } else {
    rest.water <- rbind(rest.water, fcwater.sorted[r,])
  }
}
# Combine
rest.water <- rest.water[with(as.data.frame(rest.water), order(-rest.water[,"WATER: with Nereo"]
                                                         ,-rest.water[,"WATER: with Mast"]
                                                         ,-rest.water[,"WATER: with Nereo + Mast"])), ]
fcwater.sorted2 <- rbind(fcwater.sorted[allenriched.water,]
                         , rest.water
                         , fcwater.sorted[allreduced.water,])



#### Taxonomy Names ####

# Get taxonomy names
taxonomyList.NMF <- list()
taxonomyNames.NMF <- c()
for (i in 1:nrow(fcNMF.sorted2)) {
  taxonomyList.NMF[[i]] <- unlist(strsplit(rownames(fcNMF.sorted2)[i], ";__"))
  toPaste <- paste0(taxonomyList.NMF[[i]][3], ": ", taxonomyList.NMF[[i]][5],"_",taxonomyList.NMF[[i]][6])
  if (toPaste == "NA: NA_NA") {
    toPaste <- "Unidentified"
  }
  taxonomyNames.NMF <- c(taxonomyNames.NMF, toPaste)
}
rownames(fcNMF.sorted2) <- taxonomyNames.NMF

taxonomyList.water <- list()
taxonomyNames.water <- c()
for (i in 1:nrow(fcwater.sorted2)) {
  taxonomyList.water[[i]] <- unlist(strsplit(rownames(fcwater.sorted2)[i], ";__"))
  toPaste <- paste0(taxonomyList.water[[i]][3], ": ", taxonomyList.water[[i]][5],"_",taxonomyList.water[[i]][6])
  if (toPaste == "NA: NA_NA") {
    toPaste <- "Unidentified"
  }
  taxonomyNames.water <- c(taxonomyNames.water, toPaste)
}
rownames(fcwater.sorted2) <- taxonomyNames.water


#### Filter based on >5 ####
# NMF
fcNMF.sorted3 <- fcNMF.sorted2[which(rownames(fcNMF.sorted2) %in% rownames(over5)),]
# Water
fcwater.sorted3 <- fcwater.sorted2[which(rownames(fcwater.sorted2) %in% rownames(over5)),]

# Filter out "Unidentified"
fcwater.sorted3 <- fcwater.sorted3[-grep("Unidentified", rownames(fcwater.sorted3)),]

#### Plot ####
pdf("EnrichedReducedTaxa_NMFsurface.pdf", pointsize = 14)
heatmap.2(as.matrix(fcNMF.sorted3)
          , Rowv = NA
          , cexRow = 1
          , cexCol = 1
          , labCol = forLabCol.NMF
          , Colv = NA
          , scale = "none"
          , dendrogram = "none"
          , col = orangegreen(11)

          
          , trace = "none"
          , density.info = "none"
          
          , margins = c(0,0)
          , lmat = rbind(c(0,4,5),c(2,1,6),c(0,3,7))
          , lhei = c(2,5,2)
          , lwid = c(0.5,3,4)
          
          , key.xlab = "Fold-change"
          , keysize = 0.5
          )
dev.off()

pdf("EnrichedReducedTaxa_waterColumn.pdf", pointsize = 14, width = 7, height = 14)
heatmap.2(as.matrix(fcwater.sorted3)
          , Rowv = NA
          , cexRow = 1
          , cexCol = 1
          , labCol = forLabCol.water
          , Colv = NA
          , scale = "none"
          , dendrogram = "none"
          , col = orangegreen(11)

          
          , trace = "none"
          , density.info = "none"
          
          , margins = c(0,0)
          , lmat = rbind(c(0,4,5),c(2,1,6),c(0,3,7))
          , lhei = c(2,10,2)
          , lwid = c(0.5,3,4)
          
          , key.xlab = "Fold-change"
          , keysize = 0.5
)
dev.off()


