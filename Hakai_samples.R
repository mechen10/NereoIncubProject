#!/bin/bash/R

# This looks at starfish samples to see if there is a significant difference between samples

dmFP <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/ANALYSIS_ALPHABETATAXA/beta_div/bray_curtis_dm.txt"
mfFP <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/ANALYSIS_ALPHABETATAXA/OTU_Tables_and_MP/MF_withalpha.txt"


dm <- read.delim(paste0(dmFP), row.names = 1, header = TRUE, stringsAsFactors = FALSE)
MF <- read.delim(paste0(mfFP), row.names = 1, header = TRUE, stringsAsFactors = FALSE)

# correct order
MF <- MF[rownames(dm),]

### Filter for starfisth ###
MF.star <- MF[grep("Starfish", rownames(MF)),]
dm.star <- dm[match(rownames(MF.star),rownames(dm)),match(rownames(MF.star),colnames(dm))]

### NMDS ###
NMDS <- isoMDS(d = dist(dm), k = 2)

# Filter
pos <- grep("Starfish",rownames(NMDS$points))
newNMDS <- NMDS$points[pos,]
MF.star.filt <- MF.star[rownames(newNMDS),]

quartz()
plot(newNMDS
     , pch = c(18,21)[factor(MF.star.filt$Treatment)]
     , col = c("green","blue")[factor(MF.star.filt$SubstrateType)]
)

quartz()
plot(NMDS$points
     , pch = 21
     , col = c("lightgreen","darkred","green","blue")[factor(MF$SubstrateType)]
     )
