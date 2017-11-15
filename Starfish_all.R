#!/bin/bash/R
MFFP <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/starfish/MF_star_withalpha.txt"
######## READ IN DATA #######

MF <- read.delim(paste0(MFFP), stringsAsFactors = FALSE, na.strings=c("NA","na","n/a",""), row.names = 1)


####### ALPHA #######
MF.star <- MF[which(MF[,"Project"] == "starfish"),]

metric <- "chao1"
fullmetric <- paste0(metric,"_even_2800_alpha")


quartz()
stripchart(as.numeric(MF.star[,paste0(fullmetric)])~ MF.star[,"ColRep"]
        , xlab = "Sample"
        , ylab = paste0(metric)
        , vertical = TRUE
        , pch = 21
        , col = c("darkgreen","darkblue","green","blue")
        )

levels(factor(MF.star[,"ColRep"]))
cbind(rownames(MF.star), MF.star[,paste0(fullmetric)])
