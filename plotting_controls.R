#!/bin/bash/R

# This script plots the controls against all other samples
library("MASS")
setwd("~/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/CONTROL_PLOTTING/")
filepath <- c("~/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/CONTROL_PLOTTING/jackknife_beta")
MFFP <- c("~/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/CONTROL_PLOTTING/Edited_MF.txt")


MF <- read.delim(paste0(MFFP), header=TRUE, stringsAsFactors=FALSE, row.names=1)
dm.all <- list()
nmds.all <- list()
for ( m in c("bray_curtis","unweighted_unifrac","weighted_unifrac")) {
    dm.all[[paste0(m)]] <- read.delim(paste0(filepath,"/",m, "/rare_dm/",m,"_rarefaction_160_0.txt"), header=TRUE, row.names=1)
    dm.all[[paste0(m)]] <- dm.all[[paste0(m)]][match(rownames(MF), rownames(dm.all[[paste0(m)]])),match(rownames(MF), colnames(dm.all[[paste0(m)]]))]
    nmds.all[[paste0(m)]] <- isoMDS(d = dist(dm.all[[paste0(m)]]), k = 2)
    
    pdf(paste0("negcontrols_plot",m,".pdf"))
    plot(nmds.all[[paste0(m)]][["points"]]
         , pch = ""
         , xlab="NMDS1"
         , ylab="NMDS2"
         , sub=signif(nmds.all[[paste0(m)]][["stress"]]/100,2)
         )
    points(nmds.all[[paste0(m)]][["points"]]
           , pch=19
           , col=c("yellow","white","pink","purple")[factor(MF$Nextto)]
           , cex=1
           )
    points(nmds.all[[paste0(m)]][["points"]]
        , col=c("black","green","darkred","black","black","darkgreen","blue")[factor(MF$ConPlot)]
        , bg=c("yellow",NA,NA,"pink","purple",NA,NA)[factor(MF$ConPlot)]
        , pch=c(22,21,21,22,22,21,21)[factor(MF$ConPlot)]
        , cex=c(2,1,1,2,2,1,1)[factor(MF$ConPlot)]
    )
    
    legend("topright", legend =c("Extraction Con","PCR Con 1", "PCR Con 2","","Nereo Meristem","Nereo Blade","Mastocarpus","Water")
           , pch=c(22,22,22,NA,21,21,21,21)
           , pt.bg=c("yellow","pink","purple",NA,NA,NA,NA,NA,NA)
           , col=c("black","black","black",NA,"green","darkgreen","darkred","blue")
           , bty="n"
           )
    dev.off()
    
}
