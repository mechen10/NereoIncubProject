#!/bin/bash

# For assesing pH, DO, and temp differences between treatments
# Input is response Variables (comma separated) and treatments (Column of metadata)
# If groups included, will only calculate for those specific groups.
# If groups not included, will do all.

library("optparse")
########################### OPT PARSE #################################
option_list = list(
  make_option(c("-r", "--responseVariables"), type="character",
              help="Full OTU Table"),
  make_option(c("-t", "--treatmentVariables"),
              help="Column in metadata that tells you which treatments to compare", type="character"),
  make_option(c("-f", "--filterGroups"),
              help="Column in metadata that tells you how to separate data into treatments. 
              If not included, will calculated as one data set.
              Formatting should be Column:group1,group2"
              , type="character", default = NULL),
  make_option(c("-m", "--metadata"), help="metadata", type="character"),
  make_option(c("-p","--print"), help ="Include flag if you want to save output into files", type = "logical", )
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

responseVariables = opt$responseVariables
treatmentVariables = opt$treatmentVariables
filterGroups = opt$filterGroups
metadataFP = opt$metadata
########################### FOR TESTING #################################

setwd("/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/R")
metadataFP <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/R/NMFMetadata_june232017.txt"
responseVariables <- 'pH,DO_percent_sat,Salinity_ppt,Temp_C,Growth_in_cm'
treatmentVariables <- 'Treatment'
filterGroups <- 'Type:Loneincube,Nereotest'
print <- "FALSE"

########################### LOAD DATA #################################
system("mkdir WATERCOND")

metadata <- read.delim(paste0(metadataFP)
                       , header = TRUE
                       , strip.white = TRUE
                       , na.strings = c("NA","na","")
                       , row.names = 1)
rVar <- unlist(strsplit(responseVariables, split = ","))
tVar <- unlist(strsplit(treatmentVariables, split = ","))

# If filterGroups is null, then do not filter data table first.
# Either way, make a 'vector' of metadata files to iterate through for tests
listMeta <- c()
if (is.null(filterGroups)) {
  listMeta <- c('metadata.all')
} else { # else, filter metadata each time to include only those samples
  # Get column name
  fGroups <- unlist(strsplit(filterGroups, split = ":"))
  filtCol <- fGroups[1]
  # Get groups within columns
  fGrp <- unlist(strsplit(fGroups[2], split = ","))
  for (g in fGrp) {
    assign(paste0('metadata.',g), metadata[which(metadata[,paste0(filtCol)]==g),])
    listMeta <- c(listMeta, paste0('metadata.',g))
  }
}

# Calculated ANOVAs for list of metadata files and save names of files
allResults <- c()
for (met in listMeta) {
  nameM <- gsub("metadata.","",met)
  for (rV in 1:length(rVar)) {
    nameR <- rVar[rV]
    for (tV in 1:length(tVar)) {
      nameT <- tVar[tV]
      if ( !all(is.na(get(met)[,paste0(nameR)])) ) {
          temp.met <- get(met)
          temp.lm <- lm(get(nameR) ~ get(nameT), data = get(met))
          assign(paste0("anova.",nameM,".",nameR), anova(temp.lm))
          temp.anova <- aov(temp.met[,nameR] ~ temp.met[,nameT])
          assign(paste0("tukey.",nameR,".",nameM), TukeyHSD(x = temp.anova))
          assign(paste0("pairwise.",nameR,".",nameM), pairwise.t.test(metadata.Nereotest[,nameR], metadata.Nereotest[,nameT], p.adjust.method = "fdr"))
          allResults <- c(allResults, paste0("anova.",nameM,".",nameR), paste0("pairwise.",nameR,".",nameM), paste0("tukey.temp.",nameR,".",nameM))
      }
    }
  }
}


# Make plots of all results
for (met in listMeta) {
  nameM <- gsub("metadata.","",met)
  for (rV in 1:length(rVar)) {
    nameR <- rVar[rV]
    for (tV in 1:length(tVar)) {
      nameT <- tVar[tV]
      if ( !all(is.na(get(met)[,nameR])) ) {
          toPrint <- ""
          pvalue <- get(paste0("anova.",nameM,".",nameR))[1,"Pr(>F)"]
          if ( pvalue < 0.05) {
              toPrint <- "*"
          } else if ( pvalue < 0.01) {
              toPrint <- "**"
          } else if ( pvalue < 0.001) {
              toPrint <- "***"
          }
          pdf(paste0("WATERCOND/",nameR,"by",nameT,"in",nameM,".pdf"), pointsize = 14)
          plot(get(met)[,nameR] ~ factor(get(met)[,nameT])
               , main = paste0(nameR, " in ", nameM)
               , las = 2
               , ylab = paste0(nameR)
               , xlab = paste0(nameT)
               , pch = 21
          )
          mtext(line = -3,at = c(median(1:length(levels(factor(get(met)[,nameT]))))),text = toPrint, col = "red", cex = 3.5)
          dev.off()
      }
      
    }
  }
}


# Print and capture output of files
for (r in allResults) {
  capture.output(get(r), file = paste0("WATERCOND/",r,".txt"))
}

metadata.NMF <- metadata.Nereotest[metadata.Nereotest[,"Treatment"] != "H2O",]

### Custom growth plot ###
# Two-panel growth figure
pdf("NMFgrowth_bygroup.pdf")
par( mar = c(6.1,5.1,4.1,2.1))
plot(metadata.NMF[,"Growth_in_cm"] ~ factor(metadata.NMF[,"Treatment"])
     , las = 2
     , ylab = ""
     , xlab = ""
     , pch = 21
)
title(ylab = expression("Growth in area (cm"^2*")") )
title(xlab = "Treatment", line = 5)
dev.off()

pdf("NMFgrowth_linear.pdf")
par( mar = c(6.1,5.1,4.1,2.1))
plot(metadata.NMF$AreAfter ~metadata.NMF$AreaBefore
     , xlab = ""
     , ylab = ""
)
title(ylab = expression("Area After (cm"^2*")") )
title(xlab = expression("Area Before (cm"^2*")"), line = 5)
dev.off()

# 
# beforeall <- metadata[,c("WidthBefore","LengthBefore")]
# afterall <- metadata[,c("WidthAfter","LengthAfter")]
# 
# 
# 
# quartz()
# plot(NULL,NULL
#      , xlim=c(min(beforeall, na.rm = TRUE),max(beforeall, na.rm=TRUE))
#      , ylim=c(min(afterall, na.rm=TRUE),max(afterall,na.rm=TRUE))
#      )
# points(metadata$LengthAfter ~metadata$LengthBefore, col="orange")
# points(metadata$WidthAfter ~ metadata$WidthBefore, col="blue")
# 

