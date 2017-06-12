#!/bin/Rscript
library(optparse)
################ OPT PARSE ####################

option_list = list(
  make_option(c("-m", "--mappingfile"), type="character",
              help="Mapping file with alpha diversity in it"),
  make_option(c("-a", "--alphaNames"),
              help="Comma separated list of alpha Names", type="character")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

MFPWD = opt$mappingfile
alphaNamesTMP = opt$alphaNames
alphaNames <- unlist(strsplit(alphaNamesTMP, split = ","))


#####FORTESTING ########
# setwd("/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/zz_NEREOINCUBE_16may2017/1_analysis")
# MFPWD <- "/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/zz_NEREOINCUBE_16may2017/1_analysis/ANALYSIS_ALPHABETATAXA/OTU_Tables_and_MP/MF_withalpha.txt"
# alphaNames <- c("chao1","PD_whole_tree","observed_otus")

setwd("/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis")
MFPWD <-"/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/ANALYSIS_ALPHABETATAXA/OTU_Tables_and_MP/MF_withalpha_2.txt"
alphaNames <-  'chao1_even_1000_alpha,PD_whole_tree_even_1000_alpha,observed_otus_even_1000_alpha'
alphaList <- unlist(strsplit(alphaNames, ","))

############## LOAD DATA ########
library(car)
library(xtable)
# Alpha div script

system("mkdir ALPHAPLOTS")

MF <- read.delim(paste0(MFPWD)
                 , header = TRUE
                 , row.names = 1
                 , stringsAsFactors = FALSE
                 , na.strings = c("NA","N/A","na","<NA>"))

MF.filtered <- MF[,grep("chao1|PD_whole_tree|observed_otus|ColRep|Replicate", colnames(MF))]

# Specific tests
MF.ExN <- MF.filtered[grep("ExN.Nereotest.", rownames(MF.filtered)),]
MF.ExNWater <- MF.filtered[grep("Water.Nereotest.", rownames(MF.filtered)),]
MF.LoneIncube <- MF.filtered[grep("Loneincube", rownames(MF.filtered)),]


####### PLOT ##########

# For MF.ExN
########## ***ExN*** #############
for (i in alphaList) {
  tempI <- gsub("_even_1000_alpha","", i)
  
  ######## STATS ########
  ExN.lm <- lm(MF.ExN[,paste0(i)] ~ MF.ExN$ColRep)
  anova.ExN.lm <- Anova(ExN.lm, type = 'III') # OVERALL
  
  # Separate by Experiment
  MF.ExN.ExNExN <- MF.ExN[grep("ExNExN", MF.ExN$ColRep),]
  MF.ExN.ExNNereo <- MF.ExN[grep("NereoExN", MF.ExN$ColRep),]
  MF.ExN.ExNMast <- MF.ExN[grep("MastExN", MF.ExN$ColRep),]
  MF.ExN.ExNNereoMast <- MF.ExN[grep("NereoMastExN", MF.ExN$ColRep),]
  
  TTEST.ExN.ExN.Nereo <- t.test(MF.ExN.ExNExN[,paste(i)], MF.ExN.ExNNereo[,paste0(i)])
  TTEST.ExN.ExN.Mast <- t.test(MF.ExN.ExNExN[,paste(i)], MF.ExN.ExNMast[,paste0(i)])
  TTEST.ExN.ExN.NereoMast <- t.test(MF.ExN.ExNExN[,paste(i)], MF.ExN.ExNNereoMast[,paste0(i)])
  TTEST.ExN.Nereo.Mast <- t.test(MF.ExN.ExNNereo[,paste(i)], MF.ExN.ExNMast[,paste0(i)])
  TTEST.ExN.Nereo.NereoMast <- t.test(MF.ExN.ExNNereo[,paste(i)], MF.ExN.ExNNereoMast[,paste0(i)])
  TTEST.ExN.Mast.NereoMast <- t.test(MF.ExN.ExNMast[,paste(i)], MF.ExN.ExNNereoMast[,paste0(i)])
  INDIVTtest <- c("TTEST.ExN.ExN.Nereo"
                 , "TTEST.ExN.ExN.Mast"
                 , "TTEST.ExN.ExN.NereoMast"
                 , "TTEST.ExN.Nereo.Mast"
                 , "TTEST.ExN.Nereo.NereoMast"
                 , "TTEST.ExN.Mast.NereoMast")
  allPvalues.ExN <- c()
  for (test in INDIVTtest) {
    allPvalues.ExN <- rbind(allPvalues.ExN, c(get(test)$p.value, get(test)$statistic, get(test)$parameter))
  }
  allPvalues.ExN <- cbind(allPvalues.ExN, p.adjust(allPvalues.ExN[,1], method = "fdr", n = 6))
  colnames(allPvalues.ExN) <- c("P","t","df","fdr_adj")
  rownames(allPvalues.ExN) <- c("NMFNereo"
                                , "NMFMast"
                                , "NMFNereoMast"
                                , "NereoMast"
                                , "NereoNereoMast"
                                , "MastNereoMast")
  capture.output(xtable(allPvalues.ExN, digits = 3), file = paste0("ALPHAPLOTS/allPvalues.ExN.",tempI,".txt"))
  
  
  ######## PLOTTING########
  MF.ExN.temp <- MF.ExN[,c(paste0(i), "ColRep","Replicate")]
  MF.ExN.Alpha <- reshape(MF.ExN.temp, idvar = "ColRep", timevar = "Replicate", direction = "wide")
  MF.ExN.Alpha <- data.frame(MF.ExN.Alpha
                             , row.names = 1
                             )
  MF.ExN.Alpha <- MF.ExN.Alpha[sapply(c("NereotestExNExN"
                                        ,"NereotestNereoExN"
                                        ,"NereotestMastExN"
                                        ,"NereotestNereoMastExN"), function(x) {
    grep(x, rownames(MF.ExN.Alpha))
  }),]
  
  
  ncounts <- c()
  for (ROW in 1:nrow(MF.ExN.Alpha)) {
    ncounts <- c(ncounts, sum(!is.na(MF.ExN.Alpha[ROW,])))
  }
  
  
  rownames(MF.ExN.Alpha) <- c(paste0("Meristem alone (",ncounts[1],")")
                              ,paste0("With Nereo (",ncounts[2],")")
                              ,paste0("With Mast (",ncounts[3],")")
                              ,paste0("With Both (",ncounts[4],")")
                              )
  
  pdf(paste0("ALPHAPLOTS/Alpha_div_meristemswabs_",tempI,".pdf"))
  par(mar = c(10,4,4,4))
  boxplot(t(MF.ExN.Alpha)
          , las = 2
          , col = c("gray","green","red","brown")
          , ylab = paste0("Alpha Diversity (", tempI,")")
          , main = "Alpha diversity across meristem swabs"
          )
  title(xlab = "Treatment", line = 8)
  dev.off()
  
}


############## ****ExNWater*** ############
# For MF.ExNWater
for (i in alphaList) {
  tempI <- gsub("_even_1000_alpha", "", i)
  
  ######## STATS ########
  ExNWater.lm <- lm(MF.ExNWater[,paste0(i)] ~ MF.ExNWater$ColRep)
  anova.ExNWater.lm <- Anova(ExNWater.lm, type = 'III') # OVERALL
  
  # Separate by Experiment
  MF.ExNWater.H2O <- MF.ExNWater[grep("H2O", MF.ExNWater$ColRep),]
  MF.ExNWater.ExNWater <- MF.ExNWater[grep("ExNWater", MF.ExNWater$ColRep),]
  MF.ExNWater.NereoWater <- MF.ExNWater[grep("NereoWater", MF.ExNWater$ColRep),]
  MF.ExNWater.MastWater <- MF.ExNWater[grep("MastWater", MF.ExNWater$ColRep),]
  MF.ExNWater.NereoMastWater <- MF.ExNWater[grep("NereoMastWater", MF.ExNWater$ColRep),]
  
  TTEST.ExNWater.H2O.ExNWater <- t.test(MF.ExNWater.H2O[,paste(i)], MF.ExNWater.ExNWater[,paste0(i)])
  TTEST.ExNWater.H2O.Nereo <- t.test(MF.ExNWater.H2O[,paste(i)], MF.ExNWater.NereoWater[,paste0(i)])
  TTEST.ExNWater.H2O.Mast <- t.test(MF.ExNWater.H2O[,paste(i)], MF.ExNWater.MastWater[,paste0(i)])
  TTEST.ExNWater.H2O.NereoMast <- t.test(MF.ExNWater.H2O[,paste(i)], MF.ExNWater.NereoMastWater[,paste0(i)])
  TTEST.ExNWater.ExNWater.Nereo <- t.test(MF.ExNWater.ExNWater[,paste(i)], MF.ExNWater.NereoWater[,paste0(i)])
  TTEST.ExNWater.ExNWater.Mast <- t.test(MF.ExNWater.ExNWater[,paste(i)], MF.ExNWater.MastWater[,paste0(i)])
  TTEST.ExNWater.ExNWater.NereoMast <- t.test(MF.ExNWater.ExNWater[,paste(i)], MF.ExNWater.NereoMastWater[,paste0(i)])
  TTEST.ExNWater.Nereo.Mast <- t.test(MF.ExNWater.NereoWater[,paste(i)], MF.ExNWater.MastWater[,paste0(i)])
  TTEST.ExNWater.Nereo.NereoMast <- t.test(MF.ExNWater.NereoWater[,paste(i)], MF.ExNWater.NereoMastWater[,paste0(i)])
  TTEST.ExNWater.Mast.NereoMast <- t.test(MF.ExNWater.MastWater[,paste(i)], MF.ExNWater.NereoMastWater[,paste0(i)])
  INDIVTtest <- c("TTEST.ExNWater.H2O.ExNWater"
                  , "TTEST.ExNWater.ExNWater.Nereo"
                  , "TTEST.ExNWater.ExNWater.Mast"
                  , "TTEST.ExNWater.ExNWater.NereoMast"
                  , "TTEST.ExNWater.H2O.Nereo"
                  , "TTEST.ExNWater.H2O.Mast"
                  , "TTEST.ExNWater.H2O.NereoMast"
                  , "TTEST.ExNWater.Nereo.Mast"
                  , "TTEST.ExNWater.Nereo.NereoMast"
                  , "TTEST.ExNWater.Mast.NereoMast")
  allPvalues.ExNWater <- c()
  for (test in INDIVTtest) {
    allPvalues.ExNWater <- rbind(allPvalues.ExNWater, c(get(test)$p.value, get(test)$statistic, get(test)$parameter))
  }
  allPvalues.ExNWater <- cbind(allPvalues.ExNWater, p.adjust(allPvalues.ExNWater[,1], method = "fdr", n = 10))
  colnames(allPvalues.ExNWater) <- c("P","t","df","fdr_adj")
  rownames(allPvalues.ExNWater) <- c("WaterNMF"
                                     , "NMFNereo"
                                     , "NMFMast"
                                     , "NMFNereoMast"
                                , "WaterNereo"
                                , "WaterMast"
                                , "WaterNereoMast"
                                , "NereoMast"
                                , "NereoNereoMast"
                                , "MastNereoMast")
  capture.output(xtable(allPvalues.ExNWater, digits = 3), file = paste0("ALPHAPLOTS/allPvalues.ExNWater.",tempI,".txt"))
  
  ######### PLOTTING ########
  
  MF.ExNWater.temp <- MF.ExNWater[,c(paste0(i), "ColRep","Replicate")]
  MF.ExNWater.Alpha <- reshape(MF.ExNWater.temp, idvar = "ColRep", timevar = "Replicate", direction = "wide")
  MF.ExNWater.Alpha <- data.frame(MF.ExNWater.Alpha
                             , row.names = 1
  )
  MF.ExNWater.Alpha <- MF.ExNWater.Alpha[sapply(c("NereotestH2OWater"
                                                  ,"NereotestExNWater"
                                                  ,"NereotestNereoWater"
                                                  ,"NereotestMastWater"
                                                  , "NereotestNereoMastWater"), function(x) {
    grep(x, rownames(MF.ExNWater.Alpha))
  }),]
  
  
  ncounts <- c()
  for (ROW in 1:nrow(MF.ExNWater.Alpha)) {
    ncounts <- c(ncounts, sum(!is.na(MF.ExNWater.Alpha[ROW,])))
  }
  
  
  rownames(MF.ExNWater.Alpha) <- c(paste0("Water Alone (",ncounts[1],")")
                                   ,paste0("With Meristem (",ncounts[2],")")
                                   ,paste0("With Nereo (",ncounts[3],")")
                                   ,paste0("With Mast (",ncounts[4],")")
                                   ,paste0("With Both (",ncounts[5],")")
  )
  
  pdf(paste0("ALPHAPLOTS/Alpha_div_watersamples_",tempI,".pdf"))
  par(mar = c(10,4,4,4))
  boxplot(t(MF.ExNWater.Alpha)
          , las = 2
          , col = c("blue","gray","green","red","brown")
          , ylab = paste0("Alpha Diversity (", tempI,")")
          , main = "Alpha diversity across water samples"
  )
  title(xlab = "Treatment", line = 8)
  dev.off()

  
}

###########***LONEINCUBE***###########
# For MF.LoneIncube
for (i in alphaList) {
  tempI <- gsub("_even_1000_alpha", "",i)
  
  ######## STATS ########
  LoneIncube.lm <- lm(MF.LoneIncube[,paste0(i)] ~ MF.LoneIncube$ColRep)
  anova.LoneIncube.lm <- Anova(LoneIncube.lm, type = 'III') # OVERALL
  
  # Separate by Experiment
  MF.LoneIncube.NereoNereo <- MF.LoneIncube[grep("NereoNereo", MF.LoneIncube$ColRep),]
  MF.LoneIncube.Nereowater <- MF.LoneIncube[grep("Nereowater", MF.LoneIncube$ColRep),]
  MF.LoneIncube.MastMast <- MF.LoneIncube[grep("MastMast", MF.LoneIncube$ColRep),]
  MF.LoneIncube.Mastwater <- MF.LoneIncube[grep("Mastwater", MF.LoneIncube$ColRep),]
  
  TTEST.LoneIncube.Nereo.Mast <- t.test(MF.LoneIncube.NereoNereo[,paste(i)], MF.LoneIncube.MastMast[,paste0(i)])
  TTEST.LoneIncube.Nereo.Nwater <- t.test(MF.LoneIncube.NereoNereo[,paste(i)], MF.LoneIncube.Nereowater[,paste0(i)])
  TTEST.LoneIncube.Nereo.Mwater <- t.test(MF.LoneIncube.NereoNereo[,paste(i)], MF.LoneIncube.Mastwater[,paste0(i)])
  TTEST.LoneIncube.Mast.Nwater <- t.test(MF.LoneIncube.MastMast[,paste(i)], MF.LoneIncube.Nereowater[,paste0(i)])
  TTEST.LoneIncube.Mast.Mwater <- t.test(MF.LoneIncube.MastMast[,paste(i)], MF.LoneIncube.Mastwater[,paste0(i)])
  TTEST.LoneIncube.Nwater.Mwater <- t.test(MF.LoneIncube.Nereowater[,paste(i)], MF.LoneIncube.Mastwater[,paste0(i)])
  INDIVTtest <- c("TTEST.LoneIncube.Nereo.Mast"
                  , "TTEST.LoneIncube.Nereo.Nwater"
                  , "TTEST.LoneIncube.Nereo.Mwater"
                  , "TTEST.LoneIncube.Mast.Nwater"
                  , "TTEST.LoneIncube.Mast.Mwater"
                  , "TTEST.LoneIncube.Nwater.Mwater")
  allPvalues.LoneIncube <- c()
  for (test in INDIVTtest) {
    allPvalues.LoneIncube <- rbind(allPvalues.LoneIncube, c(get(test)$p.value, get(test)$statistic, get(test)$parameter))
  }
  allPvalues.LoneIncube <- cbind(allPvalues.LoneIncube, p.adjust(allPvalues.LoneIncube[,1], method = "fdr", n = 6))
  colnames(allPvalues.LoneIncube) <- c("P","t","df","fdr_adj")
  rownames(allPvalues.LoneIncube) <- c("NereoMast"
                                , "NereoNereoWater"
                                , "NereoMastWater"
                                , "MastNereoWater"
                                , "MastMastWater"
                                , "NereoWaterMastWater")
  capture.output(xtable(allPvalues.LoneIncube, digits = 3), file = paste0("ALPHAPLOTS/allPvalues.LoneIncube.",tempI,".txt"))
  
  ###### PLOTTING #########
  MF.LoneIncube.temp <- MF.LoneIncube[,c(paste0(i), "ColRep","Replicate")]
  MF.LoneIncube.Alpha <- reshape(MF.LoneIncube.temp, idvar = "ColRep", timevar = "Replicate", direction = "wide")
  MF.LoneIncube.Alpha <- data.frame(MF.LoneIncube.Alpha
                                  , row.names = 1
  )
  MF.LoneIncube.Alpha <- MF.LoneIncube.Alpha[sapply(c("LoneincubeNereoNereo"
                                                  ,"LoneincubeMastMast"
                                                  ,"LoneincubeNereowater"
                                                  ,"LoneincubeMastwater"), function(x) {
                                                    grep(x, rownames(MF.LoneIncube.Alpha))
                                                  }),]
  
  
  ncounts <- c()
  for (ROW in 1:nrow(MF.LoneIncube.Alpha)) {
    ncounts <- c(ncounts, sum(!is.na(MF.LoneIncube.Alpha[ROW,])))
  }
  
  
  rownames(MF.LoneIncube.Alpha) <- c(paste0("Nereo Swab (",ncounts[1],")")
                                   ,paste0("Mast Swab (",ncounts[2],")")
                                   ,paste0("Nereo Water (",ncounts[3],")")
                                   ,paste0("Mast Water (",ncounts[4],")")
  )
  
  pdf(paste0("ALPHAPLOTS/Alpha_div_loneincube_",i,".pdf"))
  par(mar = c(10,4,4,4))
  boxplot(t(MF.LoneIncube.Alpha)
          , las = 2
          , col = c("green","red","blue","purple")
          , ylab = paste0("Alpha Diversity (", i,")")
          , main = "Alpha diversity"
  )
  title(xlab = "Treatment", line = 8)
  dev.off()
  
  
}

######### **NMF vs Nereo and Mast** ###########

MF.algae <- MF[grep("ExNExN|MastMast|NereoNereo|BrocktonMast|BrocktonOldNereo", MF$ColRep),]
sortedGroups <- c("NereotestExNExN"
                  , "LoneincubeNereoNereo"
                  , "LoneincubeMastMast"
                  , "EnvironmentalBrocktonOldNereo"
                  , "EnvironmentalBrocktonMast"
                  )
orderList <- c()
for (group in sortedGroups) {
  orderList <- c(orderList, grep(paste0(group), MF.algae$ColRep))
}
MF.algae <- MF.algae[orderList,]
MF.algae$ColRep <- factor(MF.algae$ColRep, levels = sortedGroups)
for (i in alphaList) {
  tempI <- gsub("_even_1000_alpha", "", i)
  MF.algae.filt <- MF.algae[,c(i, "ColRep")]
  
  ########### STATS ###########
  MFnames <- c()
  for (group in sortedGroups) {
    assign(paste0("MF.algae.",group), MF.algae[grep(group, MF.algae$ColRep),])
    MFnames <- c(MFnames, paste0("MF.algae.",group))
  }
  
  TTEST.ExN.Nereo <- t.test(get(MFnames[1])[,i],get(MFnames[2])[,i] )
  TTEST.ExN.Mast <- t.test(get(MFnames[1])[,i],get(MFnames[3])[,i] )
  TTEST.ExN.ENereo <- t.test(get(MFnames[1])[,i],get(MFnames[4])[,i] )
  TTEST.ExN.EMast <- t.test(get(MFnames[1])[,i],get(MFnames[5])[,i] )
  TTEST.Nereo.Mast <- t.test(get(MFnames[2])[,i],get(MFnames[3])[,i] )
  TTEST.Nereo.ENereo <- t.test(get(MFnames[2])[,i],get(MFnames[4])[,i] )
  TTEST.Nereo.EMast <- t.test(get(MFnames[2])[,i],get(MFnames[5])[,i] )
  TTEST.Mast.ENereo <- t.test(get(MFnames[3])[,i],get(MFnames[4])[,i] )
  TTEST.Mast.EMast <- t.test(get(MFnames[3])[,i],get(MFnames[5])[,i] )
  TTEST.ENereo.EMast <- t.test(get(MFnames[4])[,i],get(MFnames[5])[,i] )
  
  allTEST <- c("TTEST.ExN.Nereo"
  , "TTEST.ExN.Mast"
  , "TTEST.ExN.ENereo"
  , "TTEST.ExN.EMast"
  , "TTEST.Nereo.Mast"
  , "TTEST.Nereo.ENereo"
  , "TTEST.Nereo.EMast"
  , "TTEST.Mast.ENereo"
  , "TTEST.Mast.EMast"
  , "TTEST.ENereo.EMast")
  
  allPvalues.algae <- c()
  for (test in allTEST) {
    allPvalues.algae <- rbind(allPvalues.algae, c(get(test)$p.value, get(test)$statistic, get(test)$parameter))
  }
  allPvalues.algae <- cbind(allPvalues.algae, p.adjust(allPvalues.algae[,1], method = "fdr", n = 10))
  colnames(allPvalues.algae) <- c("P","t","df","fdr_adj")
  
  newRowNames <- c("NMFNereo"
                   , "NMFMast"
                   , "NMFEnvironNereo"
                   , "NMFEnvironMast"
                   , "NereoMast"
                   , "NereoEnvironNereo"
                   , "NereoEnvironMast"
                   , "MastEnvironNereo"
                   , "MastEnvironMast"
                   , "EnvironNereoEnvironMast")
  rownames(allPvalues.algae) <- newRowNames
  capture.output(xtable(allPvalues.algae, digits = 3), file = paste0("ALPHAPLOTS/allPvalues.algae.",tempI,".txt"))
   
  
  ########### PLOT ############

  pdf(paste0("ALPHAPLOTS/Alpha_algaecompare",tempI, ".pdf"), pointsize = 14)
  par(mar = c(8,5,4,4))
  plot(MF.algae.filt[,paste0(i)]~ factor(MF.algae.filt[,"ColRep"])
       , ylab = paste0("Alpha Diversity (",tempI,")")
       , xlab = NA
       , col = c("gray","green","red","darkgreen","darkred")
       , las = 2 
       , xaxt = "n")
  axis(side = 1
       , labels = c("NMF", "Lab Nereo","Lab Mast",expression(italic(In~situ)~Nereo), expression(italic(In~situ)~Mast))
       , at = c(1,2,3,4,5)
       , las = 2
  )
  title(xlab = "Seaweed Type"
        , line = 6)
  dev.off()
}

