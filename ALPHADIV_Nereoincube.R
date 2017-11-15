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

setwd("/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis")
MFPWD <-"/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/ANALYSIS_ALPHABETATAXA/OTU_Tables_and_MP/MF_withalpha.txt"
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
    ptemp <- get(test)$p.value
    ttemp <- round(get(test)$statistic,2)
    dftemp <- round(get(test)$parameter, 2)
    toPaste <- paste0("(t=",ttemp,", df=",dftemp,")")
    allPvalues.ExN <- rbind(allPvalues.ExN, c(ptemp, toPaste))
  }
  allPvalues.ExN <- cbind(signif(as.numeric(allPvalues.ExN[,1]),3), allPvalues.ExN[,2], c(signif(p.adjust(allPvalues.ExN[1:3,1], method = "fdr", n = 3),3),signif(p.adjust(allPvalues.ExN[4:6,1], method = "fdr", n = 3),3) ))
  allPvalues.ExN <- cbind(c("\\NMF NMF"
                            ,"\\NMF NMF"
                            ,"\\NMF NMF"
                            , "\\Nereo Nereo"
                            , "\\Nereo Nereo"
                            , "\\Mast Mast")
                          , c("\\Nereo Nereo"
                              , "\\Mast Mast"
                              , "\\NereoMast Nereo + Mast"
                              , "\\Mast Mast"
                              , "\\NereoMast Nereo + Mast"
                              , "\\NereoMast Nereo + Mast")
                          , allPvalues.ExN)
  colnames(allPvalues.ExN) <- c("Group 1","Group 2","p","  ","FDR adj. p")
  
  capture.output(print(xtable(allPvalues.ExN, digits = 3), include.rownames = FALSE), file = paste0("ALPHAPLOTS/allPvalues.ExN.",tempI,".txt"))
  
  # Now, the ExN vs everything else 
  ttest.ExNvsEverythingelse <- t.test(MF.ExN.ExNExN[,paste0(i)], c(MF.ExN.ExNNereo[,paste0(i)],MF.ExN.ExNMast[,paste0(i)],MF.ExN.ExNNereoMast[,paste0(i)]))
  capture.output(ttest.ExNvsEverythingelse, file = paste0("ALPHAPLOTS/ttest.ExNvsEverythingelse",tempI, ".txt"))
  
  capture.output(table(MF.ExN$ColRep), file = "ALPHAPLOTS/exn_reps.txt")
  
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
  
  
  rownames(MF.ExN.Alpha) <- c(paste0("NMF alone (",ncounts[1],")")
                              ,paste0("With Nereo (",ncounts[2],")")
                              ,paste0("With Mast (",ncounts[3],")")
                              ,paste0("With Both (",ncounts[4],")")
                              )
  
  pdf(paste0("ALPHAPLOTS/Alpha_div_meristemswabs_",tempI,".pdf"))
  par(mar = c(10,4,4,4))
  boxplot(t(MF.ExN.Alpha)
          , las = 2
          , col = c("gray","green","purple","brown")
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
    ptemp <- get(test)$p.value
    ttemp <- round(get(test)$statistic,2)
    dftemp <- round(get(test)$parameter,2)
    toPaste <- paste0("(t=",ttemp,", df=",dftemp,")")
    allPvalues.ExNWater <- rbind(allPvalues.ExNWater, c(ptemp, toPaste))
  }
  allPvalues.ExNWater <- cbind(signif(as.numeric(allPvalues.ExNWater[,1]),3), allPvalues.ExNWater[,2]
                               , c( allPvalues.ExNWater[1,1]
                                 , signif(p.adjust(allPvalues.ExNWater[2:4,1], method = "fdr", n = 3),3)
                                   , signif(p.adjust(allPvalues.ExNWater[5:7,1], method = "fdr", n = 3),3)
                                   , signif(p.adjust(allPvalues.ExNWater[8:10,1], method = "fdr", n = 3),3))
                               )
  allPvalues.ExNWater <- cbind(c("\\NMF NMF"
                                 ,"\\NMF NMF"
                                 ,"\\NMF NMF"
                                 ,"\\NMF NMF"
                                 ,"\\wateronly Water only"
                                 ,"\\wateronly Water only"
                                 ,"\\wateronly Water only"
                                 ,"\\Nereo Nereo"
                                 ,"\\Nereo Nereo"
                                 ,"\\Mast Mast")
                               , c("\\wateronly Water only"
                                   , "\\Nereo Nereo"
                                   ,"\\Mast Mast"
                                   ,"\\NereoMast Nereo + Mast"
                                   , "\\Nereo Nereo"
                                   , "\\Mast Mast"
                                   , "\\NereoMast Nereo + Mast"
                                   , "\\Mast Mast"
                                   , "\\NereoMast Nereo + Mast"
                                   , "\\NereoMast Nereo + Mast")
                               , allPvalues.ExNWater)
  
  colnames(allPvalues.ExNWater) <- c("Group 1","Group 2","p","  ","FDR adj. p")

  capture.output(print(xtable(allPvalues.ExNWater, digits = 3), include.rownames = FALSE), file = paste0("ALPHAPLOTS/allPvalues.ExNWater.",tempI,".txt"))
  
  ttest.ExNWatervsEverythingelse <- t.test(MF.ExNWater.ExNWater[,paste0(i)], c(MF.ExNWater.NereoWater[,paste0(i)],MF.ExNWater.MastWater[,paste0(i)],MF.ExNWater.NereoMastWater[,paste0(i)]))
  capture.output(ttest.ExNWatervsEverythingelse, file = paste0("ALPHAPLOTS/ttest.ExNWatervsEverythingelse",tempI, ".txt"))
  
  capture.output(table(MF.ExNWater$ColRep), file = "ALPHAPLOTS/water_reps.txt")
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
  
  
  rownames(MF.ExNWater.Alpha) <- c(paste0("Water only (",ncounts[1],")")
                                   ,paste0("NMF Alone (",ncounts[2],")")
                                   ,paste0("With Nereo (",ncounts[3],")")
                                   ,paste0("With Mast (",ncounts[4],")")
                                   ,paste0("With Both (",ncounts[5],")")
  )
  
  pdf(paste0("ALPHAPLOTS/Alpha_div_watersamples_",tempI,".pdf"))
  par(mar = c(10,4,4,4))
  boxplot(t(MF.ExNWater.Alpha)
          , las = 2
          , col = c("blue","gray","green","purple","brown")
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
    ptemp <- get(test)$p.value
    ttemp <- round(get(test)$statistic,2)
    dftemp <- round(get(test)$parameter,2)
    toPaste <- paste0("(t=",ttemp,", df=",dftemp,")")
    allPvalues.LoneIncube <- rbind(allPvalues.LoneIncube, c(ptemp, toPaste))
  }
  allPvalues.LoneIncube <- cbind(signif(as.numeric(allPvalues.LoneIncube[,1]),3), allPvalues.LoneIncube[,2], signif(p.adjust(allPvalues.LoneIncube[,1], method = "fdr", n = 6),3))
  allPvalues.LoneIncube <- cbind(c("\\Nereo Nereo"
                                   , "\\Nereo Nereo"
                                   , "\\Nereo Nereo"
                                   , "\\Mast Mast"
                                   , "\\Mast Mast"
                                   , "\\Nereowater Nereo Water")
                                 , c("\\Mast Mast"
                                     , "\\Nereowater Nereo Water"
                                     , "\\Mastwater Mast Water"
                                     , "\\Nereowater Nereo Water"
                                     , "\\Mastwater Mast Water"
                                     , "\\Mastwater Mast Water")
                                 , allPvalues.LoneIncube)
  
  colnames(allPvalues.LoneIncube) <- c("Group 1","Group 2","p","  ","FDR adj. p")
  
  capture.output(print(xtable(allPvalues.LoneIncube, digits = 3), include.rownames = FALSE), file = paste0("ALPHAPLOTS/allPvalues.LoneIncube.",tempI,".txt"))
  
  capture.output(table(MF.LoneIncube$ColRep), file = "ALPHAPLOTS/loneincube_reps.txt")
  ###### PLOTTING #########
  MF.LoneIncube.temp <- MF.LoneIncube[,c(paste0(i), "ColRep","Replicate")]
  MF.LoneIncube.Alpha <- reshape(MF.LoneIncube.temp, idvar = "ColRep", timevar = "Replicate", direction = "wide")
  MF.LoneIncube.Alpha <- data.frame(MF.LoneIncube.Alpha
                                  , row.names = 1
  )
  MF.LoneIncube.Alpha <- MF.LoneIncube.Alpha[sapply(c("NereoNereo"
                                                  ,"MastMast"
                                                  ,"Nereowater"
                                                  ,"Mastwater"), function(x) {
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
          , col = c("green","purple","lightseagreen","lightslateblue")
          , ylab = paste0("Alpha Diversity (", i,")")
          , main = "Alpha diversity"
  )
  title(xlab = "Treatment", line = 8)
  dev.off()
  
  
}

######### **NMF vs Nereo and Mast** ###########
# levels(factor(MF.algae$ColRep))
MF.algae <- MF[-grep("Starfish", MF$ColRep),]

sortedGroups <- c("NereotestExNExN"
                  , "NereotestNereoExN"
                  , "NereotestMastExN"
                  , "NereotestNereoMastExN"
                  , "EnvironmentalBrocktonYoungNereo"
                  , "LoneincubeNereoNereo"
                  , "EnvironmentalBrocktonOldNereo"
                  , "LoneincubeMastMast"
                  , "EnvironmentalBrocktonMast"
                  , "NereotestH2OWater"
                  , "NereotestExNWater"
                  , "NereotestNereoWater"
                  , "NereotestMastWater"
                  , "NereotestNereoMastWater"
                  , "LoneincubeNereowater"
                  , "LoneincubeMastwater"
                  )
newFactor <- c("Nereo"
                , "Nereo"
                , "Nereo"
                , "Nereo"
                , "Nereo"
                , "Nereo"
                , "Nereo"
                , "Mast"
                , "Mast"
                , "Water"
                , "Water"
                , "Water"
                , "Water"
                , "Water"
                , "Water"
                , "Water"
)
orderList <- c()
for (group in sortedGroups) {
  orderList <- c(orderList, grep(paste0(group), MF.algae$ColRep))
}
MF.algae <- MF.algae[orderList,]
MF.algae$ColRep <- factor(MF.algae$ColRep, levels = sortedGroups)
MF.algae$newFactor <- newFactor[factor(MF.algae$ColRep)]
capture.output(table(MF.algae$newFactor), file = "ALPHAPLOTS/algaecompare.reps.txt")

for (i in alphaList) {
  tempI <- gsub("_even_1000_alpha", "", i)
  MF.algae.filt <- MF.algae[,c(i, "ColRep", "newFactor")]
  
  ########### STATS ###########
  MFnames <- c()
  for (group in levels(factor(newFactor))) {
    assign(paste0("MF.algae.",group), MF.algae[grep(group, MF.algae$newFactor),])
    MFnames <- c(MFnames, paste0("MF.algae.",group))
  }
  
  TTEST.Nereo.Mast <- t.test(get(MFnames[1])[,i],get(MFnames[2])[,i] )
  TTEST.Nereo.Water <- t.test(get(MFnames[1])[,i],get(MFnames[3])[,i] )
  TTEST.Mast.Water <- t.test(get(MFnames[2])[,i],get(MFnames[3])[,i] )

  
  allTEST <- c("TTEST.Nereo.Mast"
  , "TTEST.Nereo.Water"
  , "TTEST.Mast.Water")
  
  allPvalues.algae <- c()
  for (test in allTEST) {
    ptemp <- get(test)$p.value
    ttemp <- round(get(test)$statistic,2)
    dftemp <- round(get(test)$parameter,2)
    toPaste <- paste0("(t=",ttemp,", df=",dftemp,")")
      
    allPvalues.algae <- rbind(allPvalues.algae, c(ptemp, toPaste))
  }
  allPvalues.algae <- cbind(signif(as.numeric(allPvalues.algae[,1]),3), allPvalues.algae[,2], signif(p.adjust(allPvalues.algae[,1], method = "fdr", n = 3),3))
  allPvalues.algae <- cbind(c("Nereo"
                              , "Nereo"
                              , "Mast")
                            , c("Mast"
                                , "Water"
                                , "Water")
                            , allPvalues.algae)
  colnames(allPvalues.algae) <- c("Group 1","Group 2","p","  ","FDR adj. p")
  
  capture.output(print(xtable(allPvalues.algae, digits = 3), include.rownames = FALSE), file = paste0("ALPHAPLOTS/allPvalues.algae.",tempI,".txt"))
   
  
  ########### PLOT ############
  pdf(paste0("ALPHAPLOTS/Alpha_algaecompare",tempI, ".pdf"), pointsize = 14)
  par(mar = c(8,5,4,4))
  plot(MF.algae.filt[,paste0(i)]~ factor(MF.algae.filt[,"ColRep"])
       , ylab = paste0("Alpha Diversity (",tempI,")")
       , xlab = NA
       , col = c(rep("green", 4),"yellowgreen","darkgreen","darkolivegreen4","purple","magenta","lightblue",rep("blue",4), rep("dodgerblue",2))
       , las = 2 
       , xaxt = "n"
       , at = c(1,2,3,4,5,6,7,11,12,16,17,18,19,20,21,22))
  axis(side = 1
       , labels = c(expression("Nereo"),expression("Mast"), "Water")
       , at = c(4,11,19)
       , las = 1
       , tick = FALSE
  )
  title(xlab = "Sample Type"
        , line = 6)
  dev.off()
}

######### **Mast vs Nereo incubated stuff** ###########

Mast.only <- MF.ExNWater[grep("Mast", MF.ExNWater$ColRep),]
notMast <- MF.ExNWater[-grep("Mast", MF.ExNWater$ColRep),]

MastorNotMast <- matrix(nrow = 3, ncol = 1)
colnames(MastorNotMast) <- c("Welch's t-Test: Mast or no Mast")
rownames(MastorNotMast) <- c("1","2","3")
for (i in 1:length(alphaList)) {
  tempI <- gsub("_even_1000_alpha","",alphaList[i])
  tempttest <- t.test(Mast.only[,alphaList[i]], notMast[,alphaList[i]])
  ptemp <- signif(tempttest$p.value,3)
  ttemp <- round(tempttest$statistic,3)
  dftemp <- round(tempttest$parameter,3)
  toPaste <- paste0(ptemp," (t=",ttemp,", df=",dftemp,")")
  
  rownames(MastorNotMast)[i] <- tempI
  MastorNotMast[paste0(tempI),1] <- toPaste
}
capture.output(xtable(MastorNotMast), file = "ALPHAPLOTS/MastorNotMast_allmetrics.txt")



