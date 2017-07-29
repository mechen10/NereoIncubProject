#!/bin/Rscript
library(optparse)
################ OPT PARSE ####################

option_list = list(
  make_option(c("-m", "--mappingfile"), type="character",
              help="Mapping file with alpha diversity in it"),
  make_option(c("-a", "--alphaNames"),
              help="Comma separated list of alpha Names", type="character"),
  make_option(c("-b", "--betaNames"),
              help="Comma separated list of alpha Names", type="character"),
  make_option(c("-B", "--BCPWD"),
              help="Bray-Curtis distance matrix"),
  make_option(c("-U", "--UWUFPWD"),
              help="Unweighted Unifrac distance matrix"),
  make_option(c("-W","--WUFPWD"),
              help="Weighted Unifrac distance matrix")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

MFPWD = opt$mappingfile
alphaNamesTMP = opt$alphaNames
alphaList <- unlist(strsplit(alphaNamesTMP, ","))
betaNamesTMP = opt$betaNames
betaList <- unlist(strsplit(betaNamesTMP, split = ","))
BCPWD = opt$BCPWD
UWUFPWD = opt$UWUFPWD
WUFPWD = opt$WUFPWD


#####FORTESTING ########
# setwd("/Users/parfreylab/Desktop/lab_member_files/melissa/ForBotanyCluster/zz_NEREOINCUBE_16may2017/1_analysis")
# MFPWD <- "./ANALYSIS_ALPHABETATAXA/OTU_Tables_and_MP/MF_withalpha.txt"
# alphaNames <- c("chao1","PD_whole_tree","observed_otus")
# 
# # setwd("/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/")
# # MFPWD <-"/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/ANALYSIS_ALPHABETATAXA/OTU_Tables_and_MP/MF_withalpha_2.txt"
# alphaNames <-  'chao1_even_1000_alpha,PD_whole_tree_even_1000_alpha,observed_otus_even_1000_alpha'
# alphaList <- unlist(strsplit(alphaNames, ","))
# betaNames <- 'BC,WUF,UWUF'
# betaList <- unlist(strsplit(betaNames,","))
# 
# BCPWD<- "./ANALYSIS_ALPHABETATAXA/beta_div/bray_curtis_dm.txt"
# WUFPWD <- "./ANALYSIS_ALPHABETATAXA/beta_div/unweighted_unifrac_dm.txt"
# UWUFPWD <- "./ANALYSIS_ALPHABETATAXA/beta_div/weighted_unifrac_dm.txt"
# # MFPWD <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/OTU_MP_filt/MF_nochlpmito_m1000.txt"


############## LOAD DATA ########
library("car")
library("xtable")
library("MASS")
library("vegan")
library("stats")
# library("multcomp")
set.seed(3)
# Alpha div script

system("mkdir ALPHAPLOTS")
system("mkdir BETAPLOTS")
system("mkdir LATEX_outputs")

dm.UWUF <- read.delim(paste0(UWUFPWD)
                      , header = TRUE
                      , row.names = 1
                      , stringsAsFactors = FALSE)

dm.WUF <- read.delim(paste0(WUFPWD)
                     , header = TRUE
                     , row.names = 1
                     , stringsAsFactors = FALSE)
dm.BC <- read.delim(paste0(BCPWD)
                    , header = TRUE
                    , row.names = 1
                    , stringsAsFactors = FALSE)

MF <- read.delim(paste0(MFPWD)
                 , header = TRUE
                 , row.names = 1
                 , stringsAsFactors = FALSE
                 , na.strings = c("NA","N/A","na","<NA>"))

MF <- MF[sapply(rownames(dm.BC), function(x) {
  grep(paste0("^",x,"$"), rownames(MF))
}),]
MF.filtered <- MF[,grep("chao1|PD_whole_tree|observed_otus|ColRep|Replicate|Treatment|Subtrate", colnames(MF))]

# Specific tests
MF.ExN <- MF.filtered[grep("ExN.Nereotest.", rownames(MF.filtered)),]
MF.ExNWater <- MF.filtered[grep("Water.Nereotest.", rownames(MF.filtered)),]
MF.LoneIncube <- MF.filtered[grep("Loneincube", rownames(MF.filtered)),]


####### PLOT ##########

# For MF.ExN
########## ***ExN*** #############
######## *ALPHA* ##########
alphaListFiles.ExN <- list()
for (i in alphaList) {
  tempI <- gsub("_even_1000_alpha","", i)
  alphaListFiles.ExN[[paste0(tempI)]] <- c()
  ######## STATS ########

  ExN.lm <- lm(MF.ExN[,paste0(i)] ~ MF.ExN$ColRep)
  anova.ExN.lm <- Anova(ExN.lm, type = 'III') # OVERALL
  assign(paste0("ANOVA.ExN.overall.",tempI,".all"), anova.ExN.lm)
  alphaListFiles.ExN[[paste0(tempI)]] <- c(alphaListFiles.ExN[[paste0(tempI)]], paste0("ANOVA.ExN.overall.",tempI,".all"))
  
  alltreatments <- c("ExN","Nereo","Mast","NereoMast")
  for (AT in alltreatments) {
    assign(paste0("MF.ExN.",AT,"ExN"), MF.ExN[grep(paste0(AT), MF.ExN$ColRep),])
  }
  
  allcomparisons <- cbind(c("ExN","ExN","ExN","Nereo","Nereo","Mast")
                          , c("Nereo","Mast","NereoMast","Mast","NereoMast","NereoMast"))
  INDIVTtest <- c()
  for (AC in 1:nrow(allcomparisons)) {
    group1 <- allcomparisons[AC,1]
    group2 <- allcomparisons[AC,2]
    assign(paste0("TTEST.ExN.",group1,".",group2,".",tempI), t.test(get(paste0("MF.ExN.",group1,"ExN"))[,paste0(i)],get(paste0("MF.ExN.",group2,"ExN"))[,paste0(i)]))
    alphaListFiles.ExN[[paste0(tempI)]] <- c(alphaListFiles.ExN[[paste0(tempI)]], paste0("TTEST.ExN.",group1,".",group2,".",tempI))
    INDIVTtest <- c(INDIVTtest, paste0("TTEST.ExN.",group1,".",group2,".",tempI))
  }
  
  allPvalues.ExN <- c()
  for (test in INDIVTtest) {
    ptemp <- get(test)$p.value
    ttemp <- round(get(test)$statistic,2)
    dftemp <- round(get(test)$parameter, 2)
    toPaste <- paste0("(t=",ttemp,", df=",dftemp,")")
    allPvalues.ExN <- rbind(allPvalues.ExN, c(ptemp, toPaste))
  }
  allPvalues.ExN <- cbind(signif(as.numeric(allPvalues.ExN[,1]),3)
                          , allPvalues.ExN[,2]
                          , c(signif(p.adjust(allPvalues.ExN[1:3,1], method = "fdr", n = 3),3)
                              ,signif(p.adjust(allPvalues.ExN[4:6,1], method = "fdr", n = 3),3) ))
  
  rownames(allPvalues.ExN) <- sapply(1:nrow(allcomparisons), function(x) paste0(allcomparisons[x,1],"vs",paste0(allcomparisons[x,2])) )
  colnames(allPvalues.ExN) <- c("p","  ","FDR adj. p")
  assign(paste0("allPvalues.ExN.",tempI), allPvalues.ExN)
  
  alphaListFiles.ExN[[paste0(tempI)]] <- c( alphaListFiles.ExN[[paste0(tempI)]], paste0("allPvalues.ExN.",tempI))
  
  #### ExN vs Everything else ##
  
  assign(paste0("TTEST.ExNvsEverythingelse.",tempI), t.test(MF.ExN.ExNExN[,paste0(i)], c(MF.ExN.NereoExN[,paste0(i)],MF.ExN.MastExN[,paste0(i)],MF.ExN.NereoMastExN[,paste0(i)])))
  alphaListFiles.ExN[[paste0(tempI)]] <- c(alphaListFiles.ExN[[paste0(tempI)]], paste0("TTEST.ExNvsEverythingelse.",tempI))
  
  system(paste0("mkdir ALPHAPLOTS/",tempI,"/"))
  
  for (a in alphaListFiles.ExN[[paste0(tempI)]]) {
    if (length(grep("ANOVA|TTEST", a)) > 0) {
      capture.output(get(a), file = paste0("ALPHAPLOTS/",tempI,"/",a,".txt"))
    }
  }
  
 
  
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
  
  pdf(paste0("ALPHAPLOTS/",tempI,"/Alpha_div_meristemswabs_",tempI,".pdf"))
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

######## *BETA* ##########

betaListFiles.ExN <- list()
for (b in betaList) {
  # Set metric
  metric <- b
  betaListFiles.ExN[[paste0(metric)]] <- c()
  dm.temp <- get(paste0("dm.",metric))
  
  # Make dm.filtered and NMDS
  assign(paste0("dm.",metric,".ExN"), dm.temp[unlist(sapply(rownames(MF.ExN), function(x) {
    grep(x, rownames(dm.temp))})),unlist(sapply(rownames(MF.ExN), function(x) {
      grep(x, colnames(dm.temp))}))])
  betaListFiles.ExN[[paste0(metric)]] <- c( betaListFiles.ExN[[paste0(metric)]], paste0("dm.",metric,".ExN"))

  assign(paste0("NMDS.",metric,".ExN")
         , isoMDS(as.matrix(get(paste0("dm.",metric,".ExN"))), y = cmdscale(as.matrix(get(paste0("dm.",metric,".ExN"))), 2)))
  betaListFiles.ExN[[paste0(metric)]] <- c( betaListFiles.ExN[[paste0(metric)]], paste0("NMDS.",metric,".ExN"))
  
  ###### STATS ##########
  # ExN
  assign(paste0("PERMANOVA.overall.",metric,".ExN")
         , adonis(get(paste0("dm.",metric,".ExN")) ~ ColRep, data = MF.ExN))
  betaListFiles.ExN[[paste0(metric)]] <- c( betaListFiles.ExN[[paste0(metric)]], paste0("PERMANOVA.overall.",metric,".ExN"))
  
  assign(paste0("betadisp.overall",metric,".ExN")
         , betadisper(dist(get(paste0("dm.",metric,".ExN"))), group = MF.ExN$ColRep))
  assign(paste0("PERMDISP.overall.",metric,".ExN")
         , anova(get(paste0("betadisp.overall",metric,".ExN"))))
  betaListFiles.ExN[[paste0(metric)]] <- c( betaListFiles.ExN[[paste0(metric)]], paste0("PERMDISP.overall.",metric,".ExN"))
  
  # Treatment types
  treatTypes <- cbind(c("ExN.ExNvsNereo","ExN.ExNvsMast","ExN.ExNvsNereoMast","ExN.NereovsMast","ExN.NereovsNereoMast","ExN.MastvsNereoMast")
                      ,c("([.]ExN[.])|([.]Nereo[.])","([.]ExN[.])|([.]Mast[.])","([.]ExN[.])|([.]NereoMast[.])","([.]Nereo[.])|([.]Mast[.])","([.]Nereo[.])|([.]NereoMast[.])","([.]Mast[.])|([.]NereoMast[.])"))
  
  for (t in 1:nrow(treatTypes)) {
    assign(paste0("dm.",treatTypes[t,1]), get(paste0("dm.",metric,".ExN"))[grep(treatTypes[t,2], rownames(get(paste0("dm.",metric,".ExN")))), grep(treatTypes[t,2], colnames(get(paste0("dm.",metric,".ExN"))))])
    assign(paste0("MF.",treatTypes[t,1]), MF.ExN[grep(treatTypes[t,2], rownames(MF.ExN)),])
    assign(paste0("PERMANOVA.",metric,".",treatTypes[t,1]), adonis(get(paste0("dm.",treatTypes[t,1])) ~ ColRep, data = get(paste0("MF.",treatTypes[t,1]))))
    
    assign(paste0("betadisp.",metric,".",treatTypes[t,1]), betadisper(dist(get(paste0("dm.",treatTypes[t,1]))), group = get(paste0("MF.",treatTypes[t,1]))[,'ColRep']))
    assign(paste0("PERMDISP",".",treatTypes[t,1]), anova(get(paste0("betadisp.",metric,".",treatTypes[t,1]))))
    
    
    betaListFiles.ExN[[paste0(metric)]] <- c(betaListFiles.ExN[[paste0(metric)]]
                                         , paste0("dm.",treatTypes[t,1])
                                         , paste0("MF.",treatTypes[t,1])
                                         , paste0("PERMANOVA.",metric,".",treatTypes[t,1])
                                         , paste0("PERMDISP",".",treatTypes[t,1]))
    }
  
  IndivListTemp <- betaListFiles.ExN[[paste0(metric)]][grep(paste0("PERMANOVA.",metric), betaListFiles.ExN[[paste0(metric)]])]
  allPValues.ExN <- matrix(nrow = length(IndivListTemp), ncol = 2)
  rownames(allPValues.ExN) <- seq(1,length(IndivListTemp), by = 1)
  colnames(allPValues.ExN) <- c("p"," ")
  for (ILT in 1:length(IndivListTemp)) {
    newName <- gsub(paste0("PERMANOVA.",metric,".ExN."), "", IndivListTemp[ILT])
    rownames(allPValues.ExN)[ILT] <- newName
    allPValues.ExN[ILT,1] <- get(IndivListTemp[ILT])$aov.tab[6]$`Pr(>F)`[1]
    allPValues.ExN[ILT,2] <- paste0("(R^2=",round(get(IndivListTemp[ILT])$aov.tab[5]$R2[1],3),", F.model=",round(get(IndivListTemp[ILT])$aov.tab[4]$F.Model[1],3),", df="
                                  ,get(IndivListTemp[ILT])$aov.tab$Df[1],","
                                  ,get(IndivListTemp[ILT])$aov.tab$Df[3],")")
    
  }
  allPValues.ExN <- cbind(allPValues.ExN
                          , c(p.adjust(allPValues.ExN[1:3,1],method = "fdr", n = 3)
                              , p.adjust(allPValues.ExN[4:6,1],method = "fdr", n = 3)))
  colnames(allPValues.ExN) <- c("p"," ","fdr_adj")
  assign(paste0("allPValues.ExN.",metric), allPValues.ExN)
  
  betaListFiles.ExN[[paste0(metric)]] <- c(betaListFiles.ExN[[paste0(metric)]],paste0("allPValues.ExN.",metric))
  
  # ExN to others
  dm.ExN.ExNvsEverything <- get(paste0("dm.",metric,".ExN"))
  MF.ExN.ExNvsEverything <- MF.ExN
  MF.ExN.ExNvsEverything$EXNCOMPARE <- ""
  for (i in 1:length(MF.ExN.ExNvsEverything$ColRep)) {
    if (MF.ExN.ExNvsEverything[i,"ColRep"] != "NereotestExNExN") {
      MF.ExN.ExNvsEverything[i,"EXNCOMPARE"] <- "OTHER"
    } else {
      MF.ExN.ExNvsEverything[i,"EXNCOMPARE"] <- "ExN"
    }
  }
  assign(paste0("PERMANOVA.",metric,".ExN.ExNvsEverything"), adonis(dm.ExN.ExNvsEverything ~ EXNCOMPARE, data = MF.ExN.ExNvsEverything))
  assign(paste0("betadisp.",metric,".ExN.ExNvsEverything"), betadisper(dist(dm.ExN.ExNvsEverything), group = MF.ExN.ExNvsEverything$EXNCOMPARE))
  assign(paste0("PERMDISP.",metric,".ExN.ExNvsEverything"), anova(get(paste0("betadisp.",metric,".ExN.ExNvsEverything"))))
  
  betaListFiles.ExN[[paste0(metric)]] <- c(betaListFiles.ExN[[paste0(metric)]],c(paste0("PERMANOVA.",metric,".ExN.ExNvsEverything"),paste0("PERMDISP.",metric,".ExN.ExNvsEverything")))
  
  # PERMANOVA of 'OTHER'
  MF.ExN.EverythingElse <- MF.ExN.ExNvsEverything[MF.ExN.ExNvsEverything$EXNCOMPARE == "OTHER",]
  dm.ExN.EverythingElse <- dm.ExN.ExNvsEverything[match(rownames(MF.ExN.EverythingElse),rownames(dm.ExN.ExNvsEverything)),match(rownames(MF.ExN.EverythingElse),colnames(dm.ExN.ExNvsEverything))]
  assign(paste0("PERMANOVA.ExN.everythingelse.",metric), adonis(dm.ExN.EverythingElse ~ Treatment,data = MF.ExN.EverythingElse))
  assign(paste0("PERMDISP.ExN.everythingelse.",metric), anova(betadisper(dist(dm.ExN.EverythingElse), group = MF.ExN.EverythingElse$Treatment)))
  
  betaListFiles.ExN[[paste0(metric)]] <- c(betaListFiles.ExN[[paste0(metric)]],c(paste0("PERMANOVA.ExN.everythingelse.",metric),paste0("PERMDISP.ExN.everythingelse.",metric)))
  
  
  # Now print everything
  system(paste0("mkdir ./BETAPLOTS/",metric,"/"))
  
  for (n in betaListFiles.ExN[[paste0(metric)]]) {
    if (length(grep("PERMANOVA|PERMDISP", paste0(n))) == 1) {
      capture.output(get(paste0(n)), file = paste0("./BETAPLOTS/",metric,"/",n,".txt"))
      
    }
  }
  
  ####### PLOT ############
  # EXN UWUF
  MF.ExN$ColRep <- factor(MF.ExN$ColRep, levels = c('NereotestExNExN','NereotestNereoExN','NereotestMastExN','NereotestNereoMastExN'))
  ExNColours <- c("darkgrey","green","purple","brown")
  
  # Make chulls
  listChulls <- c("ExN","Nereo","Mast","NereoMast")
  plot.listChulls <- list()
  for ( LCH in listChulls) {
    assign(paste0("NMDS.",metric,".ExN.",LCH), get(paste0("NMDS.",metric,".ExN"))$points[grep(paste0(".",LCH,"."), rownames(get(paste0("NMDS.",metric,".ExN"))$points), fixed = TRUE),])
    assign(paste0("NMDS.",metric,".ExN.", LCH,".chull"), chull(get(paste0("NMDS.",metric,".ExN.",LCH))))
    assign(paste0("NMDS.",metric,".ExN.", LCH,".chull"), c(get(paste0("NMDS.",metric,".ExN.", LCH,".chull")), get(paste0("NMDS.",metric,".ExN.", LCH,".chull"))[1]))
    plot.listChulls[[LCH]]<- list(NMDS = get(paste0("NMDS.",metric,".ExN.",LCH))
                                  , chulls = get(paste0("NMDS.",metric,".ExN.", LCH,".chull")))
  }

  pdf(paste0("./BETAPLOTS/",metric,"/NMDS_",metric,"_ExN.pdf"), width = 10, height = 7, pointsize = 14)
  par(fig = c(0,0.7,0,1))
  plot(get(paste0("NMDS.",metric,".ExN"))$points
       , main = "NMDS of Nereo Meristem Swabs"
       , pch = 19
       , col = ExNColours[factor(MF.ExN$ColRep)]
       , sub = paste0("Stress: ",round(get(paste0("NMDS.",metric,".ExN"))$stress/100,2))
       , xlab = "NMDS 1"
       , ylab = "NMDS 2"
       , cex = 2
  )
  for (ch in 1:length(plot.listChulls)) {
    tempCH <- names(plot.listChulls)[ch]
    lines(plot.listChulls[[tempCH]]$NMDS[plot.listChulls[[tempCH]]$chulls,]
          , col = ExNColours[ch])
  }
  
  par(fig = c(0.6,1,0,1), new = TRUE)
  plot(0,0
       , pch = ""
       , xaxt = "n"
       , yaxt = "n"
       , xlab = ""
       , ylab = ""
       , bty = "n")
  legend("top"
         , pch = 19
         , legend = c("NMF Alone","with Nereo","with Mast","with Nereo + Mast")#levels(MF.ExN$ColRep)
         , col = ExNColours
         , cex = 1)
  dev.off()

}


############## ****ExNWater*** ############
#### *ALPHA * #####
# For MF.ExNWater
alphaListFiles.ExNWater <- list()
for (i in alphaList) {
  tempI <- gsub("_even_1000_alpha", "", i)
  alphaListFiles.ExNWater[[paste0(tempI)]] <- c()
  
  ######## STATS ########
  ExNWater.lm <- lm(MF.ExNWater[,paste0(i)] ~ MF.ExNWater$ColRep)
  anova.ExNWater.lm <- Anova(ExNWater.lm, type = 'III') # OVERALL
  assign(paste0("ANOVA.ExNWater.overall.",tempI,".all"), anova.ExNWater.lm)
  alphaListFiles.ExNWater[[paste0(tempI)]] <- c(alphaListFiles.ExNWater[[paste0(tempI)]], paste0("ANOVA.ExNWater.overall.",tempI,".all"))
  
  alltreatments <- c("H2O","ExN","Nereo","Mast","NereoMast")
  for (AT in alltreatments) {
    assign(paste0("MF.ExNWater.",AT,"Water"), MF.ExNWater[grep(paste0(AT), MF.ExNWater$ColRep),])
  }
  
  allcomparisons <- cbind(c("ExN","ExN","ExN","ExN","Nereo","Nereo","Mast")
                          , c("H2O","Nereo","Mast","NereoMast","Mast","NereoMast","NereoMast"))
  INDIVTtest <- c()
  for (AC in 1:nrow(allcomparisons)) {
    group1 <- allcomparisons[AC,1]
    group2 <- allcomparisons[AC,2]
    assign(paste0("TTEST.ExNWater.",group1,".",group2,".",tempI), t.test(get(paste0("MF.ExNWater.",group1,"Water"))[,paste0(i)],get(paste0("MF.ExNWater.",group2,"Water"))[,paste0(i)]))
    alphaListFiles.ExNWater[[paste0(tempI)]] <- c(alphaListFiles.ExNWater[[paste0(tempI)]], paste0("TTEST.ExNWater.",group1,".",group2,".",tempI))
    INDIVTtest <- c(INDIVTtest, paste0("TTEST.ExNWater.",group1,".",group2,".",tempI))
  }
  
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
                                    , signif(p.adjust(allPvalues.ExNWater[5:7,1], method = "fdr", n = 3),3))
  )
  
  rownames(allPvalues.ExNWater) <- sapply(1:nrow(allcomparisons), function(x) paste0(allcomparisons[x,1],"vs",paste0(allcomparisons[x,2])) )
  colnames(allPvalues.ExNWater) <- c("p","  ","FDR adj. p")
  assign(paste0("allPvalues.ExNWater.",tempI), allPvalues.ExNWater)
  alphaListFiles.ExNWater[[tempI]] <- c(alphaListFiles.ExNWater[[tempI]], paste0("allPvalues.ExNWater.",tempI))
  
  #### ExNWater vs Everything else ##
  
  assign(paste0("TTEST.ExNWatervsEverythingelse.",tempI), t.test(MF.ExNWater.ExNWater[,paste0(i)], c(MF.ExNWater.NereoWater[,paste0(i)],MF.ExNWater.MastWater[,paste0(i)],MF.ExNWater.NereoMastWater[,paste0(i)])))
  alphaListFiles.ExNWater[[paste0(tempI)]] <- c(alphaListFiles.ExNWater[[paste0(tempI)]], paste0("TTEST.ExNWatervsEverythingelse.",tempI))
  
  for (a in alphaListFiles.ExNWater[[tempI]]) {
    if (length(grep("ANOVA|TTEST", a)) > 0) {
      capture.output(get(a), file = paste0("ALPHAPLOTS/",tempI,"/",a,".txt"))
    }
  }
  
  
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
  
  pdf(paste0("ALPHAPLOTS/",tempI,"/Alpha_div_watersamples_",tempI,".pdf"))
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

######## *BETA* ##########

betaListFiles.ExNWater <- list()
for (b in betaList) {
  # Set metric
  metric <- b
  betaListFiles.ExNWater[[paste0(metric)]] <- c()
  dm.temp <- get(paste0("dm.",metric))
  
  # Make dm.filtered and NMDS
  assign(paste0("dm.",metric,".ExNWater"), dm.temp[unlist(sapply(rownames(MF.ExNWater), function(x) {
    grep(x, rownames(dm.temp))})),unlist(sapply(rownames(MF.ExNWater), function(x) {
      grep(x, colnames(dm.temp))}))])
  betaListFiles.ExNWater[[paste0(metric)]] <- c( betaListFiles.ExNWater[[paste0(metric)]], paste0("dm.",metric,".ExNWater"))
  
  assign(paste0("NMDS.",metric,".ExNWater")
         , isoMDS(as.matrix(get(paste0("dm.",metric,".ExNWater"))), y = cmdscale(as.matrix(get(paste0("dm.",metric,".ExNWater"))), 2)))
  betaListFiles.ExNWater[[paste0(metric)]] <- c( betaListFiles.ExNWater[[paste0(metric)]], paste0("NMDS.",metric,".ExNWater"))
  
  ###### STATS ##########
  # ExNWater
  assign(paste0("PERMANOVA.overall.",metric,".ExNWater")
         , adonis(get(paste0("dm.",metric,".ExNWater")) ~ ColRep, data = MF.ExNWater))
  betaListFiles.ExNWater[[paste0(metric)]] <- c( betaListFiles.ExNWater[[paste0(metric)]], paste0("PERMANOVA.overall.",metric,".ExNWater"))
  
  assign(paste0("betadisp.overall",metric,".ExNWater")
         , betadisper(dist(get(paste0("dm.",metric,".ExNWater"))), group = MF.ExNWater$ColRep))
  assign(paste0("PERMDISP.overall.",metric,".ExNWater")
         , anova(get(paste0("betadisp.overall",metric,".ExNWater"))))
  betaListFiles.ExNWater[[paste0(metric)]] <- c( betaListFiles.ExNWater[[paste0(metric)]], paste0("PERMDISP.overall.",metric,".ExNWater"))
  
  # Treatment types
  treatTypes <- cbind(c("ExNWater.ExNvsH2O","ExNWater.ExNvsNereo","ExNWater.ExNvsMast","ExNWater.ExNvsNereoMast","ExNWater.NereovsMast","ExNWater.NereovsNereoMast","ExNWater.MastvsNereoMast")
                      ,c("([.]ExN[.])|([.]H2O[.])","([.]ExN[.])|([.]Nereo[.])","([.]ExN[.])|([.]Mast[.])","([.]ExN[.])|([.]NereoMast[.])","([.]Nereo[.])|([.]Mast[.])","([.]Nereo[.])|([.]NereoMast[.])","([.]Mast[.])|([.]NereoMast[.])"))
  
  for (t in 1:nrow(treatTypes)) {
    assign(paste0("dm.",treatTypes[t,1]), get(paste0("dm.",metric,".ExNWater"))[grep(treatTypes[t,2], rownames(get(paste0("dm.",metric,".ExNWater")))), grep(treatTypes[t,2], colnames(get(paste0("dm.",metric,".ExNWater"))))])
    assign(paste0("MF.",treatTypes[t,1]), MF.ExNWater[grep(treatTypes[t,2], rownames(MF.ExNWater)),])
    assign(paste0("PERMANOVA.",metric,".",treatTypes[t,1]), adonis(get(paste0("dm.",treatTypes[t,1])) ~ ColRep, data = get(paste0("MF.",treatTypes[t,1]))))
    betaListFiles.ExNWater[[paste0(metric)]] <- c(betaListFiles.ExNWater[[paste0(metric)]]
                                             , paste0("dm.",treatTypes[t,1])
                                             , paste0("MF.",treatTypes[t,1])
                                             , paste0("PERMANOVA.",metric,".",treatTypes[t,1]))
  }

  IndivListTemp <- betaListFiles.ExNWater[[paste0(metric)]][grep(paste0("PERMANOVA.",metric), betaListFiles.ExNWater[[paste0(metric)]])]
  allPValues.ExNWater <- matrix(nrow = length(IndivListTemp), ncol = 2)
  rownames(allPValues.ExNWater) <- seq(1,length(IndivListTemp), by = 1)
  colnames(allPValues.ExNWater) <- c("p"," ")
  for (ILT in 1:length(IndivListTemp)) {
    newName <- gsub(paste0("PERMANOVA.",metric,".ExNWater."), "", IndivListTemp[ILT])
    rownames(allPValues.ExNWater)[ILT] <- newName
    allPValues.ExNWater[ILT,1] <- get(IndivListTemp[ILT])$aov.tab[6]$`Pr(>F)`[1]
    allPValues.ExNWater[ILT,2] <- paste0("(R^2=",round(get(IndivListTemp[ILT])$aov.tab[5]$R2[1],3),", F.model=",round(get(IndivListTemp[ILT])$aov.tab[4]$F.Model[1],3),", df="
                                    ,get(IndivListTemp[ILT])$aov.tab$Df[1],","
                                    ,get(IndivListTemp[ILT])$aov.tab$Df[3],")")
    
  }
  allPValues.ExNWater <- cbind(allPValues.ExNWater
                          , c(" ",p.adjust(allPValues.ExNWater[2:4,1],method = "fdr", n = 3)
                              , p.adjust(allPValues.ExNWater[5:7,1],method = "fdr", n = 3)))
  colnames(allPValues.ExNWater) <- c("p"," ","fdr_adj")
  assign(paste0("allPValues.ExNWater.",metric), allPValues.ExNWater)
  
  betaListFiles.ExNWater[[paste0(metric)]] <- c(betaListFiles.ExNWater[[paste0(metric)]],paste0("allPValues.ExNWater.",metric))
  
  
  # ExN to others
  dm.ExNWater.ExNvsEverything <- get(paste0("dm.",metric,".ExNWater"))
  dm.ExNWater.ExNvsEverything <- dm.ExNWater.ExNvsEverything[-grep("H2O", rownames(dm.ExNWater.ExNvsEverything))
                                                             , -grep("H2O", colnames(dm.ExNWater.ExNvsEverything))]
  MF.ExNWater.ExNvsEverything <- MF.ExNWater
  MF.ExNWater.ExNvsEverything <- MF.ExNWater.ExNvsEverything[-grep("H2O", rownames(MF.ExNWater.ExNvsEverything)),]
  MF.ExNWater.ExNvsEverything$EXNCOMPARE <- ""
  for (i in 1:length(MF.ExNWater.ExNvsEverything$ColRep)) {
    if (MF.ExNWater.ExNvsEverything[i,"ColRep"] != "NereotestExNWater") {
      MF.ExNWater.ExNvsEverything[i,"EXNCOMPARE"] <- "OTHER"
    } else {
      MF.ExNWater.ExNvsEverything[i,"EXNCOMPARE"] <- "ExN"
    }
  }
  assign(paste0("PERMANOVA.",metric,".ExNWater.ExNvsEverything"), adonis(dm.ExNWater.ExNvsEverything ~ EXNCOMPARE, data = MF.ExNWater.ExNvsEverything))
  assign(paste0("betadisp.",metric,".ExNWater.ExNvsEverything"), betadisper(dist(dm.ExNWater.ExNvsEverything), group = MF.ExNWater.ExNvsEverything$EXNCOMPARE))
  assign(paste0("PERMDISP.",metric,".ExNWater.ExNvsEverything"), anova(get(paste0("betadisp.",metric,".ExNWater.ExNvsEverything"))))
  
  
  # PERMANOVA of 'OTHER'
  MF.ExNWater.EverythingElse <- MF.ExNWater.ExNvsEverything[MF.ExNWater.ExNvsEverything$EXNCOMPARE == "OTHER",]
  dm.ExNWater.EverythingElse <- dm.ExNWater.ExNvsEverything[match(rownames(MF.ExNWater.EverythingElse),rownames(dm.ExNWater.ExNvsEverything)),match(rownames(MF.ExNWater.EverythingElse),colnames(dm.ExNWater.ExNvsEverything))]
  assign(paste0("PERMANOVA.ExNWater.everythingelse.",metric), adonis(dm.ExNWater.EverythingElse ~ Treatment,data = MF.ExNWater.EverythingElse))
  assign(paste0("PERMDISP.ExNWater.everythingelse.",metric), anova(betadisper(dist(dm.ExNWater.EverythingElse), group = MF.ExNWater.EverythingElse$Treatment)))
  
  betaListFiles.ExNWater[[paste0(metric)]] <- c(betaListFiles.ExNWater[[paste0(metric)]],c(paste0("PERMANOVA.ExNWater.everythingelse.",metric),paste0("PERMDISP.ExNWater.everythingelse.",metric)))
  
  
  
  betaListFiles.ExNWater[[paste0(metric)]] <- c(betaListFiles.ExNWater[[paste0(metric)]],c(paste0("PERMANOVA.",metric,".ExNWater.ExNvsEverything"),paste0("PERMDISP.",metric,".ExN.ExNvsEverything")))
  
  
  for (n in betaListFiles.ExNWater[[paste0(metric)]]) {
    if (length(grep("PERMANOVA|PERMDISP", paste0(n))) == 1) {
      capture.output(get(paste0(n)), file = paste0("./BETAPLOTS/",metric,"/",n,".txt"))
      
    }
  }
  
  ####### PLOT ############
  # EXN 
  MF.ExNWater$ColRep <- factor(MF.ExNWater$ColRep, levels = c('NereotestH2OWater','NereotestExNWater','NereotestNereoWater','NereotestMastWater','NereotestNereoMastWater'))
  ExNWaterColours <- c("blue","darkgrey","green","purple","brown")
  
  # Make chulls
  listChulls <- c("H2O","ExN","Nereo","Mast","NereoMast")
  plot.listChulls <- list()
  for ( LCH in listChulls) {
    assign(paste0("NMDS.",metric,".ExNWater.",LCH), get(paste0("NMDS.",metric,".ExNWater"))$points[grep(paste0(".",LCH,"."), rownames(get(paste0("NMDS.",metric,".ExNWater"))$points), fixed = TRUE),])
    assign(paste0("NMDS.",metric,".ExNWater.", LCH,".chull"), chull(get(paste0("NMDS.",metric,".ExNWater.",LCH))))
    assign(paste0("NMDS.",metric,".ExNWater.", LCH,".chull"), c(get(paste0("NMDS.",metric,".ExNWater.", LCH,".chull")), get(paste0("NMDS.",metric,".ExNWater.", LCH,".chull"))[1]))
    plot.listChulls[[LCH]]<- list(NMDS = get(paste0("NMDS.",metric,".ExNWater.",LCH))
                                  , chulls = get(paste0("NMDS.",metric,".ExNWater.", LCH,".chull")))
  }
  
  pdf(paste0("./BETAPLOTS/",metric,"/NMDS_",metric,"_ExNWater.pdf"), width = 10, height = 7, pointsize = 14)
  par(fig = c(0,0.7,0,1))
  plot(get(paste0("NMDS.",metric,".ExNWater"))$points
       , main = "NMDS of Water Samples"
       , pch = 19
       , col = ExNWaterColours[factor(MF.ExNWater$ColRep)]
       , sub = paste0("Stress: ",round(get(paste0("NMDS.",metric,".ExNWater"))$stress/100,2))
       , xlab = "NMDS 1"
       , ylab = "NMDS 2"
       , cex = 2
  )
  for (ch in 1:length(plot.listChulls)) {
    tempCH <- names(plot.listChulls)[ch]
    lines(plot.listChulls[[tempCH]]$NMDS[plot.listChulls[[tempCH]]$chulls,]
          , col = ExNWaterColours[ch])
  }
  
  par(fig = c(0.6,1,0,1), new = TRUE)
  plot(0,0
       , pch = ""
       , xaxt = "n"
       , yaxt = "n"
       , xlab = ""
       , ylab = ""
       , bty = "n")
  legend("top"
         , pch = 19
         , legend = c("Water Only","NMF Alone","with Nereo","with Mast","with Nereo + Mast")#levels(MF.ExNWater$ColRep)
         , col = ExNWaterColours
         , cex = 1)
  dev.off()
  
}

###########***LONEINCUBE***###########

###### *ALPHA* #####
# For MF.LoneIncube
alphaListFiles.LoneIncube <- list()
for (i in alphaList) {
  tempI <- gsub("_even_1000_alpha", "",i)
  alphaListFiles.LoneIncube[[paste0(tempI)]] <- c()
  
  ######## STATS ########
  LoneIncube.lm <- lm(MF.LoneIncube[,paste0(i)] ~ MF.LoneIncube$ColRep)
  anova.LoneIncube.lm <- Anova(LoneIncube.lm, type = 'III') # OVERALL
  assign(paste0("ANOVA.LoneIncube.overall.",tempI,".all"), anova.LoneIncube.lm)
  alphaListFiles.LoneIncube[[paste0(tempI)]] <- c(alphaListFiles.LoneIncube[[paste0(tempI)]], paste0("ANOVA.LoneIncube.overall.",tempI,".all"))
  
  alltreatments <- c("NereoNereo","Nereowater","MastMast","Mastwater")
  for (AT in alltreatments) {
    assign(paste0("MF.LoneIncube.",AT), MF.LoneIncube[grep(paste0(AT), MF.LoneIncube$ColRep),])
  }
  
  
  allcomparisons <- cbind(c("NereoNereo","Nereowater","NereoNereo","NereoNereo","MastMast","MastMast")
                          , c("MastMast","Mastwater","Nereowater","Mastwater","Nereowater","Mastwater"))
  INDIVTtest <- c()
  for (AC in 1:nrow(allcomparisons)) {
    group1 <- allcomparisons[AC,1]
    group2 <- allcomparisons[AC,2]
    assign(paste0("TTEST.LoneIncube.",group1,".",group2,".",tempI), t.test(get(paste0("MF.LoneIncube.",group1))[,paste0(i)],get(paste0("MF.LoneIncube.",group2))[,paste0(i)]))
    alphaListFiles.LoneIncube[[paste0(tempI)]] <- c(alphaListFiles.LoneIncube[[paste0(tempI)]], paste0("TTEST.LoneIncube.",group1,".",group2,".",tempI))
    INDIVTtest <- c(INDIVTtest, paste0("TTEST.LoneIncube.",group1,".",group2,".",tempI))
  }
  
  
  allPvalues.LoneIncube <- c()
  for (test in INDIVTtest) {
    ptemp <- get(test)$p.value
    ttemp <- round(get(test)$statistic,2)
    dftemp <- round(get(test)$parameter,2)
    toPaste <- paste0("(t=",ttemp,", df=",dftemp,")")
    allPvalues.LoneIncube <- rbind(allPvalues.LoneIncube, c(ptemp, toPaste))
  }
  allPvalues.LoneIncube <- cbind(signif(as.numeric(allPvalues.LoneIncube[,1]),3), allPvalues.LoneIncube[,2], signif(p.adjust(allPvalues.LoneIncube[,1], method = "fdr", n = 6),3))
  
  rownames(allPvalues.LoneIncube) <- sapply(1:nrow(allcomparisons), function(x) paste0(allcomparisons[x,1],"vs",paste0(allcomparisons[x,2])) )
  colnames(allPvalues.LoneIncube) <- c("p","  ","FDR adj. p")
  
  assign(paste0("allPvalues.LoneIncube.",tempI), allPvalues.LoneIncube)

  alphaListFiles.LoneIncube[[paste0(tempI)]] <- c(alphaListFiles.LoneIncube[[paste0(tempI)]], paste0("allPvalues.LoneIncube.",tempI))
  
  
  for (a in alphaListFiles.LoneIncube[[paste0(tempI)]]) {
    if (length(grep("ANOVA|TTEST", a)) > 0) {
      capture.output(get(a), file = paste0("ALPHAPLOTS/",tempI,"/",a,".txt"))
    }
  }
  
  ###### PLOTTING #########
  MF.LoneIncube.temp <- MF.LoneIncube[,c(paste0(i), "ColRep","Replicate")]
  MF.LoneIncube.Alpha <- reshape(MF.LoneIncube.temp, idvar = "ColRep", timevar = "Replicate", direction = "wide")
  MF.LoneIncube.Alpha <- data.frame(MF.LoneIncube.Alpha
                                    , row.names = 1
  )
  MF.LoneIncube.Alpha <- MF.LoneIncube.Alpha[unlist(sapply(c("NereoNereo"
                                                      ,"MastMast"
                                                      ,"Nereowater"
                                                      ,"Mastwater"), function(x) {
                                                        grep(x, rownames(MF.LoneIncube.Alpha))
                                                      })),]
  
  
  ncounts <- c()
  for (ROW in 1:nrow(MF.LoneIncube.Alpha)) {
    ncounts <- c(ncounts, sum(!is.na(MF.LoneIncube.Alpha[ROW,])))
  }
  
  
  rownames(MF.LoneIncube.Alpha) <- c(paste0("Nereo Swab (",ncounts[1],")")
                                     ,paste0("Mast Swab (",ncounts[2],")")
                                     ,paste0("Nereo Water (",ncounts[3],")")
                                     ,paste0("Mast Water (",ncounts[4],")")
  )
  
  pdf(paste0("ALPHAPLOTS/",tempI,"/Alpha_div_loneincube_",i,".pdf"))
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

######## *BETA* ##########

betaListFiles.LoneIncube <- list()
for (b in betaList) {
  # Set metric
  metric <- b
  betaListFiles.LoneIncube[[paste0(metric)]] <- c()
  dm.temp <- get(paste0("dm.",metric))

  # Make dm.filtered and NMDS
  assign(paste0("dm.",metric,".LoneIncube"), dm.temp[unlist(sapply(rownames(MF.LoneIncube), function(x) {
    grep(x, rownames(dm.temp))})),unlist(sapply(rownames(MF.LoneIncube), function(x) {
      grep(x, colnames(dm.temp))}))])
  betaListFiles.LoneIncube[[paste0(metric)]] <- c( betaListFiles.LoneIncube[[paste0(metric)]], paste0("dm.",metric,".LoneIncube"))
  
  assign(paste0("NMDS.",metric,".LoneIncube")
         , isoMDS(as.matrix(get(paste0("dm.",metric,".LoneIncube"))), y = cmdscale(as.matrix(get(paste0("dm.",metric,".LoneIncube"))), 2)))
  betaListFiles.LoneIncube[[paste0(metric)]] <- c( betaListFiles.LoneIncube[[paste0(metric)]], paste0("NMDS.",metric,".LoneIncube"))
  
  ###### STATS ##########
  # LoneIncube
  assign(paste0("PERMANOVA.overall.",metric,".LoneIncube")
         , adonis(get(paste0("dm.",metric,".LoneIncube")) ~ ColRep, data = MF.LoneIncube))
  betaListFiles.LoneIncube[[paste0(metric)]] <- c( betaListFiles.LoneIncube[[paste0(metric)]], paste0("PERMANOVA.overall.",metric,".LoneIncube"))
  
  assign(paste0("betadisp.overall",metric,".LoneIncube")
         , betadisper(dist(get(paste0("dm.",metric,".LoneIncube"))), group = MF.LoneIncube$ColRep))
  assign(paste0("PERMANOVA.betadisp.overall.",metric,".LoneIncube")
         , anova(get(paste0("betadisp.overall",metric,".LoneIncube"))))
  betaListFiles.LoneIncube[[paste0(metric)]] <- c( betaListFiles.LoneIncube[[paste0(metric)]], paste0("Anova.betadisp.overall.",metric,".LoneIncube"))
  
  # Treatment types
  treatTypes <- cbind(c("LoneIncube.NereovsMast","LoneIncube.waters","LoneIncube.Nereovswater","LoneIncube.Mastvswater","LoneIncube.NereovsMwater","LoneIncube.MastvsNwater")
                      ,c("(^Nereo[.])|(^Mast[.])","^water[.]","[.]Nereo[.]","[.]Mast[.]","^Nereo[.]|^water[.].*[.]Mast[.]","^Mast[.]|^water[.].*[.]Nereo[.]"))
  
  
  for (t in 1:nrow(treatTypes)) {
    assign(paste0("dm.",treatTypes[t,1]), get(paste0("dm.",metric,".LoneIncube"))[grep(treatTypes[t,2], rownames(get(paste0("dm.",metric,".LoneIncube")))), grep(treatTypes[t,2], colnames(get(paste0("dm.",metric,".LoneIncube"))))])
    assign(paste0("MF.",treatTypes[t,1]), MF.LoneIncube[grep(treatTypes[t,2], rownames(MF.LoneIncube)),])
    assign(paste0("PERMANOVA.",metric,".",treatTypes[t,1]), adonis(get(paste0("dm.",treatTypes[t,1])) ~ ColRep, data = get(paste0("MF.",treatTypes[t,1]))))
    betaListFiles.LoneIncube[[paste0(metric)]] <- c(betaListFiles.LoneIncube[[paste0(metric)]]
                                                  , paste0("dm.",treatTypes[t,1])
                                                  , paste0("MF.",treatTypes[t,1])
                                                  , paste0("PERMANOVA.",metric,".",treatTypes[t,1]))
  }
  
  IndivListTemp <- betaListFiles.LoneIncube[[paste0(metric)]][grep(paste0("PERMANOVA.",metric), betaListFiles.LoneIncube[[paste0(metric)]])]
  allPValues.LoneIncube <- matrix(nrow = length(IndivListTemp), ncol = 2)
  rownames(allPValues.LoneIncube) <- seq(1,length(IndivListTemp), by = 1)
  colnames(allPValues.LoneIncube) <- c("p"," ")
  for (ILT in 1:length(IndivListTemp)) {
    newName <- gsub(paste0("PERMANOVA.",metric,".LoneIncube."), "", IndivListTemp[ILT])
    rownames(allPValues.LoneIncube)[ILT] <- newName
    allPValues.LoneIncube[ILT,1] <- get(IndivListTemp[ILT])$aov.tab[6]$`Pr(>F)`[1]
    allPValues.LoneIncube[ILT,2] <- paste0("(R^2=",round(get(IndivListTemp[ILT])$aov.tab[5]$R2[1],3),", F.model=",round(get(IndivListTemp[ILT])$aov.tab[4]$F.Model[1],3),", df="
                                         ,get(IndivListTemp[ILT])$aov.tab$Df[1],","
                                         ,get(IndivListTemp[ILT])$aov.tab$Df[3],")")
    
  }
  assign(paste0("allPValues.LoneIncube.",metric), allPValues.LoneIncube)
  
  betaListFiles.LoneIncube[[paste0(metric)]] <- c(betaListFiles.LoneIncube[[paste0(metric)]],paste0("allPValues.LoneIncube.",metric))
  
  for (n in betaListFiles.LoneIncube[[paste0(metric)]]) {
    if (length(grep("PERMANOVA", paste0(n))) == 1) {
      capture.output(get(paste0(n)), file = paste0("./BETAPLOTS/",metric,"/",n,".txt"))
      
    }
  }
  
  ####### PLOT ############
  # EXN UWUF
  MF.LoneIncube$ColRep <- factor(MF.LoneIncube$ColRep, levels = c('LoneincubeNereowater','LoneincubeNereoNereo','LoneincubeMastwater','LoneincubeMastMast'))
  LoneIncubeColours <- c("lightseagreen","green","lightslateblue","purple")
  
  # Make chulls
  listChulls <- c("water.Loneincube.Nereo","Nereo.Loneincube","water.Loneincube.Mast","Mast.Loneincube")
  plot.listChulls <- list()
  for ( LCH in listChulls) {
    assign(paste0("NMDS.",metric,".LoneIncube.",LCH), get(paste0("NMDS.",metric,".LoneIncube"))$points[grep(paste0(LCH), rownames(get(paste0("NMDS.",metric,".LoneIncube"))$points), fixed = TRUE),])
    assign(paste0("NMDS.",metric,".LoneIncube.", LCH,".chull"), chull(get(paste0("NMDS.",metric,".LoneIncube.",LCH))))
    assign(paste0("NMDS.",metric,".LoneIncube.", LCH,".chull"), c(get(paste0("NMDS.",metric,".LoneIncube.", LCH,".chull")), get(paste0("NMDS.",metric,".LoneIncube.", LCH,".chull"))[1]))
    plot.listChulls[[LCH]]<- list(NMDS = get(paste0("NMDS.",metric,".LoneIncube.",LCH))
                                  , chulls = get(paste0("NMDS.",metric,".LoneIncube.", LCH,".chull")))
  }
  
  pdf(paste0("./BETAPLOTS/",metric,"/NMDS_",metric,"_Loneincube.pdf"), width = 10, height = 7, pointsize = 14)
  par(fig = c(0,0.7,0,1))
  plot(get(paste0("NMDS.",metric,".LoneIncube"))$points
       , main = "NMDS of M-W Experiment"
       , pch = 19
       , col = LoneIncubeColours[factor(MF.LoneIncube$ColRep)]
       , sub = paste0("Stress: ",round(get(paste0("NMDS.",metric,".LoneIncube"))$stress/100,2))
       , xlab = "NMDS 1"
       , ylab = "NMDS 2"
       , cex = 2
  )
  for (ch in 1:length(plot.listChulls)) {
    tempCH <- names(plot.listChulls)[ch]
    lines(plot.listChulls[[tempCH]]$NMDS[plot.listChulls[[tempCH]]$chulls,]
          , col = LoneIncubeColours[ch])
  }
  
  par(fig = c(0.6,1,0,1), new = TRUE)
  plot(0,0
       , pch = ""
       , xaxt = "n"
       , yaxt = "n"
       , xlab = ""
       , ylab = ""
       , bty = "n")
  legend("top"
         , pch = 19
         , legend = c("Water-Nereo","Nereo","Water-Mast","Mast")#levels(MF.LoneIncube$ColRep)
         , col = LoneIncubeColours
         , cex = 1)
  dev.off()
  
}

######### **** ALL ALGAE **** ###########
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
MF.algae$newFactor <- factor(MF.algae$newFactor, levels = c("Nereo","Mast","Water"))
capture.output(table(MF.algae$newFactor), file = "ALPHAPLOTS/algaecompare.reps.txt")

#### *ALPHA* ####
alphaListFiles.All <- list()
for (i in alphaList) {
  tempI <- gsub("_even_1000_alpha", "", i)
  MF.algae.filt <- MF.algae[,c(i, "ColRep", "newFactor")]
  alphaListFiles.All[[paste0(tempI)]] <- c()
  
  ########### STATS ###########
  # Overall ANOVA
  assign(paste0("MF.algae.",tempI,".lm"), lm(get(i) ~ newFactor, data = MF.algae))
  assign(paste0("ANOVA.algae.",tempI), anova(get(paste0("MF.algae.",tempI,".lm"))))
  alphaListFiles.All[[paste0(tempI)]] <- c(alphaListFiles.All[[paste0(tempI)]], paste0("ANOVA.algae.",tempI))
  
  MFnames <- c()
  for (group in levels(factor(newFactor, levels = c("Nereo","Mast","Water")))) {
    assign(paste0("MF.algae.",tempI,"_",group), MF.algae[grep(group, MF.algae$newFactor),])
    MFnames <- c(MFnames, paste0("MF.algae.",tempI,"_",group))
  }
  
  allTEST <- c() 
  for (a in 1:(length(MFnames)-1)) {
    for (b in (a+1):length(MFnames)) {
      group1 <- gsub("MF.algae.observed_otus_","",MFnames[a])
      group2 <- gsub("MF.algae.observed_otus_","",MFnames[b])
      assign(paste0("TTEST.",group1,".",group2), t.test(get(MFnames[a])[,i], get(MFnames[b])[,i]))
      alphaListFiles.All[[paste0(tempI)]] <- c(alphaListFiles.All[[paste0(tempI)]], paste0("TTEST.",group1,".",group2))
      allTEST <- c(allTEST, paste0("TTEST.",group1,".",group2))
    }
  }
  
  allPvalues.algae <- c()
  for (test in allTEST) {
    ptemp <- get(test)$p.value
    ttemp <- round(get(test)$statistic,2)
    dftemp <- round(get(test)$parameter,2)
    toPaste <- paste0("(t=",ttemp,", df=",dftemp,")")
    
    allPvalues.algae <- rbind(allPvalues.algae, c(ptemp, toPaste))
  }
  allPvalues.algae <- cbind(signif(as.numeric(allPvalues.algae[,1]),3), allPvalues.algae[,2], signif(p.adjust(allPvalues.algae[,1], method = "fdr", n = 3),3))
  rownames(allPvalues.algae) <- gsub("[.]","vs",gsub("TTEST[.]","",allTEST))
  colnames(allPvalues.algae) <- c("p","  ","FDR adj. p")
  assign(paste0("allPvalues.algae.", tempI), allPvalues.algae)
  
  alphaListFiles.All[[paste0(tempI)]] <- c(alphaListFiles.All[[paste0(tempI)]], paste0("allPvalues.algae.", tempI))
  
  for (item in alphaListFiles.All[[tempI]]) {
    if (length(grep("ANOVA", item)) > 0) {
      capture.output(get(item), file = paste0("ALPHAPLOTS/",tempI,"/",item,".txt"))
    }
  }

  ########### PLOT ############
  pdf(paste0("ALPHAPLOTS/",tempI,"/Alpha_algaecompare",tempI, ".pdf"), pointsize = 14)
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

#### *BETA*#####
colorsPlot <- c( "green" # [12] "NereotestExNExN" 
                 , "green" # [17] "NereotestNereoExN"       
                 , "green" # [15] "NereotestMastExN"       
                 , "green" # [18] "NereotestNereoMastExN"  
                 , "yellowgreen" # [10] "LoneincubeNereoNereo" 
                 , "darkgreen" # [2] "EnvironmentalBrocktonOldNereo"        
                 , "darkolivegreen4" # [3] "EnvironmentalBrocktonYoungNereo"  
                 , "purple" # [8] "LoneincubeMastMast"        
                 , "magenta"# [1] "EnvironmentalBrocktonMast"            
                 , "lightblue" # [14] "NereotestH2OWater"                    
                 , "blue" # [13] "NereotestExNWater"                    
                 , "blue"# [20] "NereotestNereoWater" 
                 , "blue" # [16] "NereotestMastWater"                   
                 , "blue" # [19] "NereotestNereoMastWater"              
                 , "dodgerblue" # [11] "LoneincubeNereoWater"                 
                 , "dodgerblue" # [9] "LoneincubeMastWater"                  
)

betaListFiles.All <- list()
for (b in betaList) {
  metric <- b
  betaListFiles.All[[paste0(metric)]] <- c()
  
  assign(paste0("dm.",metric,".all"), get(paste0("dm.",metric))[-grep("Starfish", rownames(get(paste0("dm.",metric)))), -grep("Starfish", colnames(get(paste0("dm.",metric))))])
  assign(paste0("NMDS.",metric), isoMDS(as.matrix(get(paste0("dm.",metric,".all"))), y = cmdscale(as.matrix(get(paste0("dm.",metric,".all"))), 2)))
  # Reorder to make it correct order
  assign(paste0("NMDS.",metric,".points"), get(paste0("NMDS.",metric))$points[sapply(rownames(MF.algae), function(x) {
    grep(paste0("^",x,"$"), rownames(get(paste0("NMDS.",metric))$points))
  }),])
  
  ######### STATS ##############
  # Overall
  assign(paste0("PERMANOVA.algae.water.",metric), adonis(get(paste0("dm.",metric,".all")) ~ newFactor, data = MF.algae))

  assign(paste0("betadisp.",metric,".algae.water"), betadisper(dist(get(paste0("dm.",metric,".all"))), group = MF.algae$newFactor))
  assign(paste0("PERMDISP.",metric,".algae.water"), anova(get(paste0("betadisp.",metric,".algae.water"))))

  betaListFiles.All[[paste0(metric)]] <- c(betaListFiles.All[[paste0(metric)]], c(paste0("PERMANOVA.algae.water.",metric),paste0("PERMDISP.",metric,".algae.water")))  
  
  
  allPValues.algae.water <- matrix(nrow = 3, ncol = 2)
  rownames(allPValues.algae.water) <- c(1,2,3) 
  colnames(allPValues.algae.water) <- c("p"," ")
  count <- 0
  newRowNames <- c()
  for (g1 in 1:(length(unique(newFactor))-1)) {
    for (g2 in (g1+1):length(unique(newFactor))) {
      count <- count +1
      g1temp <- unique(newFactor)[g1]
      g2temp <- unique(newFactor)[g2]
      
      MF.temp <- MF.algae[grep(paste0("(",g1temp,"|",g2temp,")"), MF.algae$newFactor),]
      dm.temp <- get(paste0("dm.",metric,".all"))[sapply(rownames(MF.temp), function(x) grep(x, rownames(get(paste0("dm.",metric,".all")))))
                             , sapply(rownames(MF.temp), function(x) grep(x, colnames(get(paste0("dm.",metric,".all")))))]
      anova.temp <- adonis(dm.temp ~ newFactor, data = MF.temp)   
      ptemp <- anova.temp$aov.tab$`Pr(>F)`[1]
      rtemp <- anova.temp$aov.tab$R2[1]
      ftemp <- round(anova.temp$aov.tab$F.Model[1],3)
      dftemp <- paste0(anova.temp$aov.tab$Df[1],",",anova.temp$aov.tab$Df[3])
      toPaste <- paste0("(R^2=",round(rtemp,3),", F.model=",ftemp,", df=",dftemp,")")
      
      allPValues.algae.water[count,1] <- ptemp
      allPValues.algae.water[count,2] <- toPaste
      
      newRowNames <- rbind(newRowNames, c(g1temp, g2temp))
      assign(paste0("PERMANOVA.",metric,".",g1temp,".",g2temp), anova.temp)
      betaListFiles.All[[paste0(metric)]] <- c(betaListFiles.All[[paste0(metric)]], paste0("PERMANOVA.",metric,".",g1temp,".",g2temp))
    }
  }
  allPValues.algae.water <- cbind(signif(as.numeric(allPValues.algae.water[,1],3))
                                       , allPValues.algae.water[,2]
                                       , signif(p.adjust(allPValues.algae.water[,1], method = "fdr", n = 3),3)
  )
  rownames(allPValues.algae.water) <- paste0(newRowNames[,1],"vs",newRowNames[,2])
  colnames(allPValues.algae.water) <- c("p"," ","FDR adj. p")
  assign(paste0("allPValues.algae.water.",metric), allPValues.algae.water)
  betaListFiles.All[[paste0(metric)]] <- c(betaListFiles.All[[paste0(metric)]],paste0("allPValues.algae.water.",metric) )
  
  for (item in betaListFiles.All[[paste0(metric)]]) {
    if (length(grep("PERMANOVA|PERMDSIP", item)) > 0) {
      capture.output(get(item), file = paste0("BETAPLOTS/",metric,"/",item,".txt"))
    }
  }
  
  ###### PLOT ###########
  pdf(file = paste0("./BETAPLOTS/",metric,"/NMDS_all_",metric,".pdf"), pointsize = 14, width = 10, height = 7)
  par(fig = c(0,0.7,0,1))
  plot(get(paste0("NMDS.",metric,".points"))
       , main = paste0("NMDS plot of all samples (",metric,")")
       , sub = paste0("Stress: ",round(get(paste0("NMDS.",metric))$stress/100,2))
       , bg = colorsPlot[factor(MF.algae$ColRep)]
       , col = "black"
       , pch = 21
       # , pch = pchPlot[factor(MF.all$ColRep)]
       , cex = 2
       , xlab = "NMDS 1"
       , ylab = "NMDS 2"
  )
  par(fig = c(0.6,1,0,1), mar = c(0,0,0,0), new = TRUE)
  legend("left"
         , legend = c("Nereo Meristem (lab)"
                      , "Nereo Meristem (wild)"
                      , "Nereo blade (lab)"
                      , "Nereo blade (wild)"
                      , "Mast blade (lab)"
                      , "Mast blade (wild)"
                      , "Water alone (M-W-NMF)"
                      , "Water (M-W-NMF)"
                      , "Water (M-W)"
         )
         , pch = 21
         #, pch = c(19,18,17,18,19,18,8,11,11)
         , pt.bg = c("green", "yellowgreen", "darkgreen", "darkolivegreen4"
                     , "purple", "magenta"
                     , "lightblue", "blue", "dodgerblue" )
         , col = "black"
         , pt.cex = 2
         , cex = 1
         , y.intersp = 2
         , bty = "n")
  dev.off()
  
}




############ ****COMBINE ALL P-VALUES**** #############
# CHOSEN COMBO:
choseA <- "chao1"
choseB <- "BC"


#### For ExN ####
## Alpha ##
for (a in alphaList) {
  tempI <- gsub("_even_1000_alpha","", a)
  
  # Get labels
  factorsCompareExN <- cbind(c(" ","NMF Alone","  ","Nereo","Nereo","Mast")
                             , c("Nereo","Mast","NereoMast","Mast","NereoMast","NereoMast"))
  colnamesTop <- c("","","","Welch's t-Test","")
  colnamessubTop <- c("","","","(Richness)","")
  colnamesSecond <- c("Group 1","Group 2","p"," ","FDR adj. p")
  
  # Get ExN vs everything individual tests
  ExN.insertA <- alphaListFiles.ExN[[paste0(tempI)]][grep("ExNvsEverything", alphaListFiles.ExN[[paste0(tempI)]])]
  ptemp <- signif(get(ExN.insertA)$p.value,3)
  ttemp <- round(get(ExN.insertA)$statistic,3)
  dftemp <- round(get(ExN.insertA)$parameter,3)
  ExN.insertA.toPaste <- cbind(ptemp, paste0("(t=",ttemp,",df=",dftemp,")"), " ")
  
  # Merge everything 
  CombinedTable.ExN <- rbind(rep(" ", 3)
                             , ExN.insertA.toPaste
                             , rep(" ", 3)
                             , get(paste0("allPvalues.ExN.",tempI))[4:6,])
  
  CombinedTable.ExN <- cbind(factorsCompareExN, CombinedTable.ExN)
  CombinedTable.ExN <- rbind(colnamesTop, colnamessubTop, colnamesSecond, CombinedTable.ExN)
  
  # Add dummy row and col name so it doesn't freak out
  rownames(CombinedTable.ExN) <- seq(1,nrow(CombinedTable.ExN), by = 1)
  colnames(CombinedTable.ExN) <- seq(1, ncol(CombinedTable.ExN), by = 1)
  
  # Rename
  assign(paste0("CombinedTable.ExN.",tempI), CombinedTable.ExN)
  
  write.table(CombinedTable.ExN
              , sep = "\t"
              , row.names = FALSE
              , col.names = FALSE
              , file = paste0("LATEX_outputs/CombinedTable.ExN.",tempI,".txt"))
  
  CombinedTable.changed <- gsub("(t","($t$", CombinedTable.ExN , fixed = TRUE)
  CombinedTable.changed <- gsub("df","$df$", CombinedTable.changed , fixed = TRUE)
  CombinedTable.changed <- gsub("^p$","$p$", CombinedTable.changed , fixed = FALSE)
  CombinedTable.changed <- gsub("t-","$t$-", CombinedTable.changed , fixed = TRUE)
  CombinedTable.changed <- gsub("adj. p","adj. $p$", CombinedTable.changed , fixed = TRUE)
  
  
  
  
  capture.output(print(xtable(CombinedTable.changed)
                       , include.rownames= FALSE
                       , include.colnames = FALSE
                       , math.style.negative = TRUE
                       , math.style.exponents = TRUE
                       , sanitize.text.function = function(x) {x})
                 , file = paste0("LATEX_outputs/LATEX_CombinedTable.ExN.",tempI,".txt"))
  
}

## Beta ##
for (b in betaList){
  metric <- b
  
  # Get labels
  factorsCompareExN <- cbind(c(" ","NMF Alone","  ","Nereo","Nereo","Mast")
                             , c("Nereo","Mast","NereoMast","Mast","NereoMast","NereoMast"))
  colnamesTop <- c("","","","PERMANOVA","")
  colnamessubTop <- c("","","","(Community Dissimilarity)","")
  colnamesSecond <- c("Group 1","Group 2","p"," ","FDR adj. p")
  
  # Get ExN vs everything individual tests
  ExN.insertB <- betaListFiles.ExN[[paste0(metric)]][grep("PERMANOVA.*ExNvsEverything", betaListFiles.ExN[[paste0(metric)]])]
  ptemp <- signif(get(ExN.insertB)$aov.tab$`Pr(>F)`[1],3)
  rtemp <- round(get(ExN.insertB)$aov.tab$R2[1])
  ftemp <- round(get(ExN.insertB)$aov.tab$F.Model[1],3)
  dftemp <- paste0(round(get(ExN.insertB)$aov.tab$Df[1],3),",", round(get(ExN.insertB)$aov.tab$Df[3],3))
  ExN.insertB.toPaste <- cbind(ptemp, paste0("(R^2=",rtemp,",F.model=",ftemp,",df=",dftemp,")"), " ")
  
  # Merge everything 
  CombinedTable.ExN <- rbind(rep(" ", 3)
                             , ExN.insertB.toPaste
                             , rep(" ", 3)
                             , get(paste0("allPValues.ExN.",metric))[4:6,])
  
  CombinedTable.ExN <- cbind(factorsCompareExN, CombinedTable.ExN)
  CombinedTable.ExN <- rbind(colnamesTop,colnamessubTop, colnamesSecond, CombinedTable.ExN)
  
  # Add dummy row and col name so it doesn't freak out
  rownames(CombinedTable.ExN) <- seq(1,nrow(CombinedTable.ExN), by = 1)
  colnames(CombinedTable.ExN) <- seq(1, ncol(CombinedTable.ExN), by = 1)
  
  # Rename
  assign(paste0("CombinedTable.ExN.",metric), CombinedTable.ExN)
  
  write.table(CombinedTable.ExN
              , sep = "\t"
              , row.names = FALSE
              , col.names = FALSE
              , file = paste0("LATEX_outputs/CombinedTable.ExN.",metric,".txt"))
  
  CombinedTable.changed <- gsub("R^2","$R^2$", CombinedTable.ExN , fixed = TRUE)
  CombinedTable.changed <- gsub("df","$df$", CombinedTable.changed , fixed = TRUE)
  CombinedTable.changed <- gsub("^p$","$p$", CombinedTable.changed , fixed = FALSE)
  CombinedTable.changed <- gsub("adj. p","adj. $p$", CombinedTable.changed , fixed = TRUE)
  
  capture.output(print(xtable(CombinedTable.changed)
                       , include.rownames= FALSE
                       , include.colnames = FALSE
                       , math.style.negative = TRUE
                       , math.style.exponents = TRUE
                       , sanitize.text.function = function(x) {x})
                 , file = paste0("LATEX_outputs/LATEX_CombinedTable.ExN.",metric,".txt"))
  
  
  
  
}


#### For ExNWater ####
## Alpha##
for (a in alphaList) {
  tempI <- gsub("_even_1000_alpha", "",a)
  
  # Get labels
  factorsCompareExNWater <- cbind(c("Water Only", "","NMF Alone","  ","Nereo","Nereo","Mast")
                                  , c("NMF Alone","Nereo","Mast","NereoMast","Mast","NereoMast","NereoMast"))
  colnamesTop <- c("","","","Welch's t-Test","")
  colnamessubTop <- c("","","","(Richness)","")
  
  colnamesSecond <- c("Group 1","Group 2","p"," ","FDR adj. p")
  
  # Get ExNWater vs everything individual tests
  ExNWater.insertA <- alphaListFiles.ExNWater[[paste0(tempI)]][grep("ExNWatervsEverything", alphaListFiles.ExNWater[[paste0(tempI)]])]
  ptemp <- signif(get(ExNWater.insertA)$p.value,3)
  ttemp <- round(get(ExNWater.insertA)$statistic,3)
  dftemp <- round(get(ExNWater.insertA)$parameter,3)
  ExNWater.insertA.toPaste <- cbind(ptemp, paste0("(t=",ttemp,",df=",dftemp,")"), " ")
  
  # Merge everything 
  CombinedTable.ExNWater <- rbind( get(paste0("allPvalues.ExNWater.",tempI))[1:1,]
                                  , rep(" ", 3)
                                  , ExNWater.insertA.toPaste
                                  , rep(" ", 3)
                                  , get(paste0("allPvalues.ExNWater.",tempI))[4:6,] )
  
  CombinedTable.ExNWater <- cbind(factorsCompareExNWater, CombinedTable.ExNWater)
  CombinedTable.ExNWater <- rbind(colnamesTop, colnamessubTop, colnamesSecond, CombinedTable.ExNWater)
  
  # Add dummy row and col name so it doesn't freak out
  rownames(CombinedTable.ExNWater) <- seq(1,nrow(CombinedTable.ExNWater), by = 1)
  colnames(CombinedTable.ExNWater) <- seq(1, ncol(CombinedTable.ExNWater), by = 1)
  
  # Rename
  
  assign(paste0("CombinedTable.ExNWater.",tempI), CombinedTable.ExNWater)
  
  write.table(CombinedTable.ExNWater
              , sep = "\t"
              , row.names = FALSE
              , col.names = FALSE
              , file = paste0("LATEX_outputs/CombinedTable.ExNWater.",tempI,".txt"))
  
  CombinedTable.changed <- gsub("(t","($t$", CombinedTable.ExNWater , fixed = TRUE)
  CombinedTable.changed <- gsub("df","$df$", CombinedTable.changed , fixed = TRUE)
  CombinedTable.changed <- gsub("^p$","$p$ ", CombinedTable.changed , fixed = FALSE)
  CombinedTable.changed <- gsub("t-","$t$-", CombinedTable.changed , fixed = TRUE)
  CombinedTable.changed <- gsub("adj. p","adj. $p$", CombinedTable.changed , fixed = TRUE)
  
  capture.output(print(xtable(CombinedTable.changed)
                       , include.rownames= FALSE
                       , include.colnames = FALSE
                       , math.style.negative = TRUE
                       , math.style.exponents = TRUE
                       , sanitize.text.function = function(x) {x})
                 , file = paste0("LATEX_outputs/LATEX_CombinedTable.ExNWater.",tempI,".txt"))
  
  
  
  
  
}

## Beta ##

for (b in betaList) {
  metric <- b
  # Get labels
  factorsCompareExNWater <- cbind(c("Water Only", "","NMF Alone","  ","Nereo","Nereo","Mast")
                                  , c("NMF Alone","Nereo","Mast","NereoMast","Mast","NereoMast","NereoMast"))
  colnamesTop <- c("","","","PERMANOVA","")
  colnamessubTop <- c("","","","(Community Dissimilarity)","")
  colnamesSecond <- c("Group 1","Group 2","p"," ","FDR adj. p")
  
  # Get ExNWater vs everything individual tests
  
  ExNWater.insertB <- betaListFiles.ExNWater[[paste0(metric)]][grep("PERMANOVA.*ExNvsEverything", betaListFiles.ExNWater[[paste0(metric)]])]
  ptemp <- signif(get(ExNWater.insertB)$aov.tab$`Pr(>F)`[1],3)
  rtemp <- round(get(ExNWater.insertB)$aov.tab$R2[1])
  ftemp <- round(get(ExNWater.insertB)$aov.tab$F.Model[1],3)
  dftemp <- paste0(round(get(ExNWater.insertB)$aov.tab$Df[1],3),",", round(get(ExNWater.insertB)$aov.tab$Df[3],3))
  ExNWater.insertB.toPaste <- cbind(ptemp, paste0("(R^2=",rtemp,",F.model=",ftemp,",df=",dftemp,")"), " ")
  
  # Merge everything 
  CombinedTable.ExNWater <- rbind(get(paste0("allPValues.ExNWater.",metric))[1:1,]
                                  , rep(" ", 3)
                                  , ExNWater.insertB.toPaste
                                  , rep(" ", 3)
                                  , get(paste0("allPValues.ExNWater.",metric))[4:6,])
  
  CombinedTable.ExNWater <- cbind(factorsCompareExNWater, CombinedTable.ExNWater)
  CombinedTable.ExNWater <- rbind(colnamesTop, colnamessubTop, colnamesSecond, CombinedTable.ExNWater)
  
  # Add dummy row and col name so it doesn't freak out
  rownames(CombinedTable.ExNWater) <- seq(1,nrow(CombinedTable.ExNWater), by = 1)
  colnames(CombinedTable.ExNWater) <- seq(1, ncol(CombinedTable.ExNWater), by = 1)
  
  # Rename
  
  assign(paste0("CombinedTable.ExNWater.",metric), CombinedTable.ExNWater)
  
  write.table(CombinedTable.ExNWater
              , sep = "\t"
              , row.names = FALSE
              , col.names = FALSE
              , file = paste0("LATEX_outputs/CombinedTable.ExNWater.",metric,".txt"))
  
  CombinedTable.changed <- gsub("R^2","$R^2$", CombinedTable.ExNWater , fixed = TRUE)
  CombinedTable.changed <- gsub("df","$df$", CombinedTable.changed , fixed = TRUE)
  CombinedTable.changed <- gsub("^p$","$p$", CombinedTable.changed , fixed = FALSE)
  CombinedTable.changed <- gsub("adj. p","adj. $p$", CombinedTable.changed , fixed = TRUE)
  
  capture.output(print(xtable(CombinedTable.changed)
                       , include.rownames= FALSE
                       , include.colnames = FALSE
                       , math.style.negative = TRUE
                       , math.style.exponents = TRUE
                       , sanitize.text.function = function(x) {x})
                 , file = paste0("LATEX_outputs/LATEX_CombinedTable.ExNWater.",metric,".txt"))
  
  
  
  
}

#### For All algae ####

## Alpha##
for (a in alphaList) {
  tempI <- gsub("_even_1000_alpha", "",a)
  
  # Get labels
  factorsCompareAlgae <- cbind(c("Nereo", "Nereo","Mast")
                                  , c("Mast","Water","Water"))
  colnamesTop <- c("","","","Welch's t-Test","")
  colnamessubTop <- c("","","","(Richness)","")
  colnamesSecond <- c("Group 1","Group 2","p"," ","FDR adj. p")

  # Merge everything 
  CombinedTable.All <-  get(paste0("allPvalues.algae.",tempI))
  
  CombinedTable.All <- cbind(factorsCompareAlgae, CombinedTable.All)
  CombinedTable.All <- rbind(colnamesTop, colnamessubTop, colnamesSecond, CombinedTable.All)
  
  # Add dummy row and col name so it doesn't freak out
  rownames(CombinedTable.All) <- seq(1,nrow(CombinedTable.All), by = 1)
  colnames(CombinedTable.All) <- seq(1, ncol(CombinedTable.All), by = 1)
  
  # Rename
  
  assign(paste0("CombinedTable.All.",tempI), CombinedTable.All)
  
  write.table(CombinedTable.All
              , sep = "\t"
              , row.names = FALSE
              , col.names = FALSE
              , file = paste0("LATEX_outputs/CombinedTable.All.",tempI,".txt"))
  
  CombinedTable.changed <- gsub("(t","($t$", CombinedTable.All , fixed = TRUE)
  CombinedTable.changed <- gsub("df","$df$", CombinedTable.changed , fixed = TRUE)
  CombinedTable.changed <- gsub("^p$","$p$ ", CombinedTable.changed , fixed = FALSE)
  CombinedTable.changed <- gsub("t-","$t$-", CombinedTable.changed , fixed = TRUE)
  CombinedTable.changed <- gsub("adj. p","adj. $p$", CombinedTable.changed , fixed = TRUE)
  
  
  capture.output(print(xtable(CombinedTable.changed)
                       , include.rownames= FALSE
                       , include.colnames = FALSE
                       , math.style.negative = TRUE
                       , math.style.exponents = TRUE
                       , sanitize.text.function = function(x) {x})
                 , file = paste0("LATEX_outputs/LATEX_CombinedTable.All.",tempI,".txt"))
  
  
}

## Beta ##

for (b in betaList) {
  metric <- b
  # Get labels
  factorsCompareAlgae <- cbind(c("Nereo", "Nereo","Mast")
                               , c("Mast","Water","Water"))
  colnamesTop <- c("","","","PERMANOVA","")
  colnamessubTop <- c("","","","(Community Dissimilarity)","")
  colnamesSecond <- c("Group 1","Group 2","p"," ","FDR adj. p")
  
  # Merge everything 
  CombinedTable.All <- get(paste0("allPValues.algae.water.",metric))
  
  CombinedTable.All <- cbind(factorsCompareAlgae, CombinedTable.All)
  CombinedTable.All <- rbind(colnamesTop, colnamessubTop, colnamesSecond, CombinedTable.All)
  
  # Add dummy row and col name so it doesn't freak out
  rownames(CombinedTable.All) <- seq(1,nrow(CombinedTable.All), by = 1)
  colnames(CombinedTable.All) <- seq(1, ncol(CombinedTable.All), by = 1)
  
  # Rename
  
  assign(paste0("CombinedTable.All.",metric), CombinedTable.All)
  
  write.table(CombinedTable.All
              , sep = "\t"
              , row.names = FALSE
              , col.names = FALSE
              , file = paste0("LATEX_outputs/CombinedTable.All.",metric,".txt"))
  
  CombinedTable.changed <- gsub("R^2","$R^2$", CombinedTable.All , fixed = TRUE)
  CombinedTable.changed <- gsub("df","$df$", CombinedTable.changed , fixed = TRUE)
  CombinedTable.changed <- gsub("^p$","$p$", CombinedTable.changed , fixed = FALSE)
  CombinedTable.changed <- gsub("adj. p","adj. $p$", CombinedTable.changed , fixed = TRUE)

  capture.output(print(xtable(CombinedTable.changed)
                       , include.rownames= FALSE
                       , include.colnames = FALSE
                       , math.style.negative = TRUE
                       , math.style.exponents = TRUE
                       , sanitize.text.function = function(x) {x}
                       )
                 , file = paste0("LATEX_outputs/LATEX_CombinedTable.All.",metric,".txt"))
  
  
  
}


#### For LoneIncube ####
## Alpha##
for (a in alphaList) {
  tempI <- gsub("_even_1000_alpha", "",a)
  
  # Get labels
  factorsCompareLoneIncube <- cbind(c("Nereo", "Nereo-Water")
                                  , c("Mast","Mast-Water"))
  colnamesTop <- c("","","Welch's t-Test","")
  colnamessubTop <- c("","","(Richness)","")
  colnamesSecond <- c("Group 1","Group 2","p","")
  
  # Merge everything 
  CombinedTable.LoneIncube <- get(paste0("allPvalues.LoneIncube.",tempI))[1:2,1:2]
  
  CombinedTable.LoneIncube <- cbind(factorsCompareLoneIncube, CombinedTable.LoneIncube)
  CombinedTable.LoneIncube <- rbind(colnamesTop, colnamessubTop, colnamesSecond, CombinedTable.LoneIncube)
  
  # Add dummy row and col name so it doesn't freak out
  rownames(CombinedTable.LoneIncube) <- seq(1,nrow(CombinedTable.LoneIncube), by = 1)
  colnames(CombinedTable.LoneIncube) <- seq(1, ncol(CombinedTable.LoneIncube), by = 1)
  
  # Rename
  
  assign(paste0("CombinedTable.LoneIncube.",tempI), CombinedTable.LoneIncube)
  
  write.table(CombinedTable.LoneIncube
              , sep = "\t"
              , row.names = FALSE
              , col.names = FALSE
              , file = paste0("LATEX_outputs/CombinedTable.LoneIncube.",tempI,".txt"))
  
  CombinedTable.changed <- gsub("(t","($t$", CombinedTable.LoneIncube , fixed = TRUE)
  CombinedTable.changed <- gsub("df","$df$", CombinedTable.changed , fixed = TRUE)
  CombinedTable.changed <- gsub("^p$","$p$ ", CombinedTable.changed , fixed = FALSE)
  CombinedTable.changed <- gsub("t-","$t$-", CombinedTable.changed , fixed = TRUE)
  CombinedTable.changed <- gsub("adj. p","adj. $p$", CombinedTable.changed , fixed = TRUE)
  
  capture.output(print(xtable(CombinedTable.changed)
                       , include.rownames= FALSE
                       , include.colnames = FALSE
                       , math.style.negative = TRUE
                       , math.style.exponents = TRUE
                       , sanitize.text.function = function(x) {x})
                 , file = paste0("LATEX_outputs/LATEX_CombinedTable.LoneIncube.",tempI,".txt"))
  
  
  
  
  
}

## Beta ##

for (b in betaList) {
  metric <- b
  # Get labels
  factorsCompareLoneIncube <- cbind(c("Nereo", "Nereo-Water")
                                    , c("Mast","Mast-Water"))
  colnamesTop <- c("","","PERMANOVA","")
  colnamessubTop <- c("","","(Community Dissimilarity)","")
  colnamesSecond <- c("Group 1","Group 2","p","")
  # Get LoneIncube vs everything individual tests
  
  # Merge everything 
  CombinedTable.LoneIncube <- get(paste0("allPValues.LoneIncube.",metric))[1:2,1:2]
                                  
  CombinedTable.LoneIncube <- cbind(factorsCompareLoneIncube, CombinedTable.LoneIncube)
  CombinedTable.LoneIncube <- rbind(colnamesTop, colnamessubTop, colnamesSecond, CombinedTable.LoneIncube)
  
  # Add dummy row and col name so it doesn't freak out
  rownames(CombinedTable.LoneIncube) <- seq(1,nrow(CombinedTable.LoneIncube), by = 1)
  colnames(CombinedTable.LoneIncube) <- seq(1, ncol(CombinedTable.LoneIncube), by = 1)
  
  # Rename
  
  assign(paste0("CombinedTable.LoneIncube.",metric), CombinedTable.LoneIncube)
  
  write.table(CombinedTable.LoneIncube
              , sep = "\t"
              , row.names = FALSE
              , col.names = FALSE
              , file = paste0("LATEX_outputs/CombinedTable.LoneIncube.",metric,".txt"))
  
  CombinedTable.changed <- gsub("R^2","$R^2$", CombinedTable.LoneIncube , fixed = TRUE)
  CombinedTable.changed <- gsub("df","$df$", CombinedTable.changed , fixed = TRUE)
  CombinedTable.changed <- gsub("^p$","$p$", CombinedTable.changed , fixed = FALSE)
  CombinedTable.changed <- gsub("adj. p","adj. $p$", CombinedTable.changed , fixed = TRUE)
  
  capture.output(print(xtable(CombinedTable.changed)
                       , include.rownames= FALSE
                       , include.colnames = FALSE
                       , math.style.negative = TRUE
                       , math.style.exponents = TRUE
                       , sanitize.text.function = function(x) {x})
                 , file = paste0("LATEX_outputs/LATEX_CombinedTable.LoneIncube.",metric,".txt"))
  
  
  
  
}



###### MERGE ##########

#ExN
Table.Combined.ExN <- cbind(get(paste0("CombinedTable.ExN.",choseB)), get(paste0("CombinedTable.ExN.",choseA))[,3:5])
write.table(Table.Combined.ExN, sep = "\t", row.names = FALSE, col.names = FALSE, file = paste0("LATEX_outputs/CHOSENMETRICS_CombinedTable.ExN.txt"))

Table.Combined.ExN.Edited <- gsub("R^2","$R^2$",Table.Combined.ExN, fixed = TRUE)
Table.Combined.ExN.Edited <- gsub("^p$","$p$",Table.Combined.ExN.Edited, fixed = FALSE)
Table.Combined.ExN.Edited <- gsub("(t=","($t$=",Table.Combined.ExN.Edited, fixed = TRUE)
Table.Combined.ExN.Edited <- gsub("df","$df$",Table.Combined.ExN.Edited, fixed = TRUE)
Table.Combined.ExN.Edited <- gsub("t-","$t$-",Table.Combined.ExN.Edited, fixed = TRUE)
Table.Combined.ExN.Edited <- gsub("adj. p","adj. $p$",Table.Combined.ExN.Edited, fixed = TRUE)

Table.Combined.ExN.Edited <- apply(Table.Combined.ExN.Edited, MARGIN = c(1,2),FUN = function(x) {
  if (!is.na(as.numeric(x))) {
    if (as.numeric(x) < 0.001) {
      return("\\bf{<0.001}")
    } else if (as.numeric(x) < 0.05) {
      return(paste0("\\bf{",x,"}"))
    }else {
      return(x)
    }
  } else {
    return(x)
  }
})

capture.output(print(xtable(Table.Combined.ExN.Edited)
                     , include.colnames = FALSE
                     , include.rownames = FALSE
                     , math.style.negative = TRUE
                     , math.style.exponents = TRUE
                     , sanitize.text.function = function(x) {x}
                     )
               , file = paste0("LATEX_outputs/CHOSENMETRICS_LATEX_CombinedTable.ExN.txt"))



# ExNWater
Table.Combined.ExNWater <- cbind(get(paste0("CombinedTable.ExNWater.",choseB)), get(paste0("CombinedTable.ExNWater.",choseA))[,3:5])
write.table(Table.Combined.ExNWater, sep = "\t", row.names = FALSE, col.names = FALSE, file = paste0("LATEX_outputs/CHOSENMETRICS_CombinedTable.ExNWater.txt"))

Table.Combined.ExNWater.Edited <- gsub("R^2","$R^2$",Table.Combined.ExNWater, fixed = TRUE)
Table.Combined.ExNWater.Edited <- gsub("^p$","$p$",Table.Combined.ExNWater.Edited, fixed = FALSE)
Table.Combined.ExNWater.Edited <- gsub("(t=","($t$=",Table.Combined.ExNWater.Edited, fixed = TRUE)
Table.Combined.ExNWater.Edited <- gsub("df","$df$",Table.Combined.ExNWater.Edited, fixed = TRUE)
Table.Combined.ExNWater.Edited <- gsub("t-","$t$-",Table.Combined.ExNWater.Edited, fixed = TRUE)
Table.Combined.ExNWater.Edited <- gsub("adj. p","adj. $p$",Table.Combined.ExNWater.Edited, fixed = TRUE)

Table.Combined.ExNWater.Edited <- apply(Table.Combined.ExNWater.Edited, MARGIN = c(1,2),FUN = function(x) {
  if (!is.na(as.numeric(x))) {
    if (as.numeric(x) < 0.001) {
      return("\\bf{<0.001}")
    } else if (as.numeric(x) < 0.05) {
      return(paste0("\\bf{",x,"}"))
    }else {
      return(x)
    }
  } else {
    return(x)
  }
})

capture.output(print(xtable(Table.Combined.ExNWater.Edited)
                     , include.colnames = FALSE
                     , include.rownames = FALSE
                     , math.style.negative = TRUE
                     , math.style.exponents = TRUE
                    , sanitize.text.function = function(x) {x}
                     )
               , file = paste0("LATEX_outputs/CHOSENMETRICS_LATEX_CombinedTable.ExNWater.txt"))


# Algae ALL
Table.Combined.All <- cbind(get(paste0("CombinedTable.All.",choseB)), get(paste0("CombinedTable.All.",choseA))[,3:5])
write.table(Table.Combined.All, sep = "\t", row.names = FALSE, col.names = FALSE, file = paste0("LATEX_outputs/CHOSENMETRICS_CombinedTable.All.txt"))

Table.Combined.All.Edited <- gsub("R^2","$R^2$",Table.Combined.All, fixed = TRUE)
Table.Combined.All.Edited <- gsub("^p$","$p$",Table.Combined.All.Edited, fixed = FALSE)
Table.Combined.All.Edited <- gsub("(t=","($t$=",Table.Combined.All.Edited, fixed = TRUE)
Table.Combined.All.Edited <- gsub("df","$df$",Table.Combined.All.Edited, fixed = TRUE)
Table.Combined.All.Edited <- gsub("t-","$t$-",Table.Combined.All.Edited, fixed = TRUE)
Table.Combined.All.Edited <- gsub("adj. p","adj. $p$",Table.Combined.All.Edited, fixed = TRUE)

Table.Combined.All.Edited <- apply(Table.Combined.All.Edited, MARGIN = c(1,2),FUN = function(x) {
  if (!is.na(as.numeric(x))) {
    if (as.numeric(x) < 0.001) {
      return("\\bf{<0.001}")
    } else if (as.numeric(x) < 0.05) {
      return(paste0("\\bf{",x,"}"))
    }else {
      return(x)
    }
  } else {
    return(x)
  }
})

capture.output(print(xtable(Table.Combined.All.Edited)
                     , include.colnames = FALSE
                     , include.rownames = FALSE
                     , math.style.negative = TRUE
                     , math.style.exponents = TRUE
                     , sanitize.text.function = function(x) {x}
)
, file = paste0("LATEX_outputs/CHOSENMETRICS_LATEX_CombinedTable.All.txt"))

# LoneIncube
Table.Combined.LoneIncube <- cbind(get(paste0("CombinedTable.LoneIncube.",choseB)), get(paste0("CombinedTable.LoneIncube.",choseA))[,3:4])
write.table(Table.Combined.LoneIncube, sep = "\t", row.names = FALSE, col.names = FALSE, file = paste0("LATEX_outputs/CHOSENMETRICS_CombinedTable.LoneIncube.txt"))

Table.Combined.LoneIncube.Edited <- gsub("R^2","$R^2$",Table.Combined.LoneIncube, fixed = TRUE)
Table.Combined.LoneIncube.Edited <- gsub("^p$","$p$",Table.Combined.LoneIncube.Edited, fixed = FALSE)
Table.Combined.LoneIncube.Edited <- gsub("(t=","($t$=",Table.Combined.LoneIncube.Edited, fixed = TRUE)
Table.Combined.LoneIncube.Edited <- gsub("df","$df$",Table.Combined.LoneIncube.Edited, fixed = TRUE)
Table.Combined.LoneIncube.Edited <- gsub("t-","$t$-",Table.Combined.LoneIncube.Edited, fixed = TRUE)
Table.Combined.LoneIncube.Edited <- gsub("adj. p","adj. $p$",Table.Combined.LoneIncube.Edited, fixed = TRUE)

Table.Combined.LoneIncube.Edited <- apply(Table.Combined.LoneIncube.Edited, MARGIN = c(1,2),FUN = function(x) {
  if (!is.na(as.numeric(x))) {
    if (as.numeric(x) < 0.001) {
      return("\\bf{<0.001}")
    } else if (as.numeric(x) < 0.05) {
      return(paste0("\\bf{",x,"}"))
    }else {
      return(x)
    }
  } else {
    return(x)
  }
})

capture.output(print(xtable(Table.Combined.LoneIncube.Edited)
                     , include.colnames = FALSE
                     , include.rownames = FALSE
                     , math.style.negative = TRUE
                     , math.style.exponents = TRUE
                     , sanitize.text.function = function(x) {x}
)
, file = paste0("LATEX_outputs/CHOSENMETRICS_LATEX_CombinedTable.LoneIncube.txt"))


########### Combined Tables After Revisions ############
### EXN and EXNWATER ##
# Get single compare ExN beta
singleCompareExN.beta <- paste0("\\makecell{", Table.Combined.ExN.Edited[4:6,3]," \\\\", Table.Combined.ExN.Edited[4:6,4], "}")
# Get pairwise compare ExN beta
pairwiseCompareExN.beta <- paste0("\\makecell{", Table.Combined.ExN.Edited[7:9,5], " \\\\", Table.Combined.ExN.Edited[7:9,4],"}")

# Get control compare ExNWater beta
controlCompareExNWater.beta <- paste0("\\makecell{", Table.Combined.ExNWater.Edited[4,3]," \\\\", Table.Combined.ExNWater.Edited[4,4], "}")
# Get single compare ExNWater beta
singleCompareExNWater.beta <- paste0("\\makecell{", Table.Combined.ExNWater.Edited[5:7,3]," \\\\", Table.Combined.ExNWater.Edited[5:7,4], "}")
# Get pariwise compare ExNWater beta
pairwiseCompareExNWater.beta <- paste0("\\makecell{", Table.Combined.ExNWater.Edited[8:10,5]," \\\\", Table.Combined.ExNWater.Edited[8:10,4], "}")

# Get single compare ExN alpha
singleCompareExN.alpha <- paste0("\\makecell{", Table.Combined.ExN.Edited[4:6,6]," \\\\", Table.Combined.ExN.Edited[4:6,7], "}")
# Get pairwise compare ExN alpha
pairwiseCompareExN.alpha <- paste0("\\makecell{", Table.Combined.ExN.Edited[7:9,8], " \\\\", Table.Combined.ExN.Edited[7:9,7],"}")

# Get control compare ExNWater alpha
controlCompareExNWater.alpha <- paste0("\\makecell{", Table.Combined.ExNWater.Edited[4,6]," \\\\", Table.Combined.ExNWater.Edited[4,7], "}")
# Get single compare ExNWater alpha
singleCompareExNWater.alpha <- paste0("\\makecell{", Table.Combined.ExNWater.Edited[5:7,6]," \\\\", Table.Combined.ExNWater.Edited[5:7,7], "}")
# Get pariwise compare ExNWater alpha
pairwiseCompareExNWater.alpha <- paste0("\\makecell{", Table.Combined.ExNWater.Edited[8:10,8]," \\\\", Table.Combined.ExNWater.Edited[8:10,7], "}")


# Get first two columns for group compairsons; manual because I add in graphics
# groupLabels <- Table.Combined.ExNWater.Edited[4:10,1:2]
groupLabels <- cbind(c("\\wateronly Water Only", "","\\NMF NMF Alone", "", "\\Nereo Nereo", "\\Nereo Nereo","\\Mast Mast")
                     , c("\\NMF NMF Alone","\\Nereo Nereo", "\\Mast Mast", "\\NereoMast NereoMast", "\\Mast Mast", "\\NereoMast NereoMast","\\NereoMast NereoMast"))

# Get headers
header1 <- c("","","PERMANOVA","","Welch's $t$-test","")
header2 <- c("","","(Community Dissimilarity)","","(Richness)","")
header3 <- c("Group 1","Group 2","Water samples","NMF surface","Water samples","NMF surface")

# Get data columns
col2 <- rbind("",cbind(singleCompareExN.beta),cbind(pairwiseCompareExN.beta[1:3]))
col1 <- rbind(controlCompareExNWater.beta,cbind(singleCompareExNWater.beta),cbind(pairwiseCompareExNWater.beta))
col4 <- rbind("", cbind(singleCompareExN.alpha), cbind(pairwiseCompareExN.alpha))
col3 <- rbind(controlCompareExNWater.alpha, cbind(singleCompareExNWater.alpha), cbind(pairwiseCompareExNWater.alpha))

# combine group labels
dataCombined <- cbind(groupLabels, col1, col2, col3, col4)

# Combine headers

dataCombined <- rbind(header1,header2,header3,dataCombined)

# Print
capture.output(print(xtable(dataCombined)
                     , include.colnames = FALSE
                     , include.rownames = FALSE
                     , math.style.negative = TRUE
                     , math.style.exponents = TRUE
                     , sanitize.text.function = function(x) {x}
)
, file = paste0("LATEX_outputs/CHOSENMETRICS_LATEX_CombinedTable_postrevisions_ALL.txt"))


## ALGAE ALL ##
Table.Combined.All.Edited.postrev <- Table.Combined.All.Edited[,c(1,2,5,4,8,7)]
Table.Combined.All.Edited.postrev <- gsub("FDR adj. ","",Table.Combined.All.Edited.postrev)
# Print
capture.output(print(xtable(Table.Combined.All.Edited.postrev)
                     , include.colnames = FALSE
                     , include.rownames = FALSE
                     , math.style.negative = TRUE
                     , math.style.exponents = TRUE
                     , sanitize.text.function = function(x) {x}
)
, file = paste0("LATEX_outputs/CHOSENMETRICS_LATEX_allAlgae_postrevisions_ALLALGAE.txt"))




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



