#!/bin/RScript
set.seed(3)
library(optparse)

option_list = list(
  make_option(c("-b", "--braycurtisdm"), type="character",
              help="Distance Matrix- Braycurtis"),
  make_option(c("-u", "--unweightedunifracdm"),
              help="Distance Matrix- Unweighted Unifrac", type="character"),
  make_option(c("-w", "--weightedunifracdm"), 
              help="Distance Matrix- Weighted Unifrac", type="character"),
  make_option(c("-m", "--mappingfile"),
              help="Mappingfile", type="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

BCPWD = opt$braycurtisdm
UWUFPWD = opt$unweightedunifracdm
WUFPWD = opt$weightedunifracdm
MFPWD = opt$mappingfile
###########################
# FOR TESTING
# 
setwd("/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/zz_NEREOINCUBE_16may2017/1_analysis")
BCPWD<- "/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/zz_NEREOINCUBE_16may2017/1_analysis/ANALYSIS_ALPHABETATAXA/beta_div/bray_curtis_dm.txt"
WUFPWD <- "/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/zz_NEREOINCUBE_16may2017/1_analysis/ANALYSIS_ALPHABETATAXA/beta_div/weighted_unifrac_dm.txt"
UWUFPWD <- "/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/zz_NEREOINCUBE_16may2017/1_analysis/ANALYSIS_ALPHABETATAXA/beta_div/unweighted_unifrac_dm.txt"
MFPWD <- "/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/zz_NEREOINCUBE_16may2017/1_analysis/OTU_MP_filt/MF_nochlpmito_m1000.txt"
###########################

# MAKE BETA PLOTS
library("MASS")
library("vegan")
library("stats")

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
                 , stringsAsFactors = FALSE
                 , header = TRUE
                 , row.names= 1)
# Make sure MF has everything in it-- can just use an dm because they should be all same

MF <- MF[sapply(rownames(dm.BC), function(x) {
  grep(paste0("^",x,"$"), rownames(MF))
}),]

MF.ExN <- MF[grep("ExN.Nereotest.", rownames(MF)),]
MF.ExNWater <- MF[grep("Water.Nereotest.", rownames(MF)),]
MF.LoneIncube <- MF[grep("Loneincube", rownames(MF)),]




system("mkdir BETAPLOTS")

########## FUNCTION ############

# 
# pairwise.adonis <- function(x,factors, sim.method = 'bray', p.adjust.m ='bonferroni')
# {
#   library(vegan)
#   co = combn(unique(factors),2)
#   pairs = c()
#   F.Model =c()
#   R2 = c()
#   p.value = c()
#   
#   for(elem in 1:ncol(co)){
#     ad = adonis(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors %in% c(co[1,elem],co[2,elem])] , method =sim.method);
#     pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
#     F.Model =c(F.Model,ad$aov.tab[1,4]);
#     R2 = c(R2,ad$aov.tab[1,5]);
#     p.value = c(p.value,ad$aov.tab[1,6])
#   }
#   p.adjusted = p.adjust(p.value,method=p.adjust.m)
#   pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
#   return(pairw.res)
# } 

####### UWUF #############
dm.UWUF.ExN <- dm.UWUF[unlist(sapply(rownames(MF.ExN), function(x) {
  grep(x, rownames(dm.UWUF))})),unlist(sapply(rownames(MF.ExN), function(x) {
    grep(x, colnames(dm.UWUF))}))]

dm.UWUF.ExNWater <- dm.UWUF[unlist(sapply(rownames(MF.ExNWater), function(x) {
  grep(x, rownames(dm.UWUF))})),unlist(sapply(rownames(MF.ExNWater), function(x) {
    grep(x, colnames(dm.UWUF))}))]

dm.UWUF.LoneIncube <- dm.UWUF[unlist(sapply(rownames(MF.LoneIncube), function(x) {
  grep(x, rownames(dm.UWUF))})),unlist(sapply(rownames(MF.LoneIncube), function(x) {
    grep(x, colnames(dm.UWUF))}))]

NMDS.UWUF.ExN <- isoMDS(as.matrix(dm.UWUF.ExN), y = cmdscale(as.matrix(dm.UWUF.ExN), 2))

NMDS.UWUF.ExNWater <- isoMDS(as.matrix(dm.UWUF.ExNWater), y = cmdscale(as.matrix(dm.UWUF.ExNWater), 2))

NMDS.UWUF.LoneIncube <- isoMDS(as.matrix(dm.UWUF.LoneIncube), y = cmdscale(as.matrix(dm.UWUF.LoneIncube), 2))

###### STATS ##########
# ExN
ANOVA.UWUF.ExN <- adonis(dm.UWUF.ExN ~ ColRep, data = MF.ExN)
ANOSIM.UWUF.ExN <- anosim(dm.UWUF.ExN, grouping = MF.ExN$ColRep)
# pairwaiseAdonis.UWUF.ExN <- pairwise.adonis(dm.UWUF.ExN, factors = MF.ExN$ColRep)


# ExN Split
dm.UWUF.ExN.ExNvsNereo <- dm.UWUF.ExN[grep("([.]ExN[.])|([.]Nereo[.])", rownames(dm.UWUF.ExN)), grep("([.]ExN[.])|([.]Nereo[.])", colnames(dm.UWUF.ExN))]
MF.ExN.ExNvsNereo <- MF.ExN[grep("([.]ExN[.])|([.]Nereo[.])", rownames(MF.ExN)),]
ANOVA.UWUF.ExN.ExNvsNereo <- adonis(dm.UWUF.ExN.ExNvsNereo ~ ColRep, data = MF.ExN.ExNvsNereo)
ANOSIM.UWUF.ExN.ExNvsNereo <- anosim(dm.UWUF.ExN.ExNvsNereo, grouping = MF.ExN.ExNvsNereo$ColRep)

dm.UWUF.ExN.ExNvsMast <- dm.UWUF.ExN[grep("([.]ExN[.])|([.]Mast[.])", rownames(dm.UWUF.ExN)), grep("([.]ExN[.])|([.]Mast[.])", colnames(dm.UWUF.ExN))]
MF.ExN.ExNvsMast <- MF.ExN[grep("([.]ExN[.])|([.]Mast[.])", rownames(MF.ExN)),]
ANOVA.UWUF.ExN.ExNvsMast <- adonis(dm.UWUF.ExN.ExNvsMast ~ ColRep, data = MF.ExN.ExNvsMast)
ANOSIM.UWUF.ExN.ExNvsMast <- anosim(dm.UWUF.ExN.ExNvsMast, grouping = MF.ExN.ExNvsMast$ColRep)


dm.UWUF.ExN.ExNvsNereoMast <- dm.UWUF.ExN[grep("([.]ExN[.])|([.]NereoMast[.])", rownames(dm.UWUF.ExN)), grep("([.]ExN[.])|([.]NereoMast[.])", colnames(dm.UWUF.ExN))]
MF.ExN.ExNvsNereoMast <- MF.ExN[grep("([.]ExN[.])|([.]NereoMast[.])", rownames(MF.ExN)),]
ANOVA.UWUF.ExN.ExNvsNereoMast <- adonis(dm.UWUF.ExN.ExNvsNereoMast ~ ColRep, data = MF.ExN.ExNvsNereoMast)
ANOSIM.UWUF.ExN.ExNvsNereoMast <- anosim(dm.UWUF.ExN.ExNvsNereoMast, grouping = MF.ExN.ExNvsNereoMast$ColRep)


dm.UWUF.ExN.NereovsMast <- dm.UWUF.ExN[grep("([.]Nereo[.])|([.]Mast[.])", rownames(dm.UWUF.ExN)), grep("([.]Nereo[.])|([.]Mast[.])", colnames(dm.UWUF.ExN))]
MF.ExN.NereovsMast <- MF.ExN[grep("([.]Nereo[.])|([.]Mast[.])", rownames(MF.ExN)),]
ANOVA.UWUF.ExN.NereovsMast <- adonis(dm.UWUF.ExN.NereovsMast ~ ColRep, data = MF.ExN.NereovsMast)
ANOSIM.UWUF.ExN.NereovsMast <- anosim(dm.UWUF.ExN.NereovsMast, grouping = MF.ExN.NereovsMast$ColRep)

dm.UWUF.ExN.NereovsNereoMast <- dm.UWUF.ExN[grep("([.]Nereo[.])|([.]NereoMast[.])", rownames(dm.UWUF.ExN)), grep("([.]Nereo[.])|([.]NereoMast[.])", colnames(dm.UWUF.ExN))]
MF.ExN.NereovsNereoMast <- MF.ExN[grep("([.]Nereo[.])|([.]NereoMast[.])", rownames(MF.ExN)),]
ANOVA.UWUF.ExN.NereovsNereoMast <- adonis(dm.UWUF.ExN.NereovsNereoMast ~ ColRep, data = MF.ExN.NereovsNereoMast)
ANOSIM.UWUF.ExN.NereovsNereoMast <- anosim(dm.UWUF.ExN.NereovsNereoMast, grouping = MF.ExN.NereovsNereoMast$ColRep)


dm.UWUF.ExN.MastvsNereoMast <- dm.UWUF.ExN[grep("([.]Mast[.])|([.]NereoMast[.])", rownames(dm.UWUF.ExN)), grep("([.]Mast[.])|([.]NereoMast[.])", colnames(dm.UWUF.ExN))]
MF.ExN.MastvsNereoMast <- MF.ExN[grep("([.]Mast[.])|([.]NereoMast[.])", rownames(MF.ExN)),]
ANOVA.UWUF.ExN.MastvsNereoMast <- adonis(dm.UWUF.ExN.MastvsNereoMast ~ ColRep, data = MF.ExN.MastvsNereoMast)
ANOSIM.UWUF.ExN.MastvsNereoMast <- anosim(dm.UWUF.ExN.MastvsNereoMast, grouping = MF.ExN.MastvsNereoMast$ColRep)

# 6 comparisons; 0.05/6 = 0.00625
allPValues.ExN <- matrix(nrow= 6)
rownames(allPValues.ExN) <- c("ExNvsNereo"
                              ,"ExNvsMast"
                              ,"ExNvsNereoMast"
                              ,"NereovsMast"
                              , "NereovsNereoMast"
                              , "MastvsNereoMast")
allPValues.ExN[1,1] <- ANOVA.UWUF.ExN.ExNvsNereo$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExN[2,1] <- ANOVA.UWUF.ExN.ExNvsMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExN[3,1] <- ANOVA.UWUF.ExN.ExNvsNereoMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExN[4,1] <- ANOVA.UWUF.ExN.NereovsMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExN[5,1] <- ANOVA.UWUF.ExN.NereovsNereoMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExN[6,1] <- ANOVA.UWUF.ExN.MastvsNereoMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExN <- cbind(allPValues.ExN, p.adjust(allPValues.ExN,method = "fdr", n = 6))
colnames(allPValues.ExN) <- c("P","fdr_adj")

# ExN Water
ANOVA.UWUF.ExNWater <- adonis(dm.UWUF.ExNWater ~ ColRep, data = MF.ExNWater)
ANOSIM.UWUF.ExNWater <- anosim(dm.UWUF.ExNWater, grouping = MF.ExNWater$ColRep) # NOT WORKING
# pairwaiseAdonis.UWUF.ExNWater <- pairwise.adonis(dm.UWUF.ExNWater, factors = MF.ExNWater$ColRep)


# ExN Water split
dm.UWUF.ExNWater.ExNvsH2O <- dm.UWUF.ExNWater[grep("([.]ExN[.])|([.]H2O[.])", rownames(dm.UWUF.ExNWater)), grep("([.]ExN[.])|([.]H2O[.])", colnames(dm.UWUF.ExNWater))]
MF.ExNWater.ExNvsH2O <- MF.ExNWater[grep("([.]ExN[.])|([.]H2O[.])", rownames(MF.ExNWater)),]
ANOVA.UWUF.ExNWater.ExNvsH2O <- adonis(dm.UWUF.ExNWater.ExNvsH2O ~ ColRep, data = MF.ExNWater.ExNvsH2O)
ANOSIM.UWUF.ExNWater.ExNvsH2O <- anosim(dm.UWUF.ExNWater.ExNvsH2O, grouping = MF.ExNWater.ExNvsH2O$ColRep)

dm.UWUF.ExNWater.ExNvsNereo <- dm.UWUF.ExNWater[grep("([.]ExN[.])|([.]Nereo[.])", rownames(dm.UWUF.ExNWater)), grep("([.]ExN[.])|([.]Nereo[.])", colnames(dm.UWUF.ExNWater))]
MF.ExNWater.ExNvsNereo <- MF.ExNWater[grep("([.]ExN[.])|([.]Nereo[.])", rownames(MF.ExNWater)),]
ANOVA.UWUF.ExNWater.ExNvsNereo <- adonis(dm.UWUF.ExNWater.ExNvsNereo ~ ColRep, data = MF.ExNWater.ExNvsNereo)
ANOSIM.UWUF.ExNWater.ExNvsNereo <- anosim(dm.UWUF.ExNWater.ExNvsNereo, grouping = MF.ExNWater.ExNvsNereo$ColRep)

dm.UWUF.ExNWater.ExNvsMast <- dm.UWUF.ExNWater[grep("([.]ExN[.])|([.]Mast[.])", rownames(dm.UWUF.ExNWater)), grep("([.]ExN[.])|([.]Mast[.])", colnames(dm.UWUF.ExNWater))]
MF.ExNWater.ExNvsMast <- MF.ExNWater[grep("([.]ExN[.])|([.]Mast[.])", rownames(MF.ExNWater)),]
ANOVA.UWUF.ExNWater.ExNvsMast <- adonis(dm.UWUF.ExNWater.ExNvsMast ~ ColRep, data = MF.ExNWater.ExNvsMast)
ANOSIM.UWUF.ExNWater.ExNvsMast <- anosim(dm.UWUF.ExNWater.ExNvsMast, grouping = MF.ExNWater.ExNvsMast$ColRep)

dm.UWUF.ExNWater.ExNvsNereoMast <- dm.UWUF.ExNWater[grep("([.]ExN[.])|([.]NereoMast[.])", rownames(dm.UWUF.ExNWater)), grep("([.]ExN[.])|([.]NereoMast[.])", colnames(dm.UWUF.ExNWater))]
MF.ExNWater.ExNvsNereoMast <- MF.ExNWater[grep("([.]ExN[.])|([.]NereoMast[.])", rownames(MF.ExNWater)),]
ANOVA.UWUF.ExNWater.ExNvsNereoMast <- adonis(dm.UWUF.ExNWater.ExNvsNereoMast ~ ColRep, data = MF.ExNWater.ExNvsNereoMast)
ANOSIM.UWUF.ExNWater.ExNvsNereoMast <- anosim(dm.UWUF.ExNWater.ExNvsNereoMast, grouping = MF.ExNWater.ExNvsNereoMast$ColRep)

dm.UWUF.ExNWater.H2OvsNereo <- dm.UWUF.ExNWater[grep("([.]H2O[.])|([.]Nereo[.])", rownames(dm.UWUF.ExNWater)), grep("([.]H2O[.])|([.]Nereo[.])", colnames(dm.UWUF.ExNWater))]
MF.ExNWater.H2OvsNereo <- MF.ExNWater[grep("([.]H2O[.])|([.]Nereo[.])", rownames(MF.ExNWater)),]
ANOVA.UWUF.ExNWater.H2OvsNereo <- adonis(dm.UWUF.ExNWater.H2OvsNereo ~ ColRep, data = MF.ExNWater.H2OvsNereo)
ANOSIM.UWUF.ExNWater.H2OvsNereo <- anosim(dm.UWUF.ExNWater.H2OvsNereo, grouping = MF.ExNWater.H2OvsNereo$ColRep)

dm.UWUF.ExNWater.H2OvsMast <- dm.UWUF.ExNWater[grep("([.]H2O[.])|([.]Mast[.])", rownames(dm.UWUF.ExNWater)), grep("([.]H2O[.])|([.]Mast[.])", colnames(dm.UWUF.ExNWater))]
MF.ExNWater.H2OvsMast <- MF.ExNWater[grep("([.]H2O[.])|([.]Mast[.])", rownames(MF.ExNWater)),]
ANOVA.UWUF.ExNWater.H2OvsMast <- adonis(dm.UWUF.ExNWater.H2OvsMast ~ ColRep, data = MF.ExNWater.H2OvsMast)
ANOSIM.UWUF.ExNWater.H2OvsMast <- anosim(dm.UWUF.ExNWater.H2OvsMast, grouping = MF.ExNWater.H2OvsMast$ColRep)

dm.UWUF.ExNWater.H2OvsNereoMast <- dm.UWUF.ExNWater[grep("([.]H2O[.])|([.]NereoMast[.])", rownames(dm.UWUF.ExNWater)), grep("([.]H2O[.])|([.]Mast[.])", colnames(dm.UWUF.ExNWater))]
MF.ExNWater.H2OvsNereoMast <- MF.ExNWater[grep("([.]H2O[.])|([.]NereoMast[.])", rownames(MF.ExNWater)),]
ANOVA.UWUF.ExNWater.H2OvsNereoMast <- adonis(dm.UWUF.ExNWater.H2OvsNereoMast ~ ColRep, data = MF.ExNWater.H2OvsNereoMast)
ANOSIM.UWUF.ExNWater.H2OvsNereoMast <- anosim(dm.UWUF.ExNWater.H2OvsNereoMast, grouping = MF.ExNWater.H2OvsNereoMast$ColRep)

dm.UWUF.ExNWater.NereovsMast <- dm.UWUF.ExNWater[grep("([.]Nereo[.])|([.]Mast[.])", rownames(dm.UWUF.ExNWater)), grep("([.]Nereo[.])|([.]Mast[.])", colnames(dm.UWUF.ExNWater))]
MF.ExNWater.NereovsMast <- MF.ExNWater[grep("([.]Nereo[.])|([.]Mast[.])", rownames(MF.ExNWater)),]
ANOVA.UWUF.ExNWater.NereovsMast <- adonis(dm.UWUF.ExNWater.NereovsMast ~ ColRep, data = MF.ExNWater.NereovsMast)
ANOSIM.UWUF.ExNWater.NereovsMast <- anosim(dm.UWUF.ExNWater.NereovsMast, grouping = MF.ExNWater.NereovsMast$ColRep)

# 8 comparisons; 0.05/8 = 0.00625
allPValues.ExNWater <- matrix(nrow= 8)
rownames(allPValues.ExNWater) <- c("ExNvsH2O"
                                   ,"ExNvsNereo"
                                   ,"ExNvsMast"
                                   ,"ExNvsNereoMast"
                                   ,"H2OvsNereo"
                                   ,"H2OvsMast"
                                   ,"H2OvsNereoMast"
                                   ,"NereovsMast")
allPValues.ExNWater[1,1] <- ANOVA.UWUF.ExNWater.ExNvsH2O$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[2,1] <- ANOVA.UWUF.ExNWater.ExNvsNereo$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[3,1] <- ANOVA.UWUF.ExNWater.ExNvsMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[4,1] <- ANOVA.UWUF.ExNWater.ExNvsNereoMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[5,1] <- ANOVA.UWUF.ExNWater.H2OvsNereo$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[6,1] <- ANOVA.UWUF.ExNWater.H2OvsMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[7,1] <- ANOVA.UWUF.ExNWater.H2OvsNereoMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[8,1] <- ANOVA.UWUF.ExNWater.NereovsMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater <- cbind(allPValues.ExNWater, p.adjust(allPValues.ExNWater,method = "fdr", n = 10))
colnames(allPValues.ExNWater) <- c("P","fdr_adj")



# Lone Incube
# Add another column
MF.LoneIncube$WorS <- MF.LoneIncube$Substrate
MF.LoneIncube$WorS <- gsub("(Nereo)|(Mast)", "Seaweed", MF.LoneIncube[,'WorS'])

ANOVA.UWUF.LoneIncube <- adonis(dm.UWUF.LoneIncube ~ WorS*Treatment, data = MF.LoneIncube)
ANOSIM.UWUF.LoneIncube <- anosim(dm.UWUF.LoneIncube, grouping = MF.LoneIncube$ColRep) # NOT WORKING
# pairwaiseAdonis.UWUF.LoneIncube <- pairwise.adonis(dm.UWUF.LoneIncube, factors = MF.LoneIncube$ColRep)

# Lone Incube split
dm.UWUF.LoneIncube.NereovsMast <- dm.UWUF.LoneIncube[grep("(^Nereo[.])|(^Mast[.])", rownames(dm.UWUF.LoneIncube)), grep("(^Nereo[.])|(^Mast[.])", colnames(dm.UWUF.LoneIncube))]
MF.LoneIncube.NereovsMast <- MF.LoneIncube[grep("(^Nereo[.])|(^Mast[.])", rownames(MF.LoneIncube)),]
ANOVA.UWUF.LoneIncube.NereovsMast <- adonis(dm.UWUF.LoneIncube.NereovsMast ~ ColRep, data = MF.LoneIncube.NereovsMast)
ANOSIM.UWUF.LoneIncube.NereovsMast <- anosim(dm.UWUF.LoneIncube.NereovsMast, grouping = MF.LoneIncube.NereovsMast$ColRep)

dm.UWUF.LoneIncube.waters <- dm.UWUF.LoneIncube[grep("^water[.]", rownames(dm.UWUF.LoneIncube)), grep("^water[.]", colnames(dm.UWUF.LoneIncube))]
MF.LoneIncube.waters <- MF.LoneIncube[grep("^water[.]", rownames(MF.LoneIncube)),]
ANOVA.UWUF.LoneIncube.waters <- adonis(dm.UWUF.LoneIncube.waters ~ ColRep, data = MF.LoneIncube.waters)
ANOSIM.UWUF.LoneIncube.waters <- anosim(dm.UWUF.LoneIncube.waters, grouping = MF.LoneIncube.waters$ColRep)

dm.UWUF.LoneIncube.Nereovswater <- dm.UWUF.LoneIncube[grep("[.]Nereo[.]", rownames(dm.UWUF.LoneIncube)), grep("[.]Nereo[.]", colnames(dm.UWUF.LoneIncube))]
MF.LoneIncube.Nereovswater <- MF.LoneIncube[grep("[.]Nereo[.]", rownames(MF.LoneIncube)),]
ANOVA.UWUF.LoneIncube.Nereovswater <- adonis(dm.UWUF.LoneIncube.Nereovswater ~ ColRep, data = MF.LoneIncube.Nereovswater)
ANOSIM.UWUF.LoneIncube.Nereovswater <- anosim(dm.UWUF.LoneIncube.Nereovswater, grouping = MF.LoneIncube.Nereovswater$ColRep)

dm.UWUF.LoneIncube.Mastvswater <- dm.UWUF.LoneIncube[grep("[.]Mast[.]", rownames(dm.UWUF.LoneIncube)), grep("[.]Mast[.]", colnames(dm.UWUF.LoneIncube))]
MF.LoneIncube.Mastvswater <- MF.LoneIncube[grep("[.]Mast[.]", rownames(MF.LoneIncube)),]
ANOVA.UWUF.LoneIncube.Mastvswater <- adonis(dm.UWUF.LoneIncube.Mastvswater ~ ColRep, data = MF.LoneIncube.Mastvswater)
ANOSIM.UWUF.LoneIncube.Mastvswater <- anosim(dm.UWUF.LoneIncube.Mastvswater, grouping = MF.LoneIncube.Mastvswater$ColRep)

# 4 comparisons; 0.05/4 = 0.00625
allPValues.LoneIncube <- matrix(nrow= 4)
rownames(allPValues.LoneIncube) <- c("NereovsMast"
                                   ,"waters"
                                   ,"Nereovswater"
                                   ,"Mastvswater"
                                   )
allPValues.LoneIncube[1,1] <- ANOVA.UWUF.LoneIncube.NereovsMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.LoneIncube[2,1] <- ANOVA.UWUF.LoneIncube.waters$aov.tab[6]$`Pr(>F)`[1]
allPValues.LoneIncube[3,1] <- ANOVA.UWUF.LoneIncube.Nereovswater$aov.tab[6]$`Pr(>F)`[1]
allPValues.LoneIncube[4,1] <- ANOVA.UWUF.LoneIncube.Mastvswater$aov.tab[6]$`Pr(>F)`[1]

allPValues.LoneIncube <- cbind(allPValues.LoneIncube, p.adjust(allPValues.LoneIncube,method = "fdr", n = 4))
colnames(allPValues.LoneIncube) <- c("P","fdr_adj")



# Now print everything
system("mkdir ./BETAPLOTS/UWUF/")
capture.output(ANOVA.UWUF.ExN, file = "./BETAPLOTS/UWUF/ANOVA.UWUF.ExN.txt")
# capture.output(pairwaiseAdonis.UWUF.ExN, file = "pairwaiseAdonis.UWUF.ExN.txt")
# capture.output(allPValues.ExN, file = "allPValues.pairwise.ExN.txt")

capture.output(ANOVA.UWUF.ExNWater, file = "./BETAPLOTS/UWUF/ANOVA.UWUF.ExNWater.txt")
# capture.output(pairwaiseAdonis.UWUF.ExNWater, file = "pairwaiseAdonis.UWUF.ExNWater.txt")
capture.output(allPValues.ExNWater, file = "./BETAPLOTS/UWUF/allPValues.ExNWater.txt")
capture.output(allPValues.ExN, file = "./BETAPLOTS/UWUF/allPValues.ExN.txt")
capture.output(allPValues.LoneIncube, file = "./BETAPLOTS/UWUF/allPValues.LoneIncube.txt")


capture.output(ANOVA.UWUF.LoneIncube, file = "./BETAPLOTS/UWUF/ANOVA.UWUF.LoneIncube.txt")
# capture.output(pairwaiseAdonis.UWUF.LoneIncube, file = "pairwaiseAdonis.UWUF.LoneIncube.txt")

####### PLOT ############
# EXN UWUF
MF.ExN$ColRep <- factor(MF.ExN$ColRep, levels = c('NereotestExNExN','NereotestNereoExN','NereotestMastExN','NereotestNereoMastExN'))
ExNColours <- c("black","green","red","brown")

# Make chulls

NMDS.UWUF.ExN.ExN <- NMDS.UWUF.ExN$points[grep(".ExN.", rownames(NMDS.UWUF.ExN$points), fixed = TRUE),]
NMDS.UWUF.ExN.ExN.chull <- chull(NMDS.UWUF.ExN.ExN)
NMDS.UWUF.ExN.ExN.chull <- c(NMDS.UWUF.ExN.ExN.chull, NMDS.UWUF.ExN.ExN.chull[1])

NMDS.UWUF.ExN.Mast <- NMDS.UWUF.ExN$points[grep(".Mast.", rownames(NMDS.UWUF.ExN$points), fixed = TRUE),]
NMDS.UWUF.ExN.Mast.chull <- chull(NMDS.UWUF.ExN.Mast)
NMDS.UWUF.ExN.Mast.chull <- c(NMDS.UWUF.ExN.Mast.chull, NMDS.UWUF.ExN.Mast.chull[1])

NMDS.UWUF.ExN.Nereo <- NMDS.UWUF.ExN$points[grep(".Nereo.", rownames(NMDS.UWUF.ExN$points), fixed = TRUE),]
NMDS.UWUF.ExN.Nereo.chull <- chull(NMDS.UWUF.ExN.Nereo)
NMDS.UWUF.ExN.Nereo.chull <- c(NMDS.UWUF.ExN.Nereo.chull, NMDS.UWUF.ExN.Nereo.chull[1])

NMDS.UWUF.ExN.NereoMast <- NMDS.UWUF.ExN$points[grep(".NereoMast.", rownames(NMDS.UWUF.ExN$points), fixed = TRUE),]
NMDS.UWUF.ExN.NereoMast.chull <- chull(NMDS.UWUF.ExN.NereoMast)
NMDS.UWUF.ExN.NereoMast.chull <- c(NMDS.UWUF.ExN.NereoMast.chull, NMDS.UWUF.ExN.NereoMast.chull[1])

pdf("./BETAPLOTS/UWUF/NMDS_UWUF_ExN.pdf", width = 10, height = 7, pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(NMDS.UWUF.ExN$points
     , main = "NMDS of Nereo Meristem Swabs"
     , pch = 19
     , col = ExNColours[factor(MF.ExN$ColRep)]
     , sub = round(NMDS.UWUF.ExN$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2)
lines(NMDS.UWUF.ExN.ExN[NMDS.UWUF.ExN.ExN.chull,]
      , col = "black")
lines(NMDS.UWUF.ExN.Mast[NMDS.UWUF.ExN.Mast.chull,]
      , col = "red")
lines(NMDS.UWUF.ExN.Nereo[NMDS.UWUF.ExN.Nereo.chull,]
      , col = "green")
lines(NMDS.UWUF.ExN.NereoMast[NMDS.UWUF.ExN.NereoMast.chull,]
      , col = "brown")
par(fig = c(0.6,1,0,1), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n")
legend("center"
       , pch = 19
       , legend = c("Alone","With Nereocystis","With Mastocarpus","With Both")#levels(MF.ExN$ColRep)
       , col = ExNColours
       , cex = 1)
dev.off()


# ExN UWUF
MF.ExNWater$ColRep <- factor(MF.ExNWater$ColRep, levels = c('NereotestH2OWater','NereotestExNWater','NereotestNereoWater','NereotestMastWater','NereotestNereoMastWater'))
ExNWaterColours <- c("blue","black","green","red","brown")


# Make chulls
NMDS.UWUF.ExNWater.Water <- NMDS.UWUF.ExNWater$points[grep(".H2O.", rownames(NMDS.UWUF.ExNWater$points), fixed = TRUE),]
NMDS.UWUF.ExNWater.Water.chull <- chull(NMDS.UWUF.ExNWater.Water)
NMDS.UWUF.ExNWater.Water.chull <- c(NMDS.UWUF.ExNWater.Water.chull, NMDS.UWUF.ExNWater.Water.chull[1])

NMDS.UWUF.ExNWater.ExN <- NMDS.UWUF.ExNWater$points[grep(".ExN.", rownames(NMDS.UWUF.ExNWater$points), fixed = TRUE),]
NMDS.UWUF.ExNWater.ExN.chull <- chull(NMDS.UWUF.ExNWater.ExN)
NMDS.UWUF.ExNWater.ExN.chull <- c(NMDS.UWUF.ExNWater.ExN.chull, NMDS.UWUF.ExNWater.ExN.chull[1])

NMDS.UWUF.ExNWater.Mast <- NMDS.UWUF.ExNWater$points[grep(".Mast.", rownames(NMDS.UWUF.ExNWater$points), fixed = TRUE),]
NMDS.UWUF.ExNWater.Mast.chull <- chull(NMDS.UWUF.ExNWater.Mast)
NMDS.UWUF.ExNWater.Mast.chull <- c(NMDS.UWUF.ExNWater.Mast.chull, NMDS.UWUF.ExNWater.Mast.chull[1])

NMDS.UWUF.ExNWater.Nereo <- NMDS.UWUF.ExNWater$points[grep(".Nereo.", rownames(NMDS.UWUF.ExNWater$points), fixed = TRUE),]
NMDS.UWUF.ExNWater.Nereo.chull <- chull(NMDS.UWUF.ExNWater.Nereo)
NMDS.UWUF.ExNWater.Nereo.chull <- c(NMDS.UWUF.ExNWater.Nereo.chull, NMDS.UWUF.ExNWater.Nereo.chull[1])

NMDS.UWUF.ExNWater.NereoMast <- NMDS.UWUF.ExNWater$points[grep(".NereoMast.", rownames(NMDS.UWUF.ExNWater$points), fixed = TRUE),]
NMDS.UWUF.ExNWater.NereoMast.chull <- chull(NMDS.UWUF.ExNWater.NereoMast)
NMDS.UWUF.ExNWater.NereoMast.chull <- c(NMDS.UWUF.ExNWater.NereoMast.chull, NMDS.UWUF.ExNWater.NereoMast.chull[1])


pdf("./BETAPLOTS/UWUF/NMDS_UWUF_ExNWater.pdf", width = 10, height = 7, pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(NMDS.UWUF.ExNWater$points
     , main = "NMDS of water samples"
     , pch = 19
     , col = ExNWaterColours[factor(MF.ExNWater$ColRep)]
     , sub = round(NMDS.UWUF.ExNWater$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2)
lines(NMDS.UWUF.ExNWater.Water[NMDS.UWUF.ExNWater.Water.chull,]
      , col = "blue")
lines(NMDS.UWUF.ExNWater.ExN[NMDS.UWUF.ExNWater.ExN.chull,]
      , col = "black")
lines(NMDS.UWUF.ExNWater.Nereo[NMDS.UWUF.ExNWater.Nereo.chull,]
      , col = "green")
lines(NMDS.UWUF.ExNWater.Mast[NMDS.UWUF.ExNWater.Mast.chull,]
      , col = "red")
lines(NMDS.UWUF.ExNWater.NereoMast[NMDS.UWUF.ExNWater.NereoMast.chull,]
      , col = "brown")
par(fig = c(0.6,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n")
legend("center"
       , pch = 19
       , legend = c("Water only","Nereo Meristem Only", "Nereocystis","Mastocarpus","Nereo + Mast")
       , col = ExNWaterColours
       , cex = 1)
dev.off()

# LoneIncube UWUF
MF.LoneIncube$ColRep <- factor(MF.LoneIncube$ColRep, levels = c('LoneincubeNereowater','LoneincubeNereoNereo','LoneincubeMastwater','LoneincubeMastMast'))
LoneIncubeColours <- c("blue","green","purple","red")

NMDS.UWUF.LoneIncube.WNereo <- NMDS.UWUF.LoneIncube$points[grep("water.Loneincube.Nereo.", rownames(NMDS.UWUF.LoneIncube$points), fixed = TRUE),]
NMDS.UWUF.LoneIncube.WNereo.chull <- chull(NMDS.UWUF.LoneIncube.WNereo)
NMDS.UWUF.LoneIncube.WNereo.chull <- c(NMDS.UWUF.LoneIncube.WNereo.chull, NMDS.UWUF.LoneIncube.WNereo.chull[1])

NMDS.UWUF.LoneIncube.WMast <- NMDS.UWUF.LoneIncube$points[grep("water.Loneincube.Mast.", rownames(NMDS.UWUF.LoneIncube$points), fixed = TRUE),]
NMDS.UWUF.LoneIncube.WMast.chull <- chull(NMDS.UWUF.LoneIncube.WMast)
NMDS.UWUF.LoneIncube.WMast.chull <- c(NMDS.UWUF.LoneIncube.WMast.chull, NMDS.UWUF.LoneIncube.WMast.chull[1])

NMDS.UWUF.LoneIncube.Nereo <- NMDS.UWUF.LoneIncube$points[grep("Nereo.Loneincube.", rownames(NMDS.UWUF.LoneIncube$points), fixed = TRUE),]
NMDS.UWUF.LoneIncube.Nereo.chull <- chull(NMDS.UWUF.LoneIncube.Nereo)
NMDS.UWUF.LoneIncube.Nereo.chull <- c(NMDS.UWUF.LoneIncube.Nereo.chull, NMDS.UWUF.LoneIncube.Nereo.chull[1])

NMDS.UWUF.LoneIncube.Mast <- NMDS.UWUF.LoneIncube$points[grep("Mast.Loneincube.", rownames(NMDS.UWUF.LoneIncube$points), fixed = TRUE),]
NMDS.UWUF.LoneIncube.Mast.chull <- chull(NMDS.UWUF.LoneIncube.Mast)
NMDS.UWUF.LoneIncube.Mast.chull <- c(NMDS.UWUF.LoneIncube.Mast.chull, NMDS.UWUF.LoneIncube.Mast.chull[1])


pdf("./BETAPLOTS/UWUF/NMDS_UWUF_LoneIncube.pdf", width = 10, height = 7, pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(NMDS.UWUF.LoneIncube$points
     , main = "NMDS of Incubation Experiment"
     , pch = 19
     , col = LoneIncubeColours[factor(MF.LoneIncube$ColRep)]
     , sub = round(NMDS.UWUF.LoneIncube$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2)
lines(NMDS.UWUF.LoneIncube.WNereo[NMDS.UWUF.LoneIncube.WNereo.chull,]
      , col = "blue")
lines(NMDS.UWUF.LoneIncube.WMast[NMDS.UWUF.LoneIncube.WMast.chull,]
      , col = "purple")
lines(NMDS.UWUF.LoneIncube.Nereo[NMDS.UWUF.LoneIncube.Nereo.chull,]
      , col = "green")
lines(NMDS.UWUF.LoneIncube.Mast[NMDS.UWUF.LoneIncube.Mast.chull,]
      , col = "red")
par(fig = c(0.6,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot(0,0
     , pch = ''
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = ''
     , ylab = ''
     , bty = 'n')
legend("center"
       , pch = 19
       , legend = c("Water-Nereocystis","Nereocystis","Water-Mastocarpus","Mastocarpus")#levels(MF.LoneIncube$ColRep)
       , col = LoneIncubeColours
       , cex = 1)
dev.off()


system("cd ..")

####### WUF #############
dm.WUF.ExN <- dm.WUF[sapply(rownames(MF.ExN), function(x) {
  grep(x, rownames(dm.WUF))}),sapply(rownames(MF.ExN), function(x) {
    grep(x, colnames(dm.WUF))})]

dm.WUF.ExNWater <- dm.WUF[sapply(rownames(MF.ExNWater), function(x) {
  grep(x, rownames(dm.WUF))}),sapply(rownames(MF.ExNWater), function(x) {
    grep(x, colnames(dm.WUF))})]

dm.WUF.LoneIncube <- dm.WUF[sapply(rownames(MF.LoneIncube), function(x) {
  grep(x, rownames(dm.WUF))}),sapply(rownames(MF.LoneIncube), function(x) {
    grep(x, colnames(dm.WUF))})]

NMDS.WUF.ExN <- isoMDS(as.matrix(dm.WUF.ExN), y = cmdscale(as.matrix(dm.WUF.ExN), 2))

NMDS.WUF.ExNWater <- isoMDS(as.matrix(dm.WUF.ExNWater), y = cmdscale(as.matrix(dm.WUF.ExNWater), 2))

NMDS.WUF.LoneIncube <- isoMDS(as.matrix(dm.WUF.LoneIncube), y = cmdscale(as.matrix(dm.WUF.LoneIncube), 2))

###### STATS ##########
# ExN
ANOVA.WUF.ExN <- adonis(dm.WUF.ExN ~ ColRep, data = MF.ExN)
ANOSIM.WUF.ExN <- anosim(dm.WUF.ExN, grouping = MF.ExN$ColRep)
# pairwaiseAdonis.WUF.ExN <- pairwise.adonis(dm.WUF.ExN, factors = MF.ExN$ColRep)



# ExN Split
dm.WUF.ExN.ExNvsNereo <- dm.WUF.ExN[grep("([.]ExN[.])|([.]Nereo[.])", rownames(dm.WUF.ExN)), grep("([.]ExN[.])|([.]Nereo[.])", colnames(dm.WUF.ExN))]
MF.ExN.ExNvsNereo <- MF.ExN[grep("([.]ExN[.])|([.]Nereo[.])", rownames(MF.ExN)),]
ANOVA.WUF.ExN.ExNvsNereo <- adonis(dm.WUF.ExN.ExNvsNereo ~ ColRep, data = MF.ExN.ExNvsNereo)
ANOSIM.WUF.ExN.ExNvsNereo <- anosim(dm.WUF.ExN.ExNvsNereo, grouping = MF.ExN.ExNvsNereo$ColRep)

dm.WUF.ExN.ExNvsMast <- dm.WUF.ExN[grep("([.]ExN[.])|([.]Mast[.])", rownames(dm.WUF.ExN)), grep("([.]ExN[.])|([.]Mast[.])", colnames(dm.WUF.ExN))]
MF.ExN.ExNvsMast <- MF.ExN[grep("([.]ExN[.])|([.]Mast[.])", rownames(MF.ExN)),]
ANOVA.WUF.ExN.ExNvsMast <- adonis(dm.WUF.ExN.ExNvsMast ~ ColRep, data = MF.ExN.ExNvsMast)
ANOSIM.WUF.ExN.ExNvsMast <- anosim(dm.WUF.ExN.ExNvsMast, grouping = MF.ExN.ExNvsMast$ColRep)


dm.WUF.ExN.ExNvsNereoMast <- dm.WUF.ExN[grep("([.]ExN[.])|([.]NereoMast[.])", rownames(dm.WUF.ExN)), grep("([.]ExN[.])|([.]NereoMast[.])", colnames(dm.WUF.ExN))]
MF.ExN.ExNvsNereoMast <- MF.ExN[grep("([.]ExN[.])|([.]NereoMast[.])", rownames(MF.ExN)),]
ANOVA.WUF.ExN.ExNvsNereoMast <- adonis(dm.WUF.ExN.ExNvsNereoMast ~ ColRep, data = MF.ExN.ExNvsNereoMast)
ANOSIM.WUF.ExN.ExNvsNereoMast <- anosim(dm.WUF.ExN.ExNvsNereoMast, grouping = MF.ExN.ExNvsNereoMast$ColRep)


dm.WUF.ExN.NereovsMast <- dm.WUF.ExN[grep("([.]Nereo[.])|([.]Mast[.])", rownames(dm.WUF.ExN)), grep("([.]Nereo[.])|([.]Mast[.])", colnames(dm.WUF.ExN))]
MF.ExN.NereovsMast <- MF.ExN[grep("([.]Nereo[.])|([.]Mast[.])", rownames(MF.ExN)),]
ANOVA.WUF.ExN.NereovsMast <- adonis(dm.WUF.ExN.NereovsMast ~ ColRep, data = MF.ExN.NereovsMast)
ANOSIM.WUF.ExN.NereovsMast <- anosim(dm.WUF.ExN.NereovsMast, grouping = MF.ExN.NereovsMast$ColRep)

dm.WUF.ExN.NereovsNereoMast <- dm.WUF.ExN[grep("([.]Nereo[.])|([.]NereoMast[.])", rownames(dm.WUF.ExN)), grep("([.]Nereo[.])|([.]NereoMast[.])", colnames(dm.WUF.ExN))]
MF.ExN.NereovsNereoMast <- MF.ExN[grep("([.]Nereo[.])|([.]NereoMast[.])", rownames(MF.ExN)),]
ANOVA.WUF.ExN.NereovsNereoMast <- adonis(dm.WUF.ExN.NereovsNereoMast ~ ColRep, data = MF.ExN.NereovsNereoMast)
ANOSIM.WUF.ExN.NereovsNereoMast <- anosim(dm.WUF.ExN.NereovsNereoMast, grouping = MF.ExN.NereovsNereoMast$ColRep)


dm.WUF.ExN.MastvsNereoMast <- dm.WUF.ExN[grep("([.]Mast[.])|([.]NereoMast[.])", rownames(dm.WUF.ExN)), grep("([.]Mast[.])|([.]NereoMast[.])", colnames(dm.WUF.ExN))]
MF.ExN.MastvsNereoMast <- MF.ExN[grep("([.]Mast[.])|([.]NereoMast[.])", rownames(MF.ExN)),]
ANOVA.WUF.ExN.MastvsNereoMast <- adonis(dm.WUF.ExN.MastvsNereoMast ~ ColRep, data = MF.ExN.MastvsNereoMast)
ANOSIM.WUF.ExN.MastvsNereoMast <- anosim(dm.WUF.ExN.MastvsNereoMast, grouping = MF.ExN.MastvsNereoMast$ColRep)

# 6 comparisons; 0.05/6 = 0.00625
allPValues.ExN <- matrix(nrow= 6)
rownames(allPValues.ExN) <- c("ExNvsNereo"
                              ,"ExNvsMast"
                              ,"ExNvsNereoMast"
                              ,"NereovsMast"
                              , "NereovsNereoMast"
                              , "MastvsNereoMast")
allPValues.ExN[1,1] <- ANOVA.WUF.ExN.ExNvsNereo$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExN[2,1] <- ANOVA.WUF.ExN.ExNvsMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExN[3,1] <- ANOVA.WUF.ExN.ExNvsNereoMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExN[4,1] <- ANOVA.WUF.ExN.NereovsMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExN[5,1] <- ANOVA.WUF.ExN.NereovsNereoMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExN[6,1] <- ANOVA.WUF.ExN.MastvsNereoMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExN <- cbind(allPValues.ExN, p.adjust(allPValues.ExN,method = "fdr", n = 6))
colnames(allPValues.ExN) <- c("P","fdr_adj")


# ExN Water
ANOVA.WUF.ExNWater <- adonis(dm.WUF.ExNWater ~ ColRep, data = MF.ExNWater)
ANOSIM.WUF.ExNWater <- anosim(dm.WUF.ExNWater, grouping = MF.ExNWater$ColRep) # NOT WORKING
# pairwaiseAdonis.WUF.ExNWater <- pairwise.adonis(dm.WUF.ExNWater, factors = MF.ExNWater$ColRep)

# ExN Water split
dm.WUF.ExNWater.ExNvsH2O <- dm.WUF.ExNWater[grep("([.]ExN[.])|([.]H2O[.])", rownames(dm.WUF.ExNWater)), grep("([.]ExN[.])|([.]H2O[.])", colnames(dm.WUF.ExNWater))]
MF.ExNWater.ExNvsH2O <- MF.ExNWater[grep("([.]ExN[.])|([.]H2O[.])", rownames(MF.ExNWater)),]
ANOVA.WUF.ExNWater.ExNvsH2O <- adonis(dm.WUF.ExNWater.ExNvsH2O ~ ColRep, data = MF.ExNWater.ExNvsH2O)
ANOSIM.WUF.ExNWater.ExNvsH2O <- anosim(dm.WUF.ExNWater.ExNvsH2O, grouping = MF.ExNWater.ExNvsH2O$ColRep)

dm.WUF.ExNWater.ExNvsNereo <- dm.WUF.ExNWater[grep("([.]ExN[.])|([.]Nereo[.])", rownames(dm.WUF.ExNWater)), grep("([.]ExN[.])|([.]Nereo[.])", colnames(dm.WUF.ExNWater))]
MF.ExNWater.ExNvsNereo <- MF.ExNWater[grep("([.]ExN[.])|([.]Nereo[.])", rownames(MF.ExNWater)),]
ANOVA.WUF.ExNWater.ExNvsNereo <- adonis(dm.WUF.ExNWater.ExNvsNereo ~ ColRep, data = MF.ExNWater.ExNvsNereo)
ANOSIM.WUF.ExNWater.ExNvsNereo <- anosim(dm.WUF.ExNWater.ExNvsNereo, grouping = MF.ExNWater.ExNvsNereo$ColRep)

dm.WUF.ExNWater.ExNvsMast <- dm.WUF.ExNWater[grep("([.]ExN[.])|([.]Mast[.])", rownames(dm.WUF.ExNWater)), grep("([.]ExN[.])|([.]Mast[.])", colnames(dm.WUF.ExNWater))]
MF.ExNWater.ExNvsMast <- MF.ExNWater[grep("([.]ExN[.])|([.]Mast[.])", rownames(MF.ExNWater)),]
ANOVA.WUF.ExNWater.ExNvsMast <- adonis(dm.WUF.ExNWater.ExNvsMast ~ ColRep, data = MF.ExNWater.ExNvsMast)
ANOSIM.WUF.ExNWater.ExNvsMast <- anosim(dm.WUF.ExNWater.ExNvsMast, grouping = MF.ExNWater.ExNvsMast$ColRep)

dm.WUF.ExNWater.ExNvsNereoMast <- dm.WUF.ExNWater[grep("([.]ExN[.])|([.]NereoMast[.])", rownames(dm.WUF.ExNWater)), grep("([.]ExN[.])|([.]NereoMast[.])", colnames(dm.WUF.ExNWater))]
MF.ExNWater.ExNvsNereoMast <- MF.ExNWater[grep("([.]ExN[.])|([.]NereoMast[.])", rownames(MF.ExNWater)),]
ANOVA.WUF.ExNWater.ExNvsNereoMast <- adonis(dm.WUF.ExNWater.ExNvsNereoMast ~ ColRep, data = MF.ExNWater.ExNvsNereoMast)
ANOSIM.WUF.ExNWater.ExNvsNereoMast <- anosim(dm.WUF.ExNWater.ExNvsNereoMast, grouping = MF.ExNWater.ExNvsNereoMast$ColRep)

dm.WUF.ExNWater.H2OvsNereo <- dm.WUF.ExNWater[grep("([.]H2O[.])|([.]Nereo[.])", rownames(dm.WUF.ExNWater)), grep("([.]H2O[.])|([.]Nereo[.])", colnames(dm.WUF.ExNWater))]
MF.ExNWater.H2OvsNereo <- MF.ExNWater[grep("([.]H2O[.])|([.]Nereo[.])", rownames(MF.ExNWater)),]
ANOVA.WUF.ExNWater.H2OvsNereo <- adonis(dm.WUF.ExNWater.H2OvsNereo ~ ColRep, data = MF.ExNWater.H2OvsNereo)
ANOSIM.WUF.ExNWater.H2OvsNereo <- anosim(dm.WUF.ExNWater.H2OvsNereo, grouping = MF.ExNWater.H2OvsNereo$ColRep)

dm.WUF.ExNWater.H2OvsMast <- dm.WUF.ExNWater[grep("([.]H2O[.])|([.]Mast[.])", rownames(dm.WUF.ExNWater)), grep("([.]H2O[.])|([.]Mast[.])", colnames(dm.WUF.ExNWater))]
MF.ExNWater.H2OvsMast <- MF.ExNWater[grep("([.]H2O[.])|([.]Mast[.])", rownames(MF.ExNWater)),]
ANOVA.WUF.ExNWater.H2OvsMast <- adonis(dm.WUF.ExNWater.H2OvsMast ~ ColRep, data = MF.ExNWater.H2OvsMast)
ANOSIM.WUF.ExNWater.H2OvsMast <- anosim(dm.WUF.ExNWater.H2OvsMast, grouping = MF.ExNWater.H2OvsMast$ColRep)

dm.WUF.ExNWater.H2OvsNereoMast <- dm.WUF.ExNWater[grep("([.]H2O[.])|([.]NereoMast[.])", rownames(dm.WUF.ExNWater)), grep("([.]H2O[.])|([.]Mast[.])", colnames(dm.WUF.ExNWater))]
MF.ExNWater.H2OvsNereoMast <- MF.ExNWater[grep("([.]H2O[.])|([.]NereoMast[.])", rownames(MF.ExNWater)),]
ANOVA.WUF.ExNWater.H2OvsNereoMast <- adonis(dm.WUF.ExNWater.H2OvsNereoMast ~ ColRep, data = MF.ExNWater.H2OvsNereoMast)
ANOSIM.WUF.ExNWater.H2OvsNereoMast <- anosim(dm.WUF.ExNWater.H2OvsNereoMast, grouping = MF.ExNWater.H2OvsNereoMast$ColRep)

dm.WUF.ExNWater.NereovsMast <- dm.WUF.ExNWater[grep("([.]Nereo[.])|([.]Mast[.])", rownames(dm.WUF.ExNWater)), grep("([.]Nereo[.])|([.]Mast[.])", colnames(dm.WUF.ExNWater))]
MF.ExNWater.NereovsMast <- MF.ExNWater[grep("([.]Nereo[.])|([.]Mast[.])", rownames(MF.ExNWater)),]
ANOVA.WUF.ExNWater.NereovsMast <- adonis(dm.WUF.ExNWater.NereovsMast ~ ColRep, data = MF.ExNWater.NereovsMast)
ANOSIM.WUF.ExNWater.NereovsMast <- anosim(dm.WUF.ExNWater.NereovsMast, grouping = MF.ExNWater.NereovsMast$ColRep)

# 8 comparisons; 0.05/8 = 0.00625
allPValues.ExNWater <- matrix(nrow= 8)
rownames(allPValues.ExNWater) <- c("ExNvsH2O"
                                   ,"ExNvsNereo"
                                   ,"ExNvsMast"
                                   ,"ExNvsNereoMast"
                                   ,"H2OvsNereo"
                                   ,"H2OvsMast"
                                   ,"H2OvsNereoMast"
                                   ,"NereovsMast")
allPValues.ExNWater[1,1] <- ANOVA.WUF.ExNWater.ExNvsH2O$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[2,1] <- ANOVA.WUF.ExNWater.ExNvsNereo$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[3,1] <- ANOVA.WUF.ExNWater.ExNvsMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[4,1] <- ANOVA.WUF.ExNWater.ExNvsNereoMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[5,1] <- ANOVA.WUF.ExNWater.H2OvsNereo$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[6,1] <- ANOVA.WUF.ExNWater.H2OvsMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[7,1] <- ANOVA.WUF.ExNWater.H2OvsNereoMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[8,1] <- ANOVA.WUF.ExNWater.NereovsMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater <- cbind(allPValues.ExNWater, p.adjust(allPValues.ExNWater,method = "fdr", n = 10))
colnames(allPValues.ExNWater) <- c("P","fdr_adj")

# Lone Incube
# Add another column
MF.LoneIncube$WorS <- MF.LoneIncube$Substrate
MF.LoneIncube$WorS <- gsub("(Nereo)|(Mast)", "Seaweed", MF.LoneIncube[,'WorS'])

ANOVA.WUF.LoneIncube <- adonis(dm.WUF.LoneIncube ~ WorS*Treatment, data = MF.LoneIncube)
ANOSIM.WUF.LoneIncube <- anosim(dm.WUF.LoneIncube, grouping = MF.LoneIncube$ColRep) # NOT WORKING
# pairwaiseAdonis.WUF.LoneIncube <- pairwise.adonis(dm.WUF.LoneIncube, factors = MF.LoneIncube$ColRep)


# Lone Incube split
dm.WUF.LoneIncube.NereovsMast <- dm.WUF.LoneIncube[grep("(^Nereo[.])|(^Mast[.])", rownames(dm.WUF.LoneIncube)), grep("(^Nereo[.])|(^Mast[.])", colnames(dm.WUF.LoneIncube))]
MF.LoneIncube.NereovsMast <- MF.LoneIncube[grep("(^Nereo[.])|(^Mast[.])", rownames(MF.LoneIncube)),]
ANOVA.WUF.LoneIncube.NereovsMast <- adonis(dm.WUF.LoneIncube.NereovsMast ~ ColRep, data = MF.LoneIncube.NereovsMast)
ANOSIM.WUF.LoneIncube.NereovsMast <- anosim(dm.WUF.LoneIncube.NereovsMast, grouping = MF.LoneIncube.NereovsMast$ColRep)

dm.WUF.LoneIncube.waters <- dm.WUF.LoneIncube[grep("^water[.]", rownames(dm.WUF.LoneIncube)), grep("^water[.]", colnames(dm.WUF.LoneIncube))]
MF.LoneIncube.waters <- MF.LoneIncube[grep("^water[.]", rownames(MF.LoneIncube)),]
ANOVA.WUF.LoneIncube.waters <- adonis(dm.WUF.LoneIncube.waters ~ ColRep, data = MF.LoneIncube.waters)
ANOSIM.WUF.LoneIncube.waters <- anosim(dm.WUF.LoneIncube.waters, grouping = MF.LoneIncube.waters$ColRep)

dm.WUF.LoneIncube.Nereovswater <- dm.WUF.LoneIncube[grep("[.]Nereo[.]", rownames(dm.WUF.LoneIncube)), grep("[.]Nereo[.]", colnames(dm.WUF.LoneIncube))]
MF.LoneIncube.Nereovswater <- MF.LoneIncube[grep("[.]Nereo[.]", rownames(MF.LoneIncube)),]
ANOVA.WUF.LoneIncube.Nereovswater <- adonis(dm.WUF.LoneIncube.Nereovswater ~ ColRep, data = MF.LoneIncube.Nereovswater)
ANOSIM.WUF.LoneIncube.Nereovswater <- anosim(dm.WUF.LoneIncube.Nereovswater, grouping = MF.LoneIncube.Nereovswater$ColRep)

dm.WUF.LoneIncube.Mastvswater <- dm.WUF.LoneIncube[grep("[.]Mast[.]", rownames(dm.WUF.LoneIncube)), grep("[.]Mast[.]", colnames(dm.WUF.LoneIncube))]
MF.LoneIncube.Mastvswater <- MF.LoneIncube[grep("[.]Mast[.]", rownames(MF.LoneIncube)),]
ANOVA.WUF.LoneIncube.Mastvswater <- adonis(dm.WUF.LoneIncube.Mastvswater ~ ColRep, data = MF.LoneIncube.Mastvswater)
ANOSIM.WUF.LoneIncube.Mastvswater <- anosim(dm.WUF.LoneIncube.Mastvswater, grouping = MF.LoneIncube.Mastvswater$ColRep)

# 4 comparisons; 0.05/4 = 0.00625
allPValues.LoneIncube <- matrix(nrow= 4)
rownames(allPValues.LoneIncube) <- c("NereovsMast"
                                     ,"waters"
                                     ,"Nereovswater"
                                     ,"Mastvswater"
)
allPValues.LoneIncube[1,1] <- ANOVA.WUF.LoneIncube.NereovsMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.LoneIncube[2,1] <- ANOVA.WUF.LoneIncube.waters$aov.tab[6]$`Pr(>F)`[1]
allPValues.LoneIncube[3,1] <- ANOVA.WUF.LoneIncube.Nereovswater$aov.tab[6]$`Pr(>F)`[1]
allPValues.LoneIncube[4,1] <- ANOVA.WUF.LoneIncube.Mastvswater$aov.tab[6]$`Pr(>F)`[1]

allPValues.LoneIncube <- cbind(allPValues.LoneIncube, p.adjust(allPValues.LoneIncube,method = "fdr", n = 4))
colnames(allPValues.LoneIncube) <- c("P","fdr_adj")

# Now print everything
system("mkdir ./BETAPLOTS/WUF/")

capture.output(ANOVA.WUF.ExN, file = "./BETAPLOTS/WUF/ANOVA.WUF.ExN.txt")
# capture.output(pairwaiseAdonis.WUF.ExN, file = "pairwaiseAdonis.WUF.ExN.txt")
# capture.output(allPValues.ExN, file = "allPValues.pairwise.ExN.txt")

capture.output(ANOVA.WUF.ExNWater, file = "./BETAPLOTS/WUF/ANOVA.WUF.ExNWater.txt")
# capture.output(pairwaiseAdonis.WUF.ExNWater, file = "pairwaiseAdonis.WUF.ExNWater.txt")
capture.output(allPValues.ExNWater, file = "./BETAPLOTS/WUF/allPValues.ExNWater.txt")
capture.output(allPValues.ExN, file = "./BETAPLOTS/WUF/allPValues.ExN.txt")
capture.output(allPValues.LoneIncube, file = "./BETAPLOTS/WUF/allPValues.LoneIncube.txt")

capture.output(ANOVA.WUF.LoneIncube, file = "./BETAPLOTS/WUF/ANOVA.WUF.LoneIncube.txt")
# capture.output(pairwaiseAdonis.WUF.LoneIncube, file = "pairwaiseAdonis.WUF.LoneIncube.txt")



####### PLOT ############
# EXN WUF
MF.ExN$ColRep <- factor(MF.ExN$ColRep, levels = c('NereotestExNExN','NereotestNereoExN','NereotestMastExN','NereotestNereoMastExN'))
ExNColours <- c("black","green","red","brown")

# Make chulls

NMDS.WUF.ExN.ExN <- NMDS.WUF.ExN$points[grep(".ExN.", rownames(NMDS.WUF.ExN$points), fixed = TRUE),]
NMDS.WUF.ExN.ExN.chull <- chull(NMDS.WUF.ExN.ExN)
NMDS.WUF.ExN.ExN.chull <- c(NMDS.WUF.ExN.ExN.chull, NMDS.WUF.ExN.ExN.chull[1])

NMDS.WUF.ExN.Mast <- NMDS.WUF.ExN$points[grep(".Mast.", rownames(NMDS.WUF.ExN$points), fixed = TRUE),]
NMDS.WUF.ExN.Mast.chull <- chull(NMDS.WUF.ExN.Mast)
NMDS.WUF.ExN.Mast.chull <- c(NMDS.WUF.ExN.Mast.chull, NMDS.WUF.ExN.Mast.chull[1])

NMDS.WUF.ExN.Nereo <- NMDS.WUF.ExN$points[grep(".Nereo.", rownames(NMDS.WUF.ExN$points), fixed = TRUE),]
NMDS.WUF.ExN.Nereo.chull <- chull(NMDS.WUF.ExN.Nereo)
NMDS.WUF.ExN.Nereo.chull <- c(NMDS.WUF.ExN.Nereo.chull, NMDS.WUF.ExN.Nereo.chull[1])

NMDS.WUF.ExN.NereoMast <- NMDS.WUF.ExN$points[grep(".NereoMast.", rownames(NMDS.WUF.ExN$points), fixed = TRUE),]
NMDS.WUF.ExN.NereoMast.chull <- chull(NMDS.WUF.ExN.NereoMast)
NMDS.WUF.ExN.NereoMast.chull <- c(NMDS.WUF.ExN.NereoMast.chull, NMDS.WUF.ExN.NereoMast.chull[1])

pdf("./BETAPLOTS/WUF/NMDS_WUF_ExN.pdf", width = 10, height = 7, pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(NMDS.WUF.ExN$points
     , main = "NMDS of Nereo Meristem Swabs"
     , pch = 19
     , col = ExNColours[factor(MF.ExN$ColRep)]
     , sub = round(NMDS.WUF.ExN$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2)
lines(NMDS.WUF.ExN.ExN[NMDS.WUF.ExN.ExN.chull,]
      , col = "black")
lines(NMDS.WUF.ExN.Mast[NMDS.WUF.ExN.Mast.chull,]
      , col = "red")
lines(NMDS.WUF.ExN.Nereo[NMDS.WUF.ExN.Nereo.chull,]
      , col = "green")
lines(NMDS.WUF.ExN.NereoMast[NMDS.WUF.ExN.NereoMast.chull,]
      , col = "brown")
par(fig = c(0.6,1,0,1), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n")
legend("center"
       , pch = 19
       , legend = c("Alone","With Nereocystis","With Mastocarpus","With Both")#levels(MF.ExN$ColRep)
       , col = ExNColours
       , cex = 1)
dev.off()


# ExN WUF
MF.ExNWater$ColRep <- factor(MF.ExNWater$ColRep, levels = c('NereotestH2OWater','NereotestExNWater','NereotestNereoWater','NereotestMastWater','NereotestNereoMastWater'))
ExNWaterColours <- c("blue","black","green","red","brown")


# Make chulls
NMDS.WUF.ExNWater.Water <- NMDS.WUF.ExNWater$points[grep(".H2O.", rownames(NMDS.WUF.ExNWater$points), fixed = TRUE),]
NMDS.WUF.ExNWater.Water.chull <- chull(NMDS.WUF.ExNWater.Water)
NMDS.WUF.ExNWater.Water.chull <- c(NMDS.WUF.ExNWater.Water.chull, NMDS.WUF.ExNWater.Water.chull[1])

NMDS.WUF.ExNWater.ExN <- NMDS.WUF.ExNWater$points[grep(".ExN.", rownames(NMDS.WUF.ExNWater$points), fixed = TRUE),]
NMDS.WUF.ExNWater.ExN.chull <- chull(NMDS.WUF.ExNWater.ExN)
NMDS.WUF.ExNWater.ExN.chull <- c(NMDS.WUF.ExNWater.ExN.chull, NMDS.WUF.ExNWater.ExN.chull[1])

NMDS.WUF.ExNWater.Mast <- NMDS.WUF.ExNWater$points[grep(".Mast.", rownames(NMDS.WUF.ExNWater$points), fixed = TRUE),]
NMDS.WUF.ExNWater.Mast.chull <- chull(NMDS.WUF.ExNWater.Mast)
NMDS.WUF.ExNWater.Mast.chull <- c(NMDS.WUF.ExNWater.Mast.chull, NMDS.WUF.ExNWater.Mast.chull[1])

NMDS.WUF.ExNWater.Nereo <- NMDS.WUF.ExNWater$points[grep(".Nereo.", rownames(NMDS.WUF.ExNWater$points), fixed = TRUE),]
NMDS.WUF.ExNWater.Nereo.chull <- chull(NMDS.WUF.ExNWater.Nereo)
NMDS.WUF.ExNWater.Nereo.chull <- c(NMDS.WUF.ExNWater.Nereo.chull, NMDS.WUF.ExNWater.Nereo.chull[1])

NMDS.WUF.ExNWater.NereoMast <- NMDS.WUF.ExNWater$points[grep(".NereoMast.", rownames(NMDS.WUF.ExNWater$points), fixed = TRUE),]
NMDS.WUF.ExNWater.NereoMast.chull <- chull(NMDS.WUF.ExNWater.NereoMast)
NMDS.WUF.ExNWater.NereoMast.chull <- c(NMDS.WUF.ExNWater.NereoMast.chull, NMDS.WUF.ExNWater.NereoMast.chull[1])


pdf("./BETAPLOTS/WUF/NMDS_WUF_ExNWater.pdf", width = 10, height = 7, pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(NMDS.WUF.ExNWater$points
     , main = "NMDS of water samples"
     , pch = 19
     , col = ExNWaterColours[factor(MF.ExNWater$ColRep)]
     , sub = round(NMDS.WUF.ExNWater$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2)
lines(NMDS.WUF.ExNWater.Water[NMDS.WUF.ExNWater.Water.chull,]
      , col = "blue")
lines(NMDS.WUF.ExNWater.ExN[NMDS.WUF.ExNWater.ExN.chull,]
      , col = "black")
lines(NMDS.WUF.ExNWater.Nereo[NMDS.WUF.ExNWater.Nereo.chull,]
      , col = "green")
lines(NMDS.WUF.ExNWater.Mast[NMDS.WUF.ExNWater.Mast.chull,]
      , col = "red")
lines(NMDS.WUF.ExNWater.NereoMast[NMDS.WUF.ExNWater.NereoMast.chull,]
      , col = "brown")
par(fig = c(0.6,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n")
legend("center"
       , pch = 19
       , legend = c("Water only","Nereo Meristem Only", "Nereocystis","Mastocarpus","Nereo + Mast")
       , col = ExNWaterColours
       , cex = 1)
dev.off()

# LoneIncube WUF
MF.LoneIncube$ColRep <- factor(MF.LoneIncube$ColRep, levels = c('LoneincubeNereowater','LoneincubeNereoNereo','LoneincubeMastwater','LoneincubeMastMast'))
LoneIncubeColours <- c("blue","green","purple","red")

NMDS.WUF.LoneIncube.WNereo <- NMDS.WUF.LoneIncube$points[grep("water.Loneincube.Nereo.", rownames(NMDS.WUF.LoneIncube$points), fixed = TRUE),]
NMDS.WUF.LoneIncube.WNereo.chull <- chull(NMDS.WUF.LoneIncube.WNereo)
NMDS.WUF.LoneIncube.WNereo.chull <- c(NMDS.WUF.LoneIncube.WNereo.chull, NMDS.WUF.LoneIncube.WNereo.chull[1])

NMDS.WUF.LoneIncube.WMast <- NMDS.WUF.LoneIncube$points[grep("water.Loneincube.Mast.", rownames(NMDS.WUF.LoneIncube$points), fixed = TRUE),]
NMDS.WUF.LoneIncube.WMast.chull <- chull(NMDS.WUF.LoneIncube.WMast)
NMDS.WUF.LoneIncube.WMast.chull <- c(NMDS.WUF.LoneIncube.WMast.chull, NMDS.WUF.LoneIncube.WMast.chull[1])

NMDS.WUF.LoneIncube.Nereo <- NMDS.WUF.LoneIncube$points[grep("Nereo.Loneincube.", rownames(NMDS.WUF.LoneIncube$points), fixed = TRUE),]
NMDS.WUF.LoneIncube.Nereo.chull <- chull(NMDS.WUF.LoneIncube.Nereo)
NMDS.WUF.LoneIncube.Nereo.chull <- c(NMDS.WUF.LoneIncube.Nereo.chull, NMDS.WUF.LoneIncube.Nereo.chull[1])

NMDS.WUF.LoneIncube.Mast <- NMDS.WUF.LoneIncube$points[grep("Mast.Loneincube.", rownames(NMDS.WUF.LoneIncube$points), fixed = TRUE),]
NMDS.WUF.LoneIncube.Mast.chull <- chull(NMDS.WUF.LoneIncube.Mast)
NMDS.WUF.LoneIncube.Mast.chull <- c(NMDS.WUF.LoneIncube.Mast.chull, NMDS.WUF.LoneIncube.Mast.chull[1])


pdf("./BETAPLOTS/WUF/NMDS_WUF_LoneIncube.pdf", width = 10, height = 7, pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(NMDS.WUF.LoneIncube$points
     , main = "NMDS of Incubation Experiment"
     , pch = 19
     , col = LoneIncubeColours[factor(MF.LoneIncube$ColRep)]
     , sub = round(NMDS.WUF.LoneIncube$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2)
lines(NMDS.WUF.LoneIncube.WNereo[NMDS.WUF.LoneIncube.WNereo.chull,]
      , col = "blue")
lines(NMDS.WUF.LoneIncube.WMast[NMDS.WUF.LoneIncube.WMast.chull,]
      , col = "purple")
lines(NMDS.WUF.LoneIncube.Nereo[NMDS.WUF.LoneIncube.Nereo.chull,]
      , col = "green")
lines(NMDS.WUF.LoneIncube.Mast[NMDS.WUF.LoneIncube.Mast.chull,]
      , col = "red")
par(fig = c(0.6,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot(0,0
     , pch = ''
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = ''
     , ylab = ''
     , bty = 'n')
legend("center"
       , pch = 19
       , legend = c("Water-Nereocystis","Nereocystis","Water-Mastocarpus","Mastocarpus")#levels(MF.LoneIncube$ColRep)
       , col = LoneIncubeColours
       , cex = 1)
dev.off()



system("cd ..")



####### BRAYCURTIS #############
dm.BC.ExN <- dm.BC[sapply(rownames(MF.ExN), function(x) {
  grep(x, rownames(dm.BC))}),sapply(rownames(MF.ExN), function(x) {
    grep(x, colnames(dm.BC))})]

dm.BC.ExNWater <- dm.BC[sapply(rownames(MF.ExNWater), function(x) {
  grep(x, rownames(dm.BC))}),sapply(rownames(MF.ExNWater), function(x) {
    grep(x, colnames(dm.BC))})]

dm.BC.LoneIncube <- dm.BC[sapply(rownames(MF.LoneIncube), function(x) {
  grep(x, rownames(dm.BC))}),sapply(rownames(MF.LoneIncube), function(x) {
    grep(x, colnames(dm.BC))})]

NMDS.BC.ExN <- isoMDS(as.matrix(dm.BC.ExN), y = cmdscale(as.matrix(dm.BC.ExN), 2))

NMDS.BC.ExNWater <- isoMDS(as.matrix(dm.BC.ExNWater), y = cmdscale(as.matrix(dm.BC.ExNWater), 2))

NMDS.BC.LoneIncube <- isoMDS(as.matrix(dm.BC.LoneIncube), y = cmdscale(as.matrix(dm.BC.LoneIncube), 2))

###### STATS ##########
# ExN
ANOVA.BC.ExN <- adonis(dm.BC.ExN ~ ColRep, data = MF.ExN)
ANOSIM.BC.ExN <- anosim(dm.BC.ExN, grouping = MF.ExN$ColRep)
# pairwaiseAdonis.BC.ExN <- pairwise.adonis(dm.BC.ExN, factors = MF.ExN$ColRep)

# ExN Split
dm.BC.ExN.ExNvsNereo <- dm.BC.ExN[grep("([.]ExN[.])|([.]Nereo[.])", rownames(dm.BC.ExN)), grep("([.]ExN[.])|([.]Nereo[.])", colnames(dm.BC.ExN))]
MF.ExN.ExNvsNereo <- MF.ExN[grep("([.]ExN[.])|([.]Nereo[.])", rownames(MF.ExN)),]
ANOVA.BC.ExN.ExNvsNereo <- adonis(dm.BC.ExN.ExNvsNereo ~ ColRep, data = MF.ExN.ExNvsNereo)
ANOSIM.BC.ExN.ExNvsNereo <- anosim(dm.BC.ExN.ExNvsNereo, grouping = MF.ExN.ExNvsNereo$ColRep)

dm.BC.ExN.ExNvsMast <- dm.BC.ExN[grep("([.]ExN[.])|([.]Mast[.])", rownames(dm.BC.ExN)), grep("([.]ExN[.])|([.]Mast[.])", colnames(dm.BC.ExN))]
MF.ExN.ExNvsMast <- MF.ExN[grep("([.]ExN[.])|([.]Mast[.])", rownames(MF.ExN)),]
ANOVA.BC.ExN.ExNvsMast <- adonis(dm.BC.ExN.ExNvsMast ~ ColRep, data = MF.ExN.ExNvsMast)
ANOSIM.BC.ExN.ExNvsMast <- anosim(dm.BC.ExN.ExNvsMast, grouping = MF.ExN.ExNvsMast$ColRep)


dm.BC.ExN.ExNvsNereoMast <- dm.BC.ExN[grep("([.]ExN[.])|([.]NereoMast[.])", rownames(dm.BC.ExN)), grep("([.]ExN[.])|([.]NereoMast[.])", colnames(dm.BC.ExN))]
MF.ExN.ExNvsNereoMast <- MF.ExN[grep("([.]ExN[.])|([.]NereoMast[.])", rownames(MF.ExN)),]
ANOVA.BC.ExN.ExNvsNereoMast <- adonis(dm.BC.ExN.ExNvsNereoMast ~ ColRep, data = MF.ExN.ExNvsNereoMast)
ANOSIM.BC.ExN.ExNvsNereoMast <- anosim(dm.BC.ExN.ExNvsNereoMast, grouping = MF.ExN.ExNvsNereoMast$ColRep)


dm.BC.ExN.NereovsMast <- dm.BC.ExN[grep("([.]Nereo[.])|([.]Mast[.])", rownames(dm.BC.ExN)), grep("([.]Nereo[.])|([.]Mast[.])", colnames(dm.BC.ExN))]
MF.ExN.NereovsMast <- MF.ExN[grep("([.]Nereo[.])|([.]Mast[.])", rownames(MF.ExN)),]
ANOVA.BC.ExN.NereovsMast <- adonis(dm.BC.ExN.NereovsMast ~ ColRep, data = MF.ExN.NereovsMast)
ANOSIM.BC.ExN.NereovsMast <- anosim(dm.BC.ExN.NereovsMast, grouping = MF.ExN.NereovsMast$ColRep)

dm.BC.ExN.NereovsNereoMast <- dm.BC.ExN[grep("([.]Nereo[.])|([.]NereoMast[.])", rownames(dm.BC.ExN)), grep("([.]Nereo[.])|([.]NereoMast[.])", colnames(dm.BC.ExN))]
MF.ExN.NereovsNereoMast <- MF.ExN[grep("([.]Nereo[.])|([.]NereoMast[.])", rownames(MF.ExN)),]
ANOVA.BC.ExN.NereovsNereoMast <- adonis(dm.BC.ExN.NereovsNereoMast ~ ColRep, data = MF.ExN.NereovsNereoMast)
ANOSIM.BC.ExN.NereovsNereoMast <- anosim(dm.BC.ExN.NereovsNereoMast, grouping = MF.ExN.NereovsNereoMast$ColRep)


dm.BC.ExN.MastvsNereoMast <- dm.BC.ExN[grep("([.]Mast[.])|([.]NereoMast[.])", rownames(dm.BC.ExN)), grep("([.]Mast[.])|([.]NereoMast[.])", colnames(dm.BC.ExN))]
MF.ExN.MastvsNereoMast <- MF.ExN[grep("([.]Mast[.])|([.]NereoMast[.])", rownames(MF.ExN)),]
ANOVA.BC.ExN.MastvsNereoMast <- adonis(dm.BC.ExN.MastvsNereoMast ~ ColRep, data = MF.ExN.MastvsNereoMast)
ANOSIM.BC.ExN.MastvsNereoMast <- anosim(dm.BC.ExN.MastvsNereoMast, grouping = MF.ExN.MastvsNereoMast$ColRep)

# 6 comparisons; 0.05/6 = 0.00625
allPValues.ExN <- matrix(nrow= 6)
rownames(allPValues.ExN) <- c("ExNvsNereo"
                          ,"ExNvsMast"
                          ,"ExNvsNereoMast"
                          ,"NereovsMast"
                          , "NereovsNereoMast"
                          , "MastvsNereoMast")
allPValues.ExN[1,1] <- ANOVA.BC.ExN.ExNvsNereo$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExN[2,1] <- ANOVA.BC.ExN.ExNvsMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExN[3,1] <- ANOVA.BC.ExN.ExNvsNereoMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExN[4,1] <- ANOVA.BC.ExN.NereovsMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExN[5,1] <- ANOVA.BC.ExN.NereovsNereoMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExN[6,1] <- ANOVA.BC.ExN.MastvsNereoMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExN <- cbind(allPValues.ExN, p.adjust(allPValues.ExN,method = "fdr", n = 6))
colnames(allPValues.ExN) <- c("P","fdr_adj")


# ExN Water
ANOVA.BC.ExNWater <- adonis(dm.BC.ExNWater ~ ColRep, data = MF.ExNWater)
ANOSIM.BC.ExNWater <- anosim(dm.BC.ExNWater, grouping = MF.ExNWater$ColRep) # NOT WORKING
# pairwaiseAdonis.BC.ExNWater <- pairwise.adonis(dm.BC.ExNWater, factors = MF.ExNWater$ColRep)

# ExN Water split
dm.BC.ExNWater.ExNvsH2O <- dm.BC.ExNWater[grep("([.]ExN[.])|([.]H2O[.])", rownames(dm.BC.ExNWater)), grep("([.]ExN[.])|([.]H2O[.])", colnames(dm.BC.ExNWater))]
MF.ExNWater.ExNvsH2O <- MF.ExNWater[grep("([.]ExN[.])|([.]H2O[.])", rownames(MF.ExNWater)),]
ANOVA.BC.ExNWater.ExNvsH2O <- adonis(dm.BC.ExNWater.ExNvsH2O ~ ColRep, data = MF.ExNWater.ExNvsH2O)
ANOSIM.BC.ExNWater.ExNvsH2O <- anosim(dm.BC.ExNWater.ExNvsH2O, grouping = MF.ExNWater.ExNvsH2O$ColRep)

dm.BC.ExNWater.ExNvsNereo <- dm.BC.ExNWater[grep("([.]ExN[.])|([.]Nereo[.])", rownames(dm.BC.ExNWater)), grep("([.]ExN[.])|([.]Nereo[.])", colnames(dm.BC.ExNWater))]
MF.ExNWater.ExNvsNereo <- MF.ExNWater[grep("([.]ExN[.])|([.]Nereo[.])", rownames(MF.ExNWater)),]
ANOVA.BC.ExNWater.ExNvsNereo <- adonis(dm.BC.ExNWater.ExNvsNereo ~ ColRep, data = MF.ExNWater.ExNvsNereo)
ANOSIM.BC.ExNWater.ExNvsNereo <- anosim(dm.BC.ExNWater.ExNvsNereo, grouping = MF.ExNWater.ExNvsNereo$ColRep)

dm.BC.ExNWater.ExNvsMast <- dm.BC.ExNWater[grep("([.]ExN[.])|([.]Mast[.])", rownames(dm.BC.ExNWater)), grep("([.]ExN[.])|([.]Mast[.])", colnames(dm.BC.ExNWater))]
MF.ExNWater.ExNvsMast <- MF.ExNWater[grep("([.]ExN[.])|([.]Mast[.])", rownames(MF.ExNWater)),]
ANOVA.BC.ExNWater.ExNvsMast <- adonis(dm.BC.ExNWater.ExNvsMast ~ ColRep, data = MF.ExNWater.ExNvsMast)
ANOSIM.BC.ExNWater.ExNvsMast <- anosim(dm.BC.ExNWater.ExNvsMast, grouping = MF.ExNWater.ExNvsMast$ColRep)

dm.BC.ExNWater.ExNvsNereoMast <- dm.BC.ExNWater[grep("([.]ExN[.])|([.]NereoMast[.])", rownames(dm.BC.ExNWater)), grep("([.]ExN[.])|([.]NereoMast[.])", colnames(dm.BC.ExNWater))]
MF.ExNWater.ExNvsNereoMast <- MF.ExNWater[grep("([.]ExN[.])|([.]NereoMast[.])", rownames(MF.ExNWater)),]
ANOVA.BC.ExNWater.ExNvsNereoMast <- adonis(dm.BC.ExNWater.ExNvsNereoMast ~ ColRep, data = MF.ExNWater.ExNvsNereoMast)
ANOSIM.BC.ExNWater.ExNvsNereoMast <- anosim(dm.BC.ExNWater.ExNvsNereoMast, grouping = MF.ExNWater.ExNvsNereoMast$ColRep)

dm.BC.ExNWater.H2OvsNereo <- dm.BC.ExNWater[grep("([.]H2O[.])|([.]Nereo[.])", rownames(dm.BC.ExNWater)), grep("([.]H2O[.])|([.]Nereo[.])", colnames(dm.BC.ExNWater))]
MF.ExNWater.H2OvsNereo <- MF.ExNWater[grep("([.]H2O[.])|([.]Nereo[.])", rownames(MF.ExNWater)),]
ANOVA.BC.ExNWater.H2OvsNereo <- adonis(dm.BC.ExNWater.H2OvsNereo ~ ColRep, data = MF.ExNWater.H2OvsNereo)
ANOSIM.BC.ExNWater.H2OvsNereo <- anosim(dm.BC.ExNWater.H2OvsNereo, grouping = MF.ExNWater.H2OvsNereo$ColRep)

dm.BC.ExNWater.H2OvsMast <- dm.BC.ExNWater[grep("([.]H2O[.])|([.]Mast[.])", rownames(dm.BC.ExNWater)), grep("([.]H2O[.])|([.]Mast[.])", colnames(dm.BC.ExNWater))]
MF.ExNWater.H2OvsMast <- MF.ExNWater[grep("([.]H2O[.])|([.]Mast[.])", rownames(MF.ExNWater)),]
ANOVA.BC.ExNWater.H2OvsMast <- adonis(dm.BC.ExNWater.H2OvsMast ~ ColRep, data = MF.ExNWater.H2OvsMast)
ANOSIM.BC.ExNWater.H2OvsMast <- anosim(dm.BC.ExNWater.H2OvsMast, grouping = MF.ExNWater.H2OvsMast$ColRep)

dm.BC.ExNWater.H2OvsNereoMast <- dm.BC.ExNWater[grep("([.]H2O[.])|([.]NereoMast[.])", rownames(dm.BC.ExNWater)), grep("([.]H2O[.])|([.]Mast[.])", colnames(dm.BC.ExNWater))]
MF.ExNWater.H2OvsNereoMast <- MF.ExNWater[grep("([.]H2O[.])|([.]NereoMast[.])", rownames(MF.ExNWater)),]
ANOVA.BC.ExNWater.H2OvsNereoMast <- adonis(dm.BC.ExNWater.H2OvsNereoMast ~ ColRep, data = MF.ExNWater.H2OvsNereoMast)
ANOSIM.BC.ExNWater.H2OvsNereoMast <- anosim(dm.BC.ExNWater.H2OvsNereoMast, grouping = MF.ExNWater.H2OvsNereoMast$ColRep)

dm.BC.ExNWater.NereovsMast <- dm.BC.ExNWater[grep("([.]Nereo[.])|([.]Mast[.])", rownames(dm.BC.ExNWater)), grep("([.]Nereo[.])|([.]Mast[.])", colnames(dm.BC.ExNWater))]
MF.ExNWater.NereovsMast <- MF.ExNWater[grep("([.]Nereo[.])|([.]Mast[.])", rownames(MF.ExNWater)),]
ANOVA.BC.ExNWater.NereovsMast <- adonis(dm.BC.ExNWater.NereovsMast ~ ColRep, data = MF.ExNWater.NereovsMast)
ANOSIM.BC.ExNWater.NereovsMast <- anosim(dm.BC.ExNWater.NereovsMast, grouping = MF.ExNWater.NereovsMast$ColRep)

# 8 comparisons; 0.05/8 = 0.00625
allPValues.ExNWater <- matrix(nrow= 8)
rownames(allPValues.ExNWater) <- c("ExNvsH2O"
                          ,"ExNvsNereo"
                          ,"ExNvsMast"
                          ,"ExNvsNereoMast"
                          ,"H2OvsNereo"
                          ,"H2OvsMast"
                          ,"H2OvsNereoMast"
                          ,"NereovsMast")
allPValues.ExNWater[1,1] <- ANOVA.BC.ExNWater.ExNvsH2O$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[2,1] <- ANOVA.BC.ExNWater.ExNvsNereo$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[3,1] <- ANOVA.BC.ExNWater.ExNvsMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[4,1] <- ANOVA.BC.ExNWater.ExNvsNereoMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[5,1] <- ANOVA.BC.ExNWater.H2OvsNereo$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[6,1] <- ANOVA.BC.ExNWater.H2OvsMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[7,1] <- ANOVA.BC.ExNWater.H2OvsNereoMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater[8,1] <- ANOVA.BC.ExNWater.NereovsMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.ExNWater <- cbind(allPValues.ExNWater, p.adjust(allPValues.ExNWater,method = "fdr", n = 10))
colnames(allPValues.ExNWater) <- c("P","fdr_adj")

# Lone Incube
# Add another column
MF.LoneIncube$WorS <- MF.LoneIncube$Substrate
MF.LoneIncube$WorS <- gsub("(Nereo)|(Mast)", "Seaweed", MF.LoneIncube[,'WorS'])

ANOVA.BC.LoneIncube <- adonis(dm.BC.LoneIncube ~ WorS*Treatment, data = MF.LoneIncube)
ANOSIM.BC.LoneIncube <- anosim(dm.BC.LoneIncube, grouping = MF.LoneIncube$ColRep) # NOT WORKING
# pairwaiseAdonis.BC.LoneIncube <- pairwise.adonis(dm.BC.LoneIncube, factors = MF.LoneIncube$ColRep)


# Lone Incube split
dm.BC.LoneIncube.NereovsMast <- dm.BC.LoneIncube[grep("(^Nereo[.])|(^Mast[.])", rownames(dm.BC.LoneIncube)), grep("(^Nereo[.])|(^Mast[.])", colnames(dm.BC.LoneIncube))]
MF.LoneIncube.NereovsMast <- MF.LoneIncube[grep("(^Nereo[.])|(^Mast[.])", rownames(MF.LoneIncube)),]
ANOVA.BC.LoneIncube.NereovsMast <- adonis(dm.BC.LoneIncube.NereovsMast ~ ColRep, data = MF.LoneIncube.NereovsMast)
ANOSIM.BC.LoneIncube.NereovsMast <- anosim(dm.BC.LoneIncube.NereovsMast, grouping = MF.LoneIncube.NereovsMast$ColRep)

dm.BC.LoneIncube.waters <- dm.BC.LoneIncube[grep("^water[.]", rownames(dm.BC.LoneIncube)), grep("^water[.]", colnames(dm.BC.LoneIncube))]
MF.LoneIncube.waters <- MF.LoneIncube[grep("^water[.]", rownames(MF.LoneIncube)),]
ANOVA.BC.LoneIncube.waters <- adonis(dm.BC.LoneIncube.waters ~ ColRep, data = MF.LoneIncube.waters)
ANOSIM.BC.LoneIncube.waters <- anosim(dm.BC.LoneIncube.waters, grouping = MF.LoneIncube.waters$ColRep)

dm.BC.LoneIncube.Nereovswater <- dm.BC.LoneIncube[grep("[.]Nereo[.]", rownames(dm.BC.LoneIncube)), grep("[.]Nereo[.]", colnames(dm.BC.LoneIncube))]
MF.LoneIncube.Nereovswater <- MF.LoneIncube[grep("[.]Nereo[.]", rownames(MF.LoneIncube)),]
ANOVA.BC.LoneIncube.Nereovswater <- adonis(dm.BC.LoneIncube.Nereovswater ~ ColRep, data = MF.LoneIncube.Nereovswater)
ANOSIM.BC.LoneIncube.Nereovswater <- anosim(dm.BC.LoneIncube.Nereovswater, grouping = MF.LoneIncube.Nereovswater$ColRep)

dm.BC.LoneIncube.Mastvswater <- dm.BC.LoneIncube[grep("[.]Mast[.]", rownames(dm.BC.LoneIncube)), grep("[.]Mast[.]", colnames(dm.BC.LoneIncube))]
MF.LoneIncube.Mastvswater <- MF.LoneIncube[grep("[.]Mast[.]", rownames(MF.LoneIncube)),]
ANOVA.BC.LoneIncube.Mastvswater <- adonis(dm.BC.LoneIncube.Mastvswater ~ ColRep, data = MF.LoneIncube.Mastvswater)
ANOSIM.BC.LoneIncube.Mastvswater <- anosim(dm.BC.LoneIncube.Mastvswater, grouping = MF.LoneIncube.Mastvswater$ColRep)


# 4 comparisons; 0.05/4 = 0.00625
allPValues.LoneIncube <- matrix(nrow= 4)
rownames(allPValues.LoneIncube) <- c("NereovsMast"
                                     ,"waters"
                                     ,"Nereovswater"
                                     ,"Mastvswater"
)
allPValues.LoneIncube[1,1] <- ANOVA.BC.LoneIncube.NereovsMast$aov.tab[6]$`Pr(>F)`[1]
allPValues.LoneIncube[2,1] <- ANOVA.BC.LoneIncube.waters$aov.tab[6]$`Pr(>F)`[1]
allPValues.LoneIncube[3,1] <- ANOVA.BC.LoneIncube.Nereovswater$aov.tab[6]$`Pr(>F)`[1]
allPValues.LoneIncube[4,1] <- ANOVA.BC.LoneIncube.Mastvswater$aov.tab[6]$`Pr(>F)`[1]

allPValues.LoneIncube <- cbind(allPValues.LoneIncube, p.adjust(allPValues.LoneIncube,method = "fdr", n = 4))
colnames(allPValues.LoneIncube) <- c("P","fdr_adj")

# Now print everything
system("mkdir ./BETAPLOTS/BC/")

capture.output(ANOVA.BC.ExN, file = "./BETAPLOTS/BC/ANOVA.BC.ExN.txt")
# capture.output(pairwaiseAdonis.BC.ExN, file = "pairwaiseAdonis.BC.ExN.txt")
# capture.output(allPValues.ExN, file = "allPValues.pairwise.ExN.txt")

capture.output(ANOVA.BC.ExNWater, file = "./BETAPLOTS/BC/ANOVA.BC.ExNWater.txt")
# capture.output(pairwaiseAdonis.BC.ExNWater, file = "pairwaiseAdonis.BC.ExNWater.txt")
capture.output(allPValues.ExNWater, file = "./BETAPLOTS/BC/allPValues.ExNWater.txt")
capture.output(allPValues.ExN, file = "./BETAPLOTS/BC/allPValues.ExN.txt")
capture.output(allPValues.LoneIncube, file = "./BETAPLOTS/BC/allPValues.LoneIncube.txt")


capture.output(ANOVA.BC.LoneIncube, file = "./BETAPLOTS/BC/ANOVA.BC.LoneIncube.txt")
# capture.output(pairwaiseAdonis.BC.LoneIncube, file = "pairwaiseAdonis.BC.LoneIncube.txt")


####### PLOT ############
# EXN BC
MF.ExN$ColRep <- factor(MF.ExN$ColRep, levels = c('NereotestExNExN','NereotestNereoExN','NereotestMastExN','NereotestNereoMastExN'))
ExNColours <- c("black","green","red","brown")

# Make chulls

NMDS.BC.ExN.ExN <- NMDS.BC.ExN$points[grep(".ExN.", rownames(NMDS.BC.ExN$points), fixed = TRUE),]
NMDS.BC.ExN.ExN.chull <- chull(NMDS.BC.ExN.ExN)
NMDS.BC.ExN.ExN.chull <- c(NMDS.BC.ExN.ExN.chull, NMDS.BC.ExN.ExN.chull[1])

NMDS.BC.ExN.Mast <- NMDS.BC.ExN$points[grep(".Mast.", rownames(NMDS.BC.ExN$points), fixed = TRUE),]
NMDS.BC.ExN.Mast.chull <- chull(NMDS.BC.ExN.Mast)
NMDS.BC.ExN.Mast.chull <- c(NMDS.BC.ExN.Mast.chull, NMDS.BC.ExN.Mast.chull[1])

NMDS.BC.ExN.Nereo <- NMDS.BC.ExN$points[grep(".Nereo.", rownames(NMDS.BC.ExN$points), fixed = TRUE),]
NMDS.BC.ExN.Nereo.chull <- chull(NMDS.BC.ExN.Nereo)
NMDS.BC.ExN.Nereo.chull <- c(NMDS.BC.ExN.Nereo.chull, NMDS.BC.ExN.Nereo.chull[1])

NMDS.BC.ExN.NereoMast <- NMDS.BC.ExN$points[grep(".NereoMast.", rownames(NMDS.BC.ExN$points), fixed = TRUE),]
NMDS.BC.ExN.NereoMast.chull <- chull(NMDS.BC.ExN.NereoMast)
NMDS.BC.ExN.NereoMast.chull <- c(NMDS.BC.ExN.NereoMast.chull, NMDS.BC.ExN.NereoMast.chull[1])

pdf("./BETAPLOTS/BC/NMDS_BC_ExN.pdf", width = 10, height = 7, pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(NMDS.BC.ExN$points
     , main = "NMDS of Nereo Meristem Swabs"
     , pch = 19
     , col = ExNColours[factor(MF.ExN$ColRep)]
     , sub = round(NMDS.BC.ExN$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2)
lines(NMDS.BC.ExN.ExN[NMDS.BC.ExN.ExN.chull,]
      , col = "black")
lines(NMDS.BC.ExN.Mast[NMDS.BC.ExN.Mast.chull,]
      , col = "red")
lines(NMDS.BC.ExN.Nereo[NMDS.BC.ExN.Nereo.chull,]
      , col = "green")
lines(NMDS.BC.ExN.NereoMast[NMDS.BC.ExN.NereoMast.chull,]
      , col = "brown")
par(fig = c(0.6,1,0,1), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n")
legend("center"
       , pch = 19
       , legend = c("Alone","With Nereocystis","With Mastocarpus","With Both")#levels(MF.ExN$ColRep)
       , col = ExNColours
       , cex = 1)
dev.off()


# ExN BC
MF.ExNWater$ColRep <- factor(MF.ExNWater$ColRep, levels = c('NereotestH2OWater','NereotestExNWater','NereotestNereoWater','NereotestMastWater','NereotestNereoMastWater'))
ExNWaterColours <- c("blue","black","green","red","brown")


# Make chulls
NMDS.BC.ExNWater.Water <- NMDS.BC.ExNWater$points[grep(".H2O.", rownames(NMDS.BC.ExNWater$points), fixed = TRUE),]
NMDS.BC.ExNWater.Water.chull <- chull(NMDS.BC.ExNWater.Water)
NMDS.BC.ExNWater.Water.chull <- c(NMDS.BC.ExNWater.Water.chull, NMDS.BC.ExNWater.Water.chull[1])

NMDS.BC.ExNWater.ExN <- NMDS.BC.ExNWater$points[grep(".ExN.", rownames(NMDS.BC.ExNWater$points), fixed = TRUE),]
NMDS.BC.ExNWater.ExN.chull <- chull(NMDS.BC.ExNWater.ExN)
NMDS.BC.ExNWater.ExN.chull <- c(NMDS.BC.ExNWater.ExN.chull, NMDS.BC.ExNWater.ExN.chull[1])

NMDS.BC.ExNWater.Mast <- NMDS.BC.ExNWater$points[grep(".Mast.", rownames(NMDS.BC.ExNWater$points), fixed = TRUE),]
NMDS.BC.ExNWater.Mast.chull <- chull(NMDS.BC.ExNWater.Mast)
NMDS.BC.ExNWater.Mast.chull <- c(NMDS.BC.ExNWater.Mast.chull, NMDS.BC.ExNWater.Mast.chull[1])

NMDS.BC.ExNWater.Nereo <- NMDS.BC.ExNWater$points[grep(".Nereo.", rownames(NMDS.BC.ExNWater$points), fixed = TRUE),]
NMDS.BC.ExNWater.Nereo.chull <- chull(NMDS.BC.ExNWater.Nereo)
NMDS.BC.ExNWater.Nereo.chull <- c(NMDS.BC.ExNWater.Nereo.chull, NMDS.BC.ExNWater.Nereo.chull[1])

NMDS.BC.ExNWater.NereoMast <- NMDS.BC.ExNWater$points[grep(".NereoMast.", rownames(NMDS.BC.ExNWater$points), fixed = TRUE),]
NMDS.BC.ExNWater.NereoMast.chull <- chull(NMDS.BC.ExNWater.NereoMast)
NMDS.BC.ExNWater.NereoMast.chull <- c(NMDS.BC.ExNWater.NereoMast.chull, NMDS.BC.ExNWater.NereoMast.chull[1])


pdf("./BETAPLOTS/BC/NMDS_BC_ExNWater.pdf", width = 10, height = 7, pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(NMDS.BC.ExNWater$points
     , main = "NMDS of water samples"
     , pch = 19
     , col = ExNWaterColours[factor(MF.ExNWater$ColRep)]
     , sub = round(NMDS.BC.ExNWater$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2)
lines(NMDS.BC.ExNWater.Water[NMDS.BC.ExNWater.Water.chull,]
      , col = "blue")
lines(NMDS.BC.ExNWater.ExN[NMDS.BC.ExNWater.ExN.chull,]
      , col = "black")
lines(NMDS.BC.ExNWater.Nereo[NMDS.BC.ExNWater.Nereo.chull,]
      , col = "green")
lines(NMDS.BC.ExNWater.Mast[NMDS.BC.ExNWater.Mast.chull,]
      , col = "red")
lines(NMDS.BC.ExNWater.NereoMast[NMDS.BC.ExNWater.NereoMast.chull,]
      , col = "brown")
par(fig = c(0.6,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , xlab = ""
     , ylab = ""
     , bty = "n")
legend("center"
       , pch = 19
       , legend = c("Water only","Nereo Meristem Only", "Nereocystis","Mastocarpus","Nereo + Mast")
       , col = ExNWaterColours
       , cex = 1)
dev.off()

# LoneIncube BC
MF.LoneIncube$ColRep <- factor(MF.LoneIncube$ColRep, levels = c('LoneincubeNereowater','LoneincubeNereoNereo','LoneincubeMastwater','LoneincubeMastMast'))
LoneIncubeColours <- c("blue","green","purple","red")

NMDS.BC.LoneIncube.WNereo <- NMDS.BC.LoneIncube$points[grep("water.Loneincube.Nereo.", rownames(NMDS.BC.LoneIncube$points), fixed = TRUE),]
NMDS.BC.LoneIncube.WNereo.chull <- chull(NMDS.BC.LoneIncube.WNereo)
NMDS.BC.LoneIncube.WNereo.chull <- c(NMDS.BC.LoneIncube.WNereo.chull, NMDS.BC.LoneIncube.WNereo.chull[1])

NMDS.BC.LoneIncube.WMast <- NMDS.BC.LoneIncube$points[grep("water.Loneincube.Mast.", rownames(NMDS.BC.LoneIncube$points), fixed = TRUE),]
NMDS.BC.LoneIncube.WMast.chull <- chull(NMDS.BC.LoneIncube.WMast)
NMDS.BC.LoneIncube.WMast.chull <- c(NMDS.BC.LoneIncube.WMast.chull, NMDS.BC.LoneIncube.WMast.chull[1])

NMDS.BC.LoneIncube.Nereo <- NMDS.BC.LoneIncube$points[grep("Nereo.Loneincube.", rownames(NMDS.BC.LoneIncube$points), fixed = TRUE),]
NMDS.BC.LoneIncube.Nereo.chull <- chull(NMDS.BC.LoneIncube.Nereo)
NMDS.BC.LoneIncube.Nereo.chull <- c(NMDS.BC.LoneIncube.Nereo.chull, NMDS.BC.LoneIncube.Nereo.chull[1])

NMDS.BC.LoneIncube.Mast <- NMDS.BC.LoneIncube$points[grep("Mast.Loneincube.", rownames(NMDS.BC.LoneIncube$points), fixed = TRUE),]
NMDS.BC.LoneIncube.Mast.chull <- chull(NMDS.BC.LoneIncube.Mast)
NMDS.BC.LoneIncube.Mast.chull <- c(NMDS.BC.LoneIncube.Mast.chull, NMDS.BC.LoneIncube.Mast.chull[1])


pdf("./BETAPLOTS/BC/NMDS_BC_LoneIncube.pdf", width = 10, height = 7, pointsize = 14)
par(fig = c(0,0.7,0,1))
plot(NMDS.BC.LoneIncube$points
     , main = "NMDS of Incubation Experiment"
     , pch = 19
     , col = LoneIncubeColours[factor(MF.LoneIncube$ColRep)]
     , sub = round(NMDS.BC.LoneIncube$stress/100,2)
     , xlab = "NMDS 1"
     , ylab = "NMDS 2"
     , cex = 2)
lines(NMDS.BC.LoneIncube.WNereo[NMDS.BC.LoneIncube.WNereo.chull,]
      , col = "blue")
lines(NMDS.BC.LoneIncube.WMast[NMDS.BC.LoneIncube.WMast.chull,]
      , col = "purple")
lines(NMDS.BC.LoneIncube.Nereo[NMDS.BC.LoneIncube.Nereo.chull,]
      , col = "green")
lines(NMDS.BC.LoneIncube.Mast[NMDS.BC.LoneIncube.Mast.chull,]
      , col = "red")
par(fig = c(0.6,1,0,1), mar = c(0,0,0,0), new = TRUE)
plot(0,0
     , pch = ''
     , xaxt = 'n'
     , yaxt = 'n'
     , xlab = ''
     , ylab = ''
     , bty = 'n')
legend("center"
       , pch = 19
       , legend = c("Water-Nereocystis","Nereocystis","Water-Mastocarpus","Mastocarpus")#levels(MF.LoneIncube$ColRep)
       , col = LoneIncubeColours
       , cex = 1)
dev.off()




