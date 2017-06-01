#!/bin/Rscript
library(optparse)

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
# setwd("/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_NEREOINCUBE/ANALYSIS_ALPHABETA/R")
# MFPWD = "/Users/parfreylab/Desktop/personal_files/melissa/ForBotanyCluster/z_NEREOINCUBE/ANALYSIS_ALPHABETA/OTU_Tables_and_MP/MF_withalpha.txt"
# alphaNames = c("chao1","PD_whole_tree","observed_otus")



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
for (i in alphaNames) {
  MF.ExN.temp <- MF.ExN[,c(paste0(i,"_even_3500_alpha"), "ColRep","Replicate")]
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
  
  jpeg(paste0("ALPHAPLOTS/Alpha_div_meristemswabs_",i,".jpeg"))
  par(mar = c(10,4,4,4))
  boxplot(t(MF.ExN.Alpha)
          , las = 2
          , col = c("yellow","green","red","brown")
          , ylab = paste0("Alpha Diversity (", i,")")
          , main = "Alpha diversity across meristem swabs"
          )
  title(xlab = "Treatment", line = 8)
  dev.off()
  
}

# For MF.LoneIncube
for (i in alphaNames) {
  MF.ExNWater.temp <- MF.ExNWater[,c(paste0(i,"_even_3500_alpha"), "ColRep","Replicate")]
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
  
  jpeg(paste0("ALPHAPLOTS/Alpha_div_watersamples_",i,".jpeg"))
  par(mar = c(10,4,4,4))
  boxplot(t(MF.ExNWater.Alpha)
          , las = 2
          , col = c("yellow","green","red","brown")
          , ylab = paste0("Alpha Diversity (", i,")")
          , main = "Alpha diversity across water samples"
  )
  title(xlab = "Treatment", line = 8)
  dev.off()

  
}

# For MF.LoneIncube
for (i in alphaNames) {
  MF.LoneIncube.temp <- MF.LoneIncube[,c(paste0(i,"_even_3500_alpha"), "ColRep","Replicate")]
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
  
  jpeg(paste0("ALPHAPLOTS/Alpha_div_loneincube_",i,".jpeg"))
  par(mar = c(10,4,4,4))
  boxplot(t(MF.LoneIncube.Alpha)
          , las = 2
          , col = c("yellow","green","red","brown")
          , ylab = paste0("Alpha Diversity (", i,")")
          , main = "Alpha diversity"
  )
  title(xlab = "Treatment", line = 8)
  dev.off()
  
  
}

