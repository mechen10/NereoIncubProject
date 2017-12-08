#!/bin/bash/R

######### Comparing all DESeq results at OTU level ######
# This script compres DESeq results at the OTU level for all OTUs
# Purpose is to see how consistent patterns are for taxa with same name
# Input is output from OTUTraitsBetweenGroups.py
# We want both the fc and sig

# input for OTUTraitsBetweenGroups.py:
# /Users/melissachen/Documents/Masters/UNIVERSALCODE_git/Core_Analysis/OTUTraitsBetweenGroups.py -i /Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/DESEQ/RAW/NereotestExNWatervsNereotestNereoWater.txt,/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/DESEQ/RAW/NereotestExNWatervsNereotestMastWater.txt,/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/DESEQ/RAW/NereotestExNWatervsNereotestNereoMastWater.txt,/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/DESEQ/RAW/NereotestExNExNvsNereotestNereoExN.txt,/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/DESEQ/RAW/NereotestExNExNvsNereotestMastExN.txt,/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/DESEQ/RAW/NereotestExNExNvsNereotestNereoMastExN.txt -c log2FoldChange -o /Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/DESEQ/
# /Users/melissachen/Documents/Masters/UNIVERSALCODE_git/Core_Analysis/OTUTraitsBetweenGroups.py -i /Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/DESEQ/RAW/NereotestExNWatervsNereotestNereoWater.txt,/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/DESEQ/RAW/NereotestExNWatervsNereotestMastWater.txt,/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/DESEQ/RAW/NereotestExNWatervsNereotestNereoMastWater.txt,/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/DESEQ/RAW/NereotestExNExNvsNereotestNereoExN.txt,/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/DESEQ/RAW/NereotestExNExNvsNereotestMastExN.txt,/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/DESEQ/RAW/NereotestExNExNvsNereotestNereoMastExN.txt -c pvalue -o /Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/DESEQ/

# setwd("/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/DESEQ/")
# logChangePWD <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/DESEQ/log2FoldChange.txt"
# pValPWD <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/DESEQ/pvalue.txt"


setwd("/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/GENUSLEVEL/")
logChangePWD <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/GENUSLEVEL/Comparison_across_treatments/log2FoldChange.txt"
pValPWD <- "/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/2_analysis/GENUSLEVEL/Comparison_across_treatments/padj.txt"

#################### Load data #####################

logChange <- read.delim(paste0(logChangePWD)
                        , header = TRUE
                        , row.names = 1
                        , stringsAsFactors = FALSE)
# Get rid of last col which is just repeat of taxonomy
logChange <- logChange[-ncol(logChange)]

pVal <- read.delim(paste0(pValPWD)
                        , header = TRUE
                        , row.names = 1
                        , stringsAsFactors = FALSE)
# Get rid of last col which is just repeat of taxonomy
pVal <- pVal[-ncol(pVal)]

#################### BEGIN #####################
# First, sort so they are both alphabetical order
logChange.sorted <- logChange[order(row.names(logChange)),]
pVal.sorted <- pVal[order(row.names(pVal)),]

# Make sure names are the same
if (!all(rownames(logChange.sorted) == rownames(pVal.sorted))) {
  print("WARNING ROW NAMES NOT THE SAME")}

# Reorder column names
correctOrder <- c("NereotestExNWatervsNereotestNereoWater"
                  , "NereotestExNWatervsNereotestMastWater"  
                  , "NereotestExNWatervsNereotestNereoMastWater"
                  , "NereotestExNExNvsNereotestNereoExN" 
                  , "NereotestExNExNvsNereotestMastExN"   
                  , "NereotestExNExNvsNereotestNereoMastExN"    )
logChange.sorted <- logChange.sorted[,match(correctOrder, colnames(logChange.sorted))]
pVal.sorted <- pVal.sorted[,match(correctOrder, colnames(pVal.sorted))]


# Get max range and make red, grey, blue color palette for it
maxFC <- ceiling(max(max(logChange.sorted, na.rm = TRUE), abs(min(logChange.sorted, na.rm = TRUE))))
redbluePalette <- colorRampPalette(c("red", "grey", "blue"))
redblueColors <- redbluePalette(maxFC*2+1)

# Adjust names
nameList <- sapply(row.names(logChange.sorted)
       , function(x) {
         tempName <- strsplit(x, "..__")
       })
abbrNames <- unlist(lapply(nameList, 
       function(x) {
         if (length(x) > 1) {
           newName <- paste0(x[2], ": ", x[5],"_",x[6])
         } else {
           newName <- x[1]
         }
         return(newName)
       }
       ), use.names = FALSE)


#### Make matrix for sig ####
deseq.fc.color <- logChange.sorted

# Sort OTUs manually!
allenriched <- c()
allreduced <- c()
rest <- c()
for (r in 1:nrow(logChange.sorted)) {
  if (all(deseq.fc.color[r,] > 0, na.rm = TRUE)) {
    allenriched <- c(allenriched, r)
  } else if (all(deseq.fc.color[r,] < 0, na.rm = TRUE)) {
    allreduced <- c(allreduced,r)
  } else {
    rest <- rbind(rest, deseq.fc.color[r,])
  }
}
# Combine
# rest <- rest[with(as.data.frame(rest), order(-rest[,"Nereo.water"]
                                                               # ,-rest[,"Mast.water"]
                                                               # ,-rest[,"NereoMast.water"])), ]
enrichedTemp <- deseq.fc.color[allenriched,]
enrichSums <- order(rowSums(enrichedTemp), decreasing = TRUE)
enrichedTemp <- enrichedTemp[enrichSums,]

reducedTemp <- deseq.fc.color[allreduced,]
reduceSums <- order(rowSums(reducedTemp), decreasing = TRUE)
reducedTemp <- reducedTemp[reduceSums,]

restSums <- order(rowSums(rest), decreasing = TRUE)
rest <- rest[restSums,]

deseq.fc.color2 <- rbind(enrichedTemp, rest, reducedTemp)

deseq.sig.ALL <- pVal.sorted[match(rownames(deseq.fc.color2), rownames(pVal.sorted)),]

# Change deseq.sig into a matrix of stars for significance
deseq.sig.ALL.star <- deseq.sig.ALL
for (r in 1:nrow(deseq.sig.ALL.star)) {
  for (c in 1:ncol(deseq.sig.ALL.star)) {
    if (is.na(deseq.sig.ALL.star[r,c])) {
      deseq.sig.ALL.star[r,c] <- ""
    } else if (as.numeric(deseq.sig.ALL.star[r,c]) > 0.05) {
      deseq.sig.ALL.star[r,c] <- ""
    } else if (as.numeric(deseq.sig.ALL.star[r,c]) > 0.01) {
      deseq.sig.ALL.star[r,c] <- "*"
    } else if (as.numeric(deseq.sig.ALL.star[r,c]) > 0.001) {
      deseq.sig.ALL.star[r,c] <- "**"
    } else if (as.numeric(deseq.sig.ALL.star[r,c]) <= 0.001) {
      deseq.sig.ALL.star[r,c] <- "***"
    }
  }
}

# Now, get color
absMax <- ceiling(max(abs(deseq.fc.color2), na.rm = TRUE))
ncolors <- absMax*2+1

colorRange <- colorRampPalette(c("blue","grey","red"))

# Change names again
# Adjust names
nameList2 <- sapply(row.names(deseq.fc.color2)
                   , function(x) {
                     tempName <- strsplit(x, "..__")
                   })
abbrNames2 <- unlist(lapply(nameList2, 
                           function(x) {
                             if (length(x) > 1) {
                               newName <- paste0(x[2], ": ", x[5],"_",x[6])
                             } else {
                               newName <- x[1]
                             }
                             return(newName)
                           }
), use.names = FALSE)

rownames(deseq.fc.color2) <- make.names(abbrNames2, unique = TRUE)

######### PLOT ##############

# Make heatmap showing all OTUs in alphabetical order
pdf(file = "CompareOTUDESEq_ALL.pdf", width = 2, height = 15)
par(mar = c(5.1,1,2.1,5.5))
image(as.matrix(t(logChange.sorted))
      , col = redblueColors
      # , ylab = "OTUs"
      , axes = FALSE
      )
axis(side = 4
     , at = seq(0,1,length.out = nrow(logChange.sorted))
     , labels = abbrNames
     , las = 2
     , cex.axis = 0.2
     , tick = FALSE
     , line = -1
     )
axis(side = 1
     , at = seq(0,1, length.out = 6)
     , labels = c("WATER:Nereo", "WATER:Mast", "WATER:NereoMast"
                  , "NMF:Nereo", "NMF:Mast", "NMF:NereoMast")
     , las = 2
     , cex.axis = 0.5
     ,tick = FALSE
)
dev.off()

# Plot showing all OTUs, ordered
library("gplots")

pdf("./Heatmap_ALL_incllowabund.pdf", pointsize = 14, height = 15, width = 5)
heatmap.2(as.matrix(deseq.fc.color2)
          , Rowv = NA
          , Colv = NA
          , labCol = c("+Nereo","+Mast","+Nereo+Mast","+Nereo","+Mast","+Nereo+Mast")
          , col = colorRange(ncolors)
          , scale = "none"
          , dendrogram = "none"
          , offsetCol = 3
          , cexRow = 0.4
          
          , trace = "none"
          , density.info = "none"
          , key.xlab = "Fold-change"
          
          , na.color = "white"
          
          , margins = c(0,0)
          , lmat = rbind(c(0,4,5)
                         ,c(2,1,6)
                         ,c(0,3,7))
          , lhei = c(1,7,2)
          , lwid = c(0.5,3,4)
          , colsep = c(3)
          , sepwidth = c(0.1,0)

          , cellnote = as.matrix(deseq.sig.ALL.star)
          , notecol = "black"
          , notecex = 0.5
)
par(fig = c(0,1,0,1), mar = c(1,1,1,1), new = TRUE)
plot(0,0
     , pch = ""
     , xaxt = "n"
     , yaxt = "n"
     , bty = "n"
     , xlab = ""
     , ylab = "")
text(x = 0.05, y = c(1,0.98,0.96,0.94)
     , labels = c(" -","  *"," **","***")
     , pos = 2
     , cex = 0.8)
text(x = 0.1, y = c(1,0.98,0.96,0.94)
     , labels = c( "p = NA","p <= 0.05", "p <= 0.01","p <= 0.001")
     , pos = 4
     , cex = 0.8)
text(x = -0.80, y = c(-0.65,-0.69)
     , labels = c("_____________","WATER")
     , pos = 1
     , cex = 0.6
     , xpd = "n")
text(x = -0.33, y = c(-0.65,-0.69)
     , labels = c("_____________","NMF SURFACE")
     , pos = 1
     , cex = 0.6
     , xpd = "n")
dev.off()

