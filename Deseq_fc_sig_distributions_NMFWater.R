setwd("/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/BYOTU_notgenus/DESEQ")
data <- read.delim("/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/BYOTU_notgenus/DESEQ/deseq.fc_ALLVALUESFORHIST.txt")
data2 <- read.delim("/Users/melissachen/Documents/Masters/Project_Masters/Project_MacroalgaeSource/1_analysis/BYOTU_notgenus/DESEQ/deseq.sig.txt"
                    , row.names = 1)


data.abs <- abs(data)
pdf("Histogram_AbsFoldChange_forNMFvsWater.pdf")
hist(data.abs[,2]
     , main = "Absolute Fold-Change for NMF (red) and water (blue) samples"
     , xlab = "Fold-change"
     , col = rgb(0,0,1,0.5)
     , breaks = seq(0,10, by = 0.5)
     , ylim = c(0,100))
hist(data.abs[,1]
     , col = rgb(1,0,0,.5)
     , add = TRUE
     , breaks = seq(0,10, by= 0.5)
     , ylim = c(0,100))
dev.off()

#######
allExN_sig <- unlist(data2[,1:3])
allWater_sig <- unlist(data2[,4:6])

hist.allExN_sig <- hist(allExN_sig
                        # , breaks = seq(0,1,by = 0.025)
                        , breaks = c(0,0.001,0.01,0.05,0.1,0.5,1)
                        , plot = FALSE)
hist.allWater_sig <- hist(allWater_sig
                          # , breaks = seq(0,1,by = 0.025)
                          , breaks = c(0,0.001,0.01,0.05,0.1,0.5,1)
                          , plot = FALSE)

dat.hist.filt <- rbind(hist.allExN_sig$counts, hist.allWater_sig$counts)
colnames(dat.hist.filt) <- paste0("<", hist.allExN_sig$breaks[2:length(hist.allExN_sig$breaks)])

pdf("Histogram_sigDistr_forNMFvsWater.pdf")
barplot(dat.hist.filt
        , main = "Hist of p-value distributions for NMF(red) and Water(blue)"
        , beside = TRUE
        , col = c("red","blue")
        , border = FALSE)
dev.off()


