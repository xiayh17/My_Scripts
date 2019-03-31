require(ggplot2)
require(pheatmap)
require(reshape2)
library(RColorBrewer)
args <- commandArgs(trailing=T)
inputFile <- args[1]
outputFile <- args[2]

cc = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu"))) 

#pdf(outputFile, onefile=TRUE)
png(filename=outputFile,width = 6.0, height = 3.6, units='in',res=150)

df <- read.table(inputFile,sep='\t',header = F)
colnames(df) <- c("QualityScore","Cycle","CovariateValue","Observations")
#df$ratio <- df$CovariateValue-df$QualityScore
df$ratio <- df$QualityScore-df$CovariateValue
df2 <- df[,c('QualityScore','Cycle','ratio')]

#breaks = seq(min(df2$ratio),max(df2$ratio), length.out=100)
breaks = seq(-15,15, length.out=100)
dw <- dcast(df2, QualityScore ~ Cycle)
rownames(dw) <- dw$QualityScore
drops <- c('QualityScore')
dw2 <- dw[ , !(names(dw) %in% drops)]
pheatmap(dw2, color=cc(100), breaks=breaks,cluster_rows=F,cluster_cols = F,fontsize_row = 8,fontsize_col = 6,border_color = NA)
dev.off()
