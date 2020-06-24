require(ggplot2)
require(reshape2)
require(dplyr)
COLORS<- c("#310B68","#FC6B6B","#6BCC69","#C59E37","#FECE68","#CB3535","#036903","#FE0405","#389C66","#FC9D05",
"#073435","#386A9D","#FF9A9B","#CD0438","#04CC06","#060436","#D2686B","#373599","#37CD05","#D29BD1")
df <- read.table('snp.opt-1.barplot.txt',sep="\t",header = T)
df2 <- melt(df, id=c("gc_bin"))
text <- df2 %>% group_by(variable) %>% summarize(total = sum(value))     ## total number of each varable
p <- ggplot(df2) + geom_bar(aes(x=variable,y=value,fill=gc_bin),position = "fill",stat="identity",width=0.5) +scale_fill_manual(values =COLORS)+theme_bw()
p <- p + theme(axis.text =element_text(size=13), axis.title = element_text(size=14)) + ylab("Proportion") + ggtitle('opt-1 / opt_onestep-1 SNP FN')
p <- p + geom_text(aes(variable, c(1,1,1,1), label = total),size=4, data = tt,vjust = -0.3) ## add text on the top of the bars
p
