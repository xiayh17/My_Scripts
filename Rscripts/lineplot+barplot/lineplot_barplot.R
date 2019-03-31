require(ggplot2)
df <- read.table('five.100_bin.depth.stat',sep="\t",header = T)
ggplot(df) + geom_line(aes(x=gc,y=median,color=sample),size=0.75) + geom_bar(aes(x=gc,y=size),stat="identity",col="#FFAAAA",fill="#FFAAAA") 
p <- p + ylab("Depth") + xlab("GC Content") 
p <- p + theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=13)) +scale_x_continuous(breaks=seq(0,100,5)) + theme(legend.text = element_text(size=20))
p
