### use ggplot with plotly

require(plotly)
require(ggplot2)
require(reshape2)
df <- read.table('plotly_ggplot.txt',sep=",",header = T)
df2 <- melt(df, id.vars = c('x'))
pg <- ggplot(df2,aes(x,value,colour=variable))+geom_line()
pg <- pg + xlab("Gene body percentile (5’−>3’)") + ylab("Coverage")
pgp <- ggplotly(pg)
htmlwidgets::saveWidget(as_widget(pgp), "index.html")
