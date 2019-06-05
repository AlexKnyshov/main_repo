library(ape)
args = commandArgs(trailingOnly=TRUE)
# tab <- read.table(args[1], skip=1, stringsAsFactors=FALSE)
# print(tab[1,1:3])
# prts <- read.table(args[2], stringsAsFactors=FALSE)
# print(prts[1,4])
# print(as.numeric(unlist(strsplit(prts[1,4], "-"))))
tab <- read.csv(args[1], stringsAsFactors=FALSE)
print(tab[1,177])
start <- read.tree(args[2])
datalist <- list()
for (node in 2:start$Nnode){
	datalist[[node]] <- data.frame(char=numeric(),
                            		value=numeric(),
                            		col=numeric()) 
	tr2 <- (node-1)*2-1
	tr3 <- (node-1)*2
	for (partition in 1:length(tab)){
	#for (partition in 21:40){
		print(paste0("partition:",partition))
		print (tab[partition,2])
		print (tab[partition,tr2+1])
		print (tab[partition,tr3+1])
		pls <- tab[partition,2] - max(c(tab[partition,tr2+1],tab[partition,tr3+1]))
		print (pls)
		if (pls < 0){
			colint <- 1
		}
		else {
			colint <- 0
		}
		datalist[[node]] <- rbind(datalist[[node]], list(char=partition, value=pls, col=colint))
	}
}
# print (datalist)
library(ggtree)
library(ggplot2)
#print(datalist[[2]])
# ggplot(data=datalist[[2]], aes(x=datalist[[2]]$char, y=datalist[[2]]$value), color="black")+
#     geom_bar(stat="identity") + theme(legend.position="none",axis.title.x=element_blank(),
#                                       axis.title.y=element_blank(),
#                                       axis.text.y=element_blank(),
#                                       axis.ticks.x=element_blank(),
#                                       axis.ticks.y=element_blank(),
#                                       panel.border = element_blank(),
#                                       panel.grid.major = element_blank(),
#                                       panel.grid.minor = element_blank(),
#                                       panel.background = element_rect(fill = 'white', colour = 'black', size=0.7),
#                                       #plot.background = element_blank(),
#                                       plot.margin = margin(t=1.5, r=2, b=1.5, l=0),
#                                       axis.text.x = element_text(size=1,margin = margin(t =0.01), colour = "black"))


p1 <- ggtree(start, size=1)# + geom_tiplab(size=3) + geom_text2(aes(subset = !isTip, label=label),size=3.5 , nudge_x =0.005)#+ ggplot2::xlim(0, 30)
insets <- list()
#print ("DEBUG")
#print(datalist[[2]])
for (i in 2:90){
	#print(i)
	datatemp <- datalist[[i]]
	print (datatemp)
  	insets[[i-1]] <- ggplot(data=datatemp, aes(x=char, y=value, fill=as.character(col)))+
    geom_bar(stat="identity") + theme(legend.position="none",axis.title.x=element_blank(),
                                      axis.title.y=element_blank(),
                                      axis.text.y=element_blank(),
                                      axis.ticks.x=element_blank(),
                                      axis.ticks.y=element_blank(),
                                      panel.border = element_blank(),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.background = element_rect(fill = 'transparent', colour = 'black', size=0.7),
                                      plot.background = element_blank(),
                                      plot.margin = margin(t=1.5, r=2, b=1.5, l=0),
                                      axis.text.x = element_text(size=1,margin = margin(t =0.01), colour = "black")) +
    scale_fill_manual(values = c("red","green"))
  
}
names(insets) <- (2+Ntip(start)):(90+Ntip(start))
print(insets)
p2 <- inset(p1, insets, width=0.04, height=6,vjust=0.2,hjust=0.01)
ggsave("testplot.pdf",p2, width=11, height=8.5)