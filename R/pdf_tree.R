library (ape)
args = commandArgs(trailingOnly=TRUE)
file <- args[1]
tree <- ladderize(read.tree(file))
treename <- paste(file, ".pdf", sep="")
pdf(treename)
plot.phylo(tree, no.margin=TRUE, label.offset=0.01, cex=0.4, show.node.label = T)
dev.off()
