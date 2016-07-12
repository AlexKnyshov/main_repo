library (ape)
args = commandArgs(trailingOnly=TRUE)
file <- args[1]
file <- "BStree.tre"
tree <- read.nexus(file)
tree <- ladderize(read.tree(file))
tree
plot(tree,no.margin=TRUE, show.branch.label = T)
tree$tip.label <- gsub('.{5}$', '', tree$tip.label)
write.tree(tree, gsub('.{5}$', '', file))#"test.phy")