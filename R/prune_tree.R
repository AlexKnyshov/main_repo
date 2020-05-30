library(ape)

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  cat("Syntax: Rscript prune_tree.R [path to tree file] [path to fasta]\n")
  cat("Example: Rscript prune_tree.R tree.tre alignment.fas\n")
  quit()
} else if (length(args) == 2){
	treefilename <- args[1]
	seqfilename <- args[2]
}

tree <- read.tree(treefilename)
seqs <- read.dna(seqfilename, format="fasta", as.character = T)

phy1 <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% rownames(seqs))])

write.tree(phy1, paste0(treefilename, ".subset.tre"))