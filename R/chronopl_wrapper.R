library(ape)

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  cat("Syntax: Rscript chronopl_wrapper.R [path to tree file]\n")
  cat("Example: Rscript chronopl_wrapper.R tree.tre\n")
  quit()
} else if (length(args) == 1){
	treefilename <- args[1]
}

tree <- read.tree(treefilename)
phy1 <- chronos(tree)
write.tree(phy1, paste0(treefilename, ".ultra.tre"))