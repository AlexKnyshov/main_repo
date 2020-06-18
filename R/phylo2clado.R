library(ape)
args = commandArgs(trailingOnly=TRUE)
argsLen <- length(args);
if (argsLen == 0){
  cat("Syntax: Rscript phylo2clado.R [path to tree]\n")
  cat("Example: Rscript phylo2clado.R tree.tre\n")
  quit()
}
intree <- read.tree(args[1])
intree$edge.length <- NULL
write.tree(intree, paste0(args[1], ".edited"))