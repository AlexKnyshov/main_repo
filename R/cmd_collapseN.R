args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  cat("Syntax: Rscript cmd_collapseN.R [path to tree] [threshold]\n")
  cat("Example: Rscript cmd_collapseN.R tree.tre 10\n")
  quit()
}

library(ape)
collapse_nodes <- function(t1, threshold){
  for (i in 1:t1$Nnode){
    if (!is.na(as.numeric(t1$node.label[i])) & as.numeric(t1$node.label[i])<threshold){
      t1$edge.length[t1$edge[,1] == length(t1$tip.label)+i] <- t1$edge.length[t1$edge[,1] == length(t1$tip.label)+i] + t1$edge.length[which(t1$edge[,2] == length(t1$tip.label)+i)]
      t1$edge.length[which(t1$edge[,2] == length(t1$tip.label)+i)] <- 0
    }
  }
  t2 <- di2multi(t1)
  return(t2)
}

tree <- read.tree(args[1])
threshold <- as.numeric(args[2])

collapsed_tree <- collapse_nodes(tree, threshold)

write.tree(collapsed_tree, paste0(args[1],".newick"))