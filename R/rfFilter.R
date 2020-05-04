library(ape)

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  cat("Syntax: Rscript rfFilter.R [path to concat genetrees] [path to text file with gene names] ([path to ref tree])\n")
  cat("Example: Rscript rfFilter.R genetrees.tre genetrees.txt\n")
  cat("Example: Rscript rfFilter.R genetrees.tre genetrees.txt concattree.tre\n")
  quit()
} else if (length(args) >= 2){
	gtrees <- read.tree(args[1])
	names(gtrees) <-  readLines(args[2])
	if (length(args) == 3){
		tree <- read.tree(args[3])
	}
}

FindOutliers <- function(data) {
  lowerq = quantile(data, na.rm = T)[2]
  upperq = quantile(data, na.rm = T)[4]
  iqr = upperq - lowerq
  extreme.threshold.upper = (iqr * 1.5) + upperq
  extreme.threshold.lower = lowerq - (iqr * 1.5)
  result <- which(data > extreme.threshold.upper)
}

if (length(args) == 3){
	tree_dist <- data.frame(locus=character(), rf = numeric())
	for (t in 1:length(gtrees)){
	  print (names(gtrees[t]))
	  temptree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% gtrees[[t]]$tip.label)])
	  
	  distval = as.numeric(dist.topo(gtrees[t], temptree))
	  print (distval)
	  tree_dist <- rbind(tree_dist, data.frame(locus=names(gtrees[t]), rf=distval))
	}
} else {
	tree_dist <- data.frame(locus=character(), rf = numeric())
	for (t in 1:length(gtrees)){
	  distval <- numeric()
	  for (t2 in 1:length(gtrees)){
	    temptree1 <- drop.tip(gtrees[[t]], setdiff(gtrees[[t]]$tip.label, gtrees[[t2]]$tip.label))
	    temptree2 <- drop.tip(gtrees[[t2]], setdiff(gtrees[[t2]]$tip.label, gtrees[[t]]$tip.label))
	    distval <- c(distval, as.numeric(dist.topo(temptree1, temptree2)))
	  }
	  tree_dist <- rbind(tree_dist, data.frame(locus=names(gtrees[t]), rf=mean(distval))) 
	}
}

write(as.character(tree_dist$locus[FindOutliers(tree_dist$rf)]), "rf_loci.txt")