library(ape)
library(phangorn)
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  cat("Syntax: Rscript lbfilter.R [path to tree] [path to alignment] [seq type: dna, prot] [threshold]\n")
  cat("Example: Rscript lbfilter.R locus.tre locus.fas dna 4\n")
  quit()
}
treepath <- args[1]
seqpath <- args[2]
seqtype <- args[3]
threshval <- as.numeric(args[4])

trs <- read.tree(treepath)

splittrees <- function(treelist, threshold = 4, return_max = F){
  complete = FALSE
  i = 1
  #number of tips in each tree (initially in the input tree)
  treesize <- sapply(treelist, function(x) print(length(x$tip.label)))
  #threshold - use first elem since only 1 input tree
  treebound <- mean(treelist[[1]]$edge.length)*threshold
  print (treebound)
  while(complete == F){
    print (paste0("iteration_", i))
    i = i + 1
    templist <- list()
    complete = T
    #iterate over trees (initially 1)
    for (tr1 in treelist){
      #get ancestral and descendant nodes of the long branches
      if (length(tr1$tip.label) == 1) {
        templist <- c(templist,list(tr1))
      } else {
        nds <- tr1$edge[which(tr1$edge.length > treebound)[1],]
        print(nds)
        if (!is.na(nds[1]) && !is.na(nds[2])){
          tipward_tips <- tr1$tip.label[Descendants(tr1,nds[2])[[1]]]
          newtr1 <- drop.tip(tr1, tipward_tips)
          newtr2 <- drop.tip(tr1, tr1$tip.label[!(tr1$tip.label %in% tipward_tips)] )
          templist <- c(templist,list(newtr1, newtr2))
          complete = F
        } else {
          templist <- c(templist,list(tr1))
        }
      }
    }
    treesize <- sapply(templist, function(x) print(length(x$tip.label)) )
    treelist <- templist
  }
  if (return_max == T){
    return(treelist[[which.max(treesize)]])
  } else {
    return(treelist)
  }
}

if (seqtype == "dna"){
  locus <- read.dna(seqpath, format="fasta")
  len_adjustment <- length(locus[1,])/3
} else if (seqtype == "prot") {
  library(seqinr)
  locus <- read.alignment(seqpath, format = "fasta")
  locus <- as.matrix.alignment(locus)
  locus <- toupper(locus)
  locus <- as.AAbin.character(locus)
  len_adjustment <- length(locus[1,])
} else {
  cat ("incorrect sequence type in the command, exit\n")
  quit()
}

adjusted_thresh <- threshval / 200 * len_adjustment
trsout <- splittrees(list(trs), adjusted_thresh, T)
# trsout <- splittrees(list(trs), threshval , T)

if (length(trsout$tip.label)>0){
  locus <- locus[trimws(rownames(locus)) %in% trsout$tip.label,]
}

write.FASTA(locus, paste0(seqpath,".edited"))
