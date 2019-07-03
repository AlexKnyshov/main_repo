library(ape)

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  cat("Syntax: Rscript lbfilter.R [path to tree] [path to alignment] [threshold]\n")
  cat("Example: Rscript lbfilter.R locus.tre locus.fas 4\n")
  quit()
}
treepath <- args[1]
seqpath <- args[2]
threshval <- as.numeric(args[3])
# 
# treepath <- "RAxML_bestTree.L217.fas"
# seqpath <- "L217.fas"
# threshval <- 4
# 


trs <- read.tree(treepath)
# if (is.rooted(trs)){
#   trs <- unroot(trs)
# }

splittrees <- function(treelist, threshold = 4, return_max = F){
  complete = FALSE
  i = 1
  treesize <- sapply(treelist, function(x) print(length(x$tip.label)))
  while(complete == F){
    print (paste0("iteration_", i))
    i = i + 1
    templist <- list()
    complete = T
    #print (length(treelist))
    for (tr1 in treelist){
      # if (is.rooted(tr1)){
      #   tr1 <- unroot(tr1)
      # }
      #print (tr1)
      nds <- tr1$edge[which(tr1$edge.length > mean(tr1$edge.length)*threshold)[1],]
      #print (tr1$edge[which(tr1$edge.length > mean(tr1$edge.length)*threshold)[1],])
      #if (!is.null(nds) && !is.na(nds)){
      if (!is.na(nds)){
        if (nds[1] == length(tr1$tip.label)+1){
          newtr1 <- extract.clade(tr1,nds[2])
        } else {
          newtr1 <- extract.clade(tr1,nds[1])
        }
       # print (newtr1)
        newtr2 <- drop.tip(tr1, newtr1$tip.label)
        #print (newtr2)
        templist <- c(templist,list(newtr1, newtr2))
        complete = F
      } else {
        templist <- c(templist,list(tr1))
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
trsout <- splittrees(list(trs), threshval, T)

# plot(trs, show.tip.label = F, type="unrooted")
# edgelabels()
# nodelabels()
# newtr1 <- extract.clade(trs,134)
# plot(newtr1, show.tip.label = F, type="unrooted")
# nds <- newtr1$edge[which(newtr1$edge.length > mean(newtr1$edge.length)*4)[1],]
# if (nds[1] == length(newtr1$tip.label)+1){
#   extract.clade(newtr1,nds[2])
# }
# extract.clade(newtr1,8)
# drop.tip(newtr1, extract.clade(newtr1,nds[2])$tip.label)

locus <- read.dna(seqpath, format="fasta")
#rownames(locus)

if (length(trsout$tip.label)>0){
  seqs <- match(trsout$tip.label,trimws(rownames(locus)))
  locus <- locus[seqs,]
}

write.FASTA(locus, paste0(seqpath,".edited"))
