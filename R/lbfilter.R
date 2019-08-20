library(ape)

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
  treesize <- sapply(treelist, function(x) print(length(x$tip.label)))
  while(complete == F){
    print (paste0("iteration_", i))
    i = i + 1
    templist <- list()
    complete = T

    for (tr1 in treelist){


      nds <- tr1$edge[which(tr1$edge.length > mean(tr1$edge.length)*threshold)[1],]

      print(nds)
      if (!is.na(nds[1]) && !is.na(nds[2])){
        if (nds[1] == length(tr1$tip.label)+1){
          if (nds[2] <= length(tr1$tip.label)){
            newtr1 <- drop.tip(tr1, tr1$tip.label[nds[2]])
          } else {
            newtr1 <- extract.clade(tr1,nds[2])  
          }
        } else {

          newtr1 <- extract.clade(tr1,nds[1])
        }

        newtr2 <- drop.tip(tr1, newtr1$tip.label)

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

if (seqtype == "dna"){
  locus <- read.dna(seqpath, format="fasta")
} else if (seqtype == "prot") {
  library(seqinr)
  locus <- read.alignment(seqpath, format = "fasta")
  locus <- as.matrix.alignment(locus)
  locus <- toupper(locus)
  locus <- as.AAbin.character(locus)
} else {
  cat ("incorrect sequence type in the command, exit\n")
  quit()
}

if (length(trsout$tip.label)>0){
  locus <- locus[trimws(rownames(locus)) %in% trsout$tip.label,]
}

write.FASTA(locus, paste0(seqpath,".edited"))
