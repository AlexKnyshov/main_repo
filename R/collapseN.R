library(ape)
collapse_nodes <- function(t1, threshold){
  for (i in 1:t1$Nnode){
    #print (i)
    if (!is.na(as.numeric(t1$node.label[i])) & as.numeric(t1$node.label[i])<threshold){
      t1$edge.length[t1$edge[,1] == length(t1$tip.label)+i] <- t1$edge.length[t1$edge[,1] == length(t1$tip.label)+i] + t1$edge.length[which(t1$edge[,2] == length(t1$tip.label)+i)]
      t1$edge.length[which(t1$edge[,2] == length(t1$tip.label)+i)] <- 0
    }
  }
  t2 <- di2multi(t1)
  return(t2)
}