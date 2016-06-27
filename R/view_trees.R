#dev.off()
library (ape)
#create file list
files <- list.files(path=".", pattern="*", full.names=T, recursive=FALSE)
#read files
trees <- lapply(files, function(x)
{
  tree <- read.tree(x)
})
#read loci names separately
names <- sapply(strsplit(sapply(strsplit(files, split='bipartitions.', fixed=TRUE), function(x) (x[2])), split='.an', fixed=TRUE), function(x) (x[1]))
#set up the plot parameters
par(mfrow=c(8,8),mar=c(1,1,1,1))
#draw trees
n <- c(1:64)
for (x in n){
  print(x)
  print(trees[[x]])
  plot(trees[[x]], main=names[x])
}
