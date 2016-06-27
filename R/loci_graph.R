library (ape)
#get file list
files <- list.files(path=".", pattern="*.fas", full.names=T, recursive=FALSE)
#read files
loci <- lapply(files, function(x)
{
  locus <- read.dna(x, format="fasta", as.character = T)
  locus[order(rownames(locus)),]
  rownames(locus) <- strtrim(rownames(locus), width = 10)
  return (locus) #for normal graph
  #return (length(locus[1,])) #for histo
})

#histogram
#loci <- unlist(loci)
#hist (loci,70)

#get loci names
names <- lapply(loci, function(x){
  return (rownames(x))
})
names <- unique(unlist(names))

#set up plot parameters
max = c(1:length(names))
par(mfrow=c(length(names),1),mar=c(1,1,1,1))

#plot
for (i in max){
  name <- names[i]
  list1 <- lapply(loci, function(x){
    y <- rownames(x)
    if (name %in% y) return (TRUE)
    else return (FALSE)
  })
  list2 <- as.numeric(list1)
  name = paste(c(name,length(list2[list2==1]),round(length(list2[list2==1])/length(list2)*100, digits = 2)), collapse = " = ")
  name
  barplot(list2, main = name, axes=F, col="black", beside = T, space = 0)
}
