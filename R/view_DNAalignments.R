library (ape)
#get file list
files <- list.files(path=".", pattern="*.fas", full.names=T, recursive=FALSE)
#read files
loci <- lapply(files, function(x)
{
  locus <- read.dna(x, format="fasta", as.character = T)
  locus[order(rownames(locus)),]
  rownames(locus) <- strtrim(rownames(locus), width = 3)
  locus <- as.DNAbin(locus)
})
#set up plot parameters
par(mfrow=c(3,1), mar=c(1,1,1,1))
#plot alignments
n <- c(1:3)
for (x in n){
  image(loci[[x]], main=files[x], show.labels = T, legend = F, cex.lab=0.4, cex.axis=0.6, mgp=c(3,.1,0))
}