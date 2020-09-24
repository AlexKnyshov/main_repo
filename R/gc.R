library(ape)

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  cat("Syntax: Rscript gc.R [path to the folder with alignments]\n")
  cat("Example: Rscript gc.R ./fasta/\n")
  quit()
} else if (length(args) == 1){
	dnaseqpath <- args[1]
}

getGCcontent <- function(data, dimension){
	GCcontent <- apply(data, dimension, function (x) {
		gc <- x[x == "g" | x == "c"]
		at <- x[x == "a" | x == "t"]
		return(length(gc) / (length(gc) + length(at)))
		})
	return(GCcontent)
}

dnafiles <- list.files(path=dnaseqpath, pattern="\\.fas$", full.names=FALSE, recursive=FALSE)
locidata <- data.frame()
for (i in 1:length(dnafiles)){
	#path
	dnafile <- paste0(dnaseqpath,"/",dnafiles[i])
	dnalocus <- read.dna(dnafile, format="fasta", as.character=T, as.matrix=T)
	# pos1
	mask = c(T,F,F)
	seqbasedGC1 <- getGCcontent(dnalocus[,mask], 1) #seq based
	posbasedGC1 <- getGCcontent(dnalocus[,mask], 2) #pos based
	# pos2
	mask = c(F,T,F)
	seqbasedGC2 <- getGCcontent(dnalocus[,mask], 1) #seq based
	posbasedGC2 <- getGCcontent(dnalocus[,mask], 2) #pos based
	# pos3
	mask = c(F,F,T)
	seqbasedGC3 <- getGCcontent(dnalocus[,mask], 1) #seq based
	posbasedGC3 <- getGCcontent(dnalocus[,mask], 2) #pos based

	locidata <- rbind(locidata, data.frame(locus=dnafiles[i],
		loclen=length(dnalocus[1,]),
		taxnum=length(dnalocus[,1]),
		seqgc1mean=mean(seqbasedGC1, na.rm=T),
		seqgc1sd=sd(seqbasedGC1, na.rm=T),
		posgc1mean=mean(posbasedGC1, na.rm=T),
		posgc1sd=sd(posbasedGC1, na.rm=T),
		seqgc2mean=mean(seqbasedGC2, na.rm=T),
		seqgc2sd=sd(seqbasedGC2, na.rm=T),
		posgc2mean=mean(posbasedGC2, na.rm=T),
		posgc2sd=sd(posbasedGC2, na.rm=T),
		seqgc3mean=mean(seqbasedGC3, na.rm=T),
		seqgc3sd=sd(seqbasedGC3, na.rm=T),
		posgc3mean=mean(posbasedGC3, na.rm=T),
		posgc3sd=sd(posbasedGC3, na.rm=T)	
	))
	print (paste0(dnafiles[i], " processed"))
}
write.csv(locidata, "gc3.csv")