library(ape)
library(seqinr)
args = commandArgs(trailingOnly=TRUE)
argsLen <- length(args);
if (argsLen == 0){
  cat("Syntax: Rscript distfilter.R [path to concat dna] [path to dna folder] [path to concat protein] [path to protein folder] ([upperN] ([lowerN] ([rescale]))) \n")
  cat("Example: Rscript distfilter.R concatdna.fas ./dna/ concataa.fas ./aa/\n")
  cat("Example: Rscript distfilter.R concatdna.fas ./dna/ concataa.fas ./aa/ 1.5\n")
  cat("Example: Rscript distfilter.R concatdna.fas ./dna/ concataa.fas ./aa/ 3 3\n")
  quit()
}
if (argsLen == 4) {
	upperN = 3
	lowerN = 3
	rescale = FALSE
}
if (argsLen == 5) {
	upperN = as.numeric(args[5])
	lowerN = 3
	rescale = FALSE
}
if (argsLen == 6) {
	upperN = as.numeric(args[5])
	lowerN = as.numeric(args[6])
	rescale = FALSE
}
if (argsLen == 7) {
	upperN = as.numeric(args[5])
	lowerN = as.numeric(args[6])
	rescale = TRUE
	library(phangorn)
}
FindOutliers <- function(data) {
  lowerq = quantile(data, na.rm = T)[2]
  upperq = quantile(data, na.rm = T)[4]
  iqr = upperq - lowerq #Or use IQR(data)
  # identify the outliers
  extreme.threshold.upper = (iqr * upperN) + upperq
  extreme.threshold.lower = lowerq - (iqr * lowerN)
  return(list(upperO = which(data > extreme.threshold.upper), lowerO=which(data < extreme.threshold.lower)))
}
ComputeDistMatrix <- function(alignment) {
	alignmentmx <- as.phyDat.AAbin(alignment)
	dist_mx_dims <- c(length(alignmentmx), length(alignmentmx[[1]]))
	dist_mx <- matrix(NA, nrow =dist_mx_dims[1], ncol=dist_mx_dims[1], dimnames = list(names(alignmentmx),names(alignmentmx)) )
	logicmx <- lower.tri(dist_mx,diag = FALSE) 
	for (row1 in 1:length(dist_mx[,1])){
	  for (col1 in 1:length(dist_mx[1,])){
	    if (logicmx[row1, col1] == TRUE ){
	      dist_mx[row1, col1] <- dist.ml(alignmentmx[c(row1,col1), ], model="WAG")
	    }
	  }
	}
	dist_mx[upper.tri(dist_mx, diag = F)] <- t(dist_mx)[upper.tri(dist_mx, diag=F)]
	diag(dist_mx) <- 0
	return (dist_mx) 
}

concatdna <- args[1]
dnaseqpath <- args[2]
concataa <- args[3]
aaseqpath <- args[4]

print("check input files / folders")
dnaconcat <- read.dna(concatdna, format="fasta")
proteinconcat <- read.alignment(concataa, format = "fasta")
proteinconcat <- as.matrix.alignment(proteinconcat)
proteinconcat <- toupper(proteinconcat)
proteinconcat <- replace(proteinconcat, proteinconcat == "X", "-")
proteinconcat <- as.AAbin.character(proteinconcat)
dnafiles <- list.files(path=dnaseqpath, pattern="\\.fas$", full.names=FALSE, recursive=FALSE)
proteinfiles <- list.files(path=aaseqpath, pattern="\\.fas$", full.names=FALSE, recursive=FALSE)
commonfiles <- intersect(dnafiles, proteinfiles)
print(paste0("files to process: ", length(commonfiles)))
print("compute average taxon distance")
if (rescale == TRUE) {
	# dnadist <- ComputeDistMatrix(dnaconcat)
	dnadist <- dist.dna(dnaconcat, as.matrix = T, pairwise.deletion = T)
	proteindist <- ComputeDistMatrix(proteinconcat)
} else {
	dnadist <- dist.dna(dnaconcat, as.matrix = T, pairwise.deletion = T)
	proteindist <- dist.aa(proteinconcat, pairwise.deletion = T)
	proteindist <- as.matrix(proteindist)
	proteindist <- proteindist / length(proteinconcat[1,])
}
avgdnadist <- apply(dnadist, 1, mean)
avgproteindist <- apply(proteindist, 1, mean)
print("process files")
for (i in 1:length(commonfiles)){
	#path
	dnafile <- paste0(dnaseqpath,"/",dnafiles[i])
	proteinfile <- paste0(aaseqpath, "/", proteinfiles[i])
	#read
	dnalocus <- read.dna(dnafile, format="fasta")
	proteinlocus <- read.alignment(proteinfile, format = "fasta")
	proteinlocus <- as.matrix.alignment(proteinlocus)
	proteinlocus <- toupper(proteinlocus)
	proteinlocus <- replace(proteinlocus, proteinlocus == "X", "-")
	proteinlocus <- as.AAbin.character(proteinlocus)
	if (length(dnalocus[,1]) > 1 & length(proteinlocus[,1]) > 1){
		#get distance
		if (rescale == TRUE) {
			# dnalocusdist <- ComputeDistMatrix(dnalocus)
			dnalocusdist <- dist.dna(dnalocus, as.matrix = T, pairwise.deletion = T)
			proteinlocusdist <- ComputeDistMatrix(proteinlocus)
		} else {
			dnalocusdist <- dist.dna(dnalocus, as.matrix = T, pairwise.deletion = T)
			proteinlocusdist <- dist.aa(proteinlocus, pairwise.deletion = T)
			proteinlocusdist <- as.matrix(proteinlocusdist)
			proteinlocusdist <- proteinlocusdist / length(proteinlocus[1,])
		}
		#get average distance, correct order of all matrices
		avgdnalocusdist <- apply(dnalocusdist, 1, mean)
		nameorder <- names(avgdnalocusdist)
		avgproteinlocusdist <- apply(proteinlocusdist, 1, mean)[nameorder]
		avgdnadistsubset <- avgdnadist[nameorder]
		avgproteindistsubset <- avgproteindist[nameorder]
		#get difference
		dnadiff <- avgdnalocusdist - avgdnadistsubset
		proteindiff <- avgproteinlocusdist - avgproteindistsubset
		NdnaU <- length(FindOutliers(dnadiff)$upperO)
		NdnaL <- length(FindOutliers(dnadiff)$lowerO)
		dnaoutliers <- names(dnadiff[c(FindOutliers(dnadiff)$upperO, FindOutliers(dnadiff)$lowerO)])
		NproteinU <- length(FindOutliers(proteindiff)$upperO)
		NproteinL <- length(FindOutliers(proteindiff)$lowerO)
		proteinoutliers <- names(proteindiff[c(FindOutliers(proteindiff)$upperO, FindOutliers(proteindiff)$lowerO)])
		outliers <- union(dnaoutliers, proteinoutliers)
		if (length(outliers) / length(nameorder) < 0.5){	
			if (length(outliers)>0){
				dnalocus <- dnalocus[!(rownames(dnalocus) %in% outliers),]
				proteinlocus <- proteinlocus[!(rownames(proteinlocus) %in% outliers),]
				print(paste0("edited file: ", commonfiles[i], ", discarded ", length(outliers),
					", NTU: ", NdnaU, ", NTL: ", NdnaL, ", AAU: ", NproteinU, ", AAL: ", NproteinL))
			} else {
				print(paste0("unchanged file: ", commonfiles[i]))
			}
			write.FASTA(dnalocus, paste0(dnafile,".edited"))
			write.FASTA(proteinlocus, paste0(proteinfile,".edited"))
		} else {
			print(paste0("removed file: ", commonfiles[i], length(outliers),
				", NTU: ", NdnaU, ", NTL: ", NdnaL, ", AAU: ", NproteinU, ", AAL: ", NproteinL))
		}
	} else {
		print(paste0("removed empty file: ", commonfiles[i]))
	}

}

