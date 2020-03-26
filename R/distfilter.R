library(ape)
library(seqinr)
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  cat("Syntax: Rscript distfilter.R [path to concat dna] [path to dna folder] [path to concat protein] [path to protein folder]\n")
  cat("Example: Rscript distfilter.R concatdna.fas ./dna/ concataa.fas ./aa/\n")
  quit()
}
FindOutliers <- function(data) {
  lowerq = quantile(data, na.rm = T)[2]
  upperq = quantile(data, na.rm = T)[4]
  iqr = upperq - lowerq #Or use IQR(data)
  # we identify extreme outliers
  extreme.threshold.upper = (iqr * 3) + upperq
  extreme.threshold.lower = lowerq - (iqr * 3)
  result <- which(data > extreme.threshold.upper | data < extreme.threshold.lower)
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
dnadist <- dist.dna(dnaconcat, as.matrix = T, pairwise.deletion = T)
avgdnadist <- apply(dnadist, 1, mean)
proteindist <- dist.aa(proteinconcat, pairwise.deletion = T)
proteindist <- as.matrix(proteindist)
proteindist <- proteindist / length(proteinconcat[1,])
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
		dnalocusdist <- dist.dna(dnalocus, as.matrix = T, pairwise.deletion = T)
		proteinlocusdist <- dist.aa(proteinlocus, pairwise.deletion = T)
		proteinlocusdist <- as.matrix(proteinlocusdist)
		proteinlocusdist <- proteinlocusdist / length(proteinlocus[1,])
		#get average distance, correct order of all matrices
		avgdnalocusdist <- apply(dnalocusdist, 1, mean)
		nameorder <- names(avgdnalocusdist)
		avgproteinlocusdist <- apply(proteinlocusdist, 1, mean)[nameorder]
		avgdnadistsubset <- avgdnadist[nameorder]
		avgproteindistsubset <- avgproteindist[nameorder]
		#get difference
		dnadiff <- avgdnalocusdist - avgdnadistsubset
		proteindiff <- avgproteinlocusdist - avgproteindistsubset
		dnaoutliers <- names(dnadiff[FindOutliers(dnadiff)])
		proteinoutliers <- names(proteindiff[FindOutliers(proteindiff)])
		outliers <- union(dnaoutliers, proteinoutliers)
		if (length(outliers) / length(nameorder) < 0.5){	
			if (length(outliers)>0){
				dnalocus <- dnalocus[!(rownames(dnalocus) %in% outliers),]
				proteinlocus <- proteinlocus[!(rownames(proteinlocus) %in% outliers),]
				print(paste0("edited file: ", commonfiles[i], ", discarded ", length(outliers)))
			} else {
				print(paste0("unchanged file: ", commonfiles[i]))
			}
			write.FASTA(dnalocus, paste0(dnafile,".edited"))
			write.FASTA(proteinlocus, paste0(proteinfile,".edited"))
		} else {
			print(paste0("removed file: ", commonfiles[i]))
		}
	} else {
		print(paste0("removed empty file: ", commonfiles[i]))
	}

}

