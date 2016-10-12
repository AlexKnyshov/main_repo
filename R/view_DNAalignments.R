library (ape)
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  cat("Syntax: Rscript view_DNAalignments.R [option] [path to files] [step]\n")
  cat("Example: Rscript view_DNAalignments.R -d ./ 25\n")
  quit()
}
opt <- args[1]
filepath <- args[2]
cat("load the file list\n")
if (opt == "-d")
{
  files <- list.files(path=filepath, pattern="*.fas", full.names=T, recursive=FALSE)  
} else if (opt == "-t")
{
  files <- list.files(path=".", pattern="*bestTree*", full.names=T, recursive=FALSE)
}
#init start
s = 1
#init step
step_a <- as.numeric(args[3])
cat("read files...\n")
mainlen = 0
if (opt == "-d")
{loci <- lapply(files, function(x)
{
  locus <- read.dna(x, format="fasta", as.character = T)
  locus[order(rownames(locus)),]
  rownames(locus) <- strtrim(rownames(locus), width = 3)
  locus <- as.DNAbin(locus)
})
mainlen <<- length(loci)
} else if(opt == "-t"){
  trees <- lapply(files, function(x)
  {
    tree <- read.tree(x)
  })
  mainlen <<- length(trees)
}
#main body
if (mainlen>step_a) #if num files > step
{
  side <- sqrt(step_a)
  side <- ceiling(side)
  if (mainlen>=side*side){ #if num files > adjusted matrix
    step_a <- side*side #step = matrix
    e = step_a
  } else{ #step = num files
    step_a <- mainlen
    e = step_a
  }
}else{ #if num files <= step
e = mainlen
step_a <- e - s + 1
side <- sqrt(step_a)
side <- ceiling(side)
}
if (step_a!=args[3]){
cat(paste("step was adjusted to", step_a, "; plot dimensions are", side, "by", side, "\n"))
}else{
  cat(paste("plot dimensions are", side, "by", side, "\n"))
}
#e = step_a
# if (opt == "-d"){
# t = ceiling(length(loci)/step_a)
# } else if (opt == "-t"){
#   t = ceiling(length(trees)/step_a)
# }
t = ceiling(mainlen/step_a)
for (y in c(1:t))
{
  fname = paste("loci",y,".png", collapse = "")
  png(fname, 2048, 1335)
  par(mfrow=c(side,side), mar=c(1,1,1,1))
  cat(paste("file", y, ", starting alignment", s, ", ending alignment", e, "\n"))
  n <- c(s:e)
  pb <- txtProgressBar(min = 0, max = step_a, style = 3)
  counter <- 1
  if (opt=="-t"){
    #names <- sapply(strsplit(sapply(strsplit(files, split='bestTree.', fixed=TRUE), function(x) (x[2])), split='.fas', fixed=TRUE), function(x) (x[1]))
    names <- files
  }
  for (x in n){
    if (opt == "-d"){
    image(loci[[x]], main=files[x], show.labels = T, legend = F, cex.lab=0.4, cex.axis=0.6, mgp=c(3,.1,0))
    } else if(opt=="-t"){
      plot(trees[[x]], main=names[x], label.offset=0.01, cex=1.6, edge.width=3)
    }
    setTxtProgressBar(pb, counter)
    counter <<- counter + 1
  }
  close(pb)
  s = s+step_a
  if (e+step_a <= mainlen){
    e = e+step_a
  } else{
    e = mainlen
    step_a <- e - s + 1
    side <- sqrt(step_a)
    side <- ceiling(side)
    if (step_a != 0){
    cat(paste("step was adjusted to", step_a, "; plot dimensions are", side, "by", side, "\n"))
      }
  }
  dev.off()
}
# s = t*step_a+1
# if (opt == "-d"){
#   e = length(loci)
# } else if (opt == "-t"){
#   e = length(trees)
# }
# step_a <- e - s + 1
# side <- sqrt(step_a)
# side <- ceiling(side)
# cat(paste("step was adjusted to", step_a, "; plot dimensions are", side, "by", side, "\n"))
# fname = paste("loci",t+1,".png", collapse = "")
# png(fname, 2048, 1335)
# par(mfrow=c(side,side), mar=c(1,1,1,1))
# cat(paste("file", t+1, ", starting alignment", s, ", ending alignment", e, "\n"))
# n <- c(s:e)
# pb <- txtProgressBar(min = 1, max = step_a, style = 3)
# counter <- 1
# if (opt=="-t"){
#   #names <- sapply(strsplit(sapply(strsplit(files, split='bestTree.', fixed=TRUE), function(x) (x[2])), split='.fas', fixed=TRUE), function(x) (x[1]))
#   names <- files
# }
# for (x in n){
#   if (opt == "-d"){
#     image(loci[[x]], main=files[x], show.labels = T, legend = F, cex.lab=0.4, cex.axis=0.6, mgp=c(3,.1,0))
#   } else if(opt=="-t"){
#     plot(trees[[x]], main=names[x], label.offset=0.01, cex=1.6, edge.width=3)
#   }
#   setTxtProgressBar(pb, counter)
#   counter <<- counter + 1
# }
# close(pb)
#dev.off()
cat("done\n")