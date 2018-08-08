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
  files <- list.files(path=filepath, pattern="*.fas*", full.names=T, recursive=FALSE)  
} else if (opt == "-t")
{
  #files <- list.files(path=".", pattern="*bestTree*", full.names=T, recursive=FALSE)
  files <- list.files(path=filepath, pattern="RAxML_best*", full.names=T, recursive=FALSE)
} else if (opt == "-a")
{
  library(phangorn)
  library(seqinr) #read alignments instead of the default ape parser
  #####################################
  #reimplement image.AAbin to visualize more than three groups of AA
  #code was modified from APE 3.5
  image_alex.AAbin <- function(x, what, col, bg = "white", xlab = "", ylab = "",
                               show.labels = TRUE, cex.lab = 1, legend = TRUE, ...)
  {
    aa_l <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
              "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "X", "*")
    aa_l2 <- c(65, 67, 68, 69, 70, 71, 72, 73, 75, 76,
               77, 78, 80, 81, 82, 83, 84, 86, 87, 89, 88, 42)
    if (missing(what))
      what <- aa_l
    if (missing(col)) col <- c("blue", "pink", "purple",
                               "purple", "blue", "brown",
                               "cyan", "blue", "red",
                               "blue", "blue", "green",
                               "yellow", "green", "red",
                               "green", "green", "blue",
                               "blue", "cyan", "gray",
                               "black")
    n <- (dx <- dim(x))[1]
    s <- dx[2]
    y <- integer(N <- length(x))
    ncl <- length(what)
    col <- rep(col, length.out = ncl)
    brks <- 0.5:(ncl + 0.5)
    sm <- 0L
    for (i in ncl:1) {
      k <- aa_l2[i]#changed
      sel <- which(x == k)#location of a selected AA
      if (L <- length(sel)) {
        y[sel] <- i
        sm <- sm + L
      }
      else {
        what <- what[-i]
        col <- col[-i]
        brks <- brks[-i]
      }
    }
    dim(y) <- dx
    if (sm == N) {
      leg.co <- co <- col
      leg.txt <- what
    }
    else {
      co <- c(bg, col)
      print (co)
      leg.txt <- c(what, "Unknown")
      leg.co <- c(col, bg)
      brks <- c(-0.5, brks)
    }
    yaxt <- if (show.labels) "n" else "s"
    image.default(1:s, 1:n, t(y), col = co, xlab = xlab, ylab = ylab,
                  yaxt = yaxt, breaks = brks, ...)
    if (show.labels)
      mtext(rownames(x), side = 2, line = 0.1, at = 1:n, cex = cex.lab,
            adj = 1, las = 1)
    if (legend) {
      psr <- par("usr")
      xx <- psr[2]/2
      yy <- psr[4] * (0.5 + 0.5/par("plt")[4])
      legend(xx, yy, legend = leg.txt, pch = 22, pt.bg = leg.co,
             pt.cex = 2, bty = "n", xjust = 0.5, yjust = 0.5,
             horiz = TRUE, xpd = TRUE)
    }
  }
  #image end
  #####################################
  files <- list.files(path=filepath, pattern="*.fas*", full.names=T, recursive=FALSE)  
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
} else if(opt == "-a"){
  loci <- lapply(files, function(x)
  {
    locus <- read.alignment(x, format = "fasta")
    locus <- as.matrix.alignment(locus)
    locus <- toupper(locus)
    rownames(locus) <- strtrim(rownames(locus), width = 3)
    locus <- as.AAbin.character(locus)
  })
  mainlen <<- length(loci)
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
    image(loci[[x]], main=files[x], show.labels = T, legend = F, cex.lab=1, cex.axis=1, mgp=c(3,.1,0))
    } else if(opt=="-t"){
      plot(trees[[x]], main=names[x], label.offset=0.01, cex=1.6, edge.width=3)
    } else if(opt=="-a"){
      image_alex.AAbin(loci[[x]], main=files[x], show.labels = T, legend = F, cex.lab=1, cex.axis=1, mgp=c(3,.1,0))
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