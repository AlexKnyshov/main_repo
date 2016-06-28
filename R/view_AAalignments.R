library (ape)
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
#get file list
files <- list.files(path=".", pattern="*.fas", full.names=T, recursive=FALSE)
#read files
loci <- lapply(files, function(x)
{
  locus <- read.alignment(x, format = "fasta")
  locus <- as.matrix.alignment(locus)
  locus <- toupper(locus)
  rownames(locus) <- strtrim(rownames(locus), width = 3)
  locus <- as.AAbin.character(locus)
})
#set up plot parameters
par(mfrow=c(5,5),mar=c(1,1,1,1))
#plot alignments
n <- c(1:25)
for (x in n){
  image_alex.AAbin(loci[[x]], main=files[x], show.labels = T, legend = F, cex.lab=0.4, cex.axis=0.6, mgp=c(3,.1,0))
}
