library(ape)
library(PhyInformR)
library(MESS)
args = commandArgs(trailingOnly=TRUE)
argsLen <- length(args);
if (argsLen == 0){
  cat("Syntax: Rscript phyinfo.R [path to rates file (iqtree format)] [path to partition file] \n")
  cat("Example: Rscript phyinfo.R rates.txt partitions.txt\n")
  quit()
}

inputrates <- read.table(args[1],header=TRUE)
mxrates <- as.matrix(inputrates)
inputprts <- read.table(args[2])
sets <- within(inputprts, V2<-data.frame(do.call('rbind', strsplit(as.character(V2), '=', fixed=TRUE))))

fitAll <- data.frame(set=character(),
                     x=numeric(),
                     ypos1=numeric(),
                     ypos2=numeric(),
                     ypos3=numeric(),
                     yall=numeric())

fitAll <- data.frame(matrix(ncol = 21, nrow = length(sets[,1])))
colnames(fitAll) <- c("set", "xMax", "xMaxP1", "xMaxP2", "xMaxP3",
					"yMaxStd", "yMaxStdP1", "yMaxStdP2", "yMaxStdP3",
					"yMax", "yMaxP1", "yMaxP2", "yMaxP3",
					"auc", "aucP1", "aucP2", "aucP3",
					"aucStd", "aucStdP1", "aucStdP2", "aucStdP3")

parts <- list(c(T,T,T), c(T,F,F), c(F,T,F), c(F,F,T))

for (prt in 1:length(sets[,1])){
  prtname <- as.character(unlist(sets[prt,2][1]))
  fitAll[prt, 1] = prtname
  setData <- list()
  for (partition in 1:length(parts)){
    fit1 <- data.frame(x1=numeric(),
                       y1=numeric(),
                       y1Std = numeric())
    for (t2a in seq(0,1,0.01)){
      phyinf = PhyInformR:::site.summer(mxrates[mxrates[,1]==prt,3][parts[[partition]]],t2a)
      t2b <- data.frame(x1=t2a,
                        y1=phyinf,
                        y1Std = phyinf/length(mxrates[mxrates[,1]==prt,3][parts[[partition]]]))
      fit1 <- rbind(fit1, t2b)
    }
    setData[[partition]] <- fit1
  }

  fitAll[prt, 2] = setData[[1]]$x1[which.max(setData[[1]]$y1)]
  fitAll[prt, 3] = setData[[2]]$x1[which.max(setData[[2]]$y1)]
  fitAll[prt, 4] = setData[[3]]$x1[which.max(setData[[3]]$y1)]
  fitAll[prt, 5] = setData[[4]]$x1[which.max(setData[[4]]$y1)]
  fitAll[prt, 6] = max(setData[[1]]$y1Std)
  fitAll[prt, 7] = max(setData[[2]]$y1Std)
  fitAll[prt, 8] = max(setData[[3]]$y1Std)
  fitAll[prt, 9] = max(setData[[4]]$y1Std)
  fitAll[prt, 10] = max(setData[[1]]$y1)
  fitAll[prt, 11] = max(setData[[2]]$y1)
  fitAll[prt, 12] = max(setData[[3]]$y1)
  fitAll[prt, 13] = max(setData[[4]]$y1)
  fitAll[prt, 14] = auc(setData[[1]]$x1, setData[[1]]$y1)
  fitAll[prt, 15] = auc(setData[[2]]$x1, setData[[2]]$y1)
  fitAll[prt, 16] = auc(setData[[3]]$x1, setData[[3]]$y1)
  fitAll[prt, 17] = auc(setData[[4]]$x1, setData[[4]]$y1)
  fitAll[prt, 18] = auc(setData[[1]]$x1, setData[[1]]$y1Std)
  fitAll[prt, 19] = auc(setData[[2]]$x1, setData[[2]]$y1Std)
  fitAll[prt, 20] = auc(setData[[3]]$x1, setData[[3]]$y1Std)
  fitAll[prt, 21] = auc(setData[[4]]$x1, setData[[4]]$y1Std)

  print(paste0("gene ", prtname))
}
write.csv(fitAll, "test.csv")