dev.off()
library(ape)
trees <- read.tree("RAxML_bootstrap.an070116")
#062716
#18S.fas - 50.7%
#18S.fas.fas - 52.6%
#D2.fas.fas - 2%
#D3.fas.fas - 0%
#D3 - 0%
#total 25.8%

#070116
#total 44%
#18S 55.8%
#D2 2%
#D3 5.7%
lst <- c("Ommatides","Hypselo","Williams", "Recti", "NrHyselo", "cf_Hypselo", "Glypto")
vec <- character(0)
for (x in trees[[1]]$tip.label){
  for (y in lst){
  #  print(startsWith(x,y))
  if (startsWith(x, y) == TRUE){
    print (x)
    vec <- c(vec, x)
  }
  }
}
vec

#count:
c1 <- character(0)
for (q in 1:length(trees)){
  c1 <- c(c1, as.character(is.monophyletic(trees[[q]], vec)))
}
length(c1[c1 =="TRUE"])/length(c1)*100
#
vec1 <- c("Hypsipterigidae_sp_Thailand_4207", "Ogeriinae_CostaRica_3363")
#see
par(mfrow=c(3,4),mar=c(1,1,1,1))
for (z in 1:12){
  #plot(trees[[z]], tip.color = ifelse(trees[[z]]$tip.label %in% vec, "red","black"), main=as.character(is.monophyletic(trees[[z]], vec)), cex=0.3)
  plot(trees[[z]], tip.color = ifelse(trees[[z]]$tip.label %in% vec, "red",ifelse(trees[[z]]$tip.label %in% vec1, "green","black")), main=as.character(is.monophyletic(trees[[z]], vec)), cex=0.3)
}
