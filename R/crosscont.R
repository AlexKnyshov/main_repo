library (ape)
args = commandArgs(trailingOnly=TRUE)
filepath <- args[1]
files <- list.files(path=filepath, pattern="*.fas", full.names=T, recursive=FALSE)
write("###", "out.tab")
for (f in files){
  print("---------------------------------------------------------------")
  print(f)
  locus <- read.dna(f, format="fasta", as.character = F)
  #print(dist.dna(locus, model = "raw", pairwise.deletion = T, as.matrix = T))
	mx <- dist.dna(locus, model = "raw", pairwise.deletion = T, as.matrix = T)
	rownames(mx) <- trimws(rownames(mx))
	colnames(mx) <- trimws(colnames(mx))
	diag(mx) <- NA
	#mx[lower.tri(mx, diag = T)] <- NA
	mx2 <- which(mx < 0.01, arr.ind = TRUE)
	if (length(mx2) > 0){
	  cout1 <- character(length = 0)
	  #print(unique(rownames(mx2)))
	  for (r in 1:nrow(mx2)){
	    if(grepl("Vescia", rownames(mx)[mx2[r,1]]) == T & grepl("Vescia", rownames(mx)[mx2[r,2]])== T){
	      print("Vescia, skipped")
	    }
	    else if(grepl("Opisthacidius", rownames(mx)[mx2[r,1]]) == T & grepl("Opisthacidius", rownames(mx)[mx2[r,2]])== T){
	      print("Opisthacidius, skipped")
	    }
	    else if(grepl("Phymata", rownames(mx)[mx2[r,1]]) == T & grepl("Phpen", rownames(mx)[mx2[r,2]])== T){
	      print("Phymata, skipped")
	    }
	    else if(grepl("Phymata", rownames(mx)[mx2[r,1]]) == T & grepl("Phymata", rownames(mx)[mx2[r,2]])== T){
	      print("Phymata, skipped")
	    }
	    else if(grepl("Phpen", rownames(mx)[mx2[r,1]]) == T & grepl("Phymata", rownames(mx)[mx2[r,2]])== T){
	      print("Phymata, skipped")
	    }
	    else{
	      print(paste(c(rownames(mx)[mx2[r,1]],rownames(mx)[mx2[r,2]]), sep=',', collapse = ',')) 
	      cout1 <- c(cout1,(rownames(mx)[mx2[r,1]]))
	    }
	  }
	  if (length(cout1)>0){
	    write(f, "out.tab", append = T)
	    write(unique(cout1), "out.tab", append = T)
	    drop1 <- match(unique(cout1),rownames(mx))
	    locus <- locus[-drop1,]
	    # for (n in unique(cout1)){
	    #   print(typeof(locus))
	    #   print(which(rownames(mx)==n))
	    #   locus[which(rownames(mx)==n),, drop=T]
	    # }
	    #
	  }
	  f1 <- paste("./ccout/", unlist(strsplit(f, "//"))[2], sep='')
	  write.FASTA(locus, f1)
	}
}