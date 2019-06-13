library (ape)
args = commandArgs(trailingOnly=TRUE)
filepath <- args[1]
files <- list.files(path=filepath, pattern="*.fas", full.names=T, recursive=FALSE)
write("###", "out.tab")
for (f in files){
  print("---------------------------------------------------------------")
  print(f)
  locus <- read.dna(f, format="fasta", as.character = F)
  locus_text <- read.dna(f, format="fasta", as.character = T)
  #print(dist.dna(locus, model = "raw", pairwise.deletion = T, as.matrix = T))
	mx <- dist.dna(locus, model = "raw", pairwise.deletion = T, as.matrix = T)
	rownames(mx) <- trimws(rownames(mx))
	colnames(mx) <- trimws(colnames(mx))
	diag(mx) <- NA
  loop_OK <- TRUE
  cout1 <- character(length = 0)
  while (loop_OK == TRUE){
	#mx[lower.tri(mx, diag = T)] <- NA
	loop_OK <- FALSE
	mx2 <- which(mx < 0.01, arr.ind = TRUE)
	if (length(mx2) > 0){
	  #print(unique(rownames(mx2)))
	  for (r in 1:nrow(mx2)){
	  	name1 <- unlist(strsplit(rownames(mx)[mx2[r,1]], "\\|"))
	  	name2 <- unlist(strsplit(rownames(mx)[mx2[r,2]], "\\|"))
	  	#print (names(locus_text))
	  	if (name1[1] == name2[1]){
	  		len1 <- length(locus_text[rownames(mx)[mx2[r,1]],][locus_text[rownames(mx)[mx2[r,1]],]!="-"])
	  		len2 <- length(locus_text[rownames(mx)[mx2[r,2]],][locus_text[rownames(mx)[mx2[r,2]],]!="-"])
	  		if (len1 >= len2){
	  			# print(paste0("good:_",rownames(mx)[mx2[r,1]]))
	  			# print(paste0("excluded:_",rownames(mx)[mx2[r,2]]))
	  			cout1 <- c(cout1,(rownames(mx)[mx2[r,2]]))
	  			mx <- mx[-mx2[r,2],-mx2[r,2]]
	  			loop_OK <- TRUE
	  			break
	  		}
	  		else {
	  			# print(paste0("good:_",rownames(mx)[mx2[r,2]]))
				# print(paste0("excluded:_",rownames(mx)[mx2[r,1]]))
				cout1 <- c(cout1,(rownames(mx)[mx2[r,1]]))
				mx <- mx[-mx2[r,1],-mx2[r,1]]
				loop_OK <- TRUE
	  			break
	  		}
	  	}
	    # if(grepl("Vescia", rownames(mx)[mx2[r,1]]) == T & grepl("Vescia", rownames(mx)[mx2[r,2]])== T){
	    #   print("Vescia, skipped")
	    # }
	    # else if(grepl("Opisthacidius", rownames(mx)[mx2[r,1]]) == T & grepl("Opisthacidius", rownames(mx)[mx2[r,2]])== T){
	    #   print("Opisthacidius, skipped")
	    # }
	    # else if(grepl("Phymata", rownames(mx)[mx2[r,1]]) == T & grepl("Phpen", rownames(mx)[mx2[r,2]])== T){
	    #   print("Phymata, skipped")
	    # }
	    # else if(grepl("Phymata", rownames(mx)[mx2[r,1]]) == T & grepl("Phymata", rownames(mx)[mx2[r,2]])== T){
	    #   print("Phymata, skipped")
	    # }
	    # else if(grepl("Phpen", rownames(mx)[mx2[r,1]]) == T & grepl("Phymata", rownames(mx)[mx2[r,2]])== T){
	    #   print("Phymata, skipped")
	    # }
	    # else{
	    #   print(paste(c(rownames(mx)[mx2[r,1]],rownames(mx)[mx2[r,2]]), sep=',', collapse = ',')) 
	    #   cout1 <- c(cout1,(rownames(mx)[mx2[r,1]]))
	    # }
	  }

	}
}
	  if (length(cout1)>0){
	    write(f, "out.tab", append = T)
	    write(unique(cout1), "out.tab", append = T)
	    drop1 <- match(cout1,rownames(locus_text))
	    #print(cout1)
	    locus <- locus[-drop1,]
	    # for (n in unique(cout1)){
	    #   print(typeof(locus))
	    #   print(which(rownames(mx)==n))
	    #   locus[which(rownames(mx)==n),, drop=T]
	    # }
	    #
	  }
	  dir.create("no_iso")
	  f1 <- paste("./no_iso/", unlist(strsplit(f, "//"))[2], sep='')
	  write.FASTA(locus, f1)
}