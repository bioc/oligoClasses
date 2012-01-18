setAs("BeadStudioSet", "data.frame",
      function(from, to){
	      cn <- as.numeric(lrr(from))
	      bf <- as.numeric(baf(from))
	      x <- rep(position(from)/1e6, ncol(from))
	      ##x <- rep(position(object)[marker.index], 4)/1e6
	      is.snp <- rep(isSnp(from), ncol(from))
	      id <- rep(sampleNames(from), each=nrow(from))
	      df <- data.frame(x=x, lrr=cn, baf=bf, id=id,
			       is.snp=is.snp,
			       stringsAsFactors=FALSE)
	      df$id <- factor(df$id, ordered=TRUE, levels=unique(df$id))
	      return(df)
      })
