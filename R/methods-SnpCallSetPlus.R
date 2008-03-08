##setMethod("getA", "SnpCallSetPlus",
##          function(obj){
##            tmp <- array(NA, dim=c(nrow(antisenseThetaA(obj)),
##                               ncol(antisenseThetaA(obj)), 2),
##                         dimnames=list(rownames(antisenseThetaA(obj)),
##                           colnames(antisenseThetaA(obj)),
##                           c("antisense", "sense")))
##            tmp[,,1] <- .5*(antisenseThetaA(obj)+antisenseThetaB(obj))
##            tmp[,,2] <- .5*(senseThetaA(obj)+senseThetaB(obj))
##            return(tmp)
##          })

setMethod("calculateCopyNumber", "SnpCallSetPlus",
          function(object,
		   center.autosomes=2,
		   center.X.male=1,
		   center.X.female=2,
		   center.Y.male=1,
		   center.Y.female=0.5){
		  require(oligo) || stop("oligo package is not available")
		  annotationPackages <- unlist(strsplit(annotation(object), ","))
		  require(annotationPackages[1], character.only=TRUE) || stop(paste(annotationPackages[1], "package not available"))
		  if(length(annotationPackages) > 1){
			  require(annotationPackages[2], character.only=TRUE) || stop(paste(annotationPackages[2], "package not available"))
		  }
		  available.platforms <- affyPlatforms()		  
		  if(all(annotationPackages %in% available.platforms[c(1:4, 7,8)])){
			  A <- getA(object)
			  log2cn <- rowMeans(A, dims = 2, na.rm = TRUE)		
		  }
		  if(all(annotationPackages %in% available.platforms[c(5, 6)])){
			  A <- getA(object)
			  log2cn <- rownames(A, na.rm=TRUE)
			  ##                cm <- rowMeans(log2cn)
			  ##center the samples to have a log2 intensity of 1
			  ##(assumes on average the copy number is 2)                
			  ##                log2cn.center <- sweep(log2cn, 2, cm) + log2(2)
		  }


		  ##Homozygous and heterozygous calls have different overall means.
		  ##Seems reasonable to assume that the average copy numbers are the same		  
		  chr.matrix <- matrix(chromosome(object), ncol=ncol(log2cn), nrow=nrow(log2cn))
		  median.hom <- median(log2cn[(calls(object) == 1 | calls(object) == 3) & chr.matrix != "X"], na.rm=TRUE)
		  median.het <- median(log2cn[calls(object) == 2 & chr.matrix != "X"], na.rm=TRUE)

		  ##For each column, subtract off the overall median copy number
		  recenterByGenotype <- function(x, object, recenter.hom, recenter.het){
			  calls <- as.vector(calls(object))
			  x[calls == 1 | calls == 3] <- x[calls ==1 | calls == 3] - recenter.hom
			  x[calls == 2] <- x[calls == 2] - recenter.het
			  x
		  }
		  for(j in 1:ncol(log2cn)){
			  log2cn[, j] <- recenterByGenotype(log2cn[, j], object[, j], recenter.hom=median.hom, recenter.het=median.het)
		  }

		  ## Next, we sweep out a robust estimate of the
		  ## median from the samples (tries to put
		  ## fluorescence intensities on a similar scale for
		  ## each of the samples)
		  f <- function(x, chromosome){
			  tmp2 <- split(x, chromosome)
			  if(length(tmp2) > 15){
				  idx <- order(sapply(tmp2, "median"))
				  tmp2 <- tmp2[idx]
				  tmp3 <- tmp2[-c(1:5, (length(tmp2)-4):length(tmp2))]
				  med <- median(unlist(tmp3))
			  } else {
				  med <- median(sapply(tmp2, "median"))
			  }
			  return(med)
		  }		  
		  robust.median <- apply(log2cn, 2, f, chromosome(callset))
		  log2cn <- sweep(log2cn, 2, robust.median)
		  
		  rowSweep <- function(object, X, value, recenter, j){
			  if(length(value) == 1){
				  i <- chromosome(object) == value
			  } else {
				  i <- chromosome(object) %in% value
			  }
			  i[is.na(i)] <- FALSE
			  if(sum(i) > 1){
				  if(!missing(j)){
					  if(sum(j) < 5) warning("very few samples for calculating a robust average")
					  avg <- rowMedians(X[i, j], na.rm=TRUE)
					  X[i, j] <- sweep(X[i, j], 1, avg) + recenter
				  } else{
					  avg <- rowMedians(X[i, ], na.rm=TRUE)
					  X[i, ] <- sweep(X[i, ], 1, avg) + recenter
				  }
			  }
			  X
		  }
		  male <- object$gender == 1
		  female <- object$gender == 2
		  chromosome(object)[is.na(chromosome(object))] <- "NA"		  
		  log2cn <- rowSweep(object, log2cn,  "NA", log2(2))
		  log2cn <- rowSweep(object, log2cn,  "M", log2(2))		  
		  log2cn <- rowSweep(object, log2cn,  "X", log2(center.X.male), male)
		  log2cn <- rowSweep(object, log2cn,  "X", log2(center.X.female), female)
		  log2cn <- rowSweep(object, log2cn,  "Y", log2(center.Y.male), male)
		  log2cn <- rowSweep(object, log2cn,  "Y", log2(center.Y.female), female)		  
		  log2cn <- rowSweep(object, log2cn,  as.character(1:22), log2(center.autosomes))
		  cn <- 2^log2cn
		  return(cn)
	  })

setMethod("db", "SnpCallSetPlus",
          function(object) {
		  require(annotation(object), character.only=TRUE) || stop(paste(annotation(object), "package not available"))
		  get(annotation(object))@getdb()
	  })



