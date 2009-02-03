setMethod("db", "SnpLevelSet",
          function(object) {
            require(annotation(object), character.only=TRUE) || stop(paste(annotation(object), "package not available"))
            get(annotation(object))@getdb()
        })

setMethod("chromosome", "SnpLevelSet",
          function(object){
            if(!("chromosome" %in% fvarLabels(object))){
              require("RSQLite") || stop("RSQLite package not available")
              fs <- featureNames(object)
              sql <- "SELECT man_fsetid, chrom FROM featureSet WHERE man_fsetid LIKE 'SNP%'"              
              ##Check if two objects have been combined
              pkgs <- strsplit(annotation(object), ",")[[1]]
              if(length(pkgs) > 1){
                object2 <- object1 <- object
                annotation(object1) <- pkgs[1]
                annotation(object2) <- pkgs[2]

                tmp1 <- dbGetQuery(db(object1), sql)
                tmp2 <- dbGetQuery(db(object2), sql)
                tmp <- rbind(tmp1, tmp2)
              } else {
                tmp <- dbGetQuery(db(object), sql)
              }
              idx <- match(fs, tmp[["man_fsetid"]])
              chr <- tmp[idx, "chrom"]
            } else {
              chr <- featureData(object)$chromosome
            }
            return(chr)
          })

setReplaceMethod("chromosome", c("SnpLevelSet", "character"),
		 function(object, value){
			 chromosome(featureData(object)) <- value
			 object
		 })
setReplaceMethod("chromosome", c("AnnotatedDataFrame", "character"),
		 function(object, value){
			 pData(object)$chromosome <- value
			 object
		 })

setMethod("position", "SnpLevelSet",
          function(object){
		  if(!("position" %in% fvarLabels(object))){
			  require("RSQLite") || stop("RSQLite package not available")
			  fs <- featureNames(object)
			  sql <- "SELECT man_fsetid, physical_pos FROM featureSet WHERE man_fsetid LIKE 'SNP%'"
			  pkgs <- strsplit(annotation(object), ",")[[1]]
			  if(length(pkgs) > 1){
				  object2 <- object1 <- object
				  annotation(object1) <- pkgs[1]
				  annotation(object2) <- pkgs[2]                
              
				  tmp1 <- dbGetQuery(db(object1), sql)
				  tmp2 <- dbGetQuery(db(object2), sql)
				  tmp <- rbind(tmp1, tmp2)
			  } else{
				  tmp <- dbGetQuery(db(object), sql)
			  }
			  idx <- match(fs, tmp[["man_fsetid"]])
			  pos <- tmp[idx, "physical_pos"]
		  } else {
			  pos <- featureData(object)$position
		  }
		  return(pos)
          })

setMethod("combine", signature=signature(x="SnpLevelSet", y="SnpLevelSet"),
          function(x, y, ...){
		  ##Check that both x and y are valid objects
		  if(!validObject(x)) stop("x is not a valid object")
		  if(!validObject(y)) stop("y is not a valid object")
		  annot <- paste(sort(c(annotation(x), annotation(y))), collapse=",")
		  annotation(x) <- annotation(y) <- annot

		  if(class(x) != class(y)){
			  stop("objects must have the same class")
		  }            
		  if(storageMode(assayData(x)) != storageMode(assayData(y))){
			  stop("objects must have same storage mode for assayData")
		  }

		  fd <- combine(featureData(x), featureData(y))
		  pd <- combine(phenoData(x), phenoData(y))            
		  ad.x <- as.list(assayData(x))
		  ad.y <- as.list(assayData(y))
		  ad.xy <- mapply(rbind, ad.x, ad.y, SIMPLIFY=FALSE)
		  id.x <- match(rownames(ad.xy[[1]]), featureNames(fd))
		  ee <- combine(experimentData(x), experimentData(y))
		  assayData(x) <- ad.xy
		  storageMode(assayData(x)) <- storageMode(assayData(y))            
		  experimentData(x) <- ee
		  featureData(x) <- fd
		  phenoData(x) <- pd
		  x
          })

