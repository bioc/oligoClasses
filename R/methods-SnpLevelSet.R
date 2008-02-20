setMethod("chromosome", "SnpLevelSet",
          function(object){
            require("RSQLite") || stop("RSQLite package not available")            
            fs <- featureNames(object)
            sql <- "SELECT man_fsetid, chrom FROM featureSet WHERE man_fsetid LIKE 'SNP%'"
            tmp <- dbGetQuery(db(object), sql)
            idx <- match(fs, tmp[["man_fsetid"]])
            tmp[idx, "chrom"]
          })

setMethod("db", signature(object="DBPDInfo"),
          function(object) object@getdb())

setMethod("db", "SnpLevelSet",
          function(object) {
            require(annotation(object), character.only=TRUE) || stop(paste(annotation(object), "package not available"))
            get(annotation(object))@getdb()
        })

setMethod("position", "SnpLevelSet",
          function(object){
            require("RSQLite") || stop("RSQLite package not available")            
            fs <- featureNames(object)
            sql <- "SELECT man_fsetid, physical_pos FROM featureSet WHERE man_fsetid LIKE 'SNP%'"
            tmp <- dbGetQuery(db(object), sql)
            idx <- match(fs, tmp[["man_fsetid"]])
            tmp[idx, "physical_pos"]
          })

setMethod("combine", signature=signature(x="SnpLevelSet", y="SnpLevelSet"),
          function(x, y, ...){
            annot <- paste(sort(c(annotation(x), annotation(y))), sep=",")
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

