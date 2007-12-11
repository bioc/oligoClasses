setMethod("chromosome", "SnpLevelSet",
          function(object){
            fs <- featureNames(object)
            sql <- "SELECT man_fsetid, chrom FROM featureSet WHERE man_fsetid LIKE 'SNP%'"
            tmp <- dbGetQuery(db(object), sql)
            idx <- match(fs, tmp[["man_fsetid"]])
            tmp[idx, "chrom"]
          })

setMethod("db", "SnpLevelSet",
          function(object) {
            require(annotation(object), character.only=TRUE) || stop(paste(annotation(object), "package not available"))
            db(get(annotation(object)))
        })

setMethod("position", "SnpLevelSet",
          function(object){
            fs <- featureNames(object)
            sql <- "SELECT man_fsetid, physical_pos FROM featureSet WHERE man_fsetid LIKE 'SNP%'"
            tmp <- dbGetQuery(db(object), sql)
            idx <- match(fs, tmp[["man_fsetid"]])
            tmp[idx, "physical_pos"]
          })
