setMethod("calculateCopyNumber", "SnpCallSetPlus",
          function(object){
            require(oligo) || stop("oligo package is not available")
            available.platforms <- affyPlatforms()
            if(!(annotation(object) %in% available.platforms)){
              cn <- matrix(NA, nrow=nrow(object), ncol=ncol(object))
            } else {
              if(annotation(object) %in% available.platforms[c(1:4, 7,8)]){
                A <- getA(object)
                cn <- rowMeans(A, dims = 2, na.rm = TRUE)
                cm <- colMeans(cn)

                ##center the samples to have a log2 intensity of 1
                ##(assumes on average the copy number is 2)
                cn.center <- sweep(cn, 2, cm) + log2(2)
              }
              if(annotation(object) %in% available.platforms[c(5, 6)]){
                A <- getA(object)
                cn <- rownames(A, na.rm=TRUE)
                cm <- rowMeans(cn)

                ##center the samples to have a log2 intensity of 1
                ##(assumes on average the copy number is 2)                
                cn.center <- sweep(cn, 2, cm) + log2(2)
              }
            }
            rownames(cn) <- featureNames(object)
            colnames(cn) <- sampleNames(object)
            cn
          })
