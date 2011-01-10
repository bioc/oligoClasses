setMethod("bothStrands", "SnpFeatureSet",
          function(object){
            pkg <- annotation(object)
            set1 <- c("pd.mapping50k.xba240",
                      "pd.mapping50k.hind240",
                      "pd.mapping250k.sty",
                      "pd.mapping250k.nsp")
            set2 <- c("pd.genomewidesnp.5",
                      "pd.genomewidesnp.6")
            if (pkg %in% set1){
              return(TRUE)
            }else if (pkg %in% set2){
              return(FALSE)
            }else{
              return(NA)
            }
          })

setMethod("allele", "SnpFeatureSet",
          function(object, allele, strand){
              allele <- match.arg(allele, c("A", "B"))
              axiom <- length(grep("axiom", annotation(object))) == 1

              ## remember: select different than value below (works also
              ## on axiom)
              ac <- switch(allele, A=1, B=0)
              cc <- switch(allele, A='channel2', B='channel1')
              ccc <- setdiff(c('channel1', 'channel2'), cc)
              if (axiom){
                  sql <- paste('SELECT man_fsetid, fsetid, fid, allele, allelea, alleleb, count',
                               'FROM pmfeature',
                               'INNER JOIN featureSet USING(fsetid)',
                               'WHERE allele != ', ac)
                  tmp <- dbGetQuery(db(object), sql)
                  tmp <- tmp[order(tmp[['count']], tmp[['fsetid']], tmp[['fid']]),]
                  rownames(tmp) <- NULL

                  ## easy SNPs: [A/C] [A/G] [T/C] [T/G]
                  ##   alleleA: fid channel1
                  ##   alleleB: fid channel2
                  ## hard grn.: [A/T]
                  ##   alleleA: fid channel1
                  ##   alleleB: fid channel1 (probe2)
                  ## hard cyan: [C/G]
                  ##   alleleA: fid channel2
                  ##   alleleB: fid channel2 (probe2)
                  ## indel grn: [-/A] [-/T]
                  ##   need Angela
                  ##   alleleA: fid channel1
                  ##   alleleB: fid channel2
                  ## indel cyn: [-/C] [-/G]
                  ##   need Angela
                  ##   alleleA: fid channel2
                  ##   alleleB: fid channel1

                  easy <- c("AC", "AG", "TC", "TG")
                  hgrn <- "AT"
                  hcyn <- "CG"
                  indelg <- c("-A", "-T")
                  indelc <- c("-C", "-G")
                  types <- paste(tmp[['allelea']], tmp[['alleleb']], sep='')
                  i1 <- which(types %in% easy)
                  i2 <- which(types %in% hgrn)
                  i3 <- which(types %in% hcyn)
                  i4 <- which(types %in% indelg)
                  i5 <- which(types %in% indelc)

                  a1 <- assayDataElement(object, cc)[tmp[i1, 'fid'],,drop=FALSE]
                  a2 <- assayDataElement(object, cc)[tmp[i2, 'fid'],,drop=FALSE]
                  a3 <- assayDataElement(object, ccc)[tmp[i3, 'fid'],,drop=FALSE]
                  a4 <- assayDataElement(object, cc)[tmp[i4, 'fid'],,drop=FALSE]
                  a5 <- assayDataElement(object, ccc)[tmp[i5, 'fid'],,drop=FALSE]

                  result <- rbind(a1, a2, a3, a4, a5)
                  i <- c(i1, i2, i3, i4, i5)
                  tmp2 <- tmp[i,]
                  rownames(result) <- tmp2[, 'man_fsetid']
                  o <- order(tmp2[['fsetid']], tmp2[['fid']])
                  rm(tmp2)
                  result <- result[o,, drop=FALSE]
                  return(result)
              }else{
                  stop("Not yet implemented")
              }
              
              
          })
