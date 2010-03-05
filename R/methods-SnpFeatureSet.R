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
