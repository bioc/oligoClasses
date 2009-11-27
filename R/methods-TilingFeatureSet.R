setMethod("getM", "TilingFeatureSet",
          function(object){
            chn <- channelNames(object)
            if (!(chn %in% c("channel1", "channel2")))
              stop("Object does not have 'channel1' and 'channel2'.")
            lc1 <- log2(assayDataElement(object, "channel1"))
            lc2 <- log2(assayDataElement(object, "channel2"))
            lc1-lc2
          })

setMethod("getA", "TilingFeatureSet",
          function(object){
            chn <- channelNames(object)
            if (!(chn %in% c("channel1", "channel2")))
              stop("Object does not have 'channel1' and 'channel2'.")
            lc1 <- log2(assayDataElement(object, "channel1"))
            lc2 <- log2(assayDataElement(object, "channel2"))
            (lc1+lc2)/2
          })

setMethod("pmPosition", "TilingFeatureSet",
          function(object){
            conn <- db(object)
            tmp <- dbGetQuery(conn, "SELECT fid, position FROM pmfeature")
            tmp <- tmp[order(tmp[["fid"]]),]
            tmp[["position"]]
          })
