setMethod("db", "SnpCallSetPlus",
          function(object) {
		  require(annotation(object), character.only=TRUE) || stop(paste(annotation(object), "package not available"))
		  get(annotation(object))@getdb()
	  })



