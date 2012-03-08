integerMatrix <- function(x, scale=100) {
	if(!is(x, "matrix")) stop("argument x must be a matrix")
	dms <- dimnames(x)
	if(scale != 1){
		xx <- as.integer(x*scale)
	} else xx <- as.integer(x)
	x <- matrix(xx, nrow(x), ncol(x))
	dimnames(x) <- dms
	return(x)
}

numericMatrix <- function(x, scale=1/100) {
	return(x/scale)
}

integerArray <- function(x, scale=100){
	if(!is(x, "array")) stop("argument x must be an array")
	dims <- dim(x)
	dms <- dimnames(x)
	if(scale != 1){
		xx <- as.integer(x*scale)
	} else xx <- as.integer(x)
	x <- array(xx, dim=dims, dimnames=dms)
	return(x)
}

