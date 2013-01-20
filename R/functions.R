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

getSequenceLengths <- function(build){
	path <- system.file("extdata", package="oligoClasses")
	load(file.path(path, paste("seqlengths_", build, ".rda", sep="")))
	return(seqlengths)
}

setSequenceLengths <- function(build, names){ ## names are unique(seqnames(object))
	sl <- getSequenceLengths(build)
	sl[match(unique(names), names(sl))]
}

chromosome2integer <- function(chrom){
	chrom[chrom == "X"] <- 23; chrom[chrom == "Y"] <- 24; chrom[chrom == "XY"] <- 25; chrom[chrom=="M" | chrom == "MT" | chrom == "Mt"] <- 26
	as.integer(chrom)
}
integer2chromosome <- function(intChrom){
	charChrom <- as.character(intChrom)
	charChrom[charChrom=="23"] <- "X"
	charChrom[charChrom=="24"] <- "Y"
	charChrom[charChrom %in% c("MT", "Mt")] <- "M"
	charChrom
}

.getArm <- function(chrom, pos, genome){
	if(is.integer(chrom)) chrom <- paste("chr", integer2chromosome(chrom), sep="")
	path.gap <- system.file("extdata", package="SNPchip")
	gaps <- readRDS(list.files(path.gap, pattern=paste("gap_", genome, ".rda", sep=""), full.names=TRUE))
	centromere.starts <- start(gaps)
	centromere.ends <- end(gaps)
	names(centromere.ends) <- names(centromere.starts) <- seqnames(gaps)
	centromere.starts <- centromere.starts[chrom]
	centromere.ends <- centromere.ends[chrom]
	chr.arm <- arm <- rep(NA, length(pos))
	arm[pos <= centromere.starts] <- "p"
	arm[pos >= centromere.ends] <- "q"
	##arm <- ifelse(pos <= centromere.starts, "p", "q")
	chr.arm[!is.na(arm)] <- paste(chrom[!is.na(arm)], arm[!is.na(arm)], sep="")
	chr.arm
}
