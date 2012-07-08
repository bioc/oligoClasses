GRangesHMMList <- function(...)
{
    listData <- list(...)
    if (length(listData) == 0L) {
        unlistData <- GRangesHMM()
    } else {
        if (length(listData) == 1L && is.list(listData[[1L]]))
            listData <- listData[[1L]]
        if (!all(sapply(listData, is, "GRangesHMM")))
            stop("all elements in '...' must be GRangeHMM objects")
        unlistData <- suppressWarnings(do.call("c", unname(listData)))
    }
    end <- cumsum(elementLengths(unname(listData)))
    ans <- IRanges:::newCompressedList("GRangesHMMList",
               unlistData,
               end = end, NAMES = names(listData),
               elementMetadata = new("DataFrame", nrows = length(listData)))
    validObject(ans)
    ans
}
