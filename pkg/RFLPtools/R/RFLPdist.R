###############################################################################
## Computation of distances for RFLP data
###############################################################################

## x: data.frame with RFLP data
## distfun: function to compute distance (cf. ?dist)
RFLPdist <- function(x, distfun = dist, nrBands){
    stopifnot(is.data.frame(x))
    stopifnot(is.function(dist))
    x1 <- split(x, x$Sample)
    nrbands <- sort(unique(sapply(x1, nrow)))
    x1.bands <- sapply(x1, nrow)

    if(missing(nrBands)){
        res <- vector("list", length(nrbands))
        index <- 0
        for(i in nrbands){
            index <- index + 1
            temp <- do.call("rbind", x1[x1.bands == i])
            temp1 <- split(temp[,"MW"], factor(temp[,"Sample"]))
            res[[index]] <-  distfun(do.call("rbind", temp1))
        }
        names(res) <- nrbands
    }else{
        if(!(nrBands %in% nrbands))
            stop("No samples with given number of bands!")
        temp <- do.call("rbind", x1[x1.bands == nrBands])
        temp1 <- split(temp[,"MW"], factor(temp[,"Sample"]))
        res <-  distfun(do.call("rbind", temp1))
    }
    res
}
