###############################################################################
## Computation of distances for RFLP data - alternative approach
###############################################################################

## compute number of bands
nrBands <- function(x){
    sort(unique(sapply(split(x, x$Sample), nrow)))
}

## x: data.frame with RFLP data
## distfun: function to compute distance (e.g. ?dist)
## nrBands: number of bands 
## nrMissing: number of bands which may be nrMissing
## diag: see ?dist
## upper: see ?dist
## compares samples with number of bands in: nrBands, nrBands + 1, ..., nrBands + nrMissing
RFLPdist2 <- function(x, distfun = dist, nrBands, nrMissing, diag = FALSE, upper = FALSE){
    stopifnot(is.data.frame(x))
    x1 <- split(x, x$Sample)
    nrbands <- sort(unique(sapply(x1, nrow)))
    x1.bands <- sapply(x1, nrow)

    temp <- do.call("rbind", x1[x1.bands %in% c(nrBands:(nrBands+nrMissing))])
    temp1 <- split(temp[,"MW"], factor(temp[,"Sample"]))
    N <- length(temp1)
    temp2 <- matrix(NA, ncol = nrBands + nrMissing, nrow = N)
    rownames(temp2) <- names(temp1)
    for(i in 1:N){
        temp2[i,1:length(temp1[[i]])] <- temp1[[i]]
    }
    d <- matrix(NA, nrow = N, ncol = N)
    dfun <- function(x, y){
        m <- sum(!is.na(x))
        min(as.matrix(dist(rbind(x[1:m], t(combn(y, m)))))[-1,1])
    }
    for(i in 1:N){
        for(j in 1:i){
            if(sum(!is.na(temp2[i,])) == sum(!is.na(temp2[j,]))){
                m <- sum(!is.na(temp2[i,]))
                d[i,j] <- as.vector(distfun(rbind(temp2[i,1:m], temp2[j,1:m])))
            }
            if(sum(!is.na(temp2[i,])) > sum(!is.na(temp2[j,])))
                d[i,j] <- dfun(temp2[j,], temp2[i,])
            if(sum(!is.na(temp2[i,])) < sum(!is.na(temp2[j,])))
                d[i,j] <- dfun(temp2[i,], temp2[j,])
        }
    }
    d <- d[lower.tri(d)]
    attr(d, "Size") <- N
    attr(d, "Labels") <- dimnames(temp2)[[1]]
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- distfun
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    return(d)
}
