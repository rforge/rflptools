diffDist <- function(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2){ 
    x.diff <- t(diff(t(x))) 
    dist(x.diff, method = method, diag = diag, upper = upper, p = p) 
}
