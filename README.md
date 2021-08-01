# bigPLS
PLS regression (linear) or generalized linear regression (glm) for big data



# Cox

plsRcox
plsDR
splsDR
gplsDR
sgplsDR
sgPLS::gPLS
sgPLS::sgPLS
+ DK versions?
  
# regular PLS
  
# glm PLS

  
imputePLS






Final fit of the model.

# big data Cox (via SGD)
bigSurvSGD


# SGD cox
coxphSGD







function (x, y, strata, offset, init, control, weights, method, 
    rownames, resid = TRUE) 
{
    n <- nrow(y)
    if (is.matrix(x)) 
        nvar <- ncol(x)
    else {
        if (length(x) == 0) 
            nvar <- 0
        else nvar <- 1
    }
    time <- y[, 1]
    status <- y[, 2]
    if (length(strata) == 0) {
        sorted <- order(time)
        strata <- NULL
        newstrat <- as.integer(rep(0, n))
    }
    else {
        sorted <- order(strata, time)
        strata <- strata[sorted]
        newstrat <- as.integer(c(1 * (diff(as.numeric(strata)) != 
            0), 1))
    }
    if (missing(offset) || is.null(offset)) 
        offset <- rep(0, n)
    if (missing(weights) || is.null(weights)) 
        weights <- rep(1, n)
    else {
        if (any(weights <= 0)) 
            stop("Invalid weights, must be >0")
        weights <- weights[sorted]
    }
    stime <- as.double(time[sorted])
    sstat <- as.integer(status[sorted])
    storage.mode(weights) <- storage.mode(init) <- "double"
    if (nullmodel) {
            score <- exp(offset[sorted])
            coxres <- .C(Ccoxmart, as.integer(n), as.integer(method == 
                "efron"), stime, sstat, newstrat, as.double(score), 
                as.double(weights), resid = double(n))
            resid <- double(n)
            resid[sorted] <- coxres$resid
            names(resid) <- rownames
            list(residuals = resid)
    }
}