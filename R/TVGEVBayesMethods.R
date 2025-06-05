
##' @title Methods for the S3 Class \code{"TVGEVBayes"}
##' 
TVGEVBayesMethods <- NULL

## *************************************************************************
##' @param object A \code{TVGEVBayes} object.
##' @param type The type of coefficient.
##' @param ... Further arguments to be passed to methods.
##' 
##' @method coef TVGEVBayes
##' @export
##' @rdname TVGEVBayesMethods
coef.TVGEVBayes <- function(object, type = c("MAP", "postMean"), ...) {
    type <-  match.arg(type)
    if (type == "MAP") return(object$MAP)
    else return(object$postMean)
}

## *************************************************************************
##' @method vcov TVGEVBayes
##' @importFrom stats vcov cov
##' @export
##' @rdname TVGEVBayesMethods
##' 
vcov.TVGEVBayes <- function(object, ...) {
    MCMC <- object$MCMC
    d <- dim(MCMC)
    if (length(d) == 3) {
        MCMC <- MCMC[ , , 1:object$TVGEV$p, drop = FALSE]
        dim(MCMC) <- c(d[1] * d[2], d[3])
        dimnames(MCMC) <- list(NULL, object$TVGEV$parNames)
    } else if (length(d) == 2) {
        MCMC <- MCMC[ , 1:object$TVGEV$p, drop = FALSE]
    }
    cov(MCMC)  
}

## *************************************************************************
##' @method summary TVGEVBayes
##' @export
##' @rdname TVGEVBayesMethods
##' 
summary.TVGEVBayes <- function(object, ...) {
    o <- object
    o$coefMat <- rbind("MAP" = object$MAP,
                       "post. mean" = object$postMean,
                       "post. sd" = object$postSd,
                       "post. L"  = object$postL,
                       "post. U" = object$postU)
    class(o) <- "summary.TVGEVBayes"
    o
}

## *************************************************************************
##' @method print summary.TVGEVBayes
##' @export
print.summary.TVGEVBayes <- function(x, ...) {
    cat("o Call\n")
    print(x$call)
    cat("\no Parameters\n")
    print(t(round(x$coefMat, digits = 4)))
}

## *************************************************************************
##' @method print TVGEVBayes
##' @export
##' @rdname TVGEVBayesMethods
##' @param x A \code{TVGEVBayes} of object.
##' 
print.TVGEVBayes <- function(x, ...) {
    cat("Object with class \"TVGEVBayes\"\n")
    cat("\no Call\n")
    print(x$call)
    cat("\no MAP and posterior means\n")
    mat <- rbind("MAP" = x$MAP,
                 "post. mean" = x$postMean)
    print(t(round(mat, digits = 4)))
}

## *************************************************************************
##' Guess a default "new" time range for a \code{TVGEVBayes} object.
##'
##' Given a \code{TVGEVBayes} object we can find both the block
##' duration and the beginning and end of the last block. The default
##' new block is the next block with the same duration. Note that if
##' \code{object} is based on a block duration of \eqn{w > 2} years,
##' the underlying model does not know anything about the maxima over
##' blocks with smaller duration, hence we can not predict the maximum
##' for a block with duration \eqn{< w}.
##' 
##' @title Default "New" Time Range
##' 
##' @param object A \code{TVGEVBayes} object.
##'
##' @param timeRange A character with length one defining the new time
##' range. See the \code{\link{timeRange}} function.
##'
##' @param ... Not used yet.
##'
##' @return Alist with the two elements
##' \item{\code{timeRange} }{ A character with length one
##'     describing the time range as in \code{"2020_2019"}.
##' }
##' \item{\code{date} }{ A \code{Date} vector with length 2
##'     giving the beginning and the end of the time-range.
##' }
##'
##' @export
##' 
newTimeRange <- function(object, timeRange = NULL, ...) {
    
    if (is.null(timeRange)) {
        dnm <- object$TVGEV$date
        d <- object$data[[dnm]][nrow(object$data)]
        tR <-  seq(from = d, length.out = 2, by = object$blockDuration)
        newTimeRange <- paste(format(tR, "%Y"), collapse = "_")
        list(timeRange = newTimeRange, date = tR)
    } else {
        tR <- timeRange(timeRange)
        tR2Test <- seq(from = tR[1], length = 2, by = object$blockDuration)[2]
        if (tR[2] < tR2Test) {
            stop("'timeRange' has a duration which is shorter than ",
                 "the `object$blockDuration'")
        }
        list(timeRange = timeRange, date = tR)
    }
}

## **************************************************************************
##' Fitted values for a \code{TVGEVBayes} object. By default, the
##' fitted values for a given date/block are taken as the expectations
##' of the marginal distributions. However, marginal quantiles can
##' instead be used, for instance quantiles with probability 0.99
##' a.k.a. the 100-year return levels. In both cases the returned
##' vector contains the posterior mean of the chosen quantity and it
##' comes along with the limits of a credible interval.
##'
##' @title Fitted Values for a \code{TVGEVBayes} Object
##' 
##' @param object An object with class \code{"TVGEVBayes"}.
##'
##' @param which A location parameter for the marginal distribution,
##' either the expectation or a quantile.
##'
##' @param prob The probability to use when \code{which} is taken to
##' be \code{"quantile"}.
##'
##' @param level The credible level.
##'
##' @param ... Not used yet. 
##'
##' @return A list containing the vector \code{mean} and the two
##' vectors of bounds: \code{lower} and \code{upper}. An attribute
##' \code{"which"} is used to keep the trace of the choice made for
##' \code{which}.
##' 
##' @importFrom stats fitted quantile
##' @importFrom nieve qGEV
##' 
##' @export
##' 
##' @method fitted TVGEVBayes
##' 
fitted.TVGEVBayes <- function(object, which = c("expect", "quantile"),
                              prob = 0.99, level = 0.95,
                              ...) {
    
    which <- match.arg(which)
    
    MCMC <- object$MCMC
    d <- dim(MCMC)
    if (length(d) == 3) {
        MCMC <- MCMC[ , , 1:object$TVGEV$p, drop = FALSE]
        dim(MCMC) <- c(d[1] * d[2], d[3])
        dimnames(MCMC) <- list(NULL, object$TVGEV$parNames)
    } else if (length(d) == 2) {
        MCMC <- MCMC[ , 1:object$TVGEV$p, drop = FALSE]
    }

    theta <- array(0.0,  dim = list(nrow(MCMC), nrow(object$data), 3))
    dimnames(theta) <- list(NULL, NULL,  names(object$TVGEV$X))
    for (i in 1:3) {
        theta[ , , i] <- MCMC[ , object$TVGEV$ind[[i]], drop = FALSE] %*%
            t(object$TVGEV$X[[i]])
    }

    if (which == "expect") {
        g1 <- gamma(1 - theta[ , , "shape"])
        esp <- theta[ , , "loc"] +  (g1 - 1) * theta[ , , "scale"] /  theta[ , , "shape"]
    } else {
        esp <- nieve::qGEV(prob,
                           theta[ , , "loc"],
                           scale = theta[ , , "scale"],
                           shape = theta[ , , "shape"])
        ## print(length(esp))
        dim(esp) <- c(nrow(MCMC), nrow(object$data))
    }
        
    L <- list(mean = apply(esp, 2, mean),
              lower = apply(esp, 2, quantile, prob = (1 - level) / 2),
              upper =  apply(esp, 2, quantile, prob = 1 - (1 - level) / 2))
    attr(L, "which") <- which
    L
}


## ************************************************************************
##' Perform some computations about the predictive distribution from a
##' \code{TVGEVBayes} object. For a "new" period of time given with
##' \code{newTimeRange}, we can consider the maximum of the block maxima
##' on this period. The distribution is not a GEV distribution in general,
##' and due to the uncertainty on the parameters 
##' 
##' @title Prediction from a \code{TVGEVBayes} Object
##' 
##' @param object A \code{TVGEVBayes} object.
##'
##' @param newTimeRange A character describing the "new" time range
##' (period of time) as in \code{\link{timeRange}}.
##'
##' @param prob A vector of exceedance probabilities at which the
##' predictive distribution is to be computed.
##'
##' @param type The type of result wanted. For now only the predictive
##' Return Levels can be computed. 
##'
##' @param approx \code{Logical}. If \code{TRUE} the predictive
##' quantiles will be computed by using a linear interpolation method,
##' while a series of \code{\link{uniroot}} is used instead.
##'
##' @param trace Integer level of verbosity.
##'
##' @param ... Not used yet.
##'
##' @importFrom nieve pGEV
##' @importFrom NSGEV modelMatrices
##' @importFrom stats ppoints quantile uniroot
##' @export
##' @method predict TVGEVBayes
##' 
##' @examples  
##' newTimeRange <- "2020_2026"
##' 
predict.TVGEVBayes <- function(object,
                               newTimeRange = NULL,
                               prob,
                               type = "RL",
                               approx = FALSE,
                               trace = 0,
                               ...) {

    if (!is.null(newTimeRange)) {
        tR <- timeRange(newTimeRange)
    } else {
        nTR <- newTimeRange(object)
        newTimeRange <- nTR$timeRange
        tR <-  nTR$date
    }

    type <-  match.arg(type)
    
    if (missing(prob)) {
        prob <- c(0.900, 0.800, 0.600, 0.500, 0.250, 0.200, 0.150, 0.100,
                  0.050, 0.030, 0.025, 0.020, 0.015,
                  0.010, 0.008, 0.005, 0.004, 0.003, 0.002, 0.001)
    } else {
        if (any(prob < 0.0) || any(prob > 1.0)) {
            stop("'prob' must be a numeric vector with values between 0.0 and 1.0")
        }
        prob <- sort(prob, decreasing = TRUE)
    }
    nProb <- length(prob)

    blocks <-  blockCuts(tR, blockDuration = object$blockDuration)
    blocks <- blocks[-length(blocks)]
    mStar <- length(blocks)
    
    X <- NSGEV::modelMatrices(object = object$TVGEV, date = blocks)$X
    MCMC <- object$MCMC
    d <- dim(MCMC)
    
    if (length(d) == 3) {
        MCMC <- MCMC[ , , 1:object$TVGEV$p, drop = FALSE]
        dim(MCMC) <- c(d[1] * d[2], d[3])
        dimnames(MCMC) <- list(NULL, object$TVGEV$parNames)
    } else if (length(d) == 2) {
        MCMC <- MCMC[ , 1:object$TVGEV$p, drop = FALSE]
    }

    theta <- array(0.0,  dim = list(nrow(MCMC), nrow(X[[1]]), 3))
    dimnames(theta) <- list(NULL, NULL,  names(object$TVGEV$X))
    for (i in 1:3) {
        theta[ , , i] <- MCMC[ , object$TVGEV$ind[[i]], drop = FALSE] %*%
            t(X[[i]])
    }

    ## average over blocks
    thetaBar <- apply(theta, MARGIN = c(1, 3), mean)
    thetaMax <- apply(theta, MARGIN = 3, max)

    if (FALSE) {
        QQ <-  nieve::pGEV(q = 47,
                           loc = theta[ , , 1],
                           scale = theta[ , , 2],
                           shape = theta[ , , 3],
                           lower.tail = FALSE)
        
        dim(QQ) <- c(dim(MCMC)[1], length(blocks))
    }
    
    ## =======================================================================
    ## The posterior cumulative distribution function 
    ## =======================================================================
    FTilde <- function(q) {
        res <- rep(NA_real_, length(q))
        for (i in seq_along(res)) {
            Q <-  nieve::pGEV(q = q[i],
                              loc = theta[ , , 1],
                              scale = theta[ , , 2],
                              shape = theta[ , , 3],
                              lower.tail = TRUE)
            dim(Q) <- c(dim(MCMC)[1], length(blocks))
            res[i] <- mean(apply(Q, MARGIN = 1, FUN = prod))
        }
        res
    }
  
    quant <- rep(NA_real_, nProb)
    probTestU <- c(0.900, 0.950, 0.975, 0.980, 0.990, 0.999, 0.9999, 1 - 1e-5, 1-1e-6)
    probTestL <- c(0.300, 0.200, 0.100, 0.050, 0.010, 0.001, 0.0001, 1e-5, 1e-6)
    
    if (trace) {
        cat("Finding quantiles 'qL' and 'qU' corresponding to the greatest\n",
            " and the smallest exceedance probability.\n")
    }

    if (TRUE) {
        
        OKL <- OKU <-  FALSE
        
        for (i in seq_along(probTestU)) {
            if (!OKU) {
                yU <- quantile(nieve::qGEV(p = min(prob) / mStar,
                                           loc = thetaMax[1],   ## thetaBar[ , 1],
                                           scale = thetaMax[2], ## thetaBar[ , 2],
                                           shape = thetaMax[3], ## thetaBar[ , 3], 
                                           lower.tail = FALSE),
                               prob = probTestU[i])
            }
            if (!OKL) { 
                yL <- quantile(nieve::qGEV(p = max(prob) / mStar,
                                           loc = thetaBar[ , 1],
                                           scale = thetaBar[ , 2],
                                           shape = thetaBar[ , 3],
                                           lower.tail = FALSE),
                               prob = probTestL[i])
            }
            
            if (trace > 1) {
                cat(" yL = ", yL, " yU = ", yU, "\n")
            }
            
            if (!OKU) {
                qU <- try(uniroot(f = function(q) { FTilde(q) - 1 + min(prob) },
                                  interval = c(yL, yU)),
                          silent = TRUE)
                if (!inherits(qU, "try-error") && qU$estim.prec < 1e-4) {
                    OKU <-  TRUE
                    qU <- qU$root
                    if (trace) {
                        cat("At trial number", i, "'qU' is successfully localised: ",
                            "qU = " , qU, "\n") 
                    }
                }
                
            }
            if (!OKL) {
                qL <- try(uniroot(f = function(q) {FTilde(q) - 1 + max(prob)},
                                  interval = c(yL, yU)),
                      silent = TRUE)
                if (!inherits(qL, "try-error") && qL$estim.prec < 1e-4) {
                    OKL <-  TRUE
                    qL <- qL$root
                    if (trace) {
                        cat("At trial number", i, "'qL' is successfully localised: ",
                            "qL = ", qL, "\n")
                    }
                }  
            }
        }
        if (!OKL || !OKU) {
            stop("one of the minimal/maximal quantiles could not be localised")
        }
    } 
    
    quant[1L] <-  qL
    quant[nProb] <- qU
    
    if (nProb > 2) {
        
        ind <- (nProb - 1):2
        if (approx) {
            yGrid <- seq(from = qL, to = qU, length.out = 200)
            pGrid <- FTilde(yGrid)
            quant[ind] <- approx(x = pGrid, y = yGrid, xout = prob[ind])$y
        } else {
            
            for (i in ind) {
                q <- try(uniroot(f = function(q){FTilde(q) - 1 + prob[i]},
                                 interval = c(qL, qU)),
                         silent = TRUE)
                if (!inherits(q, "try-error") && q$estim.prec < 1e-2) {
                    OKL <-  TRUE
                    quant[i] <- q$root
                } else {
                    stop("problem in quantile determinantion. Try",
                         " 'approx = TRUE'")
                }
            }
        }
    }
    
    res <- data.frame(newTimeRange = newTimeRange, Prob = prob, Quant = quant)

    ## attr(res, "yMax") <- object$yMax
    ## attr(res, "newDuration") <- newDuration
    attr(res, "blockDuration") <- object$blockDuration
    
    class(res) <- c("predRLTVGEVBayes", "data.frame")

    res
    
}

##' *********************************************************************
##' Extract or build a prior from/for a Bayesian object.
##'
##' This method is intended to provide a description that is as close
##' as possible to that used in the prior spectification given at
##' the time of the object building.
##' 
##' @title Extract or Build a Prior from/for a Bayesian Object
##' 
##' @param object An object representing the results of a Bayesian
##' inference for some king of model.
##'
##' @param value An optional structure (usually a list) containing
##' a "default" prior for \code{object} that can be changed.
##'
##' @param ... Arguments to be passed to methods.
##'
##' @return A structure that describes a valid prior for the object
##' using the the information given in \code{value}.
##'
##' @export
##' 
## prior <- function(object, value, ...) {
##     UseMethod("prior")
## }

## ## ***************************************************************************
## ##' .. content for \description{} (no empty lines) ..
## ##'
## ##' @title Extract or Build a Prior from/for a Bayesian
## ##' \code{TVGEVBayes} Object
## ##'
## ##' @param object An object with class \code{TVGEVBayes}.
## ##'
## ##' @param value A list giving the parameters used in the prior.
## ##' 
## ##' @param type Not used yet.
## ##'
## ##' @param dist Not used yet.
## ##' 
## ##' @param ... Not used yet.
## ##' 
## ##' @return
## ##' 
## prior.TVGEVBayes <- function(object, value,
##                              type = "ind",
##                              dist = c(rep("norm", 3L)),
##                              ...) {
    
##     if (!is.list(value)) {
##         stop("'value' must be a list with suitable names")
##     }
    
##     priorType <- rep(c("mean", "cov"), times = 3L)
##     priorPar <- rep(c("mu", "sigma", "xi"), each = 2L)
##     nmsPrior <-  paste(priorType, rep("psi", 6L), priorPar, sep = "_")
##     priorList <- list()
##     pp <- fit(TVGEV$pp)
##     names(pp) <- c("mu", "sigma", "xi")
##     varProv <- c(1e7, 1e7, 10)
##     names(varProv) <- c("mu", "sigma", "xi")
    
##     if (!is.null(value)) {
        
##         if (!is.list(value) || !all(names(value) %in% nmsPrior)) {
##             msg <- paste("c(\"", paste(nmsPrior, collapse = "\", \""), "\")", sep = "")
##             stop("'prior' must be a named list with names in  \n",
##                  msg, "\n")
##         }
        
##         for (i in seq_along(nmsPrior)) {
            
##             nmi <- nm[i]                ## avoid multiple blackets
##             pi <- pp[priorPar[i]]       ## number of 'psi" parameters
##             typi <- priorType[i]        ## "mean" or "cov"?]
##             pari <- priorPar[i]         ## "mean" or "cov"?

##             if (nmi %in% names(prior)) {
##                 item <- prior[[nmi]]
##                 if (typi == "mean") {
##                     if (length(item) != pi) stop("bad length in prior for ", nmi) 
##                     else priorList[[nmi]] <- item
##                 } else if (typi == "cov") {
##                     if (pi > 1L) {
##                         if (!is.array(item) || !all.equal(dim(item), c(pi, pi)))
##                             stop("In prior", nmi, " must be a squre matrix  with ",
##                                  "size ", pi)
##                         else item <- matrix(item, nrow = 1L, ncol = 1L)
##                     }
##                     priorList[[nmi]] <- item
##                 }
##             } else {
##                 if (typi == "mean") priorList[[nmi]] <- rep(0.0, pi)
##                 else if (typi == "cov") {
##                     priorList[[nmi]] <- matrix(varProv[pari], nrow = pi, ncol = pi)
##                 }
##             }
##         }
##     }
    

## }
