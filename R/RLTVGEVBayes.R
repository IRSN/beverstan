
## ****************************************************************************
##' Compute return levels in relation with return periods.
##' 
##' @title Return Levels
##' 
##' @param object An object thar can be used to compute return levels.
##'
##' @param ... Further arguments for methods.
##'
##' @return A table of return levels, typically given as an object
##' inheriting from the \code{"data.frame"} class.
##' 
##' @export
RL <- function(object, ...) {
    UseMethod("RL")
}

## ****************************************************************************
##' Compute return levels with credible limits. The result can be used
##' to produce the classical \emph{return level plot} by using.
##' \code{\link{autoplot.TVGEVBayes}}.
##'
##' For a given period \eqn{T>1} e.g. \eqn{T = 100}, the return level
##' \eqn{\rho(T)} is defined as the quantile with exceedance
##' probability \eqn{1 / T}. In the time-varying framework, the
##' quantile is for a marginal GEV distribution hence relates to a
##' specific block in time which has to be given by using
##' \code{newTimeRange}. Hence the return level \eqn{\rho\{T;
##' \boldsymbol{\theta}(t)\}}{\rho[T; \theta(t)]} where
##' \eqn{\boldsymbol{\theta}(t)}{\theta(t)} is the vector of the thee
##' GEV parameters for the block \eqn{t}. This is deterministic
##' function of the GEV parameter and it can or should come along with
##' a credible interval.
##' 
##' @title Compute Return Levels with Credible Limits
##'
##' @param object A \code{TVGEVBayes} object.
##'
##' @param newTimeRange A time range to be passed to \code{timeRange}.
##'
##' @param period A vector of \emph{return periods}. 
##'
##' @param level Credible level to be used for the intervals.
##'
##' @param credintType The type of credible interval. \code{"HPD"}
##' corresponds to the \emph{Highest Posterior Density} interval and
##' \code{"eqtail"} corresponds to the "equal-tail" choice, where both
##' tails are given the same probability \eqn{(1 - level) / 2}.
##'
##' @param smooth Logical. If \code{TRUE}, the lower and upper credible
##' bounds will be smoothed by using \code{\link[stats]{smooth.spline}}.
##'
##' @param ... Not used yet.
##'
##' @return An object with class \code{"RL.GEVBayes"} inheriting from
##' \code{"data.frame"}.
##'
##' @note Since the class \code{"TVGEVBayes"} is devoted to time-varying models
##' it would be tedious (and potentially misleading) to derive
##' plotting positions for the observations used in the fit. So we do
##' not provide on this plot empirical points as usually shown in the
##' non time-varying (stationary) framework.
##' 
##' @seealso \code{\link[bever]{credInt}}.
##' 
##' @importFrom nieve qGEV
##' @importFrom stats median smooth.spline predict
##' @importFrom bever credInt formatLevel
##' @export
##' 
RL.TVGEVBayes <- function(object,
                          newTimeRange = NULL,
                          period = NULL,
                          level = 0.70,
                          credintType = c("HPD", "eqtail"),
                          smooth = missing(period),
                          ...) {

    newTimeRange <- newTimeRange(object, timeRange = newTimeRange)
    credintType <- match.arg(credintType)
    
    if (length(level) > 1) {
        warning("'level' can only be of length 1 for now. ",
                "Only the first element will considered")
        level <- level[1]
    }
    
    eps <-  1e-4
    fLevel <- bever::formatLevel(level)
    
    if (is.null(period)) {
        period <- c(2, 5, 10, 20, 50, 75, 100,
                    125, 150, 175, 200, 250, 300, 500, 700, 1000)
    }

    X <- modelMatrices(object = object$TVGEV,
                       date = newTimeRange$date[1L])$X
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
    for (j in 1:3) {
        theta[ , , j] <- MCMC[ , object$TVGEV$ind[[j]], drop = FALSE] %*%
            t(X[[j]])
    }
    
    if (!is.null(object$MAP)) {
        thetaMAP <- rep(NA_real_, 3)
        for (j in 1:3) {
            thetaMAP[j] <- object$MAP[object$TVGEV$ind[[j]]] %*%
                t(X[[j]])
        }
    }
            
    res <- array(NA, dim = c(length(period), 6),
                 dimnames = list(NULL,
                     c("Period", "Mode", "Median", "Mean", "L", "U")))
    
    for (i in seq_along(period)) {
        
        res[i, "Period"] <- period[i]
        
        x <- nieve::qGEV(1 - 1.0 / period[i],
                         loc = theta[ , 1L, 1L],
                         scale = theta[ , 1L, 2L],
                         shape = theta[ , 1L, 3L],
                         lower.tail = TRUE)
        
        res[i, "Mean"] <- mean(x, na.rm = TRUE)
        res[i, "Median"] <- median(x, na.rm = TRUE)
        LU <- bever::credInt(x, level = level, type = credintType)
        res[i, c("L", "U")] <- LU
          
        if (!is.null(object$MAP)) {
            res[i, "Mode"] <- nieve::qGEV(1 - 1.0 / period[i],
                                          loc = thetaMAP[1L],
                                          scale = thetaMAP[2L],
                                          shape = thetaMAP[3L],
                                          lower.tail = TRUE)
        } 
    }
    
    if (smooth) {
        for (nm in c("L", "U")) {
            ss <- smooth.spline(x = log(res[ , "Period"]), y = res[ , nm], df = 5)
            res[ , nm] <- predict(ss)$y
        }
    }
    
    res <- as.data.frame(res)
    res <- data.frame(Period = res$Period,
                      Level = fLevel,
                      Mode = res$Mode,
                      Median = res$Median,
                      Mean = res$Mean,
                      L = res$L, U = res$U)

    
    attr(res, "blockDuration") <- object[["bockDuration"]]
    attr(res, "newTimeRange") <- newTimeRange$timeRange
    attr(res, "level") <- level
    
    class(res) <- c("RL.GEVBayes", "data.frame")
    
    res
}

