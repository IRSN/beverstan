##' @description Compute confidence limits by using a sampled vector
##'     of parameter values \eqn{\psi^{[k]}}{psi[k]} (usually MCMC
##'     iterates) along with the corresponding values
##'     \eqn{\ell(\psi^{[k]})}{ell(psi[k])} of the log-likelihood. For
##'     each given confidence level \eqn{\gamma_i}, the confidence
##'     limits are the values of \eqn{\psi} for which the
##'     log-likelihood \eqn{\ell(\psi)} is equal to
##'     \eqn{\ell_{\text{max}} - \delta_i}{ell_max - delta[i]} where
##'     \eqn{\ell_{\text{max}}}{ell_max} is the maximal log-likelihood
##'     and
##'     \eqn{\delta_i:= q_{\chi^2(1)}(\gamma_i) / 2.0}{delta_i = qchisq(gamma_i, df = 1) / 2.0}.
##'
##' @details First, a threshold \eqn{u} for the log-likelihood is
##'     chosen below the smallest \eqn{\ell_{\text{max}} -
##'     \delta_i}{ell_max - delta[i]}. The points \eqn{[\psi^{[k]},
##'     \ell^{[k]}]} falling above \eqn{u} are selected and their
##'     convex hull is computed. The upper part of the convex hull
##'     provides the profile curve as a continous broken line. The two
##'     confidence limits are obtained as the interesections of the
##'     profile with the horizontal line. Inasmuch as the points have
##'     a good coverage of the region of the parameter space where the
##'     log-likelhood is high, we get a good approximation ot the
##'     profile likelihood and of the confidence limits.
##' 
##' @title Compute Confidence Limits using Sampled Vectors of
##'     Parameters and Log-Likelihood Values
##' 
##' @param param A numeric vector giving the values of the parameter
##'     of interest.
##' 
##' @param logLik A numeric vector containing the log-likelihood.
##'
##' @param level A vector of confidence levels.
##' 
##' @param plot An integer giving the "level of verbosity" of the
##'     plot. The value \code{0L} corresponds to no plot shown.
##' 
##' @return A matrix with two lines and one column by confidence
##'     level. The columns correspond to the confidence levels in
##'     increasing order so the confidence interval becomes larger as
##'     the column index increases.
##'
##' @importFrom grDevices chull adjustcolor
##' @importFrom stats qchisq approx
##' @importFrom graphics  abline points text
##' @export
##' @keywords internal
##' 
##' @examples
##' df <- within(TXMax_Dijon, Date <- as.Date(sprintf("%4d-01-01", Year)))
##' object <- TVGEVBayes(data = df, date = "Date", response = "TXMax",
##'                      design = polynomX(date = Date, degree = 1),
##'                      loc = ~ t1, scale = ~ 1, shape = ~ 1,
##'                     seed = 1234)
##' MCMC <- apply(object$MCMC, 3, rbind)
##' logLik <- -apply(MCMC, MARGIN = 1, FUN = object$TVGEV$negLogLikFun,
##'                  object = object$TVGEV)
##' ## envProf(param= MCMC[ , 1], logLik)
##' envProf(param= MCMC[ , 1], logLik, level = c(0.70, 0.95))
##' 
envProf <- function(param, logLik, level = 0.95, plot = 0) {
    
    level <- sort(level)
    nLevel <- length(level)
    n <- length(logLik)
    iMax <- which.max(logLik)
    ellMax <- logLik[iMax]
    paramMax <- param[iMax]
    delta <- qchisq(level, df = 1) / 2.0
    delta <- c(0, delta)
    dDelta <- delta[nLevel + 1] - delta[nLevel]
    delta <- c(delta, delta[nLevel + 1] + dDelta / 2.0)

    if (diff(range(param)) < 1e-10) {
        res <- matrix(param[1],
                      nrow = 2, ncol = nLevel,
                      dimnames = list(c("L", "U"), NSGEV:::formatPerc(level)))
        attr(res, "points") <- list(x = param[1], y = logLik[1])
        attr(res, "logLik") <- ellMax - delta[2:(nLevel + 1)]
        return(res) 
    }
    
    ## now delta has length nLevel + 2
    ind <- (1:n)[logLik > ellMax - delta[length(delta)]]
    indCH <- chull(param[ind], logLik[ind])
    indCH <- ind[indCH]
    nCH <- length(indCH)
    
    ## =========================================================================
    ## Since the convex hull points are provided in clockwise order,
    ## 'sdx' can take values (ignoring runs) "1, -1, 1" or "-1, 1, -1"
    ## exceptionally "1, -1" (if the numbering begins at leftmost
    ## point) and "-1, 1" (if the numbering begins at rightmost
    ## point).
    ## =========================================================================
    x <- param[indCH]
    y <- logLik[indCH]
    nx <- length(x)
    dx <- c(x[1] - x[nx], diff(x))
    sdx <- sign(dx)
    sdxNext <- c(sdx[-1], 0)
    indNS <- (1:nx)[sdx == 1 | (sdx == -1 & sdxNext == 1)]
    x <- x[indNS]
    y <- y[indNS]
    o <- order(x)
    x <- x[o]
    y <- y[o]
    if (plot > 1) {
        cat("sdx = \n")
        names(sdx) <- seq_along(sdx)
        print(rbind(sdx, sdxNext))
        print(indNS)
    }
    
    ## s <- smooth.spline(x = x, y = y)
    
    res <- matrix(NA_real_,
                  nrow = 2, ncol = nLevel,
                  dimnames = list(c("L", "U"), NSGEV:::formatPerc(level)))
    indSide <- x <= paramMax
    
    res["L", ] <- approx(x = y[indSide],
                         y = x[indSide],
                         xout = ellMax - delta[2:(nLevel + 1)])$y
    indSide <- x >= paramMax
    res["U", ] <- approx(x = y[indSide],
                         y = x[indSide],
                         xout = ellMax - delta[2:(nLevel + 1)])$y
    if (plot) {
        ind2 <- indCH[logLik[indCH] > ellMax - (delta[nLevel + 1] + dDelta )]
        param2 <- param[ind2]
        logLik2 <- logLik[ind2]
        plot(param, logLik, pch = 16,
             col = adjustcolor("SpringGreen3", alpha.f = 0.6))
        if (plot > 1) {
            points(param[indCH], logLik[indCH], pch = 16, col = "orangered")
            text(x = param[indCH], y = logLik[indCH], labels = 1:nCH, pos = 3)
            points(param2, logLik2, pch = 21, cex = 2)
            ## lines(s, col = "red", lwd = 2)
        }
        abline(h = ellMax - delta[1:(nLevel + 1)])
        abline(v = res["L", ])
        abline(v = res["U", ])
    }
    attr(res, "points") <- list(x = x, y = y)
    attr(res, "logLik") <- ellMax - delta[2:(nLevel + 1)]
    res 
}
##' 
##' @title Profile-Likelihood Confidence Interval
##'
##' @description Compute a profile likelihood confidence interval for
##'     a parameter of a \code{TVGEVBayes} object. This can be an
##'     ordinary parameter, or a function of the parameter vector such
##'     as a quantile function or a return level function.
##'
##' @note Although \code{object} represents a fitted \emph{Bayesian}
##'     model, the result is indeed a \emph{confidence} interval, not
##'     a credible interval. This interval relates to the non-Bayesian
##'     fitted model stored as \code{object$TVGEV}. The results should
##'     not depend on the prior used for the Bayes inference, provided
##'     that this prior is not too informative. Indded, the
##'     requirement is that the MCMC iterates provide a good coverage
##'     of the hig-likelihood region in the parameter space. The MCMC
##'     sampling is here used a trick to avoid complex computations
##'     such as the evaluation of the profile likelihood and more.
##' 
##' @param object A \code{TVGEVBayes} object.
##' 
##' @param fun A function of the parameter vector \code{psi} defining
##'     the quantity or "parameter" for which the profile likelihood
##'     interval will be computed. The default choice corresponds to
##'     the inference on each component of the parameter vector
##'     \eqn{\boldsymbol{\psi}}{psi}. It is a good practice to check
##'     the function e.g. by applying it on the vector \code{psi}
##'     obtained as \code{coef(object)} See \bold{Examples}.
##' 
##' @param level The confidence level.
##'
##' @param trace Integer level of verbosity.
##'
##' @param fixedParam A named numeric vector that can be used to fix
##'     the value of one or several parameteers. See \bold{Examples}.
##' 
##' @param ... Further arguments to be passed to \code{fun}.
##' 
##' @return An object with class \code{"profLik.TVGEVBayes"}
##'     inheriting from \code{"matrix"}. This is mainly a numeric
##'     matrix with \code{length(level)} columns containing the Lower
##'     and Upper confidence bounds in its two rows \code{"L"} and
##'     \code{"U"}. Some more data are attached to the object in the
##'     aim of building tools by using the \code{autoplot} method on
##'     the returned object.
##'
##' @seealso \code{\link{autoplot.profLik.TVGEVBayes}}.
##' 
##' @importFrom NSGEV profLik parNames
##' @importFrom data.table rbindlist
##' 
##' @method profLik TVGEVBayes
##' @export
##' 
##' @examples
##' df <- within(TXMax_Dijon, Date <- as.Date(sprintf("%4d-01-01", Year)))
##' fit_b <- TVGEVBayes(data = df, date = "Date", response = "TXMax",
##'                    design = polynomX(date = Date, degree = 1),
##'                    loc = ~ t1, scale = ~ 1, shape = ~ 1,
##'                    seed = 1234)
##' parNames(fit_b$TVGEV)
##' confint(fit_b$TVGEV, method = "prof", trace = 0)
##' 
##' ## It is a good practice to give the result as a named vector,
##' ## be it of length one.
##' q100 <- function(psi) {
##'    theta <- psi2theta(model = fit_b$TVGEV, psi = psi, date = "2030-01-01")
##'    c("R100" = qGEV(0.99, loc = theta[1], scale = theta[2], shape = theta[3]))
##' }
##' q100(coef(fit_b$TVGEV))
##' ## Should be fixed: one should be allowed to use `profLik` here
##' (prof_q100 <- beverstan:::profLik.TVGEVBayes(fit_b, fun = q100))
##' 
##' qfun <- function(psi) {
##'    theta <- psi2theta(model = fit_b$TVGEV, psi = psi, date = "2030-01-01")
##'    res <- qGEV(c(0.90, 0.99), loc = theta[1], scale = theta[2], shape = theta[3])
##'    names(res) <- c("RL10", "RL100")
##'    res 
##' }
##' qfun(coef(fit_b$TVGEV))
##' (prof_RLS <- beverstan:::profLik.TVGEVBayes(fit_b, fun = qfun))
##'
##' ## Check that `predict.TVGEV` compute the right profile likelihood
##' ## interval
##' (pred <- predict(fit_b$TVGEV, new = "2030-01-01", period = c(10, 100),
##'                  conf = "prof", trace = 0))
##'
##' ## What if a Gumbel distribution was used...
##' (prof_Gumb1 <-
##'     beverstan:::profLik.TVGEVBayes(fit_b, fixedParam = c("xi_0" = 0.0)))
##' (prof_Gumb2 <-
##'     beverstan:::profLik.TVGEVBayes(fit_b, fun = qfun,
##'                                    fixedParam = c("xi_0" = 0.0)))
##' 
profLik.TVGEVBayes <- function(object,
                               fun = identity,
                               level = 0.95,
                               trace = 1,
                               fixedParam = numeric(0),
                               ...) {
    
    if (!is.null(object$call$timeMAXdata)) {
        stop("This method can not be used when ",
             "`timeMAXdata` has been used in the fit")
    }

    level <- sort(level)
    nLevel <- length(level)
    fLevel <- NSGEV:::formatPerc(level)
    
    fun <- match.fun(fun)
    MCMC <- apply(object$MCMC, 3, rbind)

    if (length(fixedParam)) {
        nms <- names(fixedParam)
        if (is.null(nms) || !all(nms %in% parNames(object$TVGEV))) {
            stop("when 'fixParam' is given, it must be a named vector ",
                 "and the names must all be in `parNames(object$TVGEV)`")
        }
        for (nm in nms) {
            MCMC[ , nm] <- fixedParam[nm]
        }
    }
    
    ## ========================================================================
    ##  Compute the log-likelihood. This could be more efficiently
    ##  done at the Stan level since the log-likelihood is evaluated
    ##  for each MCMC iterate.
    ##  =======================================================================
    
    ell <- -apply(MCMC, MARGIN = 1, FUN = object$TVGEV$negLogLikFun,
                  object = object$TVGEV)

    ell <- drop(ell)
    ellMax <- max(ell)
    ellLim <- ellMax - c(0, qchisq(level, df = 1) / 2.0, Inf)
    intLev <- cut(ell, ellLim) 
    
    fTest <- fun(MCMC[1, ], ...)
    nParam <- length(fTest)
    nmsParam <- names(fTest)
    if (trace) {
        cat(sprintf("o Number of parameters to infer on: %d\n", nParam))
    }

    if (identical(fun, identity)) {
        Param <- MCMC
    } else {
        Param <- t(apply(MCMC, 1, fun, ...))
        if (nParam == 1) Param <- matrix(Param)
    }
    
    if (is.null(colnames(Param))) {
        colnames(Param) <- nmsParam <- paste0("param_", 1:nParam)
    } 
    
    CI <- array(NA_real_,
                dim = c(nParam, 2, nLevel),
                dimnames = list(nmsParam, c("L", "U"), fLevel))
    
    samples <- prof <- limits <- list()

    for (iParam in 1:nParam) {
        ci <- envProf(param = Param[ , iParam],
                      logLik = ell,
                      level = level, plot = 0)
        prof[[iParam]] <- data.frame(Name = nmsParam[iParam],
                                     x = attr(ci, "points")$x,
                                     y = attr(ci, "points")$y)
        samples[[iParam]] <- data.frame(Name = nmsParam[iParam],
                                        Param = Param[ , iParam],
                                        LogLik = ell,
                                        IntLev = intLev)
        limits[[iParam]]  <- data.frame(Name = nmsParam[iParam],
                                        L = ci["L", ],
                                        U = ci["U", ],
                                        LogLik = attr(ci, "logLik"),
                                        Level = fLevel) 
                                        
        CI[iParam, , ] <- ci
    }
    
    attr(CI, "prof") <- as.data.frame(data.table::rbindlist(prof))
    attr(CI, "samples") <- as.data.frame(data.table::rbindlist(samples))
    attr(CI, "limits") <- as.data.frame(data.table::rbindlist(limits))
    class(CI) <- c("profLik.TVGEVBayes", "array")
    CI
    
}

##' @title Print Method for \code{profLik.TVGEVBayes} Objects
##' 
##' @param x The \code{profLik.TVGEVBayes} object, created with the
##'     \code{profLik} method of the class \code{"TVGEVBayes"}.
##'
##' @param ... Not used.
##'
##' @return Nothing. Print the object.
##'
##' @seealso \code{\link{profLik.TVGEVBayes}},
##'     \code{\link{autoplot.profLik.TVGEVBayes}}.
##'
##' @export
##' @method print profLik.TVGEVBayes
##' 
print.profLik.TVGEVBayes <- function(x, ...) {
    o <- x
    for (anm in c("prof", "samples", "limits")) {
        attr(o, anm) <- NULL
    }
    class(o) <- "array"
    print.default(o)
}

##' @title Autoplot Method for \code{TVGEVBayes} Objects
##' 
##' @description This method returns a trellis graphics with one facet
##'     per "parameter". Each facet shows the profile log-likelihood
##'     obtained as the upper envelope of a scatterplot.  For a given
##'     confidence level, one horizontal line is added, corresponding
##'     to the level \eqn{\ell_{\text{max}} - \delta}{ell_max - delta}
##'     of the log-likelihood, where \eqn{\delta :=
##'     q_{\chi^2(1)}(\gamma) / 2 }{delta := qchisq(\gamma) /2.0} and
##'     \eqn{\gamma} is the confidence level. Several confidence
##'     levels can be used by passing a vector as \code{level}.
##' 
##' @details The \code{profLik.TVGEVBayes} object given in
##'     \code{object} can correspond to a function with a \emph{vector
##'     value} then defining several "parameters". The plot is then in
##'     "trellis graphics" style, with one facet per parameter.
##'
##' @param object A \code{TVGEVBayes} object.
##' 
##' @param ... Not used.
##'
##' @return An object inheriting from the \code{"ggplot"} class.
##'
##' @seealso \code{\link{TVGEVBayes}}.
##' 
##' @method autoplot profLik.TVGEVBayes
##' @export
##'
##' @examples
##' example(profLik.TVGEVBayes)
##' autoplot(prof_q100)
##' autoplot(prof_RLS)
##' autoplot(prof_Gumb1)
##' autoplot(prof_Gumb2)
autoplot.profLik.TVGEVBayes <- function(object, ...) {

    Param <- LogLik <- x <- y <- Name <- L <- U <- Level <- NULL
    
    g <- ggplot(data = attr(object, "samples")) +
        geom_point(mapping = aes(x = Param,
                                 y = LogLik,
                                 group = Param),
                   colour = "SpringGreen3", alpha = 0.2) +
        geom_line(data = attr(object, "prof"),
                  mapping = aes(x = x,
                                y = y,
                                group = Name)) +
        geom_vline(data = attr(object, "limits"),
                   mapping = aes(xintercept = L,
                                 group = Name,
                                 colour = Level,
                                 linetype = Level)) +
        geom_vline(data = attr(object, "limits"),
                   mapping = aes(xintercept = U,
                                 group = Name,
                                 colour = Level,
                                 linetype = Level)) +
        geom_hline(data = attr(object, "limits"),
                   mapping = aes(yintercept = LogLik,
                                 group = Name,
                                 colour = Level,
                                 linetype = Level)) +
        scale_colour_manual(values = c("purple", "orangered")) + 
        facet_grid( . ~ Name, scales = "free_x") +
        theme_gray()
    
    g
}
