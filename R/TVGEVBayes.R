## ***************************************************************************
##' Bayesian inference for some Time-Varying models with GEV margins.
##'
##' The Bayesian inference relies on Markov Chain Monte Carlo (MCMC)
##' sampling performed thanks to a \bold{Stan} code. The Stan code
##' needs to be compiled, which requires computing time.
##'
##' @title Bayesian Inference for some Time-Varying Models with
##' Generalized Extreme Value Margins
##'
##' @param data A data frame containing at least the two variables:
##' response and date, with arbitrary names.
##'
##' @param timeMAXdata When given, this must be an object with class
##' \code{"timeMAXdata"} specifying the largest observations for the
##' block maxima on one or several periods, usually historical
##' periods. See \code{\link{timeMAXdata}}.
##'
##' @param date,response The names of the data and the response
##' variables as in \code{\link[NSGEV]{TVGEV}}.
##'
##' @param design A \emph{call} that defines a design matrix with some
##' functions of the \code{Date} variable as its columns.
##'
##' @param loc,scale,shape  Formulas for the GEV parameters.
##'
##' @param prior Not allowed yet.
##'
##' @param blockDuration Not used yet.
##'
##' @param sampleMiss Logical. If \code{TRUE} the missing values for the
##' response are sampled. They appear in the \code{sampleFit} element of
##' the returned object, and ar e in relation with the \code{ind_miss}
##' incex vector which give their position in the \code{data} element.
##'
##' @param trace Integer level of verbosity.
##'
##' @param ... Further arguments to be passed to the
##' \code{\link[rstan]{sampling}} of the \bold{rstan} package.
##'
##' @return An object with (S3) class \code{"TVGEVBayes"}. This is
##' mainly a list with the following items
##'
##' \item{stanFit }{The object with class \code{"stanfit"} returned
##'     by \code{\link[rstan]{sampling}}. \bold{Caution}: for now
##'     the parameter names do not match the parameter names
##'     of the \code{TVGEVBayes} object.
##' }
##' \item{MCMC }{The array of MCMC iterates with the warm-up
##'     iterations discarded. The dimensions are: \emph{MCMC iterate},
##'     \emph{chain} and \emph{parameter}. The parameters have names
##'     using prefixes \code{"mu"}, \code{"sigma"} and \code{"xi"}.
##'     The remaining part of the names is given by the
##'     \code{\link[stats]{terms}} hence conforms to what would be obtained
##'     using \code{lm} with the same formula. However an exception is for
##'     the intercepts. For instance the name \code{"mu_(Intercept)"} will be
##'     replaced by \code{"mu_0"}.
##' }
##'
##' @importFrom rstan stan
##' @importFrom stats sd quantile rexp
##' @export
##'
##' @examples
##' df <- within(TXMax_Dijon, Date <- as.Date(sprintf("%4d-01-01", Year)))
##'
##' ## We now that a high temperature occured 1922-05-24, see
##' ## Infoclimat where the observations "Donnees Offcielle Meteo-France"
##' ## https://www.infoclimat.fr/climatologie/normales-records/1981-2010/dijon-longvic/valeurs/07280.html
##'
##' tMD <- timeMAXdata("1922_1922" = list("1922-05-24" = c(34.4, Inf)))
##'                    ## "1923_1923" = list("1923-09-01" = c(34.0, Inf)))
##' fit <- TVGEVBayes(data = df,
##'                   timeMAXdata = tMD,
##'                   date = "Date", response = "TXMax",
##'                   design = breaksX(date = Date, breaks = "1970-01-01", degree = 1),
##'                   loc = ~ t1 + t1_1970, scale = ~ 1, shape = ~ 1,
##'                   seed = 1234)
##'
##' ## plot the regression mean with its 95% confidence intervals, as
##' ## coputed using the \code{fitted method}.
##' autoplot(fit) + ggtitle("Kink design \"breaksX\": posterior mean")
##' ##' Return levels
##' autoplot(RL(fit)) + ggtitle("Kink design \"breaksX\": return levels")
##' autoplot(RL(fit, level = 0.95, smooth = TRUE)) +
##'     ggtitle("Kink design \"breaksX\": return levels")
##'
##' ## predictive distribution: by default for the next block
##' autoplot(predict(fit)) +
##'     ggtitle("Kink design \"breaksX\": prediction for next block")
##'
##' ## change the period.
##' autoplot(predict(fit, newTimeRange = "2020_2049")) +
##'     ggtitle("Kink design \"breaksX\": prediction for \"2020_2049\"")
##' \dontrun{
##'    library(shinystan)
##'    my_sso <- launch_shinystan(fit$stanFit)
##' }
##' fit2 <- TVGEVBayes(data = df,
##'                   timeMAXdata = tMD,
##'                   date = "Date", response = "TXMax",
##'                   design = natSplineX(date = Date, knots = "1970-01-01",
##'                                       boundaryKnots = c("1920-01-01", "2017-01-01")),
##'                   loc = ~ ns1 + ns2 + ns3 - 1, scale = ~ 1, shape = ~ 1,
##'                   seed = 2345)
##' autoplot(fit2) + ggtitle("Natural spline design: posterior mean")
##'
TVGEVBayes <- function(data,
                       timeMAXdata = NULL,
                       date, response,
                       design = NULL,
                       loc = ~ 1, scale = ~ 1, shape = ~ 1,
                       prior = NULL,
                       blockDuration = "year",
                       sampleMiss = FALSE,
                       trace = 0,
                       ...) {

    mc <- match.call()
    eDots <- list(...)

    ## Rename the response to avoid using the symbol
    colnames(data)[colnames(data) == response] <-  ".y0"

    ## manage the timeMAXdata information
    tMD <- timeMAXdata
    if (!is.null(tMD)) {
        bB <-  byBlock(tMD)
        colnames(bB) <- paste(".", colnames(bB), sep = "")

        df <- merge(x = data, y = bB, by.x = date, by.y = ".date",
                    all.x = TRUE, all.y = TRUE)
        if (any(!is.na(df$.y0) & !is.na(df$.y) & df$.y0 != df$.y) ||
            any(!is.na(df$.y0) & !is.na(df$.yL) & df$.y0 < df$.yL) ||
            any(!is.na(df$.y0) & !is.na(df$.yU) & df$.y0 > df$.yU)) {
            stop("conflict between 'timeMAXdata' and 'data'")
        }
    } else {
        df <-  data
        df$.y <- df$.yL <- df$.yU <- df$.y0
    }

    ## ========================================================================
    ## Now flag the observations in 3 categories: "observed",
    ## "missing" and "censored". For censored observations, we will
    ## use a code with values 1, 2, 3 corresponding to the finiteness
    ## of the bounds: "L", "LU", "U".
    ## ========================================================================

    code <- rep(0L, nrow(df))
    ind_all <-  1:nrow(df)
    ind_na <- ind_all[is.na(df$.y0)]

    if (!is.null(tMD)) {

        ## It can happen that some missing observation are
        ## "pseudo-missing" because they are given in MAXdata. So we first
        ## care about this.
        if (any(ind_na)) {
            ind0 <- !is.na(df$.y[ind_na])
            if (any(ind0)) {
                df$.y0[ind_na[ind0]] <- df$.y[ind_na[ind0]]
                ind_na[ind0] <- FALSE
            }
        }

        if (any(ind_na)) {
            ## censored observations
            ind1 <- is.finite(df$.yL[ind_na]) & !is.finite(df$.yU[ind_na])
            code[ind_na[ind1]] <- 1L
            ind2 <- is.finite(df$.yL[ind_na]) & is.finite(df$.yU[ind_na])
            code[ind_na[ind2]] <- 2L
            ind3 <- !is.finite(df$.yL[ind_na]) & is.finite(df$.yU[ind_na])
            code[ind_na[ind3]] <- 3L
            ind4 <- !is.finite(df$.yL[ind_na]) & !is.finite(df$.yU[ind_na])
            code[ind_na[ind4]] <- 4L
        }

        ind_obs <- ind_all[!is.na(df$.y0)] ## mind that '.y0' may have been changed
        ind_miss <- ind_all[code == 4L]
        ind_cens <- ind_all[code %in% c(1L, 2L, 3L)]
        code_cens <- code[ind_cens]

    } else {
        ind_obs <- setdiff(ind_all, ind_na)
        ind_miss <- ind_na
        ind_cens <- integer(0)
        code_cens <- numeric(0)
    }

    ## ===========================================================================
    ## We fit the TVGEV model, which will also provide us with the
    ## model matrices.  Not the rather subtle way to pass the design
    ## which must here be a call, not a character.
    ## ===========================================================================

    fit <-  do.call(NSGEV::TVGEV,
                    args = list(data = df, response = ".y0",
                        date = date,
                        design = str2lang(deparse(substitute(design),
                            width.cutoff = 500L)),
                        loc = loc, scale = scale, shape = shape))
    fit$call <- "<generated call>"

    ## ===========================================================================
    ## Prepare the data.  Mind that when a vector such as 'psi_sigma' or
    ## 'mean_psi_sigma0' turns to be of length one Stan throws an error
    ## because the object is no longer considered as a vector in the
    ## R-Stan interface. `as.array`. Here is an example
    ##     x <- 1.0; xa <- as.array(x)
    ##     dim(x)
    ##     dim(xa)
    ##     identical(x, xa)
    ## ============================================================================

    if (!sampleMiss) {
        ind_miss_bak <- ind_miss
        ind_miss <- integer(0)
    }

    data <- list(n_obs = length(ind_obs),
                 n_miss = length(ind_miss),
                 n_cens = length(ind_cens),
                 ## observations, bounds and code
                 y_obs = as.array(df$.y0[ind_obs]),
                 yL_cens = as.array(df$.yL[ind_cens]),
                 yU_cens = as.array(df$.yU[ind_cens]),
                 code_cens = as.array(code[ind_cens]),
                 ## length of the parameter blocks
                 p_mu = fit$pp["loc"],
                 p_sigma = fit$pp["scale"],
                 p_xi = fit$pp[["shape"]],
                 ## Flags indicating if the parameter is varying
                 cst_mu = fit$isCst["loc"],
                 cst_sigma = fit$isCst["scale"],
                 cst_xi = fit$isCst["shape"],
                 ## design matrices for observerations
                 X_mu_obs = fit$X[["loc"]][ind_obs, , drop = FALSE],
                 X_sigma_obs = fit$X[["scale"]][ind_obs, , drop = FALSE],
                 X_xi_obs = fit$X[["shape"]][ind_obs, , drop = FALSE],
                 ## design matrices for missing values
                 X_mu_miss = as.array(fit$X[["loc"]][ind_miss, , drop = FALSE]),
                 X_sigma_miss = as.array(fit$X[["scale"]][ind_miss, , drop = FALSE]),
                 X_xi_miss = as.array(fit$X[["shape"]][ind_miss, , drop = FALSE]),
                 ## design matrices for censored values
                 X_mu_cens = as.array(fit$X[["loc"]][ind_cens, , drop = FALSE]),
                 X_sigma_cens = as.array(fit$X[["scale"]][ind_cens, , drop = FALSE]),
                 X_xi_cens = as.array(fit$X[["shape"]][ind_cens, , drop = FALSE]),
                 ## prior mean
                 mean_psi_mu = as.array(rep(0.0, fit$pp["loc"])),
                 mean_psi_sigma = as.array(rep(1.0, fit$pp["scale"])),
                 mean_psi_xi = as.array(rep(0.0, fit$pp["shape"])),
                 ## prior covariance
                 cov_psi_mu = diag(1e7, nrow = fit$pp["loc"], ncol = fit$pp["loc"]),
                 cov_psi_sigma = diag(1e7, nrow = fit$pp["scale"], ncol = fit$pp["scale"]),
                 cov_psi_xi = diag(10.0, nrow = fit$pp["shape"], ncol = fit$pp["shape"]))

    yBar <- mean(df$.y0, na.rm = TRUE)
    yBar <- 35

    if ("chains" %in% names(eDots)) {
        nChains <- eDots$chains
    } else {
        nChains <- 4L
    }
    set.seed(123)
    if (!("inits" %in% names(eDots))) {
        inits <- list()
        for (iChain in 1:nChains) {
            ## XXX TO BE IMPROVED
            eps <- 1 + rexp(1L, rate = 20.0)
            inits[[iChain]] <-
                list(psi_mu = as.array(fit$psi[fit$ind[["loc"]]]  * eps),
                     psi_sigma = as.array(fit$psi[fit$ind[["scale"]]] * eps),
                     psi_xi = as.array(fit$psi[fit$ind[["shape"]]] * eps),
                     y_miss = as.array(rep(yBar * eps, data$n_miss) ),
                     y_cens = as.array(rep(yBar * eps, data$n_cens)))
        }
    }

    stanFit <- rstan::sampling(
        object = stanmodels$TVGEVCensor,
        data = data,              ## named list of data
        chains = nChains,         ## number of Markov chains
        ## warmup = 1000,            ## number of 'warmup' iterations per chain
        ## iter = 3000,              ## total number of iterations per chain
        ## cores = 2,                ## number of cores (could use one per chain)
        init = inits,             ## initial values
        refresh = 1,              ## progress shown?
        control = list(adapt_delta = 0.99),
        ...
        )

    A <- as.array(stanFit)
    d <-  dim(A)
    MCMC <- A[ , , 1:fit$p, drop = FALSE]
    dimnames(MCMC)[[3]] <- fit$parNames
    pn <-  dimnames(A)[[3L]]
    dim(A) <- c(d[1] * d[2], d[3])
    dimnames(A) <- list(NULL, c(fit$parName, "lp__"))
    iMAP <- which.max(A[ , "lp__"])
    MAP <- A[iMAP, 1:fit$p]
    postMean <- apply(A[ , 1:fit$p, drop = FALSE], 2L, mean)
    postSd <- apply(A[ , 1:fit$p, drop = FALSE], 2L, sd)
    postL <- apply(A[ , 1:fit$p, drop = FALSE], 2L, quantile, prob = 0.025)
    postU <- apply(A[ , 1:fit$p, drop = FALSE], 2L, quantile, prob = 0.975)

    L <- list(call = mc,
              blockDuration = blockDuration,
              data = df,
              ind_obs = ind_obs,
              ind_miss = ind_miss_bak,
              ind_cens = ind_cens,
              TVGEV = fit,
              MCMC = MCMC,
              MAP = MAP,
              postMean = postMean,
              postSd = postSd,
              postL = postL,
              postU = postU,
              stanData = data,
              stanFit = stanFit)

    class(L) <- "TVGEVBayes"

    L

}

