#' The 'beverstan' package.
#'
#' @description Bayesian inference for some Extreme-Value models such as
#' Time-Varying GEV models for block maxima.
#'
#' @docType package
#' @name beverstan-package
#' @aliases beverstan
#' @useDynLib beverstan, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom rstantools nsamples
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#'
"_PACKAGE"
