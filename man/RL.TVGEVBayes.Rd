% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RLTVGEVBayes.R
\name{RL.TVGEVBayes}
\alias{RL.TVGEVBayes}
\title{Compute Return Levels with Credible Limits}
\usage{
\method{RL}{TVGEVBayes}(
  object,
  newTimeRange = NULL,
  period = NULL,
  level = 0.7,
  credintType = c("HPD", "eqtail"),
  smooth = missing(period),
  ...
)
}
\arguments{
\item{object}{A \code{TVGEVBayes} object.}

\item{newTimeRange}{A time range to be passed to \code{timeRange}.}

\item{period}{A vector of \emph{return periods}.}

\item{level}{Credible level to be used for the intervals.}

\item{credintType}{The type of credible interval. \code{"HPD"}
corresponds to the \emph{Highest Posterior Density} interval and
\code{"eqtail"} corresponds to the "equal-tail" choice, where both
tails are given the same probability \eqn{(1 - level) / 2}.}

\item{smooth}{Logical. If \code{TRUE}, the lower and upper credible
bounds will be smoothed by using \code{\link[stats]{smooth.spline}}.}

\item{...}{Not used yet.}
}
\value{
An object with class \code{"RL.GEVBayes"} inheriting from
\code{"data.frame"}.
}
\description{
Compute return levels with credible limits. The result can be used
to produce the classical \emph{return level plot} by using.
\code{\link{autoplot.TVGEVBayes}}.
}
\details{
For a given period \eqn{T>1} e.g. \eqn{T = 100}, the return level
\eqn{\rho(T)} is defined as the quantile with exceedance
probability \eqn{1 / T}. In the time-varying framework, the
quantile is for a marginal GEV distribution hence relates to a
specific block in time which has to be given by using
\code{newTimeRange}. Hence the return level \eqn{\rho\{T;
\boldsymbol{\theta}(t)\}}{\rho[T; \theta(t)]} where
\eqn{\boldsymbol{\theta}(t)}{\theta(t)} is the vector of the thee
GEV parameters for the block \eqn{t}. This is deterministic
function of the GEV parameter and it can or should come along with
a credible interval.
}
\note{
Since the class \code{"TVGEVBayes"} is devoted to time-varying models
it would be tedious (and potentially misleading) to derive
plotting positions for the observations used in the fit. So we do
not provide on this plot empirical points as usually shown in the
non time-varying (stationary) framework.
}
\seealso{
\code{\link[bever]{credInt}}.
}
