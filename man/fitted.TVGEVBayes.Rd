% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TVGEVBayesMethods.R
\name{fitted.TVGEVBayes}
\alias{fitted.TVGEVBayes}
\title{Fitted Values for a \code{TVGEVBayes} Object}
\usage{
\method{fitted}{TVGEVBayes}(object, which = c("expect", "quantile"), prob = 0.99, level = 0.95, ...)
}
\arguments{
\item{object}{An object with class \code{"TVGEVBayes"}.}

\item{which}{A location parameter for the marginal distribution,
either the expectation or a quantile.}

\item{prob}{The probability to use when \code{which} is taken to
be \code{"quantile"}.}

\item{level}{The credible level.}

\item{...}{Not used yet.}
}
\value{
A list containing the vector \code{mean} and the two
vectors of bounds: \code{lower} and \code{upper}. An attribute
\code{"which"} is used to keep the trace of the choice made for
\code{which}.
}
\description{
Fitted values for a \code{TVGEVBayes} object. By default, the
fitted values for a given date/block are taken as the expectations
of the marginal distributions. However, marginal quantiles can
instead be used, for instance quantiles with probability 0.99
a.k.a. the 100-year return levels. In both cases the returned
vector contains the posterior mean of the chosen quantity and it
comes along with the limits of a credible interval.
}
