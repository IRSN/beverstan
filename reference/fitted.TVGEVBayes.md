# Fitted Values for a `TVGEVBayes` Object

Fitted values for a `TVGEVBayes` object. By default, the fitted values
for a given date/block are taken as the expectations of the marginal
distributions. However, marginal quantiles can instead be used, for
instance quantiles with probability 0.99 a.k.a. the 100-year return
levels. In both cases the returned vector contains the posterior mean of
the chosen quantity and it comes along with the limits of a credible
interval.

## Usage

``` r
# S3 method for class 'TVGEVBayes'
fitted(object, which = c("expect", "quantile"), prob = 0.99, level = 0.95, ...)
```

## Arguments

- object:

  An object with class `"TVGEVBayes"`.

- which:

  A location parameter for the marginal distribution, either the
  expectation or a quantile.

- prob:

  The probability to use when `which` is taken to be `"quantile"`.

- level:

  The credible level.

- ...:

  Not used yet.

## Value

A list containing the vector `mean` and the two vectors of bounds:
`lower` and `upper`. An attribute `"which"` is used to keep the trace of
the choice made for `which`.
