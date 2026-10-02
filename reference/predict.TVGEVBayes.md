# Prediction from a `TVGEVBayes` Object

Perform some computations about the predictive distribution from a
`TVGEVBayes` object. For a "new" period of time given with
`newTimeRange`, we can consider the maximum of the block maxima on this
period. The distribution is not a GEV distribution in general, and due
to the uncertainty on the parameters

## Usage

``` r
# S3 method for class 'TVGEVBayes'
predict(
  object,
  newTimeRange = NULL,
  prob,
  type = "RL",
  approx = FALSE,
  trace = 0,
  ...
)
```

## Arguments

- object:

  A `TVGEVBayes` object.

- newTimeRange:

  A character describing the "new" time range (period of time) as in
  [`timeRange`](https://irsn.github.io/beverstan/reference/timeRange.md).

- prob:

  A vector of exceedance probabilities at which the predictive
  distribution is to be computed.

- type:

  The type of result wanted. For now only the predictive Return Levels
  can be computed.

- approx:

  `Logical`. If `TRUE` the predictive quantiles will be computed by
  using a linear interpolation method, while a series of
  [`uniroot`](https://rdrr.io/r/stats/uniroot.html) is used instead.

- trace:

  Integer level of verbosity.

- ...:

  Not used yet.

## Examples

``` r
 
newTimeRange <- "2020_2026"
```
