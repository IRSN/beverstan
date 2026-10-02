# Compute Return Levels with Credible Limits

Compute return levels with credible limits. The result can be used to
produce the classical *return level plot* by using.
[`autoplot.TVGEVBayes`](https://irsn.github.io/beverstan/reference/autoplot.TVGEVBayes.md).

## Usage

``` r
# S3 method for class 'TVGEVBayes'
RL(
  object,
  newTimeRange = NULL,
  period = NULL,
  level = 0.7,
  credintType = c("HPD", "eqtail"),
  smooth = missing(period),
  ...
)
```

## Arguments

- object:

  A `TVGEVBayes` object.

- newTimeRange:

  A time range to be passed to `timeRange`.

- period:

  A vector of *return periods*.

- level:

  Credible level to be used for the intervals.

- credintType:

  The type of credible interval. `"HPD"` corresponds to the *Highest
  Posterior Density* interval and `"eqtail"` corresponds to the
  "equal-tail" choice, where both tails are given the same probability
  \\(1 - level) / 2\\.

- smooth:

  Logical. If `TRUE`, the lower and upper credible bounds will be
  smoothed by using
  [`smooth.spline`](https://rdrr.io/r/stats/smooth.spline.html).

- ...:

  Not used yet.

## Value

An object with class `"RL.GEVBayes"` inheriting from `"data.frame"`.

## Details

For a given period \\T\>1\\ e.g. \\T = 100\\, the return level
\\\rho(T)\\ is defined as the quantile with exceedance probability \\1 /
T\\. In the time-varying framework, the quantile is for a marginal GEV
distribution hence relates to a specific block in time which has to be
given by using `newTimeRange`. Hence the return level \\\rho\\T;
\boldsymbol{\theta}(t)\\\\ where \\\boldsymbol{\theta}(t)\\ is the
vector of the thee GEV parameters for the block \\t\\. This is
deterministic function of the GEV parameter and it can or should come
along with a credible interval.

## Note

Since the class `"TVGEVBayes"` is devoted to time-varying models it
would be tedious (and potentially misleading) to derive plotting
positions for the observations used in the fit. So we do not provide on
this plot empirical points as usually shown in the non time-varying
(stationary) framework.

## See also

[`credInt`](https://irsn.github.io/bever/reference/credInt.html).
