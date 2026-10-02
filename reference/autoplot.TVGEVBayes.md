# Autoplot Method for `TVGEVBayes` Objects

Autoplot a `TVGEVBayes` object, trying to show the dependence of the GEV
marginal distribution on the time/date variable.

## Usage

``` r
# S3 method for class 'TVGEVBayes'
autoplot(object, ...)
```

## Arguments

- object:

  A `TVGEVBayes` object.

- ...:

  Further arguments passed to
  [`fitted.TVGEVBayes`](https://irsn.github.io/beverstan/reference/fitted.TVGEVBayes.md),
  such as `wich` or `level`.

## Value

A graphical object inheriting from `"ggplot"`.

## See also

[`TVGEVBayes`](https://irsn.github.io/beverstan/reference/TVGEVBayes.md).
