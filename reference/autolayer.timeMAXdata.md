# Autolayer Method for `timeMAXdata` Objects

Autolayer method for `timeMAXdata` objects representing censored
observations of a same variable in time.

## Usage

``` r
# S3 method for class 'timeMAXdata'
autolayer(object, which = c("segments", "points"), ...)
```

## Arguments

- object:

  An object with class `"timeMAXdata"` representing censored
  observations of a same variables.

- which:

  Character vector giving the layer to add: `"points"` shows the
  non-censored observations as points and `"segments"` shows the
  censored observations as vertical segements.

- ...:

  Not used yet.

## Value

A graphical layer to be added to an existing `ggplot` object.
