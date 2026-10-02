# Autoplot Method for `timeMAXdata` Objects

Autoplot method for `timeMAXdata` objects representing censored
observations of a same variable in time.

## Usage

``` r
# S3 method for class 'timeMAXdata'
autoplot(object, ...)
```

## Arguments

- object:

  An object with class `"timeMAXdata"` representing censored
  observations of a same variable.

- ...:

  Not used yet.

## Value

A graphical object inheriting from `"ggplot"`.

## Examples

``` r
## block maxima style. There should be only one observation by block (year)
M <- timeMAXdata("1961_1980" = list("1961" = 7, "1973" = 40, "1979" = c(3, 30)),
                 "1930_1960" = list("1940" = c(10, Inf), "1943" = 22),
                 "1900_1929" = list("1910" = c(14, Inf), "1929" = c(20, 30))) 
g <- autoplot(M)
```
