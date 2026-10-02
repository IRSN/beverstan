# Extract Information for Block Data

Extract information for block data.

## Usage

``` r
byBlock(object, blockDuration = "year")
```

## Arguments

- object:

  An object with class `timeMAXdata`. This is mainly a list with two
  data frame elements called `MAXinfo` and `MAXdata`.

- blockDuration:

  The duration of the blocks. This is most generally the default value
  of one year, but we could aslos use `"2 years"` doe instance.

## Value

A data frame with one row by block and with the variables `y`, `yL` and
`yU` giving the value of the variable (if known), a lower bound and an
upper bound for the value of the variable.

## Examples

``` r
## block maxima style. There should be only one observation by block (year)
M <- timeMAXdata("1961_1980" = list("1961" = 7, "1973" = 40, "1979" = c(3, 30)),
                 "1930_1960" = list("1940" = c(10, Inf), "1943" = 22),
                 "1900_1929" = list("1910" = c(14, Inf), "1929" = c(20, 30))) 
byBlock <- byBlock(M)
```
