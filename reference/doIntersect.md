# Find Possible Intersections among Periods of Time

Find possible intersections between given periods of time or time
ranges.

## Usage

``` r
doIntersect(Periods, details = FALSE)
```

## Arguments

- Periods:

  Either a character vector with each of its elements
  [`timeRange`](https://irsn.github.io/beverstan/reference/timeRange.md),
  or a `"Date"` matrix with one row per period and wto colums
  corresponding to the beginning and the end.

- details:

  Logical. If `TRUE` the results are detailed and tell for all couples
  of periods if there is intersection or not.

## Value

A logical with length `1` or a data frame with one row per couple of
periods.

## Examples

``` r
per1 <- timeRange("1990_2000")
per2 <- timeRange("1970_1992")
per3 <- timeRange("1930_1960")
Periods <- rbind(per1, per2, per3)
doIntersect(Periods)
#> [1] TRUE
doIntersect(Periods, details = TRUE)
#>   Period1 Period2   Ans
#> 1    per1    per2  TRUE
#> 2    per1    per3 FALSE
#> 3    per2    per3 FALSE
doIntersect(c("1990_2000", "1970_1992", "1930_1960"), details = TRUE)
#>     Period1   Period2   Ans
#> 1 1990_2000 1970_1992  TRUE
#> 2 1990_2000 1930_1960 FALSE
#> 3 1970_1992 1930_1960 FALSE
```
