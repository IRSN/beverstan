# Define a Time Range

Define a time range from a character string in suitable format. The
format is \\start\\`_`\\end\\ where \\start\\ and \\end\\ define the
beginning and the end of the time range.

## Usage

``` r
timeRange(str)
```

## Arguments

- str:

  A character string in suitable format, see **Examples**.

## Value

An object with class `"Date"` and length 2 containing the beginning and
the end. The time-range or "period" is formed by all dates which are
`>=` the beginning and are `<` the end.

## Examples

``` r
timeRange("2020_2222")
#> [1] "2020-01-01" "2223-01-01"
## error: and underscore '_' is needed! 
try(timeRange("2020-2222"))
#> Error in timeRange("2020-2222") : 
#>   'period' must be "beginning_end" where 'beg' and 'end' are years or dates
timeRange("2020-01-01_2222")
#> [1] "2020-01-01" "2223-01-01"
timeRange("2020-01-01_2222-09-21")
#> [1] "2020-01-01" "2222-09-22"
```
