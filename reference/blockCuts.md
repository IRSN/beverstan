# Define the Time Blocks of a Given Time Range

Define the cuts (dates) for the time blocks in a given time range.

## Usage

``` r
blockCuts(timeRange, blockDuration = "year")
```

## Arguments

- timeRange:

  A time range. Can either be a character string to be passed to the
  `str` argument of the
  [`timeRange`](https://irsn.github.io/beverstan/reference/timeRange.md)
  function, or a `Date` object returned by this function.

- blockDuration:

  The duration of the blocks. Should be an integer number of years as in
  `"2 years"`.

## Value

A vector with class `"Date"` giving the cuts for the blocks. By dropping
the last element we get the beginnings of the blocks and by dropping the
first element we get the ends of the blocks.

## Examples

``` r
tB <- blockCuts("1930_2020", block = "2 years")
```
