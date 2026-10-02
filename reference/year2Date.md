# Transform a Year into a Date

Transform a year given as a string into a date corresponding either the
beginning of the year ot the begining of the following year.

## Usage

``` r
year2Date(str, type = c("start", "end"))
```

## Arguments

- str:

  A character string.

- type:

  The type of date to return. This has only an effect when only the year
  is given, as oppposed to a full date. Then, when `type` is `"start"`,
  the character `str` is expected to describe the beginning of a time
  range and the returned date will be the beginning of the year. When
  instead `type` is `"end"`, the character `str` is expected to
  describle the end of a time range and the returned date will be the
  beginning of the next year. See **Examples**.

## Value

An object of class `"Date"` representing either the beginning or the end
of a time range.

## Caution

This function is not expected to work when `str` has length `> 1`.

## Examples

``` r
year2Date("2021")
#> [1] "2021-01-01"
year2Date("2021", type = "end")
#> [1] "2022-01-01"
```
