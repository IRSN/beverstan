# Default "New" Time Range

Guess a default "new" time range for a `TVGEVBayes` object.

## Usage

``` r
newTimeRange(object, timeRange = NULL, ...)
```

## Arguments

- object:

  A `TVGEVBayes` object.

- timeRange:

  A character with length one defining the new time range. See the
  [`timeRange`](https://irsn.github.io/beverstan/reference/timeRange.md)
  function.

- ...:

  Not used yet.

## Value

Alist with the two elements

- `timeRange` :

  A character with length one describing the time range as in
  `"2020_2019"`.

- `date` :

  A `Date` vector with length 2 giving the beginning and the end of the
  time-range.

## Details

Given a `TVGEVBayes` object we can find both the block duration and the
beginning and end of the last block. The default new block is the next
block with the same duration. Note that if `object` is based on a block
duration of \\w \> 2\\ years, the underlying model does not know
anything about the maxima over blocks with smaller duration, hence we
can not predict the maximum for a block with duration \\\< w\\.
