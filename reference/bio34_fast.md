# bio34_fast: Mean Moisture of Warmest Period

Calculates the mean moisture of the specific rolling period identified
as the warmest (highest temperature).

## Usage

``` r
bio34_fast(speriod, tperiod_max_idx, cell)
```

## Arguments

- speriod:

  A numeric \*\*matrix\*\* of moisture values (means) for each rolling
  period. \*\*Rows\*\* represent spatial units (cells). \*\*Columns\*\*
  represent the rolling periods.

- tperiod_max_idx:

  An integer \*\*vector\*\* indicating the column index (1-based) of the
  warmest period for each row. Its length must be exactly equal to the
  number of rows in \`speriod\`.

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`speriod\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio34" (mean moisture of warmest
period) and "cell".
