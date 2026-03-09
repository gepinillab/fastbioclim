# bio32_fast: Mean Moisture of Most Moist Period

Calculates the mean moisture of the specific rolling period identified
as the most moist (highest moisture).

## Usage

``` r
bio32_fast(speriod, speriod_max_idx, cell)
```

## Arguments

- speriod:

  A numeric \*\*matrix\*\* of moisture values (means) for each rolling
  period. \*\*Rows\*\* represent spatial units (cells). \*\*Columns\*\*
  represent the rolling periods.

- speriod_max_idx:

  An integer \*\*vector\*\* indicating the column index (1-based) of the
  most moist period for each row. Its length must be exactly equal to
  the number of rows in \`speriod\`.

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`speriod\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio32" (mean moisture of most moist
period) and "cell".
