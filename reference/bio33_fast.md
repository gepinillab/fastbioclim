# bio33_fast: Mean Moisture of Least Moist Period

Calculates the mean moisture of the specific rolling period identified
as the least moist (lowest moisture).

## Usage

``` r
bio33_fast(speriod, speriod_min_idx, cell)
```

## Arguments

- speriod:

  A numeric \*\*matrix\*\* of moisture values (means) for each rolling
  period. \*\*Rows\*\* represent spatial units (cells). \*\*Columns\*\*
  represent the rolling periods.

- speriod_min_idx:

  An integer \*\*vector\*\* indicating the column index (1-based) of the
  least moist period for each row. Its length must be exactly equal to
  the number of rows in \`speriod\`.

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`speriod\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio33" (mean moisture of least moist
period) and "cell".
