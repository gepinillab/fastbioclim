# bio19_fast: Precipitation of Coldest Period

Calculates the total precipitation of the specific rolling period
identified as the coldest (lowest temperature).

## Usage

``` r
bio19_fast(pperiod, tperiod_min_idx, cell)
```

## Arguments

- pperiod:

  A numeric \*\*matrix\*\* of precipitation sums for each rolling
  period. \*\*Rows\*\* represent spatial units (cells). \*\*Columns\*\*
  represent the rolling periods.

- tperiod_min_idx:

  An integer \*\*vector\*\* indicating the column index (1-based) of the
  coldest period for each row. Its length must be exactly equal to the
  number of rows in \`pperiod\`.

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`pperiod\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio19" (precipitation of coldest
period) and "cell".
