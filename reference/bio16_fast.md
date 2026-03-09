# bio16_fast: Precipitation of Wettest Period

Calculates the total precipitation of the specific rolling period
identified as the wettest (highest precipitation).

## Usage

``` r
bio16_fast(pperiod, pperiod_max_idx, cell)
```

## Arguments

- pperiod:

  A numeric \*\*matrix\*\* of precipitation sums for each rolling
  period. \*\*Rows\*\* represent spatial units (cells). \*\*Columns\*\*
  represent the rolling periods.

- pperiod_max_idx:

  An integer \*\*vector\*\* indicating the column index (1-based) of the
  wettest period for each row. Its length must be exactly equal to the
  number of rows in \`pperiod\`.

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`pperiod\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio16" (precipitation of wettest
period) and "cell".
