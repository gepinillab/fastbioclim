# bio17_fast: Precipitation of Driest Period

Calculates the total precipitation of the specific rolling period
identified as the driest (lowest precipitation).

## Usage

``` r
bio17_fast(pperiod, pperiod_min_idx, cell)
```

## Arguments

- pperiod:

  A numeric \*\*matrix\*\* of precipitation sums for each rolling
  period. \*\*Rows\*\* represent spatial units (cells). \*\*Columns\*\*
  represent the rolling periods.

- pperiod_min_idx:

  An integer \*\*vector\*\* indicating the column index (1-based) of the
  driest period for each row. Its length must be exactly equal to the
  number of rows in \`pperiod\`.

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`pperiod\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio17" (precipitation of driest
period) and "cell".
