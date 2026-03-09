# bio09_fast: Mean Temperature of Driest Period

Calculates the mean temperature of the specific rolling period
identified as the driest (lowest precipitation).

## Usage

``` r
bio09_fast(tperiod, pperiod_min_idx, period_length, cell)
```

## Arguments

- tperiod:

  A numeric \*\*matrix\*\* of temperature values (means) for each
  rolling period. \*\*Rows\*\* represent spatial units (cells).
  \*\*Columns\*\* represent the rolling periods.

- pperiod_min_idx:

  An integer \*\*vector\*\* indicating the column index (1-based) of the
  driest period for each row. Its length must be exactly equal to the
  number of rows in \`tperiod\`.

- period_length:

  A single \*\*integer\*\* representing the number of units per period.

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`tperiod\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio09" (mean temperature of driest
period) and "cell".
