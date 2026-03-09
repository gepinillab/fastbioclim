# bio08_fast: Mean Temperature of Wettest Period

Calculates the mean temperature of the specific rolling period
identified as the wettest (highest precipitation).

## Usage

``` r
bio08_fast(tperiod, pperiod_max_idx, period_length, cell)
```

## Arguments

- tperiod:

  A numeric \*\*matrix\*\* of temperature values (means) for each
  rolling period. \*\*Rows\*\* represent spatial units (cells).
  \*\*Columns\*\* represent the rolling periods (typically 12).

- pperiod_max_idx:

  An integer \*\*vector\*\* indicating the column index (1-based) of the
  wettest period for each row. Its length must be exactly equal to the
  number of rows in \`tperiod\`.

- period_length:

  A single \*\*integer\*\* representing the number of units per period
  (e.g., 3 for months).

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`tperiod\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio08" (mean temperature of wettest
period) and "cell".
