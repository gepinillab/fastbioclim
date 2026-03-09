# bio24_fast: Radiation of Wettest Period

Calculates the mean solar radiation of the specific rolling period
identified as the wettest (highest precipitation).

## Usage

``` r
bio24_fast(speriod, pperiod_max_idx, cell)
```

## Arguments

- speriod:

  A numeric \*\*matrix\*\* of solar radiation values (means) for each
  rolling period. \*\*Rows\*\* represent spatial units (cells).
  \*\*Columns\*\* represent the rolling periods.

- pperiod_max_idx:

  An integer \*\*vector\*\* indicating the column index (1-based) of the
  wettest period for each row. Its length must be exactly equal to the
  number of rows in \`speriod\`.

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`speriod\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio24" (mean solar radiation of
wettest period) and "cell".
