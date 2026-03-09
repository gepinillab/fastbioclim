# bio25_fast: Radiation of Driest Period

Calculates the mean solar radiation of the specific rolling period
identified as the driest (lowest precipitation).

## Usage

``` r
bio25_fast(speriod, pperiod_min_idx, cell)
```

## Arguments

- speriod:

  A numeric \*\*matrix\*\* of solar radiation values (means) for each
  rolling period. \*\*Rows\*\* represent spatial units (cells).
  \*\*Columns\*\* represent the rolling periods.

- pperiod_min_idx:

  An integer \*\*vector\*\* indicating the column index (1-based) of the
  driest period for each row. Its length must be exactly equal to the
  number of rows in \`speriod\`.

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`speriod\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio25" (mean solar radiation of
driest period) and "cell".
