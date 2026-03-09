# bio26_fast: Radiation of Warmest Period

Calculates the mean solar radiation of the specific rolling period
identified as the warmest (highest temperature).

## Usage

``` r
bio26_fast(speriod, tperiod_max_idx, cell)
```

## Arguments

- speriod:

  A numeric \*\*matrix\*\* of solar radiation values (means) for each
  rolling period. \*\*Rows\*\* represent spatial units (cells).
  \*\*Columns\*\* represent the rolling periods.

- tperiod_max_idx:

  An integer \*\*vector\*\* indicating the column index (1-based) of the
  warmest period for each row. Its length must be exactly equal to the
  number of rows in \`speriod\`.

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`speriod\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio26" (mean solar radiation of
warmest period) and "cell".
