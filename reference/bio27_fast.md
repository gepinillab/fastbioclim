# bio27_fast: Radiation of Coldest Period

Calculates the mean solar radiation of the specific rolling period
identified as the coldest (lowest temperature).

## Usage

``` r
bio27_fast(speriod, tperiod_min_idx, cell)
```

## Arguments

- speriod:

  A numeric \*\*matrix\*\* of solar radiation values (means) for each
  rolling period. \*\*Rows\*\* represent spatial units (cells).
  \*\*Columns\*\* represent the rolling periods.

- tperiod_min_idx:

  An integer \*\*vector\*\* indicating the column index (1-based) of the
  coldest period for each row. Its length must be exactly equal to the
  number of rows in \`speriod\`.

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`speriod\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio27" (mean solar radiation of
coldest period) and "cell".
