# bio12_fast: Total Precipitation

Calculates the total precipitation (sum) across all temporal units
(e.g., 12 months).

## Usage

``` r
bio12_fast(prcp, cell)
```

## Arguments

- prcp:

  A numeric \*\*matrix\*\* of precipitation values. \*\*Rows\*\*
  represent spatial units (cells) and \*\*columns\*\* represent temporal
  units (e.g., months).

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`prcp\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio12" (total precipitation) and
"cell".
