# bio02_fast: Mean Diurnal Range

Calculates the mean of the diurnal temperature ranges (tmax - tmin)
across all temporal units.

## Usage

``` r
bio02_fast(tmin, tmax, cell)
```

## Arguments

- tmin:

  A numeric \*\*matrix\*\* of minimum temperatures. \*\*Rows\*\*
  represent spatial units (cells) and \*\*columns\*\* represent temporal
  units. Must have the exact same dimensions as \`tmax\`.

- tmax:

  A numeric \*\*matrix\*\* of maximum temperatures. \*\*Rows\*\*
  represent spatial units (cells) and \*\*columns\*\* represent temporal
  units. Must have the exact same dimensions as \`tmin\`.

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`tmin\` and \`tmax\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio02" (mean diurnal range) and
"cell".
