# bio20_fast: Mean Radiation

Calculates mean solar radiation across all temporal units (usually 12
months).

## Usage

``` r
bio20_fast(srad, cell)
```

## Arguments

- srad:

  A numeric \*\*matrix\*\* where \*\*rows\*\* represent spatial units
  (cells) and \*\*columns\*\* represent temporal units (e.g., 12
  months).

- cell:

  A numeric or character \*\*vector\*\* of original cell IDs. Its length
  must be exactly equal to the number of rows in \`srad\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio20" (the calculated mean) and
"cell" (the IDs).
