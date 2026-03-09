# bio28_fast: Mean Moisture

Calculates mean moisture across all temporal units (usually 12 months).

## Usage

``` r
bio28_fast(mois, cell)
```

## Arguments

- mois:

  A numeric \*\*matrix\*\* where \*\*rows\*\* represent spatial units
  (cells) and \*\*columns\*\* represent temporal units (e.g., 12
  months).

- cell:

  A numeric or character \*\*vector\*\* of original cell IDs. Its length
  must be exactly equal to the number of rows in \`mois\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio28" (the calculated mean) and
"cell" (the IDs).
