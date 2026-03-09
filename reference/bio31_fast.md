# bio31_fast: Moisture Seasonality (Std Dev \* 100)

Calculates Moisture Seasonality, defined as the standard deviation of
moisture values across all temporal units (e.g., 12 months), multiplied
by 100.

## Usage

``` r
bio31_fast(mois, n_units, cell)
```

## Arguments

- mois:

  A numeric \*\*matrix\*\* containing moisture values. \*\*Rows\*\*
  represent spatial units (cells) and \*\*columns\*\* represent temporal
  units.

- n_units:

  A single \*\*integer\*\* representing the number of temporal units
  (e.g., 12).

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`mois\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio31" (Moisture Seasonality) and
"cell".
