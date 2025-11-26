# bio21_fast: Highest Solar Radiation Unit

Identifies highest solar radiation unit, potentially using a static
index.

## Usage

``` r
bio21_fast(srad, cell, index_vector = NULL)
```

## Arguments

- srad:

  Matrix of solar radiation values for each unit.

- cell:

  Vector of original cell IDs.

- index_vector:

  Optional vector of unit indices (1-based). If provided, extracts Srad
  for that unit. If NULL, finds overall max Srad.

## Value

Matrix with "bio21", "cell".
