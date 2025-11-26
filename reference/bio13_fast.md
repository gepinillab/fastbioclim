# bio13_fast: Precipitation of Wettest Unit

Identifies precipitation of the wettest unit, potentially using a static
index.

## Usage

``` r
bio13_fast(prcp, cell, index_vector = NULL)
```

## Arguments

- prcp:

  Matrix of precipitation values for each unit.

- cell:

  Vector of original cell IDs.

- index_vector:

  Optional vector of unit indices (1-based). If provided, extracts Prec
  for that unit. If NULL, finds overall max Prec.

## Value

Matrix with "bio13", "cell".
