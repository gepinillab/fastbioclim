# bio14_fast: Precipitation of Driest Unit

Identifies precipitation of the driest unit, potentially using a static
index.

## Usage

``` r
bio14_fast(prcp, cell, index_vector = NULL)
```

## Arguments

- prcp:

  Matrix of precipitation values for each unit.

- cell:

  Vector of original cell IDs.

- index_vector:

  Optional vector of unit indices (1-based). If provided, extracts Prec
  for that unit. If NULL, finds overall min Prec.

## Value

Matrix with "bio14", "cell".
