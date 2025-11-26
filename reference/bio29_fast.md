# bio29_fast: Highest Moisture Unit

Identifies highest moisture unit, potentially using a static index.

## Usage

``` r
bio29_fast(mois, cell, index_vector = NULL)
```

## Arguments

- mois:

  Matrix of moisture values for each unit.

- cell:

  Vector of original cell IDs.

- index_vector:

  Optional vector of unit indices (1-based). If provided, extracts mois
  for that unit. If NULL, finds overall max mois.

## Value

Matrix with "bio29", "cell".
