# bio30_fast: Lowest Moisture Unit

Identifies lowest moisture unit, potentially using a static index.

## Usage

``` r
bio30_fast(mois, cell, index_vector = NULL)
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

Matrix with "bio30", "cell".
