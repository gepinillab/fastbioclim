# bio05_fast: Max Temperature of Warmest Unit

Identifies max temperature of the warmest unit, potentially using a
static index.

## Usage

``` r
bio05_fast(tmax, cell, index_vector = NULL)
```

## Arguments

- tmax:

  Matrix of maximum temperatures for each unit.

- cell:

  Vector of original cell IDs.

- index_vector:

  Optional vector of unit indices (1-based). If provided, extracts Tmax
  for that unit. If NULL, finds overall max Tmax.

## Value

Matrix with "bio05", "cell".
