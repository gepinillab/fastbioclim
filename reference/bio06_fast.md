# bio06_fast: Min Temperature of Coldest Unit

Identifies min temperature of the coldest unit, potentially using a
static index.

## Usage

``` r
bio06_fast(tmin, cell, index_vector = NULL)
```

## Arguments

- tmin:

  Matrix of minimum temperatures for each unit.

- cell:

  Vector of original cell IDs.

- index_vector:

  Optional vector of unit indices (1-based). If provided, extracts Tmin
  for that unit. If NULL, finds overall min Tmin.

## Value

Matrix with "bio06", "cell".
