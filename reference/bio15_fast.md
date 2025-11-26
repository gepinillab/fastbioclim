# bio15_fast: Precipitation Seasonality (CV)

Calculates coefficient of variation in precipitation across units.

## Usage

``` r
bio15_fast(prcp, bio12V, n_units, cell)
```

## Arguments

- prcp:

  Matrix containing precipitation values for each unit.

- bio12V:

  Precomputed total precipitation (BIO12 value).

- n_units:

  Integer. The total number of temporal units.

- cell:

  Vector of original cell IDs.

## Value

Matrix with "bio15", "cell".
