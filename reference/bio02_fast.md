# bio02_fast: Mean Diurnal Range

Calculates the mean of (tmax - tmin) across all temporal units.

## Usage

``` r
bio02_fast(tmin, tmax, cell)
```

## Arguments

- tmin:

  Matrix of minimum temperatures for each unit.

- tmax:

  Matrix of maximum temperatures for each unit.

- cell:

  Vector of original cell IDs.

## Value

Matrix with "bio02", "cell".
