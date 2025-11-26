# bio09_fast: Mean Temperature of Driest Period

Calculates mean temperature of the period with the lowest precipitation
sum.

## Usage

``` r
bio09_fast(tperiod, pperiod_min_idx, period_length, cell)
```

## Arguments

- tperiod:

  Matrix of temperature period sums.

- pperiod_min_idx:

  Vector indicating the index (1-based) of the driest period.

- period_length:

  Integer. Number of units per period.

- cell:

  Vector of original cell IDs.

## Value

Matrix with "bio09", "cell".
