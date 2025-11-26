# bio11_fast: Mean Temperature of Coldest Period

Calculates mean temperature of the period with the lowest temperature
sum.

## Usage

``` r
bio11_fast(tperiod, tperiod_min_idx, period_length, cell)
```

## Arguments

- tperiod:

  Matrix of temperature period sums.

- tperiod_min_idx:

  Vector indicating the index (1-based) of the coldest period.

- period_length:

  Integer. Number of units per period.

- cell:

  Vector of original cell IDs.

## Value

Matrix with "bio11", "cell".
