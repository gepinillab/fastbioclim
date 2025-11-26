# bio10_fast: Mean Temperature of Warmest Period

Calculates mean temperature of the period with the highest temperature
sum.

## Usage

``` r
bio10_fast(tperiod, tperiod_max_idx, period_length, cell)
```

## Arguments

- tperiod:

  Matrix of temperature period sums.

- tperiod_max_idx:

  Vector indicating the index (1-based) of the warmest period.

- period_length:

  Integer. Number of units per period.

- cell:

  Vector of original cell IDs.

## Value

Matrix with "bio10", "cell".
