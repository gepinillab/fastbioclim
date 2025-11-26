# bio18_fast: Precipitation of Warmest Period

Calculates precipitation sum of the period with the highest temperature
sum.

## Usage

``` r
bio18_fast(pperiod, tperiod_max_idx, cell)
```

## Arguments

- pperiod:

  Matrix of precipitation period sums.

- tperiod_max_idx:

  Vector indicating the index (1-based) of the warmest period.

- cell:

  Vector of original cell IDs.

## Value

Matrix with "bio18", "cell".
