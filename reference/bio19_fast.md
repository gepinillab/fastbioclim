# bio19_fast: Precipitation of Coldest Period

Calculates precipitation sum of the period with the lowest temperature
sum.

## Usage

``` r
bio19_fast(pperiod, tperiod_min_idx, cell)
```

## Arguments

- pperiod:

  Matrix of precipitation period sums.

- tperiod_min_idx:

  Vector indicating the index (1-based) of the coldest period.

- cell:

  Vector of original cell IDs.

## Value

Matrix with "bio19", "cell".
