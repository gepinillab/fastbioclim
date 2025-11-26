# bio17_fast: Precipitation of Driest Period

Calculates precipitation sum of the period with the lowest precipitation
sum.

## Usage

``` r
bio17_fast(pperiod, pperiod_min_idx, cell)
```

## Arguments

- pperiod:

  Matrix of precipitation period sums.

- pperiod_min_idx:

  Vector indicating the index (1-based) of the driest period.

- cell:

  Vector of original cell IDs.

## Value

Matrix with "bio17", "cell".
