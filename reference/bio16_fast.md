# bio16_fast: Precipitation of Wettest Period

Calculates precipitation sum of the period with the highest
precipitation sum.

## Usage

``` r
bio16_fast(pperiod, pperiod_max_idx, cell)
```

## Arguments

- pperiod:

  Matrix of precipitation period sums (output from \`var_periods\`).

- pperiod_max_idx:

  Vector indicating the index (1-based) of the wettest period.

- cell:

  Vector of original cell IDs.

## Value

Matrix with "bio16", "cell".
