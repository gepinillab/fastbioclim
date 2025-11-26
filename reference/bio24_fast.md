# bio24_fast: Solar Radiation of Wettest Period

Calculates solar radiation mean of the period with the highest
precipitation sum.

## Usage

``` r
bio24_fast(speriod, pperiod_max_idx, cell)
```

## Arguments

- speriod:

  Matrix of solar radiation period means (output from \`var_periods\`).

- pperiod_max_idx:

  Vector indicating the index (1-based) of the wettest period.

- cell:

  Vector of original cell IDs.

## Value

Matrix with "bio24", "cell".
