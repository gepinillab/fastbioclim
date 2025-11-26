# bio25_fast: Solar Radiation of Driest Period

Calculates solar radiation mean of the period with the highest
precipitation sum.

## Usage

``` r
bio25_fast(speriod, pperiod_min_idx, cell)
```

## Arguments

- speriod:

  Matrix of solar radiation period means (output from \`var_periods\`).

- pperiod_min_idx:

  Vector indicating the index (1-based) of the driest period.

- cell:

  Vector of original cell IDs.

## Value

Matrix with "bio25", "cell".
