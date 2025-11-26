# bio26_fast: Solar Radiation of Warmest Period

Calculates solar radiation mean of the period with the highest
temperature mean.

## Usage

``` r
bio26_fast(speriod, tperiod_max_idx, cell)
```

## Arguments

- speriod:

  Matrix of solar radiation period means (output from \`var_periods\`).

- tperiod_max_idx:

  Vector indicating the index (1-based) of the warmest period.

- cell:

  Vector of original cell IDs.

## Value

Matrix with "bio26", "cell".
