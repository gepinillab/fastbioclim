# bio27_fast: Solar Radiation of Coldest Period

Calculates solar radiation mean of the period with the lowest
temperature mean.

## Usage

``` r
bio27_fast(speriod, tperiod_min_idx, cell)
```

## Arguments

- speriod:

  Matrix of solar radiation period means (output from \`var_periods\`).

- tperiod_min_idx:

  Vector indicating the index (1-based) of the coldest period.

- cell:

  Vector of original cell IDs.

## Value

Matrix with "bio27", "cell".
