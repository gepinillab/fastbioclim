# bio32_fast: Moisture of the Most Moist Period

Calculates moisture mean of the most moist period.

## Usage

``` r
bio32_fast(speriod, speriod_max_idx, cell)
```

## Arguments

- speriod:

  Matrix of moisture period means (output from \`var_periods\`).

- speriod_max_idx:

  Vector indicating the index (1-based) of the most moist period.

- cell:

  Vector of original cell IDs.

## Value

Matrix with "bio32", "cell".
