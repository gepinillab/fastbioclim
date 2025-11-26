# bio33_fast: Moisture of the Least Moist Period

Calculates moisture mean of the least moist period.

## Usage

``` r
bio33_fast(speriod, speriod_min_idx, cell)
```

## Arguments

- speriod:

  Matrix of moisture period means (output from \`var_periods\`).

- speriod_min_idx:

  Vector indicating the index (1-based) of the least moist period.

- cell:

  Vector of original cell IDs.

## Value

Matrix with "bio33", "cell".
