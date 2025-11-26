# bio35_fast: Moisture of Coldest Period

Calculates moisture mean of the period with the lowest temperature mean.

## Usage

``` r
bio35_fast(speriod, tperiod_min_idx, cell)
```

## Arguments

- speriod:

  Matrix of moisture period means (output from \`var_periods\`).

- tperiod_min_idx:

  Vector indicating the index (1-based) of the coldest period.

- cell:

  Vector of original cell IDs.

## Value

Matrix with "bio35", "cell".
