# bio34_fast: Moisture of Warmest Period

Calculates moisture mean of the period with the highest temperature
mean.

## Usage

``` r
bio34_fast(speriod, tperiod_max_idx, cell)
```

## Arguments

- speriod:

  Matrix of moisture period means (output from \`var_periods\`).

- tperiod_max_idx:

  Vector indicating the index (1-based) of the warmest period.

- cell:

  Vector of original cell IDs.

## Value

Matrix with "bio34", "cell".
