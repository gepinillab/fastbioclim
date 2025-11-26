# bio08_fast: Mean Temperature of Wettest Period

Calculates mean temperature of the period with the highest precipitation
sum.

## Usage

``` r
bio08_fast(tperiod, pperiod_max_idx, period_length, cell)
```

## Arguments

- tperiod:

  Matrix of temperature period sums (output from \`var_periods\`).

- pperiod_max_idx:

  Vector indicating the index (1-based) of the wettest period for each
  row.

- period_length:

  Integer. Number of units per period.

- cell:

  Vector of original cell IDs.

## Value

Matrix with "bio08", "cell".
