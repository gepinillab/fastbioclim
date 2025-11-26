# Calculate Sliding Periods for Temporal Analysis

Defines sliding periods of a specified length for temporal calculations,
handling wrap-around cases for circular data.

## Usage

``` r
compute_periods(n_units, period_length, circular = TRUE)
```

## Arguments

- n_units:

  Integer. The total number of temporal units (layers).

- period_length:

  Integer. The number of consecutive units in each period.

- circular:

  Logical. If TRUE, allows periods to wrap around.

## Value

List where each element contains \`period_length\` consecutive unit
indices.
