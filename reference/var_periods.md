# Calculate Temporal Period Aggregates

Computes aggregates (sums or mean) over defined temporal periods and
identifies min/max periods.

## Usage

``` r
var_periods(variable, periodos, n_units, period_length, stat)
```

## Arguments

- variable:

  Matrix of temporal values (rows=pixels/cells, cols=units).

- periodos:

  List of unit groupings (output from compute_periods).

- n_units:

  Integer. Total number of temporal units.

- period_length:

  Integer. Length of each period.

- stat:

  Character. Either \`"mean"\` or \`"sum"\`, specifying whether to
  calculate the average or total

## Value

Matrix with period sums (\`num_periods\` columns named "period_X"),
"min_idx", and "max_idx" columns.
