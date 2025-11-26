# Tiled, Out-of-Core Custom Variable Summarization

Internal function to calculate custom summary statistics for very large
datasets by processing them in tiles.

## Usage

``` r
stats_fast(
  variable_path,
  n_units,
  stats = c("mean", "max", "min"),
  period_length = 3,
  period_stats = "mean",
  circular = TRUE,
  inter_variable_path = NULL,
  inter_stats = NULL,
  max_unit_path = NULL,
  min_unit_path = NULL,
  max_period_path = NULL,
  min_period_path = NULL,
  max_interactive_path = NULL,
  min_interactive_path = NULL,
  prefix_variable = "var",
  suffix_inter_max = "inter_high",
  suffix_inter_min = "inter_low",
  user_region = NULL,
  tile_degrees = 5,
  output_dir = tempdir(),
  write_raw_vars = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- variable_path:

  Path to primary variable rasters.

- n_units:

  Integer, number of layers per variable.

- stats:

  Character vector of stats to compute.

- prefix_variable:

  Character, prefix for output files.

- ...:

  Other arguments including inter_variable_path, period_length,
  circular, static index paths, etc.

## Value

Character string: Path to the temporary directory containing
intermediate \`.qs2\` files.

## See also

The user-facing wrapper function \`derive_statistics()\`.
