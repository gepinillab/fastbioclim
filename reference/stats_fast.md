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
  output_prefix = "var",
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

- period_length:

  Integer, length of a calculation period.

- period_stats:

  Character, either \`"mean"\` or \`"sum"\`; statistic used to summarize
  each period.

- circular:

  Logical, whether to wrap periods.

- inter_variable_path:

  (Optional) Path to the interactive variable rasters, one path per
  layer.

- inter_stats:

  (Optional) Character vector of interactive statistics to compute
  (\`"max_inter"\`, \`"min_inter"\`).

- max_unit_path, min_unit_path, max_period_path, min_period_path,
  max_interactive_path, min_interactive_path:

  (Optional) Path to a single-layer raster with a static index.

- output_prefix:

  Character, prefix for output files.

- suffix_inter_max, suffix_inter_min:

  Character, suffixes for the \`"max_inter"\` and \`"min_inter"\`
  statistic names.

- user_region:

  (Optional) An \`sf\` or \`terra::SpatVector\` object. If \`NULL\`, the
  full extent of the input rasters is used.

- tile_degrees:

  Numeric, size of processing tiles.

- output_dir:

  Character, path for temporary files.

- write_raw_vars:

  Logical. Currently ignored: writing raw variable tiles is controlled
  by the \`BIOCLIM_DEBUG_RAW_VARS\` environment variable.

- verbose:

  Logical, If \`TRUE\`, prints messages.

- ...:

  Other arguments including inter_variable_path, period_length,
  circular, static index paths, etc.

## Value

Character string: Path to the temporary directory containing
intermediate \`.qs2\` files.

## See also

The user-facing wrapper function \`derive_statistics()\`.
