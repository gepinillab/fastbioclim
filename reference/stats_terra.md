# In-Memory Custom Variable Summarization

Internal function to calculate custom summary statistics using \`terra\`
functions. Designed for datasets that fit into RAM.

## Usage

``` r
stats_terra(
  variable,
  stats = c("mean", "max", "min", "cv_cli", "max_period", "min_period"),
  period_length = 3,
  period_stats = "mean",
  circular = TRUE,
  inter_variable = NULL,
  inter_stats = c("max_inter", "min_inter"),
  max_unit = NULL,
  min_unit = NULL,
  max_period = NULL,
  min_period = NULL,
  max_interactive = NULL,
  min_interactive = NULL,
  output_prefix = "var",
  suffix_inter_max = "inter_high",
  suffix_inter_min = "inter_low",
  gdal_opt = c("COMPRESS=DEFLATE", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS"),
  overwrite = FALSE,
  output_dir = tempdir(),
  verbose = TRUE,
  ...
)
```

## Arguments

- variable:

  A SpatRaster for the primary variable.

- stats:

  Character vector of stats to compute.

- period_length:

  Integer, length of a calculation period.

- period_stats:

  Character, either \`"mean"\` or \`"sum"\`; statistic used to summarize
  each period.

- circular:

  Logical, whether to wrap periods.

- inter_variable:

  (Optional) A \`terra::SpatRaster\` for an interactive variable.

- inter_stats:

  Character vector of interactive statistics to compute
  (\`"max_inter"\`, \`"min_inter"\`).

- max_unit, min_unit, max_period, min_period, max_interactive,
  min_interactive:

  (Optional) \`terra::SpatRaster\` with a static index.

- output_prefix:

  Character, prefix for output layer names.

- suffix_inter_max, suffix_inter_min:

  Character, suffixes for the \`"max_inter"\` and \`"min_inter"\`
  statistic names.

- gdal_opt:

  Character vector of GDAL options for writing.

- overwrite:

  Logical, whether to overwrite existing files.

- output_dir:

  Character, path to save final rasters.

- verbose:

  Logical, If \`TRUE\`, prints messages.

- ...:

  Other arguments including inter_variable, period_length, circular,
  static index SpatRasters, etc.

## Value

A \`terra::SpatRaster\` object containing the calculated summary layers.

## See also

The user-facing wrapper function \`derived_statistics()\`.
