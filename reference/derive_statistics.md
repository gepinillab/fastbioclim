# Derive Custom Summary Statistics from Climate Variables

Calculates a wide range of custom summary statistics for a primary
climate variable, with options for interactions with a second variable.
This function serves as a smart wrapper that automatically selects the
most efficient processing workflow (in-memory vs. tiled).

## Usage

``` r
derive_statistics(
  variable,
  stats = c("mean", "max", "min"),
  inter_variable = NULL,
  inter_stats = NULL,
  prefix_variable = "var",
  suffix_inter_max = "inter_high",
  suffix_inter_min = "inter_low",
  output_dir = tempdir(),
  period_length = 3,
  period_stats = "mean",
  circular = TRUE,
  user_region = NULL,
  method = c("auto", "tiled", "terra"),
  tile_degrees = 5,
  gdal_opt = c("COMPRESS=DEFLATE", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS"),
  overwrite = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- variable:

  A \`terra::SpatRaster\` object for the primary variable.

- stats:

  A character vector of statistics to compute for the primary variable.
  Supported: \`"mean"\`, \`"max"\`, \`"min"\`, \`"sum"\`, \`"stdev"\`,
  \`"cv_cli"\`, \`"max_period"\`, \`"min_period"\`.

- inter_variable:

  (Optional) A \`terra::SpatRaster\` for an interactive variable.

- inter_stats:

  (Optional) A character vector of interactive statistics to compute.
  Requires \`inter_variable\`. Supported: \`"max_inter"\`,
  \`"min_inter"\`.

- prefix_variable:

  A character string used as the prefix for all output file names (e.g.,
  \`prefix_variable = "wind"\` results in "wind_mean.tif",
  "wind_max.tif").

- suffix_inter_max:

  Character. Suffix for the "max_inter" statistic name. Default:
  "inter_high".

- suffix_inter_min:

  Character. Suffix for the "min_inter" statistic name. Default:
  "inter_low".

- output_dir:

  The directory where the final summary rasters will be saved.

- period_length:

  Integer. The number of temporal units defining a "period". Default: 3.

- period_stats:

  Character. The statistic ("mean" or "sum") to summarize data over each
  period. Default: "mean".

- circular:

  Logical. If \`TRUE\` (the default), period calculations wrap around.

- user_region:

  (Optional) An \`sf\` or \`terra::SpatVector\` object defining the area
  of interest.

- method:

  The processing method. See Details for more information.

- tile_degrees:

  (Tiled method only) The approximate size of processing tiles.

- gdal_opt:

  (Optional) A character vector of GDAL creation options for the output
  GeoTIFF files.

- overwrite:

  (Optional) Logical. If \`FALSE\` (the default), the function will stop
  if output files already exist.

- verbose:

  Logical, If \`TRUE\`, prints messages.

- ...:

  Additional arguments, primarily for passing static index
  \`SpatRaster\` objects. See the "Static Indices" section.

## Value

A \`terra::SpatRaster\` object pointing to the newly created summary
rasters.

## Details

This function provides a flexible alternative to \`derive_bioclim()\`
for any multi-layer climate variable (e.g., wind speed, humidity). It
unifies two processing backends, controlled by the \`method\` argument:

- \`"auto"\`: (Default) Intelligently chooses between "terra" and
  "tiled".

- \`"terra"\`: Forces the fast, in-memory workflow.

- \`"tiled"\`: Forces the memory-safe, out-of-core workflow. Requires
  that all input SpatRasters point to files on disk.

## Static Indices

For advanced control, provide pre-calculated index rasters as named
\`SpatRaster\` objects via the \`...\` argument (e.g., \`max_unit =
max_idx_rast\`). Supported indices: \`max_unit\`, \`min_unit\`,
\`max_period\`, \`min_period\`, \`max_interactive\`,
\`min_interactive\`.
