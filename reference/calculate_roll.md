# Calculate Rolling Temporal Averages for SpatRasters

Calculates temporal summaries for each time unit over a moving window of
cycles. This function is designed for time series where fundamental time
\*\*units\*\* (e.g., months) are grouped into repeating \*\*cycles\*\*
(e.g., years).

## Usage

``` r
calculate_roll(
  x,
  window_size,
  freq = 12,
  step = 1,
  fun = "mean",
  name_template = "{prefix}_w{start_window}-{end_window}_u{idx_unit}",
  output_prefix = "output",
  output_dir = tempdir(),
  user_region = NULL,
  method = c("auto", "tiled", "terra"),
  tile_degrees = 5,
  gdal_opt = c("COMPRESS=DEFLATE", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS"),
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- x:

  A \`terra::SpatRaster\` object where each layer represents a time
  \*\*unit\*\*.

- window_size:

  Integer. The size of the moving window, measured in the number of
  \*\*cycles\*\*. For example, if the data cycle is annual (\`freq =
  12\`), a \`window_size\` of 20 represents a 20-year window.

- freq:

  Integer. The number of time \*\*units\*\* (layers) that constitute one
  complete \*\*cycle\*\*. Common examples: 12 for monthly units in a
  yearly cycle, or 24 for hourly units in a daily cycle.

- step:

  Integer. The number of \*\*cycles\*\* to slide the window by for each
  iteration. Default is 1.

- fun:

  Character. The name of the summary function (e.g., "mean"). Default is
  "mean".

- name_template:

  A character string defining the template for output filenames, using
  \`glue\` syntax. Default:
  \`"{prefix}\_w{start_window}-{end_window}\_u{idx_unit}"\`. Available
  placeholders are:

  - \`{prefix}\`: The value from \`output_prefix\`.

  - \`{start_window}\`: The starting \*\*cycle\*\* index of the window.

  - \`{end_window}\`: The ending \*\*cycle\*\* index of the window.

  - \`{idx_unit}\`: The index of the time \*\*unit\*\* within the cycle
    (e.g., the month number).

- output_prefix:

  A character string for output filenames. Default is "output".

- output_dir:

  Directory to save the final GeoTIFF files.

- user_region:

  (Optional) An \`sf\` or \`terra::SpatVector\` for clipping.

- method:

  Processing method: "auto", "tiled", or "terra".

- tile_degrees:

  (Tiled method only) Approximate size of processing tiles.

- gdal_opt:

  (Optional) GDAL creation options for GeoTIFFs.

- overwrite:

  Logical. If \`FALSE\` (default), stops if output files exist.

- verbose:

  Logical, If \`TRUE\`, prints messages.

## Value

A \`terra::SpatRaster\` object pointing to the newly created files.
