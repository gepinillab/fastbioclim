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

  (Optional) An \`sf\` or \`terra::SpatVector\` object. If provided, the
  input raster \`x\` is clipped and masked to this region before
  processing. The output raster's extent is the same of the
  \`user_region\`.

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

A \`terra::SpatRaster\` object pointing to the newly created files, with
the following characteristics:

- \*\*Number of layers:\*\* The number of layers is determined by the
  number of rolling windows processed (controlled by \`window_size\` and
  \`step\`) multiplied by the cycle frequency (\`freq\`).

- \*\*Layer names:\*\* Layer names are constructed based on the
  \`name_template\` argument, incorporating the window range and unit
  index (e.g., 'output_w01-20_u01').

- \*\*Extent:\*\* If \`user_region\` is provided, the extent of the
  output raster will be clipped to match that region. Otherwise, the
  extent will be the same as the input raster \`x\`.

## Examples

``` r
# \donttest{
# The example raster "prcp.tif" is included in the package's `inst/extdata` directory.
# Load example data from Lesotho (Montlhy time series from 2016-01 to 2020-12)
raster_path <- system.file("extdata", "prcp.tif", package = "fastbioclim")
# Load the SpatRaster from the file
prcp_ts <- terra::rast(raster_path)
# The data has 60 layers (5 years of monthly data).
# We want to calculate a 3-year rolling average on monthly data.
# Therefore, the window size is 3 (number of years) 
# and the lenght of the cycle is 12 (number of months).
# We also want a moving window of one month, so the number of steps is 1.
n_years <- 3
n_months <- 12
n_steps <- 1
# Define a temporary directory for the output files
output_dir <- file.path(tempdir(), "roll_prcp_avg")
# Run the calculate_average function
roll_avg <- calculate_roll(
  x = prcp_ts,
  window_size = n_years,
  freq = n_months,
  step = n_steps,
  fun = "mean",
  output_prefix = "prcp_roll_avg",
  output_dir = output_dir,
  overwrite = TRUE,
  verbose = FALSE
)
# Print the resulting SpatRaster summary of 36 layesr (n_year * n_months)
print(roll_avg)
#> class       : SpatRaster 
#> size        : 49, 57, 36  (nrow, ncol, nlyr)
#> resolution  : 0.04166673, 0.04166674  (x, y)
#> extent      : 26.95833, 29.33334, -30.66667, -28.625  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=longlat +ellps=WGS84 +no_defs 
#> sources     : prcp_roll_avg_w1-3_u01.tif  
#>               prcp_roll_avg_w1-3_u02.tif  
#>               prcp_roll_avg_w1-3_u03.tif  
#>               ... and 33 more sources
#> names       : prcp_~3_u01, prcp_~3_u02, prcp_~3_u03, prcp_~3_u04, prcp_~3_u05, prcp_~3_u06, ... 
#> min values  :        81.1,        63.0,    56.13334,    26.33333,    10.73333,         4.8, ... 
#> max values  :       194.5,       151.5,   155.56667,    71.26666,    46.73333,        20.1, ... 
# Clean up the created files
unlink(output_dir, recursive = TRUE)
# }
```
