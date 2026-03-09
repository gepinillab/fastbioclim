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
  output_prefix = "var",
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

- output_prefix:

  A character string used as the prefix for all output file names (e.g.,
  \`output_prefix = "wind"\` results in "wind_mean.tif",
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

  Integer. The number of temporal units (e.g., months) that define a
  "period". This is used for all period-based statistics, such as
  \`"max_period"\` and interactive statistics like \`"max_inter"\`. The
  same period length is applied to both the primary \`variable\` and the
  \`inter_variable\` when calculating these statistics. Default: 3.

- period_stats:

  Character. The statistic ("mean" or "sum") to summarize data over each
  period. Default: "mean".

- circular:

  Logical. If \`TRUE\` (the default), period calculations wrap around.

- user_region:

  (Optional) An \`sf\` or \`terra::SpatVector\` object. If provided, the
  input raster \`x\` is clipped and masked to this region before
  processing. The output raster's extent is the same of the
  \`user_region\`.

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
rasters, with the following characteristics:

- \*\*Number of layers:\*\* The number of layers is equal to the total
  number of statistics requested in the \`stats\` and \`inter_stats\`
  arguments.

- \*\*Layer names:\*\* Layer names are constructed by combining the
  \`output_prefix\` with the name of each statistic (e.g.,
  'prefix_mean', 'prefix_max'). For interactive statistics, the names
  are 'prefix_suffix_inter_high' and 'prefix_suffix_inter_low', where
  the suffixes are controlled by the \`suffix_inter_max\` and
  \`suffix_inter_min\` arguments.

- \*\*Extent:\*\* If a \`user_region\` is provided, the extent of the
  output raster will be clipped to match that region. Otherwise, the
  extent will be the same as the input \`variable\` raster.

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

## Examples

``` r
# \donttest{
# The example raster "prcp.tif" is included in the package's `inst/extdata` directory.
# Load example data from Lesotho (Montlhy time series from 2016-01 to 2020-12)
raster_path <- system.file("extdata", "prcp.tif", package = "fastbioclim")
# Load the SpatRaster from the file
prcp_ts <- terra::rast(raster_path)
# The data has 60 layers (5 years of monthly data), so we create an
# index to group layers by month (1 to 12).
monthly_index <- rep(1:12, times = 5)
# Define a temporary directory for the output files
output_dir <- file.path(tempdir(), "prcp_stats")
# Run the calculate_average function
monthly_avg <- calculate_average(
  x = prcp_ts,
  index = monthly_index,
  output_names = "prcp_avg",
  output_dir = output_dir,
  overwrite = TRUE,
  verbose = FALSE
)
# Once the monthly averaged is obtained, we can use it to derive statictics of precipitation
# Here we will obtain four summary statitics: mean, maximum, minimum, and standard deviation
prcp_stats <- derive_statistics(
  variable = monthly_avg,
  stats = c("mean", "max", "min", "stdev"),
  prefix_variable = "prcp",
  output_dir = output_dir,
  overwrite = TRUE
)
#> Using 'auto' method to select workflow...
#> Full rasters appear to fit in memory. Selecting 'terra' workflow.
#> SpatRasters have same extent, number of rows and columns, projection, resolution, and origin
#> Writing GeoTIFFs...
#> Processing complete. Final rasters are in: /tmp/RtmpWglfG0/prcp_stats
# Print the resulting SpatRaster summary with the four requested layers
print(prcp_stats)
#> class       : SpatRaster 
#> size        : 49, 57, 4  (nrow, ncol, nlyr)
#> resolution  : 0.04166673, 0.04166674  (x, y)
#> extent      : 26.95833, 29.33334, -30.66667, -28.625  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=longlat +ellps=WGS84 +no_defs 
#> sources     : var_mean.tif  
#>               var_max.tif  
#>               var_min.tif  
#>               var_stdev.tif  
#> names       : var_mean, var_max, var_min, var_stdev 
#> min values  : 37.43666,   75.44,    3.24,  27.61666 
#> max values  : 82.51000,  183.60,   13.70,  63.43156 
# Clean up the created files
unlink(output_dir, recursive = TRUE)
# }
```
