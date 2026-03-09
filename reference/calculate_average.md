# Calculate Averages for SpatRasters

Calculates temporal averages for a multi-layer SpatRaster. This function
serves as a smart wrapper, automatically selecting between an in-memory
(\`terra\`) or out-of-core (\`tiled\`) workflow based on data size.

## Usage

``` r
calculate_average(
  x,
  index,
  output_names = NULL,
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

  A \`terra::SpatRaster\` object with multiple layers representing a
  time series.

- index:

  A numeric or integer vector defining the grouping for aggregation. Its
  length must equal the number of layers in \`x\`. For example, to
  average 360 monthly layers into 12 monthly means, \`index\` would be
  \`rep(1:12, 30)\`.

- output_names:

  A character vector of names for the output layers. Its length must
  equal the number of unique groups in \`index\`. If \`NULL\`, names
  like "avg_unit_1" are generated.

- output_dir:

  The directory where the final averaged raster layers will be saved as
  GeoTIFF files.

- user_region:

  (Optional) An \`sf\` or \`terra::SpatVector\` object. If provided, the
  input raster \`x\` is clipped and masked to this region before
  processing. The output raster's extent is the same of the
  \`user_region\`.

- method:

  The processing method: "auto", "tiled", or "terra".

- tile_degrees:

  (Tiled method only) The approximate size of processing tiles.

- gdal_opt:

  (Optional) GDAL creation options for the output GeoTIFFs.

- overwrite:

  Logical. If \`FALSE\` (default), stops if output files exist.

- verbose:

  Logical, If \`TRUE\`, prints messages.

## Value

A \`terra::SpatRaster\` object pointing to the newly created files, with
the following characteristics:

- \*\*Number of layers:\*\* The number of layers will be equal to the
  number of unique values in the \`index\` argument.

- \*\*Layer names:\*\* Layer names are determined by the
  \`output_names\` argument. If \`NULL\`, they will be generated
  automatically (e.g., 'avg_unit_01', 'avg_unit_02', etc.).

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
# The data has 60 layers (5 years of monthly data), so we create an
# index to group layers by month (1 to 12).
monthly_index <- rep(1:12, times = 5)
# Define a temporary directory for the output files
output_dir <- file.path(tempdir(), "monthly_prcp_avg")
# Run the calculate_average function
monthly_avg <- calculate_average(
  x = prcp_ts,
  index = monthly_index,
  output_names = "prcp_avg",
  output_dir = output_dir,
  overwrite = TRUE,
  verbose = FALSE
)
# Print the resulting SpatRaster summary
print(monthly_avg)
#> class       : SpatRaster 
#> size        : 49, 57, 12  (nrow, ncol, nlyr)
#> resolution  : 0.04166673, 0.04166674  (x, y)
#> extent      : 26.95833, 29.33334, -30.66667, -28.625  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=longlat +ellps=WGS84 +no_defs 
#> sources     : prcp_avg_01.tif  
#>               prcp_avg_02.tif  
#>               prcp_avg_03.tif  
#>               ... and 9 more sources
#> names       : prcp_avg_01, prcp_avg_02, prcp_avg_03, prcp_avg_04, prcp_avg_05, prcp_avg_06, ... 
#> min values  :       74.18,        69.3,       54.36,       36.42,        7.36,        3.38, ... 
#> max values  :      183.60,       172.6,      148.84,      100.18,       32.78,       13.70, ... 
# Clean up the created files
unlink(output_dir, recursive = TRUE)
# }
```
