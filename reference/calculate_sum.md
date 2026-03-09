# Calculate Sums for SpatRasters

Calculates temporal sums for a multi-layer SpatRaster. This function
serves as a smart wrapper, automatically selecting between an in-memory
(\`terra\`) or out-of-core (\`tiled\`) workflow based on data size.

## Usage

``` r
calculate_sum(
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
  length must equal the number of layers in \`x\`. For example, to sum
  365 daily layers into 12 monthly totals, \`index\` would group the
  days by month.

- output_names:

  A character vector of names for the output layers. Its length must
  equal the number of unique groups in \`index\`. If \`NULL\`, names
  like "sum_unit_1" are generated.

- output_dir:

  The directory where the final summed raster layers will be saved as
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
  automatically (e.g., 'sum_unit_01', 'sum_unit_02', etc.).

- \*\*Extent:\*\* If \`user_region\` is provided, the extent of the
  output raster will be clipped to match that region. Otherwise, the
  extent will be the same as the input raster \`x\`.

## Examples

``` r
# \donttest{
# The example raster "prcp.tif" is included in the package's `inst/extdata` directory.
# Load example data from Lesotho (Monthly time series from 2016-01 to 2020-12)
raster_path <- system.file("extdata", "prcp.tif", package = "fastbioclim")
# Load the SpatRaster from the file
prcp_ts <- terra::rast(raster_path)
# To calculate the total annual precipitation, we group the 60 monthly layers by year.
annual_index <- rep(2016:2020, each = 12)
# Define a temporary directory for the output files
output_dir <- file.path(tempdir(), "annual_prcp_sum")
# Run the calculate_sum function
annual_sum <- calculate_sum(
  x = prcp_ts,
  index = annual_index,
  output_names = "prcp_sum",
  output_dir = output_dir,
  overwrite = TRUE,
  verbose = FALSE
)
# Print the resulting SpatRaster summary (should have 5 layers)
print(annual_sum)
#> class       : SpatRaster 
#> size        : 49, 57, 5  (nrow, ncol, nlyr)
#> resolution  : 0.04166673, 0.04166674  (x, y)
#> extent      : 26.95833, 29.33334, -30.66667, -28.625  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=longlat +ellps=WGS84 +no_defs 
#> sources     : prcp_sum_2016.tif  
#>               prcp_sum_2017.tif  
#>               prcp_sum_2018.tif  
#>               ... and 2 more sources
#> names       : prcp_sum_2016, prcp_sum_2017, prcp_sum_2018, prcp_sum_2019, prcp_sum_2020 
#> min values  :         430.0,         466.8,         388.8,         392.7,         551.2 
#> max values  :         954.6,        1004.5,         939.2,         986.9,        1096.0 
# Clean up the created files
unlink(output_dir, recursive = TRUE)
# }
```
