# Derive Comprehensive Bioclimatic Variables

Calculates up to 35 bioclimatic variables from average monthly climate
SpatRasters (or other temporal units). This function serves as a smart
wrapper that automatically selects the most efficient processing
workflow (in-memory vs. tiled) based on data size and user-defined
region of interest.

## Usage

``` r
derive_bioclim(
  bios,
  tmin = NULL,
  tmax = NULL,
  prcp = NULL,
  tavg = NULL,
  srad = NULL,
  mois = NULL,
  output_dir = tempdir(),
  period_length = 3,
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

- bios:

  Numeric vector specifying which bioclimatic variables (1-35) to
  compute.

- tmin, tmax, prcp, tavg, srad, mois:

  (Optional) \`terra::SpatRaster\` objects containing the climate data
  for each temporal unit (e.g., 12 monthly layers). All provided rasters
  must have the same geometry and number of layers.

- output_dir:

  The directory where the final bioclimatic variable rasters will be
  saved. The directory will be created if it does not exist. The default
  is temporal directory created by \`tempdir\`.

- period_length:

  Integer. The number of temporal units (e.g., months) that define a
  "period" for calculating summary variables like BIO8 (Mean Temp of
  Wettest Quarter). Defaults to 3, representing quarters for monthly
  data.

- circular:

  Logical. If \`TRUE\` (the default), period calculations will wrap
  around the beginning and end of the time series (e.g., for monthly
  data, Dec-Jan-Feb is considered a valid period).

- user_region:

  (Optional) An \`sf\` or \`terra::SpatVector\` object defining the area
  of interest. If provided, all calculations will be clipped to this
  region.

- method:

  The processing method. See Details for more information.

- tile_degrees:

  (Tiled method only) The approximate size of processing tiles in
  degrees. Ignored if the 'terra' workflow is used.

- gdal_opt:

  (Optional) A character vector of GDAL creation options for the output
  GeoTIFF files. Controls compression, threading, etc.

- overwrite:

  (Optional) Logical. If \`FALSE\` (the default), the function will stop
  immediately if any target output files already exist.

- verbose:

  Logical, If \`TRUE\`, prints messages.

- ...:

  Additional arguments, primarily for passing static index rasters. See
  the "Static Indices" section for details.

## Value

A \`terra::SpatRaster\` object pointing to the newly created bioclimatic
variable rasters, with the following characteristics:

- \*\*Number of layers:\*\* The number of layers is equal to the number
  of unique variables requested in the \`bios\` argument.

- \*\*Layer names:\*\* Layer names are standardized (e.g., 'bio01',
  'bio12') corresponding to the requested variable numbers.

- \*\*Extent:\*\* If a \`user_region\` is provided, the extent of the
  output raster will be clipped to match that region. Otherwise, the
  extent will be the same as the input rasters.

The returned object contains the following calculated variables:

- bio01:

  Mean Temperature

- bio02:

  Mean Diurnal Range

- bio03:

  Isothermality

- bio04:

  Temperature Seasonality

- bio05:

  Max Temperature of Warmest Unit

- bio06:

  Min Temperature of Coldest Unit

- bio07:

  Temperature Range

- bio08:

  Mean Temperature of Wettest Period

- bio09:

  Mean Temperature of Driest Period

- bio10:

  Mean Temperature of Warmest Period

- bio11:

  Mean Temperature of Coldest Period

- bio12:

  Precipitation Sum

- bio13:

  Precipitation of Wettest Unit

- bio14:

  Precipitation of Driest Unit

- bio15:

  Precipitation Seasonality

- bio16:

  Precipitation of Wettest Period

- bio17:

  Precipitation of Driest Period

- bio18:

  Precipitation of Warmest Period

- bio19:

  Precipitation of Coldest Period

- bio20:

  Mean Radiation

- bio21:

  Highest Radiation Unit

- bio22:

  Lowest Radiation Unit

- bio23:

  Radiation Seasonality

- bio24:

  Radiation of Wettest Period

- bio25:

  Radiation of Driest Period

- bio26:

  Radiation of Warmest Period

- bio27:

  Radiation of Coldest Period

- bio28\*:

  Mean Moisture

- bio29\*:

  Highest Moisture Unit

- bio30\*:

  Lowest Moisture Unit

- bio31\*:

  Moisture Seasonality

- bio32\*:

  Mean Moisture of Most Moist Period

- bio33\*:

  Mean Moisture of Least Moist Period

- bio34\*:

  Mean Moisture of Warmest Period

- bio35\*:

  Mean Moisture of Coldest Period

## Details

This function unifies two processing backends. The \`method\` argument
controls which is used:

- \`"auto"\`: (Default) Intelligently chooses between "terra" and
  "tiled" based on estimated memory requirements.

- \`"terra"\`: Forces the fast, in-memory workflow. May fail on very
  large datasets.

- \`"tiled"\`: Forces the memory-safe, out-of-core workflow. Ideal for
  very large datasets. Requires that the input SpatRasters point to
  files on disk.

Period-based variables (e.g., BIO8, BIO10) are calculated using a moving
window defined by \`period_length\`.

## Note

\*The original moisture variables proposed in the ANUCLIM manual are
based on the Moisture Index (MI). However, this function allows users to
calculate moisture-based bioclimatic variables using other units of
moisture as inputs, offering greater flexibility in input data usage.

## Static Indices

For advanced use cases, such as time-series analysis or defining
specific seasons, you can provide pre-calculated index rasters to
override the dynamic calculations. These are passed as named
\`SpatRaster\` objects via the \`...\` argument (e.g., \`warmest_period
= my_warmest_idx_rast\`). The wrapper function automatically handles
passing them to the appropriate workflow.

When using the "tiled" workflow, these static index rasters \*\*must\*\*
be file-backed (i.e., not held entirely in memory). Supported static
indices include:

- \`warmest_unit\`, \`coldest_unit\`, \`wettest_unit\`, \`driest_unit\`

- \`high_rad_unit\`, \`low_rad_unit\`, \`high_mois_unit\`,
  \`low_mois_unit\`

- \`warmest_period\`, \`coldest_period\`, \`wettest_period\`,
  \`driest_period\`

- \`high_mois_period\`, \`low_mois_period\`

## References

O’Donnell, M. S., & Ignizio, D. A. (2012). Bioclimatic predictors for
supporting ecological applications in the conterminous United States.
ANUCLIM 6.1 User Guide. Centre for Resource and Environmental Studies,
The Australian National University.

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
output_dir <- file.path(tempdir(), "prcp_bios")
# Run the calculate_average function
monthly_avg <- calculate_average(
  x = prcp_ts,
  index = monthly_index,
  output_names = "prcp_avg",
  output_dir = output_dir,
  overwrite = TRUE,
  verbose = FALSE
)
# Once the monthly averaged is obtained, we can use it to obtain bioclimatic variables based
# just in precipitation (bios 12, 13, 14, 15, 16, 17).
prcp_bios <- derive_bioclim(
  bios = 12:17,
  prcp = monthly_avg,
  output_dir = output_dir,
  overwrite = TRUE,
  verbose = FALSE
)
# Print the resulting SpatRaster summary with the four requested layers
print(prcp_bios)
#> class       : SpatRaster 
#> size        : 49, 57, 6  (nrow, ncol, nlyr)
#> resolution  : 0.04166673, 0.04166674  (x, y)
#> extent      : 26.95833, 29.33334, -30.66667, -28.625  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=longlat +ellps=WGS84 +no_defs 
#> sources     : bio12.tif  
#>               bio13.tif  
#>               bio14.tif  
#>               ... and 3 more sources
#> names       :  bio12,  bio13, bio14,    bio15, bio16, bio17 
#> min values  : 449.24,  75.44,  3.24, 63.89158, 220.8, 12.82 
#> max values  : 990.12, 183.60, 13.70, 81.78966, 504.6, 52.54 
# Clean up the created files
unlink(output_dir, recursive = TRUE)
# }
```
