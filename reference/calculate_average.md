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

  (Optional) An \`sf\` or \`terra::SpatVector\` object for clipping.

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

A \`terra::SpatRaster\` object pointing to the newly created files.
