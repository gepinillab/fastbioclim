# In-Memory Bioclimatic Variable Calculation

Internal function to calculate bioclimatic variables using \`terra\`
functions. It is designed for datasets that can fit into RAM.

## Usage

``` r
bioclim_terra(
  bios,
  tmin = NULL,
  tmax = NULL,
  tavg = NULL,
  prcp = NULL,
  srad = NULL,
  mois = NULL,
  period_length = 3,
  circular = TRUE,
  gdal_opt = c("COMPRESS=DEFLATE", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS"),
  overwrite = FALSE,
  output_dir = tempdir(),
  verbose = TRUE,
  ...
)
```

## Arguments

- bios:

  Numeric vector of variables to compute.

- period_length:

  Integer, length of a calculation period.

- circular:

  Logical, whether to wrap periods.

- gdal_opt:

  Character vector of GDAL options for writing.

- overwrite:

  Logical, whether to overwrite existing files.

- output_dir:

  Character, path to save final rasters.

- verbose:

  Logical, If \`TRUE\`, prints messages.

- ...:

  \`SpatRaster\` objects for climate variables (e.g., \`tmin\`,
  \`tmax\`) and static indices (e.g., \`warmest_period\`).

## Value

A \`terra::SpatRaster\` object pointing to the newly created files.

## See also

The user-facing wrapper function \`derive_bioclim()\`.
