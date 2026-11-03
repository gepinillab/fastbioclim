# In-Memory Time Series Aggregation (Internal)

In-Memory Time Series Aggregation (Internal)

## Usage

``` r
aggregate_terra(
  x,
  index,
  output_names,
  output_dir,
  gdal_opt,
  overwrite,
  verbose,
  aggregation_type = "mean"
)
```

## Arguments

- x:

  A \`terra::SpatRaster\` with the time series to aggregate.

- index:

  Vector giving the group each layer belongs to (see \`terra::tapp()\`).

- output_names:

  Character vector with the names of the aggregated output layers.

- output_dir:

  Character, path to save final rasters.

- gdal_opt:

  Character vector of GDAL options for writing.

- overwrite:

  Logical, whether to overwrite existing files.

- verbose:

  Logical, If \`TRUE\`, prints messages.

- aggregation_type:

  Character, either \`"mean"\` or \`"sum"\`.
