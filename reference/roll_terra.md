# In-Memory Rolling Time Series Averaging (Internal)

In-Memory Rolling Time Series Averaging (Internal)

## Usage

``` r
roll_terra(
  x,
  window_size,
  freq,
  step,
  fun,
  output_names_list,
  output_dir,
  gdal_opt,
  overwrite,
  verbose
)
```

## Arguments

- x:

  A \`terra::SpatRaster\` where each layer is a time unit.

- window_size:

  Integer, size of the moving window in number of cycles.

- freq:

  Integer, number of layers (units) that make up one cycle.

- step:

  Integer, number of cycles the window slides at each step.

- fun:

  Character, name of the summary function (e.g., \`"mean"\`).

- output_names_list:

  List of character vectors with the output layer names, one element per
  window.

- output_dir:

  Character, path to save final rasters.

- gdal_opt:

  Character vector of GDAL options for writing.

- overwrite:

  Logical, whether to overwrite existing files.

- verbose:

  Logical, If \`TRUE\`, prints messages.
