# Tiled, Out-of-Core Rolling Time Series Averaging (Internal)

Tiled, Out-of-Core Rolling Time Series Averaging (Internal)

## Usage

``` r
roll_fast(
  paths,
  window_size,
  freq,
  step,
  fun,
  output_names_list,
  user_region,
  tile_degrees,
  output_dir,
  verbose
)
```

## Arguments

- paths:

  Character vector of file paths of the input rasters (e.g., from
  \`terra::sources()\`).

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

- user_region:

  (Optional) An \`sf\` or \`terra::SpatVector\` object. If \`NULL\`, the
  full extent of the input rasters is used.

- tile_degrees:

  Numeric, size of processing tiles.

- output_dir:

  Character, path for temporary files.

- verbose:

  Logical, If \`TRUE\`, prints messages.
