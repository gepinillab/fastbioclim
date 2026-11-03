# Tiled, Out-of-Core Time Series Aggregation (Internal)

Tiled, Out-of-Core Time Series Aggregation (Internal)

## Usage

``` r
aggregate_fast(
  paths,
  index,
  output_names,
  user_region,
  tile_degrees,
  output_dir,
  verbose,
  aggregation_type = "mean"
)
```

## Arguments

- paths:

  Character vector of file paths of the input rasters (e.g., from
  \`terra::sources()\`).

- index:

  Vector giving the group each layer belongs to (see \`terra::tapp()\`).

- output_names:

  Character vector with the names of the aggregated output layers.

- user_region:

  (Optional) An \`sf\` or \`terra::SpatVector\` object. If \`NULL\`, the
  full extent of the input rasters is used.

- tile_degrees:

  Numeric, size of processing tiles.

- output_dir:

  Character, path for temporary files.

- verbose:

  Logical, If \`TRUE\`, prints messages.

- aggregation_type:

  Character, either \`"mean"\` or \`"sum"\`.
