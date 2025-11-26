# Tiled, Out-of-Core Bioclimatic Variable Calculation

Internal function to calculate bioclimatic variables for very large
datasets by processing them in tiles. It reads data from file paths
using \`exactextractr\` and performs calculations with \`Rfast\`.

## Usage

``` r
bioclim_fast(
  bios,
  n_units,
  tmin_path = NULL,
  tmax_path = NULL,
  prcp_path = NULL,
  tavg_path = NULL,
  srad_path = NULL,
  mois_path = NULL,
  period_length = 3,
  circular = TRUE,
  user_region = NULL,
  tile_degrees = 5,
  output_dir = tempdir(),
  verbose = TRUE,
  ...
)
```

## Arguments

- bios:

  Numeric vector of variables to compute.

- n_units:

  Integer, number of layers per variable.

- period_length:

  Integer, length of a calculation period.

- circular:

  Logical, whether to wrap periods.

- user_region:

  An \`sf\` or \`SpatVector\` object for the area of interest.

- tile_degrees:

  Numeric, size of processing tiles.

- output_dir:

  Character, path for temporary files.

- verbose:

  Logical, If \`TRUE\`, prints messages.

- ...:

  File paths for climate variables (e.g., \`tmin_path\`) and static
  indices (e.g., \`warmest_period_path\`).

## Value

Character string: Path to the temporary directory containing
intermediate \`.qs2\` files, to be used by an assembly function.

## See also

The user-facing wrapper function \`derive_bioclim()\`.
