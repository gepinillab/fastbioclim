# Create an Inverse Cell ID Translation Function

Generates a function to translate cell IDs from a source raster grid to
a target raster grid, considering potential offsets and different
dimensions. Handles cases where the target grid is a subset (e.g.,
cropped/masked) of the source grid.

## Usage

``` r
define_translate(ncol_src, ncol_tgt, row_offset, col_offset)
```

## Arguments

- ncol_src:

  Integer. Number of columns in the source raster.

- ncol_tgt:

  Integer. Number of columns in the target raster.

- row_offset:

  Integer. Row offset of the target grid's top-left corner relative to
  the source grid's top-left corner (0-based).

- col_offset:

  Integer. Column offset of the target grid's top-left corner relative
  to the source grid's top-left corner (0-based).

## Value

A function that takes a vector of source cell IDs (\`cell_src\`) and
returns a vector of corresponding target cell IDs. Cells falling outside
the target grid bounds will have \`NA_integer\_\` as their target ID.
