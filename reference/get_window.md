# Create Moving Window Summaries of a SpatRaster

This internal helper function calculates moving window summaries
(periods) across the layers of a SpatRaster. It generates a new
SpatRaster where each layer represents the sum of the layers from the
original raster that fall within a specific window. The \`circular\`
argument controls whether the windows "wrap around" from the end of the
time series to the beginning.

## Usage

``` r
get_window(x, period, circular)
```

## Arguments

- x:

  A \`terra::SpatRaster\` object where each layer represents a temporal
  unit.

- period:

  An integer specifying the size (number of layers) of the moving
  window. For monthly data (\`nlyr(x) = 12\`), a \`period\` of 3
  corresponds to a quarter.

- circular:

  A logical. If \`TRUE\`, the window wraps from the last layer to the
  first. If \`FALSE\`, windows are only created where they fit
  completely within the original sequence of layers.

## Value

A \`terra::SpatRaster\` object where each layer is the pixel-wise sum
for a corresponding period.

- If \`circular = TRUE\`, the output has the same number of layers as
  the input \`x\`.

- If \`circular = FALSE\`, the output has \`nlyr(x) - period + 1\`
  layers.
