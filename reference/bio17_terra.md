# bio17_terra: Precipitation of Driest Period

Calculates the total precipitation of the specific rolling period
identified as the driest.

## Usage

``` r
bio17_terra(wet, driest_period)
```

## Arguments

- wet:

  A \`SpatRaster\` object where each layer is the precipitation sum for
  a rolling period.

- driest_period:

  A single-layer \`SpatRaster\` where cell values are integers
  indicating the layer index (1-based) of the driest period.

## Value

A single-layer \`SpatRaster\` with the precipitation of the driest
period, named "bio17".
