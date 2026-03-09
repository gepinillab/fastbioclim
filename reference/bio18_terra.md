# bio18_terra: Precipitation of Warmest Period

Calculates the total precipitation of the specific rolling period
identified as the warmest.

## Usage

``` r
bio18_terra(wet, warmest_period)
```

## Arguments

- wet:

  A \`SpatRaster\` object where each layer is the precipitation sum for
  a rolling period.

- warmest_period:

  A single-layer \`SpatRaster\` where cell values are integers
  indicating the layer index (1-based) of the warmest period.

## Value

A single-layer \`SpatRaster\` with the precipitation of the warmest
period, named "bio18".
