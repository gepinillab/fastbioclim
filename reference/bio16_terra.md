# bio16_terra: Precipitation of Wettest Period

Calculates the total precipitation of the specific rolling period
identified as the wettest.

## Usage

``` r
bio16_terra(wet, wettest_period)
```

## Arguments

- wet:

  A \`SpatRaster\` object where each layer is the precipitation sum for
  a rolling period.

- wettest_period:

  A single-layer \`SpatRaster\` where cell values are integers
  indicating the layer index (1-based) of the wettest period.

## Value

A single-layer \`SpatRaster\` with the precipitation of the wettest
period, named "bio16".
