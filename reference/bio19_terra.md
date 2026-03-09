# bio19_terra: Precipitation of Coldest Period

Calculates the total precipitation of the specific rolling period
identified as the coldest.

## Usage

``` r
bio19_terra(wet, coldest_period)
```

## Arguments

- wet:

  A \`SpatRaster\` object where each layer is the precipitation sum for
  a rolling period.

- coldest_period:

  A single-layer \`SpatRaster\` where cell values are integers
  indicating the layer index (1-based) of the coldest period.

## Value

A single-layer \`SpatRaster\` with the precipitation of the coldest
period, named "bio19".
