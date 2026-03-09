# bio11_terra: Mean Temperature of Coldest Period

Calculates the mean temperature of the specific rolling period
identified as the coldest.

## Usage

``` r
bio11_terra(tmp, coldest_period)
```

## Arguments

- tmp:

  A \`SpatRaster\` object where each layer represents the mean
  temperature for a rolling period.

- coldest_period:

  A single-layer \`SpatRaster\` where cell values are integers
  indicating the layer index (1-based) of the coldest period.

## Value

A single-layer \`SpatRaster\` containing the mean temperature of the
coldest period, named "bio11".
