# bio10_terra: Mean Temperature of Warmest Period

Calculates the mean temperature of the specific rolling period
identified as the warmest.

## Usage

``` r
bio10_terra(tmp, warmest_period)
```

## Arguments

- tmp:

  A \`SpatRaster\` object where each layer represents the mean
  temperature for a rolling period.

- warmest_period:

  A single-layer \`SpatRaster\` where cell values are integers
  indicating the layer index (1-based) of the warmest period.

## Value

A single-layer \`SpatRaster\` containing the mean temperature of the
warmest period, named "bio10".
