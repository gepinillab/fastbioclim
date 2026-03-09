# bio08_terra: Mean Temperature of Wettest Period

Calculates the mean temperature of the specific rolling period
identified as the wettest.

## Usage

``` r
bio08_terra(tmp, wettest_period)
```

## Arguments

- tmp:

  A \`SpatRaster\` object where each layer represents the mean
  temperature for a rolling period.

- wettest_period:

  A single-layer \`SpatRaster\` where cell values are integers
  indicating the layer index (1-based) of the wettest period.

## Value

A single-layer \`SpatRaster\` containing the mean temperature of the
wettest period, named "bio08".
