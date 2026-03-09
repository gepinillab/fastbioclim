# bio34_terra: Mean Moisture of Warmest Period

Calculates the mean moisture of the specific rolling period identified
as the warmest.

## Usage

``` r
bio34_terra(pmois, warmest_period)
```

## Arguments

- pmois:

  A \`SpatRaster\` object where each layer is the mean moisture for a
  rolling period.

- warmest_period:

  A single-layer \`SpatRaster\` where cell values are integers
  indicating the layer index (1-based) of the warmest period.

## Value

A single-layer \`SpatRaster\` with the moisture of the warmest period,
named "bio34".
