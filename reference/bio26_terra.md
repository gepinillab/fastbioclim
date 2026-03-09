# bio26_terra: Radiation of Warmest Period

Calculates the mean solar radiation of the specific rolling period
identified as the warmest.

## Usage

``` r
bio26_terra(prad, warmest_period)
```

## Arguments

- prad:

  A \`SpatRaster\` object where each layer is the mean solar radiation
  for a rolling period.

- warmest_period:

  A single-layer \`SpatRaster\` where cell values are integers
  indicating the layer index (1-based) of the warmest period.

## Value

A single-layer \`SpatRaster\` with the solar radiation of the warmest
period, named "bio26".
