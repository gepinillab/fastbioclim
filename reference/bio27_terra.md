# bio27_terra: Radiation of Coldest Period

Calculates the mean solar radiation of the specific rolling period
identified as the coldest.

## Usage

``` r
bio27_terra(prad, coldest_period)
```

## Arguments

- prad:

  A \`SpatRaster\` object where each layer is the mean solar radiation
  for a rolling period.

- coldest_period:

  A single-layer \`SpatRaster\` where cell values are integers
  indicating the layer index (1-based) of the coldest period.

## Value

A single-layer \`SpatRaster\` with the solar radiation of the coldest
period, named "bio27".
