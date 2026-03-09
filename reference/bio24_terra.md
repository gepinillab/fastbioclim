# bio24_terra: Radiation of Wettest Period

Calculates the mean solar radiation of the specific rolling period
identified as the wettest.

## Usage

``` r
bio24_terra(prad, wettest_period)
```

## Arguments

- prad:

  A \`SpatRaster\` object where each layer is the mean solar radiation
  for a rolling period.

- wettest_period:

  A single-layer \`SpatRaster\` where cell values are integers
  indicating the layer index (1-based) of the wettest period.

## Value

A single-layer \`SpatRaster\` with the solar radiation of the wettest
period, named "bio24".
