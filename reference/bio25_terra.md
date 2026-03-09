# bio25_terra: Radiation of Driest Period

Calculates the mean solar radiation of the specific rolling period
identified as the driest.

## Usage

``` r
bio25_terra(prad, driest_period)
```

## Arguments

- prad:

  A \`SpatRaster\` object where each layer is the mean solar radiation
  for a rolling period.

- driest_period:

  A single-layer \`SpatRaster\` where cell values are integers
  indicating the layer index (1-based) of the driest period.

## Value

A single-layer \`SpatRaster\` with the solar radiation of the driest
period, named "bio25".
