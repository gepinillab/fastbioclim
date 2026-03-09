# bio09_terra: Mean Temperature of Driest Period

Calculates the mean temperature of the specific rolling period
identified as the driest.

## Usage

``` r
bio09_terra(tmp, driest_period)
```

## Arguments

- tmp:

  A \`SpatRaster\` object where each layer represents the mean
  temperature for a rolling period.

- driest_period:

  A single-layer \`SpatRaster\` where cell values are integers
  indicating the layer index (1-based) of the driest period.

## Value

A single-layer \`SpatRaster\` containing the mean temperature of the
driest period, named "bio09".
