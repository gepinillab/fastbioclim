# bio32_terra: Mean Moisture of Most Moist Period

Calculates the mean moisture of the specific rolling period identified
as the most moist.

## Usage

``` r
bio32_terra(pmois, high_mois_period)
```

## Arguments

- pmois:

  A \`SpatRaster\` object where each layer is the mean moisture for a
  rolling period.

- high_mois_period:

  A single-layer \`SpatRaster\` where cell values are integers
  indicating the layer index (1-based) of the most moist period.

## Value

A single-layer \`SpatRaster\` with the moisture of the most moist
period, named "bio32".
