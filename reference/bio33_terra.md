# bio33_terra: Mean Moisture of Least Moist Period

Calculates the mean moisture of the specific rolling period identified
as the least moist.

## Usage

``` r
bio33_terra(pmois, low_mois_period)
```

## Arguments

- pmois:

  A \`SpatRaster\` object where each layer is the mean moisture for a
  rolling period.

- low_mois_period:

  A single-layer \`SpatRaster\` where cell values are integers
  indicating the layer index (1-based) of the least moist period.

## Value

A single-layer \`SpatRaster\` with the moisture of the least moist
period, named "bio33".
