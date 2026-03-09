# bio35_terra: Mean Moisture of Coldest Period

Calculates the mean moisture of the specific rolling period identified
as the coldest.

## Usage

``` r
bio35_terra(pmois, coldest_period)
```

## Arguments

- pmois:

  A \`SpatRaster\` object where each layer is the mean moisture for a
  rolling period.

- coldest_period:

  A single-layer \`SpatRaster\` where cell values are integers
  indicating the layer index (1-based) of the coldest period.

## Value

A single-layer \`SpatRaster\` with the moisture of the coldest period,
named "bio35".
