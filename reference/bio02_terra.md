# bio02_terra: Mean Diurnal Range

Calculates the mean of the diurnal temperature range (tmax - tmin)
across all temporal units (layers).

## Usage

``` r
bio02_terra(tmin, tmax)
```

## Arguments

- tmin:

  A \`SpatRaster\` object of minimum temperatures, where each layer is a
  temporal unit. Must have the same dimensions and number of layers as
  \`tmax\`.

- tmax:

  A \`SpatRaster\` object of maximum temperatures, where each layer is a
  temporal unit. Must have the same dimensions and number of layers as
  \`tmin\`.

## Value

A single-layer \`SpatRaster\` with the mean diurnal range, named
"bio02".
