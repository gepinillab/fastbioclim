# bio15_terra: Precipitation Seasonality (CV)

Calculates the Coefficient of Variation (CV) of precipitation across all
layers.

## Usage

``` r
bio15_terra(prcp)
```

## Arguments

- prcp:

  A \`SpatRaster\` object where each layer represents precipitation for
  a temporal unit.

## Value

A single-layer \`SpatRaster\` with the precipitation seasonality, named
"bio15".

## Note

The formula adds 1 to the mean to avoid division by zero in arid areas.
