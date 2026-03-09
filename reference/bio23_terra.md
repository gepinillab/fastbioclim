# bio23_terra: Radiation Seasonality (CV)

Calculates the Coefficient of Variation (CV) of solar radiation across
all layers.

## Usage

``` r
bio23_terra(srad)
```

## Arguments

- srad:

  A \`SpatRaster\` object where each layer represents solar radiation
  for a temporal unit.

## Value

A single-layer \`SpatRaster\` with the solar radiation seasonality,
named "bio23".
