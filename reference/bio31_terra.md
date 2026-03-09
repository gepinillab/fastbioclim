# bio31_terra: Moisture Seasonality (Std Dev \* 100)

Calculates the standard deviation of moisture across all layers,
multiplied by 100.

## Usage

``` r
bio31_terra(mois)
```

## Arguments

- mois:

  A \`SpatRaster\` object where each layer represents moisture for a
  temporal unit.

## Value

A single-layer \`SpatRaster\` with the moisture seasonality, named
"bio31".
