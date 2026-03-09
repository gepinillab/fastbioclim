# bio04_terra: Temperature Seasonality (Std Dev \* 100)

Calculates the standard deviation of average temperatures across all
layers, multiplied by 100.

## Usage

``` r
bio04_terra(tavg)
```

## Arguments

- tavg:

  A \`SpatRaster\` object where each layer represents average
  temperature for a temporal unit.

## Value

A single-layer \`SpatRaster\` with the temperature seasonality, named
"bio04".
