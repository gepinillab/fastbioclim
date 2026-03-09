# bio01_terra: Mean Temperature

Calculates mean temperature across all temporal units (layers).

## Usage

``` r
bio01_terra(tavg)
```

## Arguments

- tavg:

  A \`SpatRaster\` object where each layer represents average
  temperature for a temporal unit (e.g., 12 months).

## Value

A single-layer \`SpatRaster\` with the calculated mean temperature,
named "bio01".
