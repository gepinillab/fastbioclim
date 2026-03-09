# bio21_terra: Highest Moisture Unit

Identifies the highest solar radiation of the unit with the highest
value.

## Usage

``` r
bio21_terra(srad, high_rad_unit = NULL)
```

## Arguments

- srad:

  A \`SpatRaster\` object where each layer represents solar radiation
  for a temporal unit.

- high_rad_unit:

  (Optional) A single-layer \`SpatRaster\` where cell values are
  integers indicating a static layer index (1-based). If \`NULL\`, the
  overall maximum across all layers is calculated.

## Value

A single-layer \`SpatRaster\` with the maximum solar radiation, named
"bio21".
