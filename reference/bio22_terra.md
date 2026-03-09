# bio22_terra: Lowest Radiation Unit

Identifies the lowest solar radiation of the unit with the lowest value.

## Usage

``` r
bio22_terra(srad, low_rad_unit = NULL)
```

## Arguments

- srad:

  A \`SpatRaster\` object where each layer represents solar radiation
  for a temporal unit.

- low_rad_unit:

  (Optional) A single-layer \`SpatRaster\` where cell values are
  integers indicating a static layer index (1-based). If \`NULL\`, the
  overall minimum across all layers is calculated.

## Value

A single-layer \`SpatRaster\` with the minimum solar radiation, named
"bio22".
