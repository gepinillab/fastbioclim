# bio30_terra: Lowest Moisture Unit

Identifies the lowest moisture of the unit with the lowest value.

## Usage

``` r
bio30_terra(mois, low_mois_unit = NULL)
```

## Arguments

- mois:

  A \`SpatRaster\` object where each layer represents moisture for a
  temporal unit.

- low_mois_unit:

  (Optional) A single-layer \`SpatRaster\` where cell values are
  integers indicating a static layer index (1-based). If \`NULL\`, the
  overall minimum across all layers is calculated.

## Value

A single-layer \`SpatRaster\` with the minimum moisture, named "bio30".
