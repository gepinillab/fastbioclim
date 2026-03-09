# bio29_terra: Highest Moisture Unit

Identifies the highest moisture of the unit with the highest value.

## Usage

``` r
bio29_terra(mois, high_mois_unit = NULL)
```

## Arguments

- mois:

  A \`SpatRaster\` object where each layer represents moisture for a
  temporal unit.

- high_mois_unit:

  (Optional) A single-layer \`SpatRaster\` where cell values are
  integers indicating a static layer index (1-based). If \`NULL\`, the
  overall maximum across all layers is calculated.

## Value

A single-layer \`SpatRaster\` with the maximum moisture, named "bio29".
