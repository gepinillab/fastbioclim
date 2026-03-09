# bio14_terra: Precipitation of Driest Unit

Identifies the precipitation of the driest temporal unit (layer).

## Usage

``` r
bio14_terra(prcp, driest_unit = NULL)
```

## Arguments

- prcp:

  A \`SpatRaster\` object where each layer represents precipitation for
  a temporal unit.

- driest_unit:

  (Optional) A single-layer \`SpatRaster\` where cell values are
  integers indicating a static layer index (1-based) from which to
  extract the value. If \`NULL\` (the default), the overall minimum
  across all layers is calculated.

## Value

A single-layer \`SpatRaster\` with the precipitation of the driest unit,
named "bio14".
