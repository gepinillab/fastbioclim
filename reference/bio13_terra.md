# bio13_terra: Precipitation of Wettest Unit

Identifies the precipitation of the wettest temporal unit (layer).

## Usage

``` r
bio13_terra(prcp, wettest_unit = NULL)
```

## Arguments

- prcp:

  A \`SpatRaster\` object where each layer represents precipitation for
  a temporal unit.

- wettest_unit:

  (Optional) A single-layer \`SpatRaster\` where cell values are
  integers indicating a static layer index (1-based) from which to
  extract the value. If \`NULL\` (the default), the overall maximum
  across all layers is calculated.

## Value

A single-layer \`SpatRaster\` with the precipitation of the wettest
unit, named "bio13".
