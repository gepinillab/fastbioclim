# bio07_terra: Temperature Range (bio05 - bio06)

Calculates the difference between the Maximum Temperature of the Warmest
Unit (bio05) and the Minimum Temperature of the Coldest Unit (bio06).

## Usage

``` r
bio07_terra(bio05, bio06)
```

## Arguments

- bio05:

  A single-layer \`SpatRaster\` of Bio05 values.

- bio06:

  A single-layer \`SpatRaster\` of Bio06 values.

## Value

A single-layer \`SpatRaster\` with the temperature annual range, named
"bio07".
