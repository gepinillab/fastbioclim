# bio03_terra: Isothermality

Calculates Isothermality, defined as (bio02 / bio07) \* 100.

## Usage

``` r
bio03_terra(bio02, bio07)
```

## Arguments

- bio02:

  A single-layer \`SpatRaster\` of Mean Diurnal Range (Bio02).

- bio07:

  A single-layer \`SpatRaster\` of Temperature Range (Bio07).

## Value

A single-layer \`SpatRaster\` with the calculated isothermality, named
"bio03".
