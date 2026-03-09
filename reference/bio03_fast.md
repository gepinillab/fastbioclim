# bio03_fast: Isothermality

Calculates Isothermality representing the ratio of mean diurnal range to
temperature range: (bio02 / bio07) \* 100.

## Usage

``` r
bio03_fast(bio02V, bio07V, cell)
```

## Arguments

- bio02V:

  A numeric \*\*vector\*\* or \*\*single-column matrix\*\* of Bio02
  values (Mean Diurnal Range). Length (or number of rows) must match
  \`bio07V\` and \`cell\`.

- bio07V:

  A numeric \*\*vector\*\* or \*\*single-column matrix\*\* of Bio07
  values (Temperature Range). Length (or number of rows) must match
  \`bio02V\` and \`cell\`.

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  length/rows of \`bio02V\` and \`bio07V\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio03" and "cell".
