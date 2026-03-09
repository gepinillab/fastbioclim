# bio07_fast: Temperature Range (bio05 - bio06)

Calculates the Temperature Range, defined as the difference between the
Maximum Temperature of the Warmest Unit (bio05) and the Minimum
Temperature of the Coldest Unit (bio06).

## Usage

``` r
bio07_fast(bio05V, bio06V, cell)
```

## Arguments

- bio05V:

  A numeric \*\*vector\*\* or \*\*single-column matrix\*\* of Bio05
  values. Length (or number of rows) must match \`bio06V\` and \`cell\`.

- bio06V:

  A numeric \*\*vector\*\* or \*\*single-column matrix\*\* of Bio06
  values. Length (or number of rows) must match \`bio05V\` and \`cell\`.

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  length/rows of \`bio05V\` and \`bio06V\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio07" (annual range) and "cell".
