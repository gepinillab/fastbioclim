# \#' @title bio21_fast: Highest Radiation Unit

Identifies the highest solar radiation of the temporal unit with the
highest value. If \`index_vector\` is \`NULL\`, it calculates the
row-wise maximum. If \`index_vector\` is provided, it extracts the value
from the specific column index for each row.

## Usage

``` r
bio21_fast(srad, cell, index_vector = NULL)
```

## Arguments

- srad:

  A numeric \*\*matrix\*\* of solar radiation values. \*\*Rows\*\*
  represent spatial units (cells) and \*\*columns\*\* represent temporal
  units (e.g., months).

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`srad\`.

- index_vector:

  (Optional) An integer \*\*vector\*\* of column indices (1-based). If
  provided, its length must be exactly equal to the number of rows in
  \`srad\`. Values must be between 1 and \`ncol(srad)\`. This is
  typically used to extract the solar radiation of a specific unit
  identified by another metric.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio21" (the maximum solar radiation)
and "cell".
