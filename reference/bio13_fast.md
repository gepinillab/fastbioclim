# bio13_fast: Precipitation of Wettest Unit

Identifies the precipitation of the wettest temporal unit (e.g., month).
If \`index_vector\` is \`NULL\`, it calculates the row-wise maximum. If
\`index_vector\` is provided, it extracts the value from the specific
column index for each row.

## Usage

``` r
bio13_fast(prcp, cell, index_vector = NULL)
```

## Arguments

- prcp:

  A numeric \*\*matrix\*\* of precipitation values. \*\*Rows\*\*
  represent spatial units (cells) and \*\*columns\*\* represent temporal
  units (e.g., months).

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`prcp\`.

- index_vector:

  (Optional) An integer \*\*vector\*\* of column indices (1-based). If
  provided, its length must be exactly equal to the number of rows in
  \`prcp\`. Values must be between 1 and \`ncol(prcp)\`. This is
  typically used to extract the precipitation of the specific month
  identified as the wettest by another metric.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio13" (precipitation of wettest
unit) and "cell".
