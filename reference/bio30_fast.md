# bio30_fast: Lowest Moisture Unit

Identifies the lowest moisture of the temporal unit with the lowest
value. If \`index_vector\` is \`NULL\`, it calculates the row-wise
minimum. If \`index_vector\` is provided, it extracts the value from the
specific column index for each row.

## Usage

``` r
bio30_fast(mois, cell, index_vector = NULL)
```

## Arguments

- mois:

  A numeric \*\*matrix\*\* of moisture values. \*\*Rows\*\* represent
  spatial units (cells) and \*\*columns\*\* represent temporal units
  (e.g., months).

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`mois\`.

- index_vector:

  (Optional) An integer \*\*vector\*\* of column indices (1-based). If
  provided, its length must be exactly equal to the number of rows in
  \`mois\`. Values must be between 1 and \`ncol(mois)\`. This is
  typically used to extract the moisture of a specific unit identified
  by another metric.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio30" (the minimum moisture) and
"cell".
