# bio15_fast: Precipitation Seasonality (CV)

Calculates the Coefficient of Variation (CV) of precipitation. The
formula used is: \`(StandardDeviation / (Mean + 1)) \* 100\`. (The "+1"
is added to the mean to avoid division by zero in completely arid
areas).

## Usage

``` r
bio15_fast(prcp, bio12V, n_units, cell)
```

## Arguments

- prcp:

  A numeric \*\*matrix\*\* of precipitation values. \*\*Rows\*\*
  represent spatial units (cells) and \*\*columns\*\* represent temporal
  units (e.g., 12 months).

- bio12V:

  A numeric \*\*vector\*\* or \*\*single-column matrix\*\* of Total
  Annual Precipitation (Bio12). Its length (or number of rows) must be
  exactly equal to the number of rows in \`prcp\`.

- n_units:

  A single \*\*integer\*\* representing the number of temporal units
  (e.g., 12). This is used to calculate the mean precipitation from the
  total \`bio12V\`.

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`prcp\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio15" (Precipitation Seasonality)
and "cell".
