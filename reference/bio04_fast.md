# bio04_fast: Temperature Seasonality (Std Dev \* 100)

Calculates Temperature Seasonality, defined as the standard deviation of
average temperatures across all temporal units (e.g., 12 months),
multiplied by 100.

## Usage

``` r
bio04_fast(tavg, cell)
```

## Arguments

- tavg:

  A numeric \*\*matrix\*\* of average temperatures. \*\*Rows\*\*
  represent spatial units (cells) and \*\*columns\*\* represent temporal
  units.

- cell:

  A vector of original cell IDs. Its length must be exactly equal to the
  number of rows in \`tavg\`.

## Value

A \*\*matrix\*\* with dimensions \`c(N, 2)\`, where N is the number of
input cells. The columns are named "bio04" (seasonality) and "cell".
