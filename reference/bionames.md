# Print Bioclimatic Variable Names

This function prints the names of bioclimatic variables based on the
specified indices.

## Usage

``` r
bionames(bios = 1:35)
```

## Arguments

- bios:

  Numeric vector indicating the indices of bioclimatic variables to
  print. Default is 1:35, which prints all variable names.

## Value

None. Prints the names of the selected bioclimatic variables to the
console.

## Examples

``` r
bionames()           # Print all bioclimatic variable names
#> bio01: Mean Temperature
#> bio02: Mean Diurnal Range
#> bio03: Isothermality
#> bio04: Temperature Seasonality
#> bio05: Max Temperature of Warmest Unit
#> bio06: Min Temperature of Coldest Unit
#> bio07: Temperature Range
#> bio08: Mean Temperature of Wettest Period
#> bio09: Mean Temperature of Driest Period
#> bio10: Mean Temperature of Warmest Period
#> bio11: Mean Temperature of Coldest Period
#> bio12: Total Precipitation
#> bio13: Precipitation of Wettest Unit
#> bio14: Precipitation of Driest Unit
#> bio15: Precipitation Seasonality
#> bio16: Precipitation of Wettest Period
#> bio17: Precipitation of Driest Period
#> bio18: Precipitation of Warmest Period
#> bio19: Precipitation of Coldest Period
#> bio20: Mean Radiation
#> bio21: Highest Radiation Unit
#> bio22: Lowest Radiation Unit
#> bio23: Radiation Seasonality
#> bio24: Radiation of Wettest Period
#> bio25: Radiation of Driest Period
#> bio26: Radiation of Warmest Period
#> bio27: Radiation of Coldest Period
#> bio28: Mean Moisture
#> bio29: Highest Moisture Unit
#> bio30: Lowest Moisture Unit
#> bio31: Moisture Seasonality
#> bio32: Mean Moisture of Most Moist Period
#> bio33: Mean Moisture of Least Moist Period
#> bio34: Mean Moisture of Warmest Period
#> bio35: Mean Moisture of Coldest Period
bionames(c(1, 5, 12)) # Print names for variables 1, 5, and 12
#> bio01: Mean Temperature
#> bio05: Max Temperature of Warmest Unit
#> bio12: Total Precipitation
```
