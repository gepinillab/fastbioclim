[![CRAN status](https://www.r-pkg.org/badges/version/fastbioclim)](https://CRAN.R-project.org/package=fastbioclim)
[![downloads](https://cranlogs.r-pkg.org:443/badges/grand-total/fastbioclim?color=orange)](https://cranlogs.r-pkg.org:443/badges/grand-total/fastbioclim?color=orange)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check.yaml](https://github.com/gepinillab/fastbioclim/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/gepinillab/fastbioclim/actions/workflows/R-CMD-check.yaml)

# fastbioclim <a href="https://gepinillab.github.io/fastbioclim/"><img src="man/figures/logo.png" align="right" height="123" alt="fastbioclim website" /></a>

`fastbioclim` is an R package for creating custom-time bioclimatic and derived environmental summary variables from supplied raster data. It is designed to overcome computational bottlenecks and methodological inflexibility by automatically switching between processing frameworks to handle large-scale extents on standard hardware.

## Overview

Working with large climate datasets often presents a major challenge: the data is too large to fit into memory. `fastbioclim` addresses this gap by providing an efficient, unified interface that automatically switches between two frameworks:

1.  **In-Memory (`"terra"`)**: For smaller datasets, `fastbioclim` uses the `terra` package to maximize speed.
2.  **On-Disk Tiling (`"tiled"`)**: For rasters that exceed available RAM, it automatically switches to an on-disk tiling framework. This approach leverages the high-performance `Rfast` and `exactextractr` packages to process data in chunks, ensuring scalability on standard personal computers without requiring the end-user to manage these decisions.

The core of the package is a set of flexible functions—`derive_bioclim()` and `derive_statistics()`—that automatically select the optimal processing framework, providing a seamless experience for generating temporally-matched environmental variables.

## Key Features

*   **Smart `auto` Method**: Automatically manages memory by switching between the in-memory ("terra") and on-disk ("tiled") frameworks based on raster size and available RAM.
*   **Scalability**: Accelerates large-scale raster processing for variable creation, from smaller regional extents to global datasets.
*   **Extended Bioclimatic Variables**: Calculates an extended set of 35 bioclimatic variables, incorporating standard indices as well as variables based on moisture and solar radiation (bio20–35).
*   **Derived Summary Statistics**: Provides a flexible `derive_statistics()` function to compute summary statistics (mean, max, min, standard deviation, etc.) for any other time series data available (e.g., wind speed, evapotranspiration, cloud cover).
*   **Custom Timeframes**: Easily define time periods beyond standard quarters, enabling the summarization of data from finer time units such as weeks or days.
*   **Advanced Control**: Supports the use of a fixed temporal index (static indices) for single units or entire periods, which is critical for analyzing trends due to the temporal shifting of the indices.

## Installation

You can install the development version of `fastbioclim` from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("gepinillab/fastbioclim")
# Install to get the package example data 
remotes::install_github("gepinillab/egdata.fastbioclim")
```

## Usage

The package provides two primary core functions for variable calculation:

*   **`derive_bioclim()`**: For calculating the standard and extended set of 35 bioclimatic variables.
*   **`derive_statistics()`**: For deriving summary statistics from any other environmental variable.

*Note: The package also includes aggregation functions like calculate_average(), calculate_roll(), and calculate_sum() to easily prepare your time-series data*

### Example: A Quick Demonstration

This example demonstrates the core functionality using simple, self-contained mock data.

```r
library(fastbioclim)
library(terra)
library(future.apply)
library(progressr)

tmin_ecu <- system.file("extdata/ecuador/", package = "egdata.fastbioclim") |>
  list.files("tmin", full.names = TRUE) |> rast()
tmax_ecu <- system.file("extdata/ecuador/", package = "egdata.fastbioclim") |>
  list.files("tmax", full.names = TRUE) |> rast()
prcp_ecu <- system.file("extdata/ecuador/", package = "egdata.fastbioclim") |>
  list.files("prcp", full.names = TRUE) |> rast()

# The function will automatically use the fast in-memory "terra" method for this small dataset
output_dir_bioclim <- file.path(tempdir(), "bioclim_ecuador")

bioclim_vars <- derive_bioclim(
  bios = 1:19,
  tmin = tmin_ecu,
  tmax = tmax_ecu,
  prcp = prcp_ecu,
  output_dir = output_dir_bioclim,
  overwrite = TRUE
)

plot(bioclim_vars[[c("bio01", "bio12")]])
```

```r
# Derive environmental summary variables for a different factor (e.g., wind speed)
wind_rast <- system.file("extdata/ecuador/", package = "egdata.fastbioclim") |>
  list.files("wind", full.names = TRUE) |> rast()
output_dir_custom <- file.path(tempdir(), "wind_ecuador")

custom_stats <- derive_statistics(
  variable = wind_rast,
  stats = c("mean", "max", "stdev"),
  output_prefix = "wind",
  output_dir = output_dir_custom,
  overwrite = TRUE
)

plot(custom_stats)
```
### Handling Large-Scale Data (The Tiled Workflow)

The real power of `fastbioclim` shines with large datasets. The `method = "auto"` setting in `derive_bioclim()` and `derive_statistics()` handles this automatically.

When the wrapper function detects that the input rasters are too large to fit in memory, it seamlessly switches to the tiled workflow.

**Important Requirement:** For the tiled workflow to function, your input `SpatRaster` objects must be pointing to files on disk, not held entirely in memory.

```r
# Conceptual example for large, file-based rasters
tmin_neo <- system.file("extdata/neotropics/", package = "egdata.fastbioclim") |>
  list.files("tmin", full.names = TRUE) |> rast()
tmax_neo <- system.file("extdata/neotropics/", package = "egdata.fastbioclim") |>
  list.files("tmax", full.names = TRUE) |> rast()
prcp_neo <- system.file("extdata/neotropics/", package = "egdata.fastbioclim") |>
  list.files("prcp", full.names = TRUE) |> rast()
output_dir_bios <- file.path(tempdir(), "bioclim_neotropics")

# Optional: ACTIVATE PROGRESS BAR
progressr::handlers(global = TRUE)
# Optional: DEFINE PARALLEL PLAN FOR EVEN FASTER PROCESSING
future::plan("multisession", workers = 4)

# The call is identical. `derive_bioclim` will detect the large file size
# and automatically use the memory-safe tiled method.
large_scale_vars <- derive_bioclim(
  bios = 1:19,
  tmin = tmin_neo,
  tmax = tmax_neo,
  prcp = prcp_neo,
  output_dir = output_dir_bios,
  tile_degrees = 20,
  overwrite = TRUE
)
print(large_scale_vars)
plot(large_scale_vars[["bio11"]])
```

## Under Active Development

This R package is currently under active development. While it is functional, it may contain bugs or undergo changes to the API.

Contributions, bug reports, and feature requests are highly encouraged. Please [open an issue](https://github.com/gepinillab/fastbioclim/issues) on our GitHub repository to provide feedback.
