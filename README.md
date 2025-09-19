# `fastbioclim`: Scalable Derivation of Climate Variables

`fastbioclim` is an R package for efficiently deriving standard bioclimatic and custom summary variables from large-scale climate raster data. It is designed to overcome the memory limitations of traditional approaches by intelligently switching between processing backends.

## Overview

Working with large climate datasets often presents a major challenge: the data is too large to fit into memory. `fastbioclim` solves this problem by providing a powerful and unified interface with a dual-backend architecture:

1.  **In-Memory (`"terra"`)**: For datasets that fit comfortably in RAM, `fastbioclim` uses the highly optimized `terra` package for maximum speed.
2.  **In-Disk (`"tiled"`)**: For datasets that exceed available RAM, it automatically switches to a memory-safe tiled workflow. This backend uses `exactextractr` and `Rfast` to process the data chunk by chunk, ensuring that even continent-scale analyses can run on a standard computer.

The core of the package is a set of "smart" wrapper functions—`derive_bioclim()` and `derive_statistics()`—that automatically select the best backend, providing a seamless experience for the user.

## Key Features

*   **Smart `auto` Method**: Automatically detects data size and chooses the optimal processing workflow (in-memory vs. tiled).
*   **Scalability**: Robustly handles raster datasets of any size, from small study areas to global extents.
*   **Comprehensive Bioclimatic Variables**: Calculates the full set of 35 standard bioclimatic variables, including those based on temperature, precipitation, solar radiation, and moisture.
*   **Custom Statistics**: Provides a flexible `derive_statistics()` function to compute custom summaries (mean, max, min, standard deviation, etc.) for any climate variable.
*   **Flexible Time Periods**: Easily define custom periods (e.g., weeks, semesters) instead of being limited to standard quarters.
*   **Advanced Control**: Supports static index rasters for advanced time-series analysis and defining specific seasonal periods.

## Installation

You can install the development version of `fastbioclim` from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("gepinillab/fastbioclim")
# Install to get the package example data 
remotes::install_github("gepinillab/egdata.fastbioclim")
```

## Usage

The package provides two primary, user-facing functions:

*   **`derive_bioclim()`**: For calculating the standard set of 19-35 bioclimatic variables.
*   **`derive_statistics()`**: For calculating custom summary statistics on any variable.

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

# The function will automatically use the fast "terra" method for this small dataset
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
# 3. Derive custom summary statistics for a different variable (e.g., wind speed)
wind_rast <- system.file("extdata/ecuador/", package = "egdata.fastbioclim") |>
  list.files("wind", full.names = TRUE) |> rast()
output_dir_custom <- file.path(tempdir(), "wind_ecuador")

custom_stats <- derive_statistics(
  variable = wind_rast,
  stats = c("mean", "max", "stdev"),
  prefix_variable = "wind",
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
