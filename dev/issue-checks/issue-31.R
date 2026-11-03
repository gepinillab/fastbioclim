# Check for issue #31: tiled outputs lose the CRS identity (EPSG code dropped).
#
# Run with: Rscript dev/issue-checks/issue-31.R

suppressPackageStartupMessages({
  library(fastbioclim)
  library(terra)
})
cat("Issue #31: Tiled outputs lose the CRS identity (EPSG code dropped)\n")
cat("fastbioclim version:", as.character(packageVersion("fastbioclim")), "\n")

old_s2 <- suppressMessages(sf::sf_use_s2())
set.seed(31)
dir_in <- tempfile("in_")
dir.create(dir_in)
# 4 x 5 cells of 10 km in EPSG:3035 (ETRS89 / LAEA Europe).
template <- rast(nrows = 4, ncols = 5, xmin = 4000000, xmax = 4050000,
                 ymin = 3000000, ymax = 3040000, crs = "EPSG:3035")
files <- file.path(dir_in, sprintf("prcp_%02d.tif", 1:12))
for (i in 1:12) {
  r <- rast(template)
  values(r) <- runif(ncell(r), 0, 100)
  writeRaster(r, files[i], datatype = "FLT8S", overwrite = TRUE)
}
p <- rast(files)

run <- function(method) {
  res <- suppressWarnings(suppressMessages(
    derive_bioclim(bios = 12, prcp = p, method = method, tile_degrees = 1e6,
                   output_dir = tempfile(), verbose = FALSE)
  ))
  suppressMessages(sf::sf_use_s2(old_s2))
  res
}
terra_out <- run("terra")
tiled_out <- run("tiled")
code <- function(x) crs(x, describe = TRUE)$code

cat("EXPECTED: EPSG code 3035 on both outputs; same.crs(input, tiled) TRUE\n")
cat("GOT:      terra code", code(terra_out), "| tiled code", code(tiled_out),
    "| same.crs(input, tiled):", same.crs(p, tiled_out), "\n")
ok <- identical(code(terra_out), "3035") && identical(code(tiled_out), "3035") &&
  isTRUE(same.crs(p, tiled_out))
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
