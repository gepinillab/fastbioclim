# Check for issue #22: derive_statistics(stats = "stdev") gives the population SD
# on terra and the sample SD on tiled.
#
# Run with: Rscript dev/issue-checks/issue-22.R

suppressPackageStartupMessages({
  library(fastbioclim)
  library(terra)
})
cat("Issue #22: derive_statistics(stats = \"stdev\"): population SD on terra, sample SD on tiled\n")
cat("fastbioclim version:", as.character(packageVersion("fastbioclim")), "\n")

old_s2 <- suppressMessages(sf::sf_use_s2())
set.seed(22)
dir_in <- tempfile("in_")
dir.create(dir_in)
template <- rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3, crs = "EPSG:4326")
vals <- matrix(rnorm(ncell(template) * 12, 10, 5), ncell(template))
files <- file.path(dir_in, sprintf("tmin_%02d.tif", 1:12))
for (i in 1:12) {
  r <- rast(template)
  values(r) <- vals[, i]
  writeRaster(r, files[i], datatype = "FLT8S", overwrite = TRUE)
}
tmin <- rast(files)

run <- function(method) {
  res <- suppressWarnings(suppressMessages(
    derive_statistics(tmin, stats = "stdev", method = method, output_prefix = "tmin",
                      output_dir = tempfile(), verbose = FALSE)
  ))
  suppressMessages(sf::sf_use_s2(old_s2))
  res
}
terra_sd <- values(run("terra"))[, 1]
tiled_sd <- values(run("tiled"))[, 1]
sample_sd <- apply(vals, 1, sd)

cat("EXPECTED: sample SD (sd()) on both backends; cell 1 =", signif(sample_sd[1], 7), "\n")
cat("GOT:      cell 1 terra =", signif(terra_sd[1], 7), "tiled =", signif(tiled_sd[1], 7), "\n")
ok <- isTRUE(all.equal(terra_sd, sample_sd, tolerance = 1e-5)) &&
  isTRUE(all.equal(tiled_sd, sample_sd, tolerance = 1e-5))
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
