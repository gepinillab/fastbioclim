# Check for issue #23: stats_fast() treats a region whose values are all <= 0 as
# empty and uses the full extent.
#
# Run with: Rscript dev/issue-checks/issue-23.R

suppressPackageStartupMessages({
  library(fastbioclim)
  library(terra)
})
cat("Issue #23: stats_fast() treats a region whose values are all <= 0 as empty and uses the full extent\n")
cat("fastbioclim version:", as.character(packageVersion("fastbioclim")), "\n")

old_s2 <- suppressMessages(sf::sf_use_s2())
set.seed(23)
dir_in <- tempfile("in_")
dir.create(dir_in)
template <- rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10, crs = "EPSG:4326")
files <- file.path(dir_in, sprintf("neg_%02d.tif", 1:12))
for (i in 1:12) {
  r <- rast(template)
  values(r) <- runif(ncell(r), -20, -1)   # all values negative
  writeRaster(r, files[i], datatype = "FLT8S", overwrite = TRUE)
}
neg <- rast(files)
region <- vect(ext(2, 6, 2, 6), crs = "EPSG:4326")

warns <- character()
tiled <- withCallingHandlers(
  suppressMessages(derive_statistics(neg, stats = "mean", user_region = region, method = "tiled",
                                     output_dir = tempfile(), verbose = FALSE)),
  warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
suppressMessages(sf::sf_use_s2(old_s2))
terra_out <- suppressMessages(derive_statistics(neg, stats = "mean", user_region = region,
                                                method = "terra", output_dir = tempfile(),
                                                verbose = FALSE))

got_ext <- as.vector(ext(tiled))
empty_warn <- any(grepl("empty", warns))
cat("EXPECTED: tiled extent 2 6 2 6, no 'empty' warning, values equal to terra\n")
cat("GOT:      tiled extent", paste(got_ext, collapse = " "),
    "| 'empty' warning:", empty_warn, "\n")
ok <- isTRUE(all.equal(unname(got_ext), c(2, 6, 2, 6))) && !empty_warn &&
  compareGeom(tiled, terra_out, stopOnError = FALSE, messages = FALSE) &&
  isTRUE(all.equal(values(tiled)[, 1], values(terra_out)[, 1], tolerance = 1e-5))
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
