# Check for issue #27: tiled backends read the wrong layers for layer subsets and
# reject multi-band files.
#
# Run with: Rscript dev/issue-checks/issue-27.R

suppressPackageStartupMessages({
  library(fastbioclim)
  library(terra)
})
cat("Issue #27: Tiled backends read the wrong layers for layer subsets and reject multi-band files\n")
cat("fastbioclim version:", as.character(packageVersion("fastbioclim")), "\n")

old_s2 <- suppressMessages(sf::sf_use_s2())
set.seed(27)
template <- rast(nrows = 4, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 4, crs = "EPSG:4326")
# One 24-band file. Bands 13-24 differ clearly from bands 1-12.
mb_file <- tempfile(fileext = ".tif")
mb <- rast(template, nlyrs = 24)
values(mb) <- cbind(matrix(runif(ncell(template) * 12, 0, 10), ncol = 12),
                    matrix(runif(ncell(template) * 12, 100, 200), ncol = 12))
writeRaster(mb, mb_file, datatype = "FLT8S", overwrite = TRUE)
mb <- rast(mb_file)
sub <- mb[[13:24]]
region <- vect(ext(mb), crs = crs(mb))   # full extent

try_run <- function(expr) {
  res <- tryCatch(suppressWarnings(suppressMessages(expr)), error = function(e) e)
  suppressMessages(sf::sf_use_s2(old_s2))
  res
}
same <- function(a, b) {
  !inherits(a, "error") && !inherits(b, "error") &&
    isTRUE(all.equal(unname(values(a)), unname(values(b)), tolerance = 1e-5))
}
show <- function(a, b) {
  if (inherits(a, "error")) return(paste("error:", conditionMessage(a)))
  paste0("max |tiled - terra| = ", signif(max(abs(values(a) - values(b))), 4))
}

# derive_bioclim(): a 12-layer subset of a multi-band file.
bio_terra <- try_run(derive_bioclim(bios = 12, prcp = sub, method = "terra",
                                    output_dir = tempfile(), verbose = FALSE))
bio_tiled <- try_run(derive_bioclim(bios = 12, prcp = sub, method = "tiled",
                                    output_dir = tempfile(), verbose = FALSE))
# calculate_average(): the same subset. A full-extent user_region avoids issue #12.
avg_terra <- try_run(calculate_average(sub, index = 1:12, output_names = "s", method = "terra",
                                       output_dir = tempfile(), verbose = FALSE))
avg_tiled <- try_run(calculate_average(sub, index = 1:12, output_names = "s", method = "tiled",
                                       user_region = region, output_dir = tempfile(),
                                       verbose = FALSE))

ok_bio <- same(bio_tiled, bio_terra)
ok_avg <- same(avg_tiled, avg_terra)
cat("EXPECTED: tiled == terra for derive_bioclim(bios = 12) and calculate_average() on mb[[13:24]]\n")
cat("GOT:      derive_bioclim:", if (ok_bio) "equal" else show(bio_tiled, bio_terra),
    "| calculate_average:", if (ok_avg) "equal" else show(avg_tiled, avg_terra), "\n")
cat("RESULT:", if (ok_bio && ok_avg) "PASS" else "FAIL", "\n")
