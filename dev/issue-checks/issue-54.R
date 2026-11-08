# Check for issue #54: tiled calculate_*() fail with "No intermediate '.qs2' files
# found" when the input files share a band name.
#
# Run with: Rscript dev/issue-checks/issue-54.R

suppressPackageStartupMessages({
  library(fastbioclim)
  library(terra)
})
cat("Issue #54: Tiled calculate_*() fail with \"No intermediate '.qs2' files found\" when input files share a band name\n")
cat("fastbioclim version:", as.character(packageVersion("fastbioclim")), "\n")

old_s2 <- suppressMessages(sf::sf_use_s2())
set.seed(54)
dir_in <- tempfile("in_")
dir.create(dir_in)
template <- rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3, crs = "EPSG:4326")
vals <- matrix(rnorm(ncell(template) * 6, 10, 5), ncell(template))
files <- file.path(dir_in, sprintf("x_%02d.tif", 1:6))
for (i in 1:6) {
  r <- rast(template)
  values(r) <- vals[, i]
  names(r) <- "lyr.1"          # every file stores the same band name
  writeRaster(r, files[i], datatype = "FLT8S", overwrite = TRUE)
}
x <- rast(files)
index <- rep(1:2, each = 3)
# A full-extent user_region, because the tiled calculate_*() fail without one (#12, not
# part of this check).
region <- vect(ext(x), crs = crs(x))

run <- function(fun) {
  tryCatch({
    res <- suppressWarnings(suppressMessages(
      fun(x, index = index, user_region = region, method = "tiled", output_dir = tempfile(),
          verbose = FALSE)
    ))
    values(res)
  }, error = function(e) conditionMessage(e), finally = suppressMessages(sf::sf_use_s2(old_s2)))
}
avg <- run(calculate_average)
sm <- run(calculate_sum)
exp_avg <- cbind(rowMeans(vals[, 1:3]), rowMeans(vals[, 4:6]))
exp_sum <- cbind(rowSums(vals[, 1:3]), rowSums(vals[, 4:6]))
matches <- function(v, e) is.matrix(v) && isTRUE(all.equal(unname(v), unname(e), tolerance = 1e-5))
describe <- function(v) if (is.matrix(v)) paste("values, cell 1 =", paste(signif(v[1, ], 6), collapse = ", ")) else paste("error:", v)

cat("EXPECTED: group means and sums on the tiled backend; cell 1 mean =",
    paste(signif(exp_avg[1, ], 6), collapse = ", "), "\n")
cat("GOT:      calculate_average", describe(avg), "\n")
cat("          calculate_sum    ", describe(sm), "\n")
cat("RESULT:", if (matches(avg, exp_avg) && matches(sm, exp_sum)) "PASS" else "FAIL", "\n")
