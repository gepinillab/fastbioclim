# Check for issue #55: stats_fast() ignores write_raw_vars because
# BIOCLIM_DEBUG_RAW_VARS overrides it.
#
# Run with: Rscript dev/issue-checks/issue-55.R

suppressPackageStartupMessages({
  library(fastbioclim)
  library(terra)
})
cat("Issue #55: stats_fast(): write_raw_vars argument is ignored because BIOCLIM_DEBUG_RAW_VARS overrides it\n")
cat("fastbioclim version:", as.character(packageVersion("fastbioclim")), "\n")

old_s2 <- suppressMessages(sf::sf_use_s2())
set.seed(55)
dir_in <- tempfile("in_")
dir.create(dir_in)
template <- rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3, crs = "EPSG:4326")
files <- file.path(dir_in, sprintf("tmin_%02d.tif", 1:12))
for (i in 1:12) {
  r <- rast(template)
  values(r) <- rnorm(ncell(template), 10, 5)
  writeRaster(r, files[i], datatype = "FLT8S", overwrite = TRUE)
}

# Number of raw variable tiles written by stats_fast().
n_raw <- function(...) {
  dir_out <- suppressWarnings(suppressMessages(
    fastbioclim:::stats_fast(variable_path = files, n_units = 12, stats = "mean",
                             output_prefix = "tmin", output_dir = tempfile(), verbose = FALSE, ...)
  ))
  suppressMessages(sf::sf_use_s2(old_s2))
  length(list.files(dir_out, pattern = "^raw_"))
}

Sys.unsetenv("BIOCLIM_DEBUG_RAW_VARS")
arg_true <- n_raw(write_raw_vars = TRUE)
Sys.setenv(BIOCLIM_DEBUG_RAW_VARS = "TRUE")
arg_false <- n_raw(write_raw_vars = FALSE)
no_arg <- n_raw()
Sys.unsetenv("BIOCLIM_DEBUG_RAW_VARS")

cat("EXPECTED: raw tiles written with write_raw_vars = TRUE (> 0), none with FALSE even if the env var is TRUE (0), and the env var used when the argument is missing (> 0)\n")
cat("GOT:      TRUE:", arg_true, "| FALSE with env var:", arg_false, "| missing with env var:", no_arg, "\n")
ok <- arg_true > 0 && arg_false == 0 && no_arg > 0
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
