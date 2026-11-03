# Check for issue #11: tiled + user_region writes pixels from outside the region
# into wrong cells, and keeps edge cells that the terra backend masks.
#
# Run with: Rscript dev/issue-checks/issue-11.R
# Uses the installed fastbioclim (see compare.R to run it against two versions).

suppressPackageStartupMessages({
  library(fastbioclim)
  library(terra)
  library(sf)
})
cat("Issue #11: Tiled + user_region writes pixels from outside the region into wrong cells\n")
cat("fastbioclim version:", as.character(packageVersion("fastbioclim")), "\n")

old_s2 <- suppressMessages(sf::sf_use_s2())
dir_in <- tempfile("in_")
dir.create(dir_in)

# 10 x 10 raster, 12 single-layer files; every layer holds its source cell id.
template <- rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10,
                 crs = "EPSG:4326")
files <- file.path(dir_in, sprintf("id_%02d.tif", 1:12))
for (f in files) {
  r <- rast(template)
  values(r) <- 1:100
  writeRaster(r, f, datatype = "FLT8S", overwrite = TRUE)
}
id <- rast(files)

# L-shaped region whose edges cut through cells.
reg <- st_as_sf(st_union(
  st_as_sfc(st_bbox(c(xmin = 2.3, ymin = 4, xmax = 6.3, ymax = 6.3), crs = st_crs(4326))),
  st_as_sfc(st_bbox(c(xmin = 4, ymin = 2, xmax = 6.3, ymax = 4), crs = st_crs(4326)))
))

run <- function(method) {
  res <- tryCatch(
    suppressWarnings(suppressMessages(
      derive_statistics(id, stats = "mean", user_region = reg, method = method,
                        tile_degrees = 2, output_dir = tempfile(), verbose = FALSE)
    )),
    error = function(e) e
  )
  suppressMessages(sf::sf_use_s2(old_s2))
  res
}
terra_res <- run("terra")
tiled_res <- run("tiled")

if (inherits(terra_res, "error") || inherits(tiled_res, "error")) {
  msg <- if (inherits(tiled_res, "error")) conditionMessage(tiled_res) else conditionMessage(terra_res)
  cat("EXPECTED: tiled and terra outputs identical (same extent, NA cells and values)\n")
  cat("GOT:      error:", msg, "\n")
  cat("RESULT: FAIL\n")
} else {
  tv <- values(terra_res)[, 1]
  dv <- values(tiled_res)[, 1]
  tv[is.nan(tv)] <- NA
  dv[is.nan(dv)] <- NA
  same_geom <- compareGeom(terra_res, tiled_res, stopOnError = FALSE, messages = FALSE)
  same_na <- identical(which(is.na(tv)), which(is.na(dv)))
  same_val <- isTRUE(all.equal(tv[!is.na(tv)], dv[!is.na(dv)]))
  diff_cells <- which(xor(is.na(tv), is.na(dv)) | (!is.na(tv) & !is.na(dv) & abs(tv - dv) > 1e-6))
  cat("EXPECTED: tiled and terra outputs identical; terra non-NA cells:",
      paste(which(!is.na(tv)), collapse = ","), "\n")
  cat("GOT:      tiled non-NA cells:", paste(which(!is.na(dv)), collapse = ","), "\n")
  if (length(diff_cells) > 0) {
    cat("          cells that differ:", paste(diff_cells, collapse = ","),
        "| terra:", paste(tv[diff_cells], collapse = ","),
        "| tiled:", paste(dv[diff_cells], collapse = ","), "\n")
  }
  cat("RESULT:", if (same_geom && same_na && same_val) "PASS" else "FAIL", "\n")
}
