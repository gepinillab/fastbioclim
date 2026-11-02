make_tiled_series <- function(path, name) {
  template <- terra::rast(
    nrows = 2,
    ncols = 2,
    nlyr = 1,
    xmin = 0,
    xmax = 2,
    ymin = 0,
    ymax = 2,
    crs = "EPSG:4326"
  )
  files <- file.path(path, paste0(name, "_", seq_len(12), ".tif"))

  for (i in seq_along(files)) {
    terra::writeRaster(
      terra::setValues(template, rep(i, terra::ncell(template))),
      files[[i]],
      overwrite = TRUE
    )
  }

  terra::rast(files)
}

expect_s2_restored <- function(code) {
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(TRUE)

  suppressMessages(withCallingHandlers(
    force(code),
    warning = function(w) {
      known_warnings <- c(
        "x is already of type POLYGON.",
        "No CRS specified for polygons; assuming they have the same CRS as the raster."
      )
      if (conditionMessage(w) %in% known_warnings) {
        invokeRestart("muffleWarning")
      }
    }
  ))

  expect_true(sf::sf_use_s2())
}

test_that("tiled backends restore the previous sf s2 setting", {
  skip_on_cran()
  path <- withr::local_tempdir()
  prcp <- make_tiled_series(path, "prcp")
  variable <- make_tiled_series(path, "variable")
  region <- terra::vect("POLYGON((0 0, 2 0, 2 2, 0 2, 0 0))", crs = "EPSG:4326")
  mockery::stub(calculate_average, "write_layers", function(...) "unused.tif")
  mockery::stub(calculate_average, "terra::rast", function(...) variable)
  mockery::stub(calculate_roll, "write_layers", function(...) "unused.tif")
  mockery::stub(calculate_roll, "terra::rast", function(...) variable)

  expect_s2_restored(derive_bioclim(
    bios = 12,
    prcp = prcp,
    method = "tiled",
    user_region = region,
    output_dir = file.path(path, "bioclim"),
    overwrite = TRUE,
    verbose = FALSE
  ))

  expect_s2_restored(derive_statistics(
    variable = variable,
    stats = "mean",
    method = "tiled",
    user_region = region,
    output_dir = file.path(path, "statistics"),
    overwrite = TRUE,
    verbose = FALSE
  ))

  expect_s2_restored(calculate_average(
    x = variable,
    index = rep(1:3, each = 4),
    method = "tiled",
    user_region = region,
    output_dir = file.path(path, "average"),
    overwrite = TRUE,
    verbose = FALSE
  ))

  expect_s2_restored(calculate_roll(
    x = variable,
    window_size = 1,
    method = "tiled",
    user_region = region,
    output_dir = file.path(path, "roll"),
    overwrite = TRUE,
    verbose = FALSE
  ))
})
