# The base code for this test were created with Gemini Pro 2.5 (Jan 2025) using 
# as an input the function and asking to generated test. It was later modified and checked manually by
# Gonzalo E. Pinilla-Buitrago.

# tests/testthat/test-derive_statistics.R
# --- Test Helper Functions ---
# A helper to create dummy file-backed rasters for statistics tests.
create_dummy_stat_rasters <- function(dir) {
  r <- rast(nrows = 10, ncols = 10, nlyr = 12, crs = "EPSG:4326",
            xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  
  # Create a file-backed primary variable
  var_vals <- runif(ncell(r) * nlyr(r), 50, 150)
  var_rast <- setValues(r, var_vals)
  names(var_rast) <- paste0("var_", 1:12)
  writeRaster(var_rast, file.path(dir, "variable.tif"), overwrite = TRUE)

  # Create a file-backed interactive variable
  inter_vals <- runif(ncell(r) * nlyr(r), 0, 1)
  inter_rast <- setValues(r, inter_vals)
  names(inter_rast) <- paste0("inter_", 1:12)
  writeRaster(inter_rast, file.path(dir, "inter_variable.tif"), overwrite = TRUE)
  
  return(list(
    variable = rast(file.path(dir, "variable.tif")),
    inter_variable = rast(file.path(dir, "inter_variable.tif"))
  ))
}

# --- Tests ---

test_that("Input validation and initial checks work correctly", {
  local_dir <- withr::local_tempdir()
  dummy_rasters <- create_dummy_stat_rasters(local_dir)
  
  # 1. Stops if the main `variable` is missing
  expect_error(
    derive_statistics(variable = NULL, stats = "mean"),
    regexp = "The 'variable' SpatRaster must be provided"
  )
  
  # 2. Stops if static indices are not SpatRaster objects
  expect_error(
    derive_statistics(variable = dummy_rasters$variable, max_unit = "not_a_raster"),
    regexp = "All static indices provided via '...' must be SpatRaster objects"
  )
  
  # 3. Stops if `user_region` does not overlap
  non_overlapping_poly <- vect("POLYGON((100 100, 101 100, 101 101, 100 101, 100 100))", crs="EPSG:4326")
  expect_error(
    derive_statistics(variable = dummy_rasters$variable, user_region = non_overlapping_poly),
    regexp = "does not overlap with the input rasters"
  )
})


test_that("File `overwrite` logic correctly handles dynamic filenames", {
  local_out_dir <- withr::local_tempdir()
  dummy_rasters <- create_dummy_stat_rasters(withr::local_tempdir())
  
  # Create a fake existing output file based on prefix and stat
  prefix <- "test_var"
  existing_file_path <- file.path(local_out_dir, paste0(prefix, "_mean.tif"))
  file.create(existing_file_path)

  # Mock the downstream function to avoid actual work
  stub(derive_statistics, 'stats_terra', function(...) rast())
  
  # 1. Stops if a file exists and overwrite is FALSE
  expect_error(
    derive_statistics(variable = dummy_rasters$variable, stats = "mean",
                    output_prefix = prefix, output_dir = local_out_dir),
    regexp = "Potential output files already exist"
  )
  
  # 2. Does NOT stop if overwrite is TRUE (the function itself has no warning here)
  expect_no_error(
     derive_statistics(variable = dummy_rasters$variable, stats = "mean",
                    output_prefix = prefix, output_dir = local_out_dir,
                    overwrite = TRUE, verbose = FALSE)
  )
})


test_that("Workflow dispatch logic selects the correct method", {
  local_dir <- withr::local_tempdir()
  dummy_rasters <- create_dummy_stat_rasters(local_dir)
  
  # Mock downstream functions to prevent actual computation
  stub(derive_statistics, 'write_layers', function(...) "called_write_layers")
  stub(derive_statistics, 'terra::rast', function(x) x) # Pass through results
  
  # --- Case 1: User forces 'terra' workflow ---
  stub(derive_statistics, 'stats_terra', function(...) "called_terra_workflow")
  res_terra <- derive_statistics(variable = dummy_rasters$variable,
                                 method = "terra", verbose = FALSE)
  expect_equal(res_terra, "called_terra_workflow")
  
  # --- Case 2: User forces 'tiled' workflow ---
  stub(derive_statistics, 'stats_terra', function(...) stop("Terra should not be called!"))
  stub(derive_statistics, 'stats_fast', function(...) "called_tiled_workflow")
  res_tiled <- derive_statistics(variable = dummy_rasters$variable,
                                 method = "tiled", verbose = FALSE)
  expect_equal(res_tiled, "called_write_layers")
  
  # --- Case 3: 'auto' selects 'terra' for small data ---
  stub(derive_statistics, 'terra::mem_info', function(...) c(fits_mem = 1))
  stub(derive_statistics, 'stats_terra', function(...) "called_terra_auto")
  stub(derive_statistics, 'stats_fast', function(...) stop("Tiled should not be called!"))
  res_auto_terra <- derive_statistics(variable = dummy_rasters$variable,
                                      method = "auto", verbose = FALSE)
  expect_equal(res_auto_terra, "called_terra_auto")

  # --- Case 4: 'auto' selects 'tiled' for large data ---
  stub(derive_statistics, 'terra::mem_info', function(...) c(fits_mem = 0))
  stub(derive_statistics, 'stats_terra', function(...) stop("Terra should not be called!"))
  stub(derive_statistics, 'stats_fast', function(...) "called_tiled_auto")
  res_auto_tiled <- derive_statistics(variable = dummy_rasters$variable,
                                      method = "auto", verbose = FALSE)
  expect_equal(res_auto_tiled, "called_write_layers")
})


test_that("Arguments are passed correctly to the chosen workflow", {
  local_dir <- withr::local_tempdir()
  dummy_rasters <- create_dummy_stat_rasters(local_dir)
  static_rast <- rast(dummy_rasters$variable, vals = 1, nlyr=1)
  
  # FIX: Stub the final terra::rast call to prevent file system errors.
  # This will intercept the "ok" string from the write_layers stub.
  stub(derive_statistics, 'terra::rast', function(x) x)

  # --- Test 1: Arguments passed to 'terra' workflow ---
  terra_args <- NULL
  stub(derive_statistics, 'stats_terra', function(...) {
    terra_args <<- list(...)
    # This return is fine because the terra workflow doesn't call terra::rast at the end
    return(rast()) 
  })
  
  derive_statistics(variable = dummy_rasters$variable, 
                    inter_variable = dummy_rasters$inter_variable,
                    stats = c("mean", "max"), output_prefix = "wind",
                    method = "terra", max_unit = static_rast, verbose = FALSE)
  
  expect_true("variable" %in% names(terra_args))
  expect_s4_class(terra_args$variable, "SpatRaster")
  expect_true("inter_variable" %in% names(terra_args))
  expect_true("max_unit" %in% names(terra_args))
  expect_equal(terra_args$output_prefix, "wind")
  
  # --- Test 2: Arguments passed to 'tiled' workflow ---
  tiled_args <- NULL
  stub(derive_statistics, 'stats_fast', function(...) {
    tiled_args <<- list(...)
    return("ok")
  })
  stub(derive_statistics, 'write_layers', function(...) "ok")
  
  writeRaster(static_rast, file.path(local_dir, "static.tif"), overwrite = TRUE)
  
  # This call will now succeed
  derive_statistics(variable = dummy_rasters$variable, 
                    inter_variable = dummy_rasters$inter_variable,
                    stats = c("mean", "max"), output_prefix = "wind",
                    method = "tiled", 
                    max_unit = rast(file.path(local_dir, "static.tif")), 
                    verbose = FALSE)
  
  expect_true("variable" %in% names(tiled_args))
  expect_type(tiled_args$variable, "character")
  expect_match(tiled_args$variable, "variable.tif")
  
  expect_true("inter_variable" %in% names(tiled_args))
  expect_match(tiled_args$inter_variable, "inter_variable.tif")
  
  expect_true("max_unit_path" %in% names(tiled_args))
  expect_match(tiled_args$max_unit_path, "static.tif")
})