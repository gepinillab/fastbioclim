# The base code for this test were created with Gemini Pro 2.5 (Jan 2025) using 
# as an input the function and asking to generated test. It was later modified and checked manually by
# Gonzalo E. Pinilla-Buitrago.
# tests/testthat/test-derive_bioclim.R

# --- Test Helper Functions ---
# It's good practice to create helpers for complex setup. This function creates
# small, file-backed dummy rasters for our tests.
create_dummy_rasters <- function(dir) {
  # Create a template raster
  r <- rast(nrows = 10, ncols = 10, nlyr = 12, crs = "EPSG:4326",
            xmin = 0, xmax = 10, ymin = 0, ymax = 10)
  
  # Create file-backed tmin, tmax, prcp rasters
  tmin_vals <- runif(ncell(r) * nlyr(r), 0, 10)
  tmin_rast <- setValues(r, tmin_vals)
  names(tmin_rast) <- paste0("tmin_", 1:12)
  writeRaster(tmin_rast, file.path(dir, "tmin.tif"), overwrite = TRUE)

  tmax_vals <- tmin_vals + runif(ncell(r) * nlyr(r), 5, 10)
  tmax_rast <- setValues(r, tmax_vals)
  names(tmax_rast) <- paste0("tmax_", 1:12)
  writeRaster(tmax_rast, file.path(dir, "tmax.tif"), overwrite = TRUE)
  
  prcp_vals <- runif(ncell(r) * nlyr(r), 20, 100)
  prcp_rast <- setValues(r, prcp_vals)
  names(prcp_rast) <- paste0("prcp_", 1:12)
  writeRaster(prcp_rast, file.path(dir, "prcp.tif"), overwrite = TRUE)
  
  return(list(
    tmin = rast(file.path(dir, "tmin.tif")),
    tmax = rast(file.path(dir, "tmax.tif")),
    prcp = rast(file.path(dir, "prcp.tif"))
  ))
}

# --- Tests ---

test_that("Input validation and initial checks work correctly", {
  # Create a temporary directory for outputs that is cleaned up automatically
  local_dir <- withr::local_tempdir()
  
  # 1. Stops if no raster inputs are provided at all
  expect_error(
    derive_bioclim(bios = 1:2),
    regexp = "No input SpatRaster objects were provided"
  )
  
  # 2. Stops if static indices are not SpatRaster objects
  dummy_rasters <- create_dummy_rasters(withr::local_tempdir())
  expect_error(
    derive_bioclim(bios = 1, tmin = dummy_rasters$tmin, warmest_unit = "not_a_raster"),
    regexp = "All static indices provided via '...' must be SpatRaster objects"
  )
  
  # 3. Stops if `user_region` does not overlap with input rasters
  non_overlapping_poly <- vect("POLYGON((100 100, 101 100, 101 101, 100 101, 100 100))", crs="EPSG:4326")
  expect_error(
    derive_bioclim(bios = 1, tmin = dummy_rasters$tmin, user_region = non_overlapping_poly),
    regexp = "does not overlap with the input rasters"
  )
})


test_that("File `overwrite` logic is correctly handled", {
  local_out_dir <- withr::local_tempdir()
  dummy_rasters <- create_dummy_rasters(withr::local_tempdir())
  
  # Create a fake existing output file
  existing_file_path <- file.path(local_out_dir, "bio01.tif")
  file.create(existing_file_path)
  
  # 1. Stops if a file exists and overwrite is FALSE (default)
  expect_error(
    derive_bioclim(bios = 1, tmin = dummy_rasters$tmin, tmax = dummy_rasters$tmax, output_dir = local_out_dir,
                 # Mock the downstream function to avoid actual work
                 bioclim_terra = function(...) "ok"),
    regexp = "Output files already exist"
  )
  
  # 2. Warns if a file exists and overwrite is TRUE
  expect_warning(
    derive_bioclim(bios = 1, tmin = dummy_rasters$tmin, tmax = dummy_rasters$tmax, output_dir = local_out_dir,
                 overwrite = TRUE, verbose = FALSE,
                 # Mock the downstream function
                 bioclim_terra = function(...) "ok"),
    regexp = "Overwriting existing files"
  )
})


test_that("Workflow dispatch logic selects the correct method", {
  local_dir <- withr::local_tempdir()
  dummy_rasters <- create_dummy_rasters(local_dir)

  stub(derive_bioclim, 'write_layers', function(...) "final_files")
  
  # FIX: Stub the final terra::rast call for ALL cases in this block.
  # This identity function will simply return whatever it's given,
  # preventing the file system error.
  stub(derive_bioclim, 'terra::rast', function(x) x)

  # --- Case 1: User forces 'terra' workflow ---
  stub(derive_bioclim, 'bioclim_terra', function(...) "called_terra_workflow")
  
  res_terra <- derive_bioclim(bios = 1, tmin = dummy_rasters$tmin,
                              method = "terra", verbose = FALSE)
  # The call no longer errors, and res_terra now gets the correct string.
  expect_equal(res_terra, "called_terra_workflow")
  
  # --- Case 2: User forces 'tiled' workflow ---
  stub(derive_bioclim, 'bioclim_terra', function(...) stop("Terra should not be called!"))
  stub(derive_bioclim, 'bioclim_fast', function(...) "called_tiled_workflow")

  res_tiled <- derive_bioclim(bios = 1, tmin = dummy_rasters$tmin,
                              method = "tiled", verbose = FALSE)
  expect_equal(res_tiled, "final_files")
  
  # --- Case 3: 'auto' method selects 'terra' for small data ---
  stub(derive_bioclim, 'terra::mem_info', function(...) c(fits_mem = 1))
  stub(derive_bioclim, 'bioclim_terra', function(...) "called_terra_auto")
  stub(derive_bioclim, 'bioclim_fast', function(...) stop("Tiled should not be called!"))

  res_auto_terra <- derive_bioclim(bios = 1, tmin = dummy_rasters$tmin,
                                   method = "auto", verbose = FALSE)
  expect_equal(res_auto_terra, "called_terra_auto")
  
  # --- Case 4: 'auto' method selects 'tiled' for large data ---
  stub(derive_bioclim, 'terra::mem_info', function(...) c(fits_mem = 0))
  stub(derive_bioclim, 'bioclim_terra', function(...) stop("Terra should not be called!"))
  stub(derive_bioclim, 'bioclim_fast', function(...) "called_tiled_auto")
  
  res_auto_tiled <- derive_bioclim(bios = 1, tmin = dummy_rasters$tmin,
                                   method = "auto", verbose = FALSE)
  expect_equal(res_auto_tiled, "final_files")
})


test_that("Static indices are passed correctly to the chosen workflow", {
  local_dir <- withr::local_tempdir()
  dummy_rasters <- create_dummy_rasters(local_dir)
  static_rast <- rast(dummy_rasters$tmin, vals = 2, nlyr=1)
  names(static_rast) <- "static"
  
  # Stub for the final terra::rast call (from previous fix)
  stub(derive_bioclim, 'terra::rast', function(x) x)
  
  # FIX: Stub `write_layers` to prevent the real function from being called
  # in the tiled workflow. This completes the isolation of derive_bioclim.
  stub(derive_bioclim, 'write_layers', function(...) { 
    return("path/to/final_files") 
  })
  
  # --- Test 1: Passed correctly to 'terra' workflow ---
  terra_args <- NULL
  stub(derive_bioclim, 'bioclim_terra', function(...) {
    terra_args <<- list(...)
    return("ok") # This is now harmless
  })
  
  derive_bioclim(bios = 1, tmin = dummy_rasters$tmin, method = "terra",
                 warmest_period = static_rast, verbose = FALSE)
  
  expect_true("warmest_period" %in% names(terra_args))
  expect_s4_class(terra_args$warmest_period, "SpatRaster")
  
  # --- Test 2: Passed correctly to 'tiled' workflow ---
  tiled_args <- NULL
  stub(derive_bioclim, 'bioclim_fast', function(...) {
    tiled_args <<- list(...)
    return("ok") # This is now harmless because it's passed to the STUBBED write_layers
  })
  writeRaster(static_rast, file.path(local_dir, "static.tif"), overwrite = TRUE)
  
  # This call will now succeed because the call to write_layers is intercepted
  derive_bioclim(bios = 1, tmin = dummy_rasters$tmin, method = "tiled",
                 warmest_period = rast(file.path(local_dir, "static.tif")), 
                 verbose = FALSE)
  
  expect_true("warmest_period_path" %in% names(tiled_args))
  expect_type(tiled_args$warmest_period_path, "character")
  expect_match(tiled_args$warmest_period_path, "static.tif")
})