#' Tiled, Out-of-Core Custom Variable Summarization
#'
#' Internal function to calculate custom summary statistics for very large datasets
#' by processing them in tiles.
#'
#' @param variable_path Path to primary variable rasters.
#' @param n_units Integer, number of layers per variable.
#' @param stats Character vector of stats to compute.
#' @param prefix_variable Character, prefix for output files.
#' @param ... Other arguments including inter_variable_path, period_length, circular,
#'   static index paths, etc.
#' @return Character string: Path to the temporary directory containing
#'   intermediate `.qs` files.
#' @keywords internal
#' @seealso The user-facing wrapper function `derive_statistics()`.
stats_fast <- function(variable_path,
  n_units,
  stats = c("mean", "max", "min"),
  period_length = 3,
  period_stats = "mean",
  circular = TRUE,
  inter_variable_path = NULL,
  inter_stats = NULL, # Default to NULL if inter_variable_path not given
  max_unit_path = NULL,
  min_unit_path = NULL,
  max_period_path = NULL,
  min_period_path = NULL,
  max_interactive_path = NULL,
  min_interactive_path = NULL,
  prefix_variable = "var",
  suffix_inter_max = "inter_high",
  suffix_inter_min = "inter_low",
  user_region = NULL,
  tile_degrees = 5,
  output_dir = tempdir(),
  write_raw_vars = FALSE,
  ...) {

  # --- 0. Input Validation, Dependency Mapping, Static Index Parsing ---
  dot_args <- list(...) # Currently not used but good practice
  static_index_paths <- list()

  # Validate core inputs
  if (missing(variable_path) || !is.character(variable_path) || length(variable_path) == 0) stop("'variable_path' is required.")
  if (missing(n_units) || !is.numeric(n_units) || n_units <= 0 || length(variable_path) != n_units) {
    stop(sprintf("'n_units' required and must match length of 'variable_path' (%d).", n_units))
  }
  if (!is.numeric(period_length) || period_length <= 0 || period_length > n_units) stop("'period_length' invalid.")
  if (!(period_stats %in% c("mean", "sum"))) stop("'period_stats' must be 'mean' or 'sum'.")

  valid_stats_opts <- c("mean", "max", "min", "sum", "stdev", "cv_cli", "max_period", "min_period")
  invalid_stats <- setdiff(stats, valid_stats_opts)
  if (length(invalid_stats) > 0) {
    stop(paste("Invalid 'stats' provided:", paste(invalid_stats, collapse = ", ")))
  }

  valid_inter_stats_opts <- c("max_inter", "min_inter")
  if (!is.null(inter_stats)) {
    invalid_inter_stats <- setdiff(inter_stats, valid_inter_stats_opts)
    if (length(invalid_inter_stats) > 0) {
      stop(paste("Invalid 'inter_stats' provided:", paste(invalid_inter_stats, collapse = ", ")))
    }
  } else {
    inter_stats <- c() # Ensure it's an empty vector if NULL
  }

  if (length(c(stats, inter_stats)) == 0) stop("No statistics requested to calculate.")


  # Parse and validate static index paths
  static_inputs_map <- list(
    max_unit = max_unit_path, 
    min_unit = min_unit_path,
    max_period = max_period_path, 
    min_period = min_period_path,
    max_interactive = max_interactive_path, 
    min_interactive = min_interactive_path
  )
  for (idx_name in names(static_inputs_map)) {
    idx_path <- static_inputs_map[[idx_name]]
    if (!is.null(idx_path)) {
      if (is.character(idx_path) && length(idx_path) == 1 && file.exists(idx_path)) {
        static_index_paths[[idx_name]] <- idx_path
        message("Using static index for: ", idx_name)
      } else {
        warning("Static index path for '", idx_name, "' (", idx_path,") is invalid or file not found. Ignoring.")
      }
    }
  }

  # Determine required inputs based on requested stats
  req_variable_data <- TRUE # Always needed
  req_inter_variable_data <- length(inter_stats) > 0 && !is.null(inter_variable_path)
  req_period_calculation_for_var <- any(c("max_period", "min_period") %in% stats) || 
                  (req_inter_variable_data && length(inter_stats) > 0) # inter stats need var periods
  req_period_calculation_for_inter_var <- req_inter_variable_data && length(inter_stats) > 0

  if (req_inter_variable_data && (is.null(inter_variable_path) || length(inter_variable_path) != n_units)) {
    stop(sprintf("'inter_variable_path' must be provided with %d layers if 'inter_stats' are requested.", n_units))
  }
  if (!req_inter_variable_data && length(inter_stats) > 0) {
    warning("'inter_stats' requested, but 'inter_variable_path' not provided or invalid. Interactive stats will not be calculated.")
    inter_stats <- c() # Clear inter_stats if no inter_variable_path
    req_period_calculation_for_inter_var <- FALSE
  }
  if (req_inter_variable_data && is.null(inter_variable_path)) {
    # This case should be caught above, but as a safeguard
    warning("'inter_variable_path' is NULL but interactive stats might be implied. Disabling interactive stats.")
    inter_stats <- c()
    req_inter_variable_data <- FALSE
    req_period_calculation_for_inter_var <- FALSE
  }


  # --- 1. Setup Input Paths ---
  paths <- c()
  path_variables_map <- list() # To map extracted columns to their roles

  names(variable_path) <- paste0("var_", seq_len(n_units))
  paths <- c(paths, variable_path)
  path_variables_map$variable <- names(variable_path)

  if (req_inter_variable_data && !is.null(inter_variable_path)) {
    names(inter_variable_path) <- paste0("intervar_", seq_len(n_units))
    paths <- c(paths, inter_variable_path)
    path_variables_map$inter_variable <- names(inter_variable_path)
  }

  static_paths_vec <- c()
  if (length(static_index_paths) > 0) {
    path_variables_map$static <- names(static_index_paths)
    static_paths_vec <- unlist(static_index_paths)
    # Prefix static index names for clarity in merged raster stack
    names(static_paths_vec) <- paste0("idx_", names(static_index_paths))
    paths <- c(paths, static_paths_vec)
  }

  if (length(paths) == 0) stop("No relevant input data paths identified.")

  # --- Check Geometry Consistency ---
  message("Checking geometry of input rasters...")
  ref_rast_geom_check <- NULL # Define outside tryCatch
  tryCatch({
    ref_rast_geom_check <- terra::rast(paths[1]) # Load first raster for reference
    ref_crs_check <- terra::crs(ref_rast_geom_check)
    ref_ext_check <- terra::ext(ref_rast_geom_check)
    if (length(paths) > 1) {
      for (i in 2:length(paths)) {
        current_rast_info <- try(terra::rast(paths[i]), silent = TRUE)
        if (inherits(current_rast_info, "try-error")) {
          stop("Failed to read header for: ", paths[i])
        }
        # Compare geometry with the reference
        terra::compareGeom(ref_rast_geom_check, current_rast_info, stopOnError = TRUE, messages = FALSE)
      }
    }
  }, error = function(e) {
    stop("Input rasters (including static indices) do not have the same geometry. ", e$message)
  })
  # Use the already loaded ref_rast_geom_check for template info later to avoid reloading
  # Store its CRS and Extent for use in region definition
  ref_crs <- terra::crs(ref_rast_geom_check) 
  ref_ext <- terra::ext(ref_rast_geom_check)

  # --- 2. Create Temporary Directory ---
  stats_qs_dir <- file.path(output_dir, paste0("stats_vars_qs_", basename(tempfile(pattern = ""))))
  if (!dir.exists(stats_qs_dir)) dir.create(stats_qs_dir, recursive = TRUE)

  # --- 3. Define Processing Region & Template Geometry ---
  base_map <- NULL
  if (!is.null(user_region)) {
    message("Using user-provided region.")
    if (inherits(user_region, "SpatVector")) base_map <- sf::st_as_sf(user_region)
    else if (inherits(user_region, "sf") || inherits(user_region, "sfc")) base_map <- sf::st_as_sf(sf::st_geometry(user_region))
    else stop("'user_region' must be an sf object or a terra SpatVector.")
    if (!all(sf::st_is_valid(base_map))) base_map <- sf::st_make_valid(base_map)
    if (sf::st_crs(base_map) != ref_crs) {
      base_map <- sf::st_transform(base_map, ref_crs)
    }
    region_bbox <- sf::st_bbox(base_map)
    raster_bbox_sf <- sf::st_bbox(c(xmin= ref_ext[1] |> as.numeric(), 
    ymin=ref_ext[3] |> as.numeric(), xmax=ref_ext[2] |> as.numeric(), ymax=ref_ext[4] |> as.numeric()), 
                                  crs = sf::st_crs(ref_crs))
    if (region_bbox$xmax < raster_bbox_sf$xmin || region_bbox$xmin > raster_bbox_sf$xmax ||
          region_bbox$ymax < raster_bbox_sf$ymin || region_bbox$ymin > raster_bbox_sf$ymax) {
      stop("Provided user_region does not overlap with the extent of the input rasters.")
    }
  } else {
    message("No user_region provided. Using the full extent of input rasters.")
    base_map <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(ref_ext), crs = ref_crs))
  }

  original_extent_vec <- as.vector(terra::ext(ref_rast_geom_check))
  original_dims_vec <- c(terra::nrow(ref_rast_geom_check), terra::ncol(ref_rast_geom_check))
  original_crs_txt <- terra::crs(ref_rast_geom_check, proj = TRUE)
  original_res_vec <- terra::res(ref_rast_geom_check) # Store original resolution
  original_ncol <- original_dims_vec[2]

  target_extent_vec <- original_extent_vec
  target_dims_vec <- original_dims_vec
  target_crs_txt <- original_crs_txt
  target_res_vec <- original_res_vec # Initialize target_res with original
  target_ncol <- original_ncol


  if (!is.null(user_region)) {
    target_template_rast <- NULL
    tryCatch({
      # Crop using the already loaded ref_rast_geom_check
      target_template_rast <- terra::crop(ref_rast_geom_check, terra::vect(base_map), mask = TRUE)
      if (terra::ncell(target_template_rast) > 0 && sum(terra::values(target_template_rast,na.rm=TRUE)) > 0 ) { # Check if not all NA
        target_extent_vec <- as.vector(terra::ext(target_template_rast))
        target_dims_vec <- c(terra::nrow(target_template_rast), terra::ncol(target_template_rast))
        target_res_vec <- terra::res(target_template_rast) # Get resolution of cropped raster
        target_ncol <- target_dims_vec[2]
      } else {
      warning("Cropping resulted in empty or all-NA raster, using original geometry.")
      }
    rm(target_template_rast)
    }, error = function(e) {
      warning("Could not derive exact cropped geometry, using full raster extent info. Error: ", e$message)
    })
  }

  template_info <- list(
    original_geom = list(
      extent = original_extent_vec,
      dimensions = original_dims_vec,
      crs = original_crs_txt,
      res = original_res_vec
    ),
    target_geom = list(
      extent = target_extent_vec,
      dimensions = target_dims_vec,
      crs = target_crs_txt,
      res = target_res_vec
    )
  )

  translate_cell_fun <- NULL # Initialize to NULL
  if (!is.null(user_region)) {
    # Use original raster (ref_rast_geom_check) for offset calculation relative to its top-left
    # Target coordinates are the top-left of the *target extent*
    # Ensure res from template_info$target_geom$res is used if different from original

    target_x_coord_for_offset <- template_info$target_geom$extent[1] + template_info$target_geom$res[1] / 2
    target_y_coord_for_offset <- template_info$target_geom$extent[4] - template_info$target_geom$res[2] / 2

    # colFromX/rowFromY expect coordinates of cell *centers*
    # The ref_rast_geom_check is the *original* full raster
    col_offset <- tryCatch(
      terra::colFromX(ref_rast_geom_check, target_x_coord_for_offset) - 1L,
      error = function(e) {
        warning("Could not precisely determine column offset. Assuming 0 offset. Error: ", e$message); 0L
    })
    row_offset <- tryCatch(
      terra::rowFromY(ref_rast_geom_check, target_y_coord_for_offset) - 1L,
      error = function(e) {
        warning("Could not precisely determine row offset. Assuming 0 offset. Error: ", e$message); 0L
    })

    col_offset <- max(0L, col_offset)
    row_offset <- max(0L, row_offset)

    translate_cell_fun <- define_translate(
      ncol_src = original_ncol,
      ncol_tgt = target_ncol,
      row_offset = row_offset,
      col_offset = col_offset
    )
  }

  template_info_file <- file.path(stats_qs_dir, "template_info.qs")
  tryCatch({
    qs::qsave(template_info, template_info_file)
  }, error = function(e) {
    stop("Failed to save template geometry information: ", e$message)
  })
  rm(ref_rast_geom_check) # Clean up the reference raster object

  # --- 4. Create Spatial Tiles ---
  sf::sf_use_s2(FALSE) # Recommended for st_intersection with planar data
  grid_bbox <- sf::st_bbox(base_map)
  # Ensure tile_degrees is positive
  if(tile_degrees <= 0) stop("'tile_degrees' must be positive.")
  rtt_grid <- sf::st_make_grid(grid_bbox, cellsize = tile_degrees,
            what = "polygons", crs = sf::st_crs(base_map))
  rtt <- sf::st_intersection(base_map, rtt_grid) # Intersect with the actual processing area
  rtt <- sf::st_collection_extract(rtt, "POLYGON") # Ensure only polygons
  rtt <- sf::st_as_sf(rtt) # Convert to sf object
  rtt <- rtt[!sf::st_is_empty(rtt), ] # Remove empty geometries
  ntiles <- nrow(rtt)
  if (ntiles == 0) stop("No overlapping tiles found for the processing area.")
  message("Processing area divided into ", ntiles, " tiles.")

  # --- 5. Define Intermediate File Paths (.qs) ---
  all_stat_names_to_calc <- c()
  if ("mean" %in% stats) all_stat_names_to_calc <- c(all_stat_names_to_calc, paste0(prefix_variable, "_mean"))
  if ("max" %in% stats) all_stat_names_to_calc <- c(all_stat_names_to_calc, paste0(prefix_variable, "_max"))
  if ("min" %in% stats) all_stat_names_to_calc <- c(all_stat_names_to_calc, paste0(prefix_variable, "_min"))
  if ("sum" %in% stats) all_stat_names_to_calc <- c(all_stat_names_to_calc, paste0(prefix_variable, "_sum"))
  if ("stdev" %in% stats) all_stat_names_to_calc <- c(all_stat_names_to_calc, paste0(prefix_variable, "_stdev"))
  if ("cv_cli" %in% stats) all_stat_names_to_calc <- c(all_stat_names_to_calc, paste0(prefix_variable, "_cv"))
  if ("max_period" %in% stats) all_stat_names_to_calc <- c(all_stat_names_to_calc, paste0(prefix_variable, "_max_period"))
  if ("min_period" %in% stats) all_stat_names_to_calc <- c(all_stat_names_to_calc, paste0(prefix_variable, "_min_period"))
  if ("max_inter" %in% inter_stats) all_stat_names_to_calc <- c(all_stat_names_to_calc, paste0(prefix_variable, "_", suffix_inter_max))
  if ("min_inter" %in% inter_stats) all_stat_names_to_calc <- c(all_stat_names_to_calc, paste0(prefix_variable, "_", suffix_inter_min))

  stats_output_qs_paths <- list()
  for (stat_name_full in all_stat_names_to_calc) {
    stats_output_qs_paths[[stat_name_full]] <- file.path(stats_qs_dir, paste0(stat_name_full, "_", seq_len(ntiles), ".qs"))
  }
    
  raw_paths_list <- list()

  # Enable the debug mode: Sys.setenv(BIOCLIM_DEBUG_RAW_VARS = "TRUE")
  write_raw_vars <- identical(toupper(Sys.getenv("BIOCLIM_DEBUG_RAW_VARS")), "TRUE")
  if (write_raw_vars) {
    message("DEBUG MODE: Writing raw variable tiles because BIOCLIM_DEBUG_RAW_VARS is set to TRUE.")
  }
  if (write_raw_vars) {
    raw_paths_list$variable <- file.path(stats_qs_dir, 
      paste0("raw_", prefix_variable, "_", seq_len(ntiles), ".qs"))
    if (req_inter_variable_data && !is.null(inter_variable_path)) {
      raw_paths_list$inter_variable <- file.path(stats_qs_dir, 
        paste0("raw_intervar_", prefix_variable, "_", seq_len(ntiles), ".qs"))
    }
  }

  # --- 6. Parallel Processing Loop ---
  p <- progressr::progressor(steps = ntiles)

  # Variables to export to parallel workers
  # Note: define_translate is a closure, so its environment (containing offsets, ncols) is captured.
  # The specific function `translate_cell_fun` generated by it needs to be passed if user_region is active.
  export_vars <- c("paths", "path_variables_map", "rtt", "ntiles", "stats", "inter_stats",
                  "n_units", "period", "period_stats", "circular",
                  "prefix_variable", "suffix_inter_max", "suffix_inter_min",
                  "req_inter_variable_data", "req_period_calculation_for_var", "req_period_calculation_for_inter_var",
                  "check_evar", "compute_periods", "var_periods",
                  "stats_output_qs_paths", "write_raw_vars", "stats_qs_dir", "user_region"
  )
  if (!is.null(translate_cell_fun)) export_vars <- c(export_vars, "translate_cell_fun")
  if (exists("raw_paths_list") && length(raw_paths_list) > 0) export_vars <- c(export_vars, "raw_paths_list")


  future.apply::future_lapply(seq_len(ntiles), function(tile_idx) {
  p(message = sprintf("Processing tile %d of %d", tile_idx, ntiles))

  tile_data_prepared <- list() # To store matrices like variable_matrix, inter_variable_matrix
  static_indices_for_tile <- list() # To store extracted static index values

  # Load all rasters for the current tile
  # `paths` contains all variable, inter_variable, and static_index paths
  raster_stack_for_tile <- tryCatch({ terra::rast(paths) }, error = function(e) { NULL })
  if (is.null(raster_stack_for_tile)) {
    warning(sprintf("Tile %d: Failed to load raster stack. Skipping.", tile_idx))
    return(NULL)
  }
  names(raster_stack_for_tile) <- names(paths) # Critical for exct_fun_combined

  current_tile_geom <- sf::st_geometry(rtt[tile_idx, ])

  # Define extraction function for this specific set of variables
  # This needs to be defined *inside* future_lapply or passed correctly with its environment
  exct_fun_combined_stats <- function(df, coverage_fractions) {
    # Determine a reference column from the primary variable to find non-NA cells
    # Use the first layer of the primary variable
    ref_col_name_for_na <- path_variables_map$variable[1]

    # If df is NULL or ref_col_name_for_na doesn't exist, return NULL early
    if (is.null(df) || !(ref_col_name_for_na %in% names(df))) return(NULL)

    nonaID <- which(!is.na(df[[ref_col_name_for_na]]))
    if (length(nonaID) == 0L) return(NULL)

    all_colnames_in_df <- names(df)

    # Extract cell IDs
    cell_id_col_name <- "cell" # Default from exact_extract
    if (! (cell_id_col_name %in% all_colnames_in_df) ) {
      warning("Tile ", tile_idx, ": 'cell' column not found in exact_extract output.")
      return(NULL)
    }
    cell_ids_raw_extracted <- df[[cell_id_col_name]][nonaID]

    # Translate cell IDs if a user_region was provided and cropping occurred
    cell_ids_final <- cell_ids_raw_extracted
    row_indices_to_keep <- seq_along(nonaID) # Keep all rows by default

    if (!is.null(user_region) && !is.null(translate_cell_fun)) {
      cell_ids_translated <- translate_cell_fun(cell_ids_raw_extracted)
      valid_translation_mask <- !is.na(cell_ids_translated)

      if (sum(valid_translation_mask) == 0) {
        warning("Tile ", tile_idx, ": No cells mapped to target region after translation.")
        return(NULL)
      }
      cell_ids_final <- cell_ids_translated[valid_translation_mask]
      row_indices_to_keep <- which(valid_translation_mask) # Indices relative to nonaID
    }

    # --- Extract Primary Variable Data ---
    var_cols_to_extract <- path_variables_map$variable
    # Ensure these columns actually exist in df
    var_cols_present <- intersect(var_cols_to_extract, all_colnames_in_df)
    if (length(var_cols_present) != length(var_cols_to_extract)){
      warning("Tile ", tile_idx, ": Not all expected primary variable columns found in extracted data.")
      # Decide if this is fatal or if we proceed with what's found
    }
    variable_df_subset <- df[nonaID[row_indices_to_keep], var_cols_present, drop = FALSE]
    variable_mat <- Rfast::data.frame.to_matrix(variable_df_subset)
    colnames(variable_mat) <- var_cols_present # Preserve original raster layer names

    # --- Extract Interactive Variable Data (if applicable) ---
    inter_variable_mat <- NULL
    if (req_inter_variable_data && !is.null(path_variables_map$inter_variable)) {
      inter_var_cols_to_extract <- path_variables_map$inter_variable
      inter_var_cols_present <- intersect(inter_var_cols_to_extract, all_colnames_in_df)
      if (length(inter_var_cols_present) == length(inter_var_cols_to_extract)) {
        inter_variable_df_subset <- df[nonaID[row_indices_to_keep], inter_var_cols_present, drop = FALSE]
        inter_variable_mat <- Rfast::data.frame.to_matrix(inter_variable_df_subset)
        colnames(inter_variable_mat) <- inter_var_cols_present
      } else {
      warning("Tile ", tile_idx, ": Not all expected interactive variable columns found.")
      }
    }

    # --- Extract Static Index Data (if applicable) ---
    static_idx_values_list <- list()
    if (!is.null(path_variables_map$static) && length(path_variables_map$static) > 0) {
      # Recall static paths were named "idx_originalname" when added to `paths`
      # So, `all_colnames_in_df` will contain these `idx_*` names
      # `path_variables_map$static` contains the original names ("max_unit", etc.)
      for (original_idx_name in path_variables_map$static) {
        # Construct the expected column name in df (e.g., "idx_max_unit")
        colname_in_df_for_static_idx <- paste0("idx_", original_idx_name)
        if (colname_in_df_for_static_idx %in% all_colnames_in_df) {
          static_idx_values_list[[original_idx_name]] <- df[[colname_in_df_for_static_idx]][nonaID[row_indices_to_keep]]
        } else {
        # This case should ideally not happen if paths were set up correctly
          warning("Tile ", tile_idx, ": Expected static index column '", colname_in_df_for_static_idx, "' not found.")
        }
      }
    }

    return(list(
      variable_matrix = variable_mat,
      inter_variable_matrix = inter_variable_mat,
      static_indices = static_idx_values_list,
      cell_ids = cell_ids_final))
  }

  # Perform Extraction
  extracted_data_list_raw <- tryCatch({
    exactextractr::exact_extract(
      raster_stack_for_tile, current_tile_geom,
      fun = exct_fun_combined_stats,
      include_cell = TRUE, # Cell IDs are crucial
      stack_apply = FALSE) # Process rasters together for exct_fun
    }, error = function(e) {
      warning(sprintf("Tile %d: exact_extract failed: %s", tile_idx, e$message))
      NULL
  })

  if (is.null(extracted_data_list_raw) || length(extracted_data_list_raw) == 0 || is.null(extracted_data_list_raw[[1]])) {
    warning(sprintf("Tile %d: No data extracted or extraction function returned NULL.", tile_idx))
    return(NULL)
  }
  # Assuming single polygon feature per tile after st_intersection and st_collection_extract
  extracted_data_for_tile <- extracted_data_list_raw[, 1]
  if (!is.list(extracted_data_for_tile) || is.null(extracted_data_for_tile$variable_matrix) || nrow(extracted_data_for_tile$variable_matrix) == 0) {
    warning(sprintf("Tile %d: Extracted data is invalid or empty.", tile_idx))
    return(NULL)
  }

  tile_data_prepared$variable_matrix <- check_evar(extracted_data_for_tile$variable_matrix)
  tile_data_prepared$cell_ids <- extracted_data_for_tile$cell_ids

  if (write_raw_vars && "variable" %in% names(raw_paths_list)) {
    rio::export(cbind(tile_data_prepared$variable_matrix, cell = tile_data_prepared$cell_ids),
                      raw_paths_list$variable[tile_idx])
  }

  if (req_inter_variable_data && !is.null(extracted_data_for_tile$inter_variable_matrix)) {
    tile_data_prepared$inter_variable_matrix <- check_evar(extracted_data_for_tile$inter_variable_matrix)
    if (write_raw_vars && "inter_variable" %in% names(raw_paths_list)) {
      rio::export(cbind(tile_data_prepared$inter_variable_matrix, cell = tile_data_prepared$cell_ids),
                        raw_paths_list$inter_variable[tile_idx])
    }
  } else if (req_inter_variable_data) {
    warning(sprintf("Tile %d: Interactive variable data expected but not found after extraction.", tile_idx))
    # Downgrade requirements if inter_variable matrix is missing
    # req_inter_variable_data_tile <- FALSE # This would need to be passed to stat calculation logic
    # For now, assume subsequent checks for matrix availability will handle it
  }

  static_indices_for_tile <- extracted_data_for_tile$static_indices

  # --- Period Calculations (if needed) ---
  defined_periods <- NULL
  if (req_period_calculation_for_var || req_period_calculation_for_inter_var) {
    defined_periods <- compute_periods(n_units = n_units, period_length = period_length, circular = circular)
  }

  if (req_period_calculation_for_var && !is.null(defined_periods)) {
    if (!is.null(tile_data_prepared$variable_matrix)) {
    tile_data_prepared$variable_period_summary <- var_periods(
      variable = tile_data_prepared$variable_matrix,
      periodos = defined_periods,
      stat = period_stats)
    } else {
      warning(sprintf("Tile %d: Variable matrix not available for period calculation.", tile_idx))
      req_period_calculation_for_var <- FALSE
    }
  }

  if (req_period_calculation_for_inter_var && !is.null(defined_periods)) {
    if (!is.null(tile_data_prepared$inter_variable_matrix)) {
      tile_data_prepared$inter_variable_period_summary <- var_periods(
      variable = tile_data_prepared$inter_variable_matrix,
      periodos = defined_periods,
      stat = period_stats)
    } else {
      warning(sprintf("Tile %d: Interactive variable matrix not available for period calculation.", tile_idx))
      req_period_calculation_for_inter_var <- FALSE
    }
  }

  # --- Calculate Requested Statistics ---
  # Combine stats and inter_stats for iteration, but handle their specific needs
  all_stats_to_attempt <- unique(c(stats, inter_stats))

  for (current_stat_type in all_stats_to_attempt) {
    result_value <- NULL
    output_stat_name_full <- NULL # To match keys in stats_output_qs_paths

    # Check prerequisites for the current stat type for *this tile*
    can_calculate_stat <- TRUE
    if (is.null(tile_data_prepared$variable_matrix)) {
      can_calculate_stat <- FALSE # Basic requirement
    } else if (current_stat_type %in% c("max_period", "min_period") && 
        (is.null(tile_data_prepared$variable_period_summary))) {
    can_calculate_stat <- FALSE
    } else if (current_stat_type %in% c("max_inter", "min_inter") &&
        (is.null(tile_data_prepared$variable_period_summary) || is.null(tile_data_prepared$inter_variable_period_summary))) {
    can_calculate_stat <- FALSE
    }

    if (!can_calculate_stat) {
      warning(sprintf("Tile %d: Skipping stat '%s' due to missing prerequisite data for this tile.", tile_idx, current_stat_type))
      next # Skip to the next statistic
    }

    # Calculate Value
    if (current_stat_type == "mean") {
      result_value <- Rfast::rowmeans(tile_data_prepared$variable_matrix)
      output_stat_name_full <- paste0(prefix_variable, "_mean")
    } else if (current_stat_type == "max") {
      idx_vec <- static_indices_for_tile$max_unit
      if (!is.null(idx_vec)) {
        # Ensure idx_vec is valid for matrix indexing
        idx_vec <- as.integer(idx_vec)
        valid_idx_mask <- !is.na(idx_vec) & idx_vec >= 1 & idx_vec <= ncol(tile_data_prepared$variable_matrix)
        result_value <- rep(NA_real_, nrow(tile_data_prepared$variable_matrix))
        if (any(valid_idx_mask)) {
          result_value[valid_idx_mask] <- tile_data_prepared$variable_matrix[cbind(which(valid_idx_mask), idx_vec[valid_idx_mask])]
        }
      } else {
        result_value <- Rfast::rowMaxs(tile_data_prepared$variable_matrix, value = TRUE)
      }
      output_stat_name_full <- paste0(prefix_variable, "_max")
    } else if (current_stat_type == "min") {
      idx_vec <- static_indices_for_tile$min_unit
      if (!is.null(idx_vec)) {
        idx_vec <- as.integer(idx_vec)
        valid_idx_mask <- !is.na(idx_vec) & idx_vec >= 1 & idx_vec <= ncol(tile_data_prepared$variable_matrix)
        result_value <- rep(NA_real_, nrow(tile_data_prepared$variable_matrix))
        if (any(valid_idx_mask)) {
          result_value[valid_idx_mask] <- tile_data_prepared$variable_matrix[cbind(which(valid_idx_mask), idx_vec[valid_idx_mask])]
        }
      } else {
        result_value <- Rfast::rowMins(tile_data_prepared$variable_matrix, value = TRUE)
      }
      output_stat_name_full <- paste0(prefix_variable, "_min")
    } else if (current_stat_type == "sum") {
      result_value <- Rfast::rowsums(tile_data_prepared$variable_matrix)
      output_stat_name_full <- paste0(prefix_variable, "_sum")
    } else if (current_stat_type == "stdev") {
      result_value <- Rfast::rowVars(tile_data_prepared$variable_matrix, std = TRUE)
      output_stat_name_full <- paste0(prefix_variable, "_stdev")
    } else if (current_stat_type == "cv_cli") {
      # (stdev / (1 + mean)) * 100 ; if mean is 0, CV is 0 or NA depending on stdev.
      # Let's make it 0 if mean is 0 to avoid NaN/Inf issues.
      means <- Rfast::rowmeans(tile_data_prepared$variable_matrix)
      stdevs <- Rfast::rowVars(tile_data_prepared$variable_matrix, std = TRUE)
      # Add a small epsilon if mean is exactly zero but stdev is not, or handle as 0.
      # The `1 + means` handles negative means okay if data can be negative.
      # If means can be 0, `stdevs / (means + (means==0 & stdevs !=0) * 1e-9)` is one way.
      # A simpler approach like bio15_fun:
      result_value <- ifelse( (1 + means) == 0, 0, (stdevs / (1 + means)) * 100)
      result_value[is.na(means) | is.na(stdevs)] <- NA # Ensure NAs propagate
      output_stat_name_full <- paste0(prefix_variable, "_cv")
    } else if (current_stat_type == "max_period") {
      # Use the period summary of the primary variable
      # `variable_period_summary` has columns: period_1, period_2, ..., min_idx, max_idx
      num_actual_period_cols <- ncol(tile_data_prepared$variable_period_summary) - 2
      idx_vec <- static_indices_for_tile$max_period %||% tile_data_prepared$variable_period_summary[, "max_idx"]

      idx_vec <- as.integer(idx_vec)
      valid_idx_mask <- !is.na(idx_vec) & idx_vec >= 1 & idx_vec <= num_actual_period_cols
      result_value <- rep(NA_real_, nrow(tile_data_prepared$variable_period_summary))
      if (any(valid_idx_mask)) {
        result_value[valid_idx_mask] <- tile_data_prepared$variable_period_summary[cbind(which(valid_idx_mask), idx_vec[valid_idx_mask])]
      }
      output_stat_name_full <- paste0(prefix_variable, "_max_period")
    } else if (current_stat_type == "min_period") {
      num_actual_period_cols <- ncol(tile_data_prepared$variable_period_summary) - 2
      idx_vec <- static_indices_for_tile$min_period %||% tile_data_prepared$variable_period_summary[, "min_idx"]

      idx_vec <- as.integer(idx_vec)
      valid_idx_mask <- !is.na(idx_vec) & idx_vec >= 1 & idx_vec <= num_actual_period_cols
      result_value <- rep(NA_real_, nrow(tile_data_prepared$variable_period_summary))
      if (any(valid_idx_mask)) {
      result_value[valid_idx_mask] <- tile_data_prepared$variable_period_summary[cbind(which(valid_idx_mask), idx_vec[valid_idx_mask])]
      }
      output_stat_name_full <- paste0(prefix_variable, "_min_period")
    } else if (current_stat_type == "max_inter") {
      # Value of *variable's* period, for the period where *inter_variable* is max
      num_actual_period_cols_var <- ncol(tile_data_prepared$variable_period_summary) - 2
      # Index comes from inter_variable_period_summary or static_interactive index
      idx_period_of_inter_max <- static_indices_for_tile$max_interactive %||% tile_data_prepared$inter_variable_period_summary[, "max_idx"]

      idx_period_of_inter_max <- as.integer(idx_period_of_inter_max)
      valid_idx_mask <- !is.na(idx_period_of_inter_max) & idx_period_of_inter_max >= 1 & idx_period_of_inter_max <= num_actual_period_cols_var
      result_value <- rep(NA_real_, nrow(tile_data_prepared$variable_period_summary))

      if (any(valid_idx_mask)) {
        result_value[valid_idx_mask] <- tile_data_prepared$variable_period_summary[cbind(which(valid_idx_mask), idx_period_of_inter_max[valid_idx_mask])]
      }
      output_stat_name_full <- paste0(prefix_variable, "_", suffix_inter_max)
    } else if (current_stat_type == "min_inter") {
      num_actual_period_cols_var <- ncol(tile_data_prepared$variable_period_summary) - 2
      idx_period_of_inter_min <- static_indices_for_tile$min_interactive %||% tile_data_prepared$inter_variable_period_summary[, "min_idx"]

      idx_period_of_inter_min <- as.integer(idx_period_of_inter_min)
      valid_idx_mask <- !is.na(idx_period_of_inter_min) & idx_period_of_inter_min >= 1 & idx_period_of_inter_min <= num_actual_period_cols_var
      result_value <- rep(NA_real_, nrow(tile_data_prepared$variable_period_summary))

      if (any(valid_idx_mask)) {
        result_value[valid_idx_mask] <- tile_data_prepared$variable_period_summary[cbind(which(valid_idx_mask), idx_period_of_inter_min[valid_idx_mask])]
      }
      output_stat_name_full <- paste0(prefix_variable, "_", suffix_inter_min)
    }

    # Save the calculated statistic if successful
    if (!is.null(result_value) && !is.null(output_stat_name_full) && output_stat_name_full %in% names(stats_output_qs_paths)) {
      rio::export(cbind(value = result_value, cell = tile_data_prepared$cell_ids),
      stats_output_qs_paths[[output_stat_name_full]][tile_idx])
    } else if (is.null(output_stat_name_full) && current_stat_type %in% all_stat_names_to_calc) {
      # This means a stat was requested, loop entered, but not handled by if/else if. Should not happen.
      warning(sprintf("Tile %d: Stat '%s' was attempted but no calculation logic matched or output name not resolved.", tile_idx, current_stat_type))
    } else if (!is.null(output_stat_name_full) && !(output_stat_name_full %in% names(stats_output_qs_paths)) ) {
      warning(sprintf("Tile %d: Output name '%s' for stat '%s' does not match any predefined Q_S path.", tile_idx, output_stat_name_full, current_stat_type))
    }
  }

  # --- Cleanup for the tile ---
    rm(tile_data_prepared, 
       static_indices_for_tile, 
       extracted_data_for_tile, 
       raster_stack_for_tile, 
       current_tile_geom, 
       extracted_data_list_raw)
    gc() # Garbage collect
    return(NULL) 
  }, 
  future.seed = TRUE,
  future.globals = structure(export_vars, names = export_vars),
  future.packages = c("sf", "terra", "exactextractr", "Rfast", "rio", "purrr")
  )

  message("Tiled computation finished.")
  return(stats_qs_dir)
}