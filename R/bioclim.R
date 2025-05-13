#' Compute Specified Bioclimatic Variables
#'
#' Calculates specified bioclimatic variables (1-35) from climate rasters
#' with user-defined temporal units, optionally using static indices. AOI is defined
#' by `user_region` or defaults to the full raster extent. Uses parallel processing,
#' spatial tiling, and `exactextractr`.
#'
#' @param bios Numeric vector specifying which bioclimatic variables (1-35) to compute.
#' @param n_units Integer. The number of temporal units (layers) per input variable.
#' @param tmin_path Character vector of `n_units` paths to minimum temperature rasters.
#' @param tmax_path Character vector of `n_units` paths to maximum temperature rasters.
#' @param prec_path Character vector of `n_units` paths to precipitation rasters.
#' @param tavg_path Character vector of `n_units` paths to optional average temperature rasters.
#' @param srad_path Character vector of `n_units` paths to solar radiation rasters.
#' @param mois_path Character vector of `n_units` paths to moisture rasters.
#' @param period_length Integer. The number of units defining a "period". Default: 3.
#' @param circular Logical. Calculate periods wrapping around the cycle? Default: TRUE.
#' @param user_region Optional. An `sf` or `terra::SpatVector` object defining the
#'   processing area. If `NULL`, the full extent of the input rasters is used.
#' @param tile_degrees Numeric. Approximate size of processing tiles in degrees. Default: 5.
#' @param temp_dir Character. Path for temporary tile files. Default: `tempdir()`.
#' @param write_raw_vars Logical. Save intermediate extracted climate data. Default: FALSE.
#' @param ... Additional arguments, including optional paths to static index rasters
#'   (e.g., `warmest_period_path`).
#'
#' @return Character string: Path to the temporary directory containing intermediate `.qs` files.
#'
#' @details Calculates BIOs 1-35. The Area of Interest (AOI) is determined by `user_region`.
#'   If `user_region` is NULL, the AOI is the full spatial extent of the input climate rasters.
#'   Static indices provided via `...` override dynamic calculations.
#'   Uses `exactextractr` for extraction and `Rfast` for matrix calculations.
#'   Assemble results using a compatible `write_layers` function.
#' @export
bioclim_vars <- function(bios,
  n_units,
  tmin_path = NULL,
  tmax_path = NULL,
  prec_path = NULL,
  tavg_path = NULL,
  srad_path = NULL,
  mois_path = NULL,
  period_length = 3,
  circular = TRUE,
  user_region = NULL,
  tile_degrees = 5,
  temp_dir = tempdir(),
  write_raw_vars = FALSE,
  ...) {
    
    # --- 0. Input Validation, Dependency Mapping, Static Index Parsing ---
    dot_args <- list(...)
    static_index_paths <- list()
    valid_static_indices <- c("warmest_unit", "coldest_unit", "wettest_unit", "driest_unit",
                              "high_rad_unit", "low_rad_unit", "high_mois_unit", "low_mois_unit", 
                              "warmest_period", "coldest_period", "wettest_period", "driest_period",
                              "high_mois_period", "low_mois_period")
    for (arg_name in names(dot_args)) {
      path_suffix <- "_path"
      if (endsWith(arg_name, path_suffix)) { 
        index_name <- sub(path_suffix, "", arg_name)
        if (index_name %in% valid_static_indices) {
          if(is.character(dot_args[[arg_name]]) && length(dot_args[[arg_name]]) == 1 && file.exists(dot_args[[arg_name]])) { 
            static_index_paths[[index_name]] <- dot_args[[arg_name]]
            message("Using static index for: ", index_name)
          } else { 
            warning("Static index path '", arg_name, "' invalid/not found. Ignoring.")
          }
        }
      }
    }
    if (missing(n_units) || !is.numeric(n_units) || n_units <= 0) stop("'n_units' required.")
    if (!is.numeric(period_length) || period_length <= 0 || period_length > n_units) stop("'period_length' invalid.")
    if (!is.numeric(bios) || any(bios < 1) || any(bios > 35)) stop("'bios' must be 1-35.")
    bios <- sort(unique(bios))
    needs <- list(tmin = c(2, 3, 6, 7), 
                  tmax = c(2, 3, 5, 7), 
                  tavg = c(1, 4, 5, 6, 8, 9, 10, 11, 18, 19, 26, 27, 34, 35), 
                  prec = c(8, 9, 12, 13, 14, 15, 16, 17, 18, 19, 24, 25),
                  srad = c(20, 21, 22, 23, 24, 25, 26, 27),
                  mois = c(28, 29, 30, 31, 32, 33, 34, 35),
                  bio02 = c(3), 
                  bio05 = c(7, 3), 
                  bio06 = c(7, 3), 
                  bio07 = c(3), 
                  bio12 = c(15), 
                  temp_p = c(8, 9, 10, 11, 18, 19, 26, 27, 34, 35), 
                  prec_p = c(8, 9, 16, 17, 18, 19, 24, 25),
                  srad_p = c(24, 25, 26, 27),
                  mois_p = c(32, 33, 34, 35))
    bios_to_calculate <- bios
    check_again <- TRUE
    while(check_again) { 
      check_again <- FALSE
      if (3 %in% bios_to_calculate && !(7 %in% bios_to_calculate)) {
        bios_to_calculate <- c(bios_to_calculate, 7)
        check_again <- TRUE
      }
      if (3 %in% bios_to_calculate && !(2 %in% bios_to_calculate)) {
        bios_to_calculate <- c(bios_to_calculate, 2)
        check_again <- TRUE
      }
      if (7 %in% bios_to_calculate && !(5 %in% bios_to_calculate)) {
        bios_to_calculate <- c(bios_to_calculate, 5)
        check_again <- TRUE
      }
      if(7 %in% bios_to_calculate && !(6 %in% bios_to_calculate)) {
        bios_to_calculate <- c(bios_to_calculate, 6)
        check_again <- TRUE}
      if(15 %in% bios_to_calculate && !(12 %in% bios_to_calculate)) {
        bios_to_calculate <- c(bios_to_calculate, 12)
        check_again <- TRUE
      }
      bios_to_calculate <- sort(unique(bios_to_calculate))
    }
    req_tmin_direct <- any(needs$tmin %in% bios_to_calculate)
    req_tmax_direct <- any(needs$tmax %in% bios_to_calculate)
    req_prec_direct <- any(needs$prec %in% bios_to_calculate)
    req_srad_direct <- any(needs$srad %in% bios_to_calculate)
    req_mois_direct <- any(needs$mois %in% bios_to_calculate)
    req_tavg_value <- any(needs$tavg %in% bios_to_calculate) || any(needs$temp_p %in% bios_to_calculate)
    req_tavg_load <- req_tavg_value && !is.null(tavg_path)
    req_tavg_calc <- req_tavg_value && is.null(tavg_path)
    req_tmin_path <- req_tmin_direct || req_tavg_calc
    req_tmax_path <- req_tmax_direct || req_tavg_calc
    req_prec_path <- req_prec_direct || any(needs$prec_p %in% bios_to_calculate)
    req_srad_path <- req_srad_direct || any(needs$srad_p %in% bios_to_calculate)
    req_mois_path <- req_mois_direct || any(needs$mois_p %in% bios_to_calculate)
    if (req_tmin_path && is.null(tmin_path)) stop("tmin_path required.")
    if (req_tmax_path && is.null(tmax_path)) stop("tmax_path required.")
    if (req_prec_path && is.null(prec_path)) stop("prec_path required.")
    if (req_srad_path && is.null(srad_path)) stop("srad_path required.")
    if (req_mois_path && is.null(mois_path)) stop("mois_path required.")  
    if (req_tmin_path && length(tmin_path) != n_units) stop(sprintf("tmin_path length error (need %d).", n_units))
    if (req_tmax_path && length(tmax_path) != n_units) stop(sprintf("tmax_path length error (need %d).", n_units))
    if (req_prec_path && length(prec_path) != n_units) stop(sprintf("prec_path length error (need %d).", n_units))
    if (req_tavg_load && length(tavg_path) != n_units) stop(sprintf("tavg_path length error (need %d).", n_units))
    if (req_srad_path && length(srad_path) != n_units) stop(sprintf("srad_path length error (need %d).", n_units))
    if (req_mois_path && length(mois_path) != n_units) stop(sprintf("mois_path length error (need %d).", n_units))
    # --- 1. Setup Input Paths ---
    paths <- c()
    path_variables <- list()
    climate_paths <- c()
    if (req_tmin_path) {
      names(tmin_path) <- paste0("tmin_", seq_len(n_units))
      climate_paths <- c(climate_paths, tmin_path)
    }
    if (req_tmax_path) {
      names(tmax_path) <- paste0("tmax_", seq_len(n_units))
      climate_paths <- c(climate_paths, tmax_path)
    }
    if (req_prec_path) {
      names(prec_path) <- paste0("prec_", seq_len(n_units))
      climate_paths <- c(climate_paths, prec_path)
    }
    if (req_tavg_load) {
      names(tavg_path) <- paste0("tavg_", seq_len(n_units))
      climate_paths <- c(climate_paths, tavg_path)
    }
    if (req_srad_path) {
      names(srad_path) <- paste0("srad_", seq_len(n_units))
      climate_paths <- c(climate_paths, srad_path)
    }
    if (req_mois_path) {
      names(mois_path) <- paste0("mois_", seq_len(n_units))
      climate_paths <- c(climate_paths, mois_path)
    }
    paths <- climate_paths
    path_variables$climate <- names(paths)
    static_paths_vec <- c()
    if (length(static_index_paths) > 0) {
      path_variables$static <- names(static_index_paths)
      static_paths_vec <- unlist(static_index_paths)
      names(static_paths_vec) <- paste0("idx_", names(static_index_paths))
      paths <- c(paths, static_paths_vec)
    }
    if (length(paths) == 0) stop("No relevant input data paths identified.")
    
    # --- Check Geometry Consistency ---
    message("Checking geometry of input rasters...")
    tryCatch({
      ref_rast <- terra::rast(paths[1])
      ref_crs <- terra::crs(ref_rast)
      ref_ext <- terra::ext(ref_rast)
      if (length(paths) > 1) {
        for(i in 2:length(paths)) {
          current_rast_info <- try(terra::rast(paths[i]), silent = TRUE)
          if (inherits(current_rast_info, "try-error")) {
            stop("Failed to read header for: ", paths[i])
          }
          terra::compareGeom(ref_rast, current_rast_info, stopOnError = TRUE, messages = FALSE)
        }
      }
      message("Input raster geometries appear consistent.")
    }, error = function(e) {
      stop("Input rasters (including static indices) do not have the same geometry. ", e$message)
    })
    rm(ref_rast)
  
    
    
    # --- 2. Create Temporary Directory ---
    bioclima_dir <- file.path(temp_dir, paste0("bioclima_qs_", basename(tempfile(pattern = ""))))
    if (!dir.exists(bioclima_dir)) dir.create(bioclima_dir, recursive = TRUE)
    message("Intermediate .qs files will be stored in: ", bioclima_dir)
    
    # --- 3. Define Processing Region ---
    base_map <- NULL
    if (!is.null(user_region)) {
      message("Using user-provided region.")
      if (inherits(user_region, "SpatVector")) base_map <- sf::st_as_sf(user_region)
      else if (inherits(user_region, "sf") || inherits(user_region, "sfc")) base_map <- sf::st_as_sf(sf::st_geometry(user_region))
      else stop("'user_region' must be an sf object or a terra SpatVector.")
      if (!all(sf::st_is_valid(base_map))) base_map <- sf::st_make_valid(base_map)
      # Transform user_region CRS to match raster CRS if needed
      if(sf::st_crs(base_map) != ref_crs) {
        message("Transforming user_region CRS to match raster CRS.")
        base_map <- sf::st_transform(base_map, ref_crs)
      }
      # Check if user region overlaps with raster extent
      region_bbox <- sf::st_bbox(base_map)
      raster_bbox <- c(xmin = ref_ext[1], ymin = ref_ext[3], xmax = ref_ext[2], ymax = ref_ext[4])
      # Simple bounding box overlap check
      if(region_bbox$xmax < raster_bbox["xmin.xmin"] || region_bbox$xmin > raster_bbox["xmax.xmax"] ||
        region_bbox$ymax < raster_bbox["ymin.ymin"] || region_bbox$ymin > raster_bbox["ymax.ymax"]) {
          stop("Provided user_region does not overlap with the extent of the input rasters.")
        }
        
      } else {
        message("No user_region provided. Using the full extent of input rasters.")
        # Create an sf polygon from the raster extent obtained earlier
        base_map <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(ref_ext), crs = ref_crs))
      }
        
      # --- Extract and Save Template Geometry Info ---
      message("Extracting template geometry information...")
      ref_rast_geom <- terra::rast(paths[1])
      original_extent_vec <- as.vector(terra::ext(ref_rast_geom))
      original_dims_vec <- c(terra::nrow(ref_rast_geom), terra::ncol(ref_rast_geom))
      original_crs_txt <- terra::crs(ref_rast_geom, proj = TRUE)
      original_ncol <- original_dims_vec[2]

      # Determine target geometry
      target_extent_vec <- original_extent_vec
      target_dims_vec <- original_dims_vec
      target_crs_txt <- original_crs_txt
      target_ncol <- original_ncol

      if (!is.null(user_region)) {
          message("Deriving target geometry from user_region.")
          target_template_rast <- NULL
          tryCatch({
              target_template_rast <- terra::crop(ref_rast_geom, terra::vect(base_map), mask = TRUE)
              if(terra::ncell(target_template_rast) > 0){
                  target_extent_vec <- as.vector(terra::ext(target_template_rast))
                  target_dims_vec <- c(terra::nrow(target_template_rast), terra::ncol(target_template_rast))
                  target_ncol <- target_dims_vec[2]
                  target_res <- terra::res(target_template_rast)
                  message("Target template geometry based on cropped reference raster.")
              } else {
                  warning("Cropping resulted in empty raster, using original geometry.")
              }
              rm(target_template_rast)
          }, error = function(e){
              warning("Could not derive exact cropped geometry, using full raster extent info. Error: ", e$message)
          })
      } else {
          message("Target template geometry based on full input raster extent.")
      }


      template_info <- list(
          original_geom = list(
              extent = original_extent_vec,
              dimensions = original_dims_vec,
              crs = original_crs_txt
          ),
          target_geom = list(
              extent = target_extent_vec,
              dimensions = target_dims_vec,
              crs = target_crs_txt,
              res = target_res
          )
      )
  
      # --- Calculate Translation Parameters ---
      if (!is.null(user_region)) {
        message("Calculating cell ID translation parameters...")
      # Use the original raster *header* info (ref_rast_geom) for offset calculation
      target_xmin <- template_info$target_geom$extent[1] + template_info$target_geom$res[1] / 2
      target_ymax <- template_info$target_geom$extent[4] - template_info$target_geom$res[2] / 2
      # Handle potential floating point inaccuracies by finding the *nearest* cell
      col_offset <- tryCatch(
          terra::colFromX(ref_rast_geom, target_xmin) - 1L,
          error = function(e) {
              warning("Could not precisely determine column offset using terra::colFromX, possibly due to edge alignment. Assuming 0 offset. Error: ", e$message)
              0L # Default to 0 if calculation fails
          }
      )
      row_offset <- tryCatch(
          terra::rowFromY(ref_rast_geom, target_ymax) - 1L,
          error = function(e) {
              warning("Could not precisely determine row offset using terra::rowFromY, possibly due to edge alignment. Assuming 0 offset. Error: ", e$message)
              0L # Default to 0 if calculation fails
          }
      )

      # Ensure offsets are non-negative (should be if target is within original)
      col_offset <- max(0L, col_offset)
      row_offset <- max(0L, row_offset)

      message(sprintf("Translation params: ncol_src=%d, ncol_tgt=%d, row_offset=%d, col_offset=%d",
                      original_ncol, target_ncol, row_offset, col_offset))

      # --- Create the specific translation function ---
      translate_cell <- define_translate(
          ncol_src = original_ncol,
          ncol_tgt = target_ncol,
          row_offset = row_offset,
          col_offset = col_offset
      )
      }
  
      # Save the template info within that directory
      template_info_file <- file.path(bioclima_dir, "template_info.qs")
      tryCatch({
          qs::qsave(template_info, template_info_file)
          message("Template geometry information saved to: ", template_info_file)
      }, error = function(e){
          stop("Failed to save template geometry information: ", e$message)
      })
      rm(ref_rast_geom)
      # --- 4. Create Spatial Tiles ---
      sf::sf_use_s2(FALSE)
      grid_bbox <- sf::st_bbox(base_map)
      rtt_grid <- sf::st_make_grid(grid_bbox, cellsize = tile_degrees, 
                                   what = "polygons", crs = sf::st_crs(base_map))
      rtt <- sf::st_intersection(base_map, rtt_grid)
      rtt <- sf::st_collection_extract(rtt, c("POLYGON"))
      rtt <- sf::st_as_sf(rtt)
      rtt <- rtt[!sf::st_is_empty(rtt),]
      ntiles <- nrow(rtt)
      if (ntiles == 0) stop("No overlapping tiles found for the processing area.")
      message("Rasters divided into ", ntiles, " tiles for processing.")
      
      
      # --- 5. Define Intermediate File Paths (.qs) ---
      bios_qs_paths <- list()
      for (bio_num in bios) { 
        bio_name <- paste0("bio", sprintf("%02d", bio_num))
        bios_qs_paths[[bio_name]] <- file.path(bioclima_dir, paste0(bio_name, "_", seq_len(ntiles), ".qs"))
      }
      if (write_raw_vars) {
        raw_paths_list <- list()
        if (req_tmin_path) raw_paths_list$tmin <- file.path(bioclima_dir, paste0("raw_tmin_", seq_len(ntiles), ".qs"))
        if (req_tmax_path) raw_paths_list$tmax <- file.path(bioclima_dir, paste0("raw_tmax_", seq_len(ntiles), ".qs"))
        if (req_prec_path) raw_paths_list$prec <- file.path(bioclima_dir, paste0("raw_prec_", seq_len(ntiles), ".qs"))
        if (req_tavg_value) raw_paths_list$tavg <- file.path(bioclima_dir, paste0("raw_tavg_", seq_len(ntiles), ".qs"))
        if (req_srad_path) raw_paths_list$srad <- file.path(bioclima_dir, paste0("raw_srad_", seq_len(ntiles), ".qs"))
        if (req_mois_path) raw_paths_list$mois <- file.path(bioclima_dir, paste0("raw_mois_", seq_len(ntiles), ".qs"))
      }
      # --- 6. Parallel Processing Loop ---
      # progressr::with_progress({
      p <- progressr::progressor(steps = ntiles)
      worker_req_tavg <- any(needs$tavg %in% bios_to_calculate) || any(needs$temp_p %in% bios_to_calculate)
      worker_req_temp_p <- any(needs$temp_p %in% bios_to_calculate)
      worker_req_prec_p <- any(needs$prec_p %in% bios_to_calculate)
      worker_req_srad_p <- any(needs$srad_p %in% bios_to_calculate)
      worker_req_mois_p <- any(needs$mois_p %in% bios_to_calculate)
      export_vars <- c("paths", "path_variables", "rtt", "ntiles", "bios", "bios_to_calculate", 
                        "n_units", "period_length", "circular", "req_tmin_path", "req_tmax_path", 
                        "req_prec_path", "req_tavg_load", "req_tavg_calc", "req_srad_path", "req_mois_path", 
                        "worker_req_tavg", "worker_req_temp_p", "worker_req_prec_p", 
                        "worker_req_srad_p", "worker_req_mois_p", 
                        "check_evar", "compute_periods", "var_periods", 
                        paste0("bio", 1:19, "_fun"), "bios_qs_paths", "write_raw_vars", "bioclima_dir")
      if (exists("raw_paths_list")) export_vars <- c(export_vars, "raw_paths_list")
      vals <- future.apply::future_lapply(seq_len(ntiles), function(x) {
        p(message = sprintf("Processing tile %d of %d", x, ntiles))
        # p(message = "HOLA")
        tile_results <- list()
        static_indices_tile <- list()
        evars_stack_tile <- tryCatch({ terra::rast(paths) }, error = function(e) { NULL })
        if (is.null(evars_stack_tile)) { 
          warning(sprintf("Tile %d: Failed load",x))
          return(NULL)
        }
        names(evars_stack_tile) <- names(paths)
        rt0_geom <- sf::st_geometry(rtt[x, ])

        exct_fun_combined <- function(df, coverage_fractions) {
          ref_col_name <- path_variables$climate[1]
          nonaID <- which(!is.na(df[[ref_col_name]]))
          if (length(nonaID) == 0L) return(NULL)
          all_colnames <- names(df)
          climate_colnames <- path_variables$climate
          static_colnames_expected <- names(paths)[grepl("^idx_", names(paths))]
          static_colnames_present <- intersect(static_colnames_expected, all_colnames)
          cell_id_col_name <- "cell"
          cell_ids_extracted <- df[[cell_id_col_name]][nonaID]
          if (!is.null(user_region)) {
            cell_ids_extracted <- translate_cell(cell_ids_extracted)
            noMap <- !is.na(cell_ids_extracted)
            cell_ids_extracted <- cell_ids_extracted[noMap]
            climate_df_subset <- df[nonaID[noMap], climate_colnames, drop = FALSE]
          } else {
            climate_df_subset <- df[nonaID, climate_colnames, drop = FALSE]
          }
          climate_mat <- Rfast::data.frame.to_matrix(climate_df_subset)
          colnames(climate_mat) <- climate_colnames
          static_idx_list <- list()
          if (length(static_colnames_present) > 0) {
            for (sc_name in static_colnames_present) { 
              original_idx_name <- sub("^idx_", "", sc_name)
              if (!is.null(user_region)) {
                static_idx_list[[original_idx_name]] <- df[[sc_name]][nonaID[noMap]]
              } else {
                static_idx_list[[original_idx_name]] <- df[[sc_name]][nonaID]
              }
            }
          }
          return(list(climate_matrix = climate_mat, 
                      static_indices = static_idx_list,
                      cell_ids = cell_ids_extracted))
        }
        
        # Perform Extraction
        extracted_data_list <- tryCatch({ 
          exactextractr::exact_extract(evars_stack_tile, rt0_geom, fun = exct_fun_combined, include_cell = TRUE) 
        }, error = function(e){ NULL })
        if (is.null(extracted_data_list) || length(extracted_data_list) == 0 || is.null(extracted_data_list[[1]])) return(NULL)
        extracted_data <- extracted_data_list[, 1]
        if (!is.list(extracted_data) || is.null(extracted_data$climate_matrix) || nrow(extracted_data$climate_matrix) == 0) return(NULL)
        climate_matrix <- extracted_data$climate_matrix
        static_indices_tile <- extracted_data$static_indices
        cell_ids <- extracted_data$cell_ids
        
        # Conditional Data Prep
        if (req_tmin_path) {
          tile_results$tmin <- check_evar(climate_matrix[, grep(pattern = "^tmin_", path_variables$climate), drop = FALSE])
          if (write_raw_vars && "tmin" %in% names(raw_paths_list)) {
            rio::export(cbind(tile_results$tmin, cell = cell_ids), raw_paths_list$tmin[x])
          }
        }
        if (req_tmax_path) {
          tile_results$tmax <- check_evar(climate_matrix[, grep(pattern = "^tmax_", path_variables$climate), drop = FALSE])
          if (write_raw_vars && "tmax" %in% names(raw_paths_list)) {
            rio::export(cbind(tile_results$tmax, cell = cell_ids), raw_paths_list$tmax[x])
          }
        }
        if (req_prec_path) {
          tile_results$prec_vals <- check_evar(climate_matrix[, grep(pattern = "^prec_", path_variables$climate), drop = FALSE])
          if (write_raw_vars && "prec" %in% names(raw_paths_list)) {
            rio::export(cbind(tile_results$prec_vals, cell = cell_ids), raw_paths_list$prec[x])
          }
        }
        if (req_srad_path) {
          tile_results$srad_vals <- check_evar(climate_matrix[, grep(pattern = "^srad_", path_variables$climate), drop = FALSE])
          if (write_raw_vars && "srad" %in% names(raw_paths_list)) {
            rio::export(cbind(tile_results$srad_vals, cell = cell_ids), raw_paths_list$srad[x])
          }
        }
        if (req_mois_path) {
          tile_results$mois_vals <- check_evar(climate_matrix[, grep(pattern = "^mois_", path_variables$climate), drop = FALSE])
          if (write_raw_vars && "mois" %in% names(raw_paths_list)) {
            rio::export(cbind(tile_results$mois_vals, cell = cell_ids), raw_paths_list$mois[x])
          }
        }
        tavg_available <- FALSE
        if (worker_req_tavg) { 
          if (req_tavg_load) {
            tavg_cols <- grep(pattern = "^tavg_", path_variables$climate)
            if (all(tavg_cols %in% colnames(climate_matrix))) {
              tile_results$temperature_avg <- check_evar(climate_matrix[, tavg_cols, drop = FALSE])
              tavg_available < TRUE
              if (write_raw_vars && "tavg" %in% names(raw_paths_list)) {
                rio::export(cbind(tile_results$temperature_avg, cell = cell_ids), raw_paths_list$tavg[x])
              }
            } else { 
              warning("Tile ", x, ": Failed Tavg load.")
            }
          } else if (req_tavg_calc) { 
            if (!is.null(tile_results$tmin) && !is.null(tile_results$tmax)) { 
              tile_results$temperature_avg <- (tile_results$tmax + tile_results$tmin) / 2
              tavg_available <- TRUE
              if (is.vector(tile_results$temperature_avg)) {
                tile_results$temperature_avg <- t(tile_results$temperature_avg)
                colnames(tile_results$temperature_avg) <- paste0("tavg_",seq_len(n_units))
              }
              if (write_raw_vars && ("tavg" %in% names(raw_paths_list))) {
                rio::export(cbind(tile_results$temperature_avg, cell = cell_ids),raw_paths_list$tavg[x])
                
              }
            } else { 
              warning("Tile ", x, ": Cannot calc Tavg.") 
            } 
          } 
        }
        temp_p_available <- FALSE
        prec_p_available <- FALSE
        srad_p_available <- FALSE
        mois_p_available <- FALSE
        if (worker_req_temp_p || worker_req_prec_p || worker_req_srad_p || worker_req_mois_p) { 
          periodos <- compute_periods(n_units = n_units, period_length = period_length, circular = circular)
          if (worker_req_temp_p && tavg_available) { 
            tile_results$tempr_periods <- var_periods(variable = tile_results$temperature_avg, 
                                                      periodos = periodos,
                                                      n_units = n_units,
                                                      period_length = period_length,
                                                      stat = "mean")
            temp_p_available <- TRUE
          } else if (worker_req_temp_p) {
            warning("Tile ", x, ": Cannot calc Temp Periods.")
          } 
          if (worker_req_prec_p && !is.null(tile_results$prec_vals)) {
            tile_results$preci_periods <- var_periods(variable = tile_results$prec_vals, 
                                                      periodos = periodos, 
                                                      n_units = n_units, 
                                                      period_length = period_length,
                                                      stat = "sum")
            prec_p_available <- TRUE
          } else if (worker_req_prec_p) {
            warning("Tile ", x, ": Cannot calc Prec Periods.")
          }
          if (worker_req_srad_p && !is.null(tile_results$srad_vals)) {
            tile_results$srad_periods <- var_periods(variable = tile_results$srad_vals, 
                                                      periodos = periodos, 
                                                      n_units = n_units, 
                                                      period_length = period_length,
                                                      stat = "mean")
            srad_p_available <- TRUE
          } else if (worker_req_srad_p) {
            warning("Tile ", x, ": Cannot calc Srad Periods.")
          }
          if (worker_req_mois_p && !is.null(tile_results$mois_vals)) {
            tile_results$mois_periods <- var_periods(variable = tile_results$mois_vals, 
                                                      periodos = periodos, 
                                                      n_units = n_units, 
                                                      period_length = period_length,
                                                      stat = "mean")
            mois_p_available <- TRUE
          } else if (worker_req_mois_p) {
            warning("Tile ", x, ": Cannot calc Mois Periods.")
          }
        }
        
        # Calculate BIOs
        calculated_bios_in_tile <- list()
        # Order bios
        bios_order <- c(1, 2, 4, 5, 6, 7, 3, 12, 13, 14, 15, 16, 17, 10, 11, 20, 21, 22, 23, 28, 29, 30, 31, 32, 33, 
                        8, 9, 18, 19, 24, 25, 26, 27, 34, 35)
        bios_to_calculate <- bios_order[bios_order %in% bios_to_calculate]
        for (bio_num in bios_to_calculate) {
          can_calculate <- TRUE
          if (bio_num %in% c(1, 4) && !tavg_available) {can_calculate <- FALSE}
          if (bio_num %in% c(5, 6) && 
            is.null(static_indices_tile[[paste0(ifelse(bio_num == 5, "warmest", "coldest"), "_unit")]]) && 
            !tavg_available) {
            can_calculate <- FALSE
          }
          if (bio_num == 5 && is.null(static_indices_tile[["warmest_unit"]]) && is.null(tile_results$tmax)) {can_calculate <- FALSE}
          if (bio_num == 6 && is.null(static_indices_tile[["coldest_unit"]]) && is.null(tile_results$tmin)) {can_calculate <- FALSE}
          if (bio_num %in% c(2, 3, 7) && (is.null(tile_results$tmin) || is.null(tile_results$tmax))) {can_calculate <- FALSE}
          if (bio_num %in% c(12, 13, 14) && is.null(tile_results$prec_vals)) {can_calculate <- FALSE}
          if (bio_num %in% c(20, 21, 22, 23) && is.null(tile_results$srad_vals)) {can_calculate <- FALSE}
          if (bio_num %in% c(28, 29, 30, 31) && is.null(tile_results$mois_vals)) {can_calculate <- FALSE}
          if (bio_num %in% c(8, 9, 10, 11) && (!temp_p_available)) {can_calculate <- FALSE}
          if (bio_num %in% c(16, 17, 18, 19) && (!prec_p_available)) {can_calculate <- FALSE}
          if (bio_num %in% c(24, 25, 26, 27) && (!srad_p_available)) {can_calculate <- FALSE}
          if (bio_num %in% c(32, 33, 34, 35) && (!mois_p_available)) {can_calculate <- FALSE}
          if (bio_num %in% c(8, 24) && is.null(static_indices_tile[["wettest_period"]]) && !prec_p_available) {can_calculate <- FALSE}
          if (bio_num %in% c(9, 25) && is.null(static_indices_tile[["driest_period"]]) && !prec_p_available) {can_calculate <- FALSE}
          if (bio_num %in% c(18, 26, 34) && is.null(static_indices_tile[["warmest_period"]]) && !temp_p_available) {can_calculate <- FALSE}
          if (bio_num %in% c(19, 27, 35) && is.null(static_indices_tile[["coldest_period"]]) && !temp_p_available) {can_calculate <- FALSE}
          if (bio_num %in% c(26, 27) && (!srad_p_available || !temp_p_available)) {can_calculate <- FALSE}
          if (bio_num %in% c(34, 35) && (!mois_p_available || !temp_p_available)) {can_calculate <- FALSE}
          if (bio_num %in% c(8, 9, 18, 19) && (!temp_p_available || !prec_p_available)) {can_calculate <- FALSE}
          if (bio_num == 7 && (is.null(calculated_bios_in_tile$bio05) || is.null(calculated_bios_in_tile$bio06))) {can_calculate <- FALSE}
          if (bio_num == 3 && (is.null(calculated_bios_in_tile$bio02) || is.null(calculated_bios_in_tile$bio07))) {can_calculate <- FALSE}
          if (bio_num == 15 && (is.null(tile_results$prec_vals) || is.null(calculated_bios_in_tile$bio12))) {can_calculate <- FALSE}
          if (!can_calculate) { 
            if(bio_num %in% bios) warning("Tile ", x, ": Skip BIO", bio_num, " prereqs.")
            next
          }
          idx_vec_unit <- NULL
          idx_vec_period <- NULL
          if (bio_num == 5) idx_vec_unit <- static_indices_tile$warmest_unit
          if (bio_num == 6) idx_vec_unit <- static_indices_tile$coldest_unit
          if (bio_num == 13) idx_vec_unit <- static_indices_tile$wettest_unit
          if (bio_num == 14) idx_vec_unit <- static_indices_tile$driest_unit
          if (bio_num == 8) idx_vec_period <- static_indices_tile$wettest_period %||% tile_results$preci_periods[, "max_idx"]
          if (bio_num == 9) idx_vec_period <- static_indices_tile$driest_period %||% tile_results$preci_periods[, "min_idx"]
          if (bio_num == 10) idx_vec_period <- static_indices_tile$warmest_period %||% tile_results$tempr_periods[, "max_idx"]
          if (bio_num == 11) idx_vec_period <- static_indices_tile$coldest_period %||% tile_results$tempr_periods[, "min_idx"]
          if (bio_num == 16) idx_vec_period <- static_indices_tile$wettest_period %||% tile_results$preci_periods[, "max_idx"]
          if (bio_num == 17) idx_vec_period <- static_indices_tile$driest_period %||% tile_results$preci_periods[, "min_idx"]
          if (bio_num == 18) idx_vec_period <- static_indices_tile$warmest_period %||% tile_results$tempr_periods[, "max_idx"]
          if (bio_num == 19) idx_vec_period <- static_indices_tile$coldest_period %||% tile_results$tempr_periods[, "min_idx"]
          if (bio_num == 21) idx_vec_unit <- static_indices_tile$high_rad_unit
          if (bio_num == 22) idx_vec_unit <- static_indices_tile$low_rad_unit
          if (bio_num == 29) idx_vec_unit <- static_indices_tile$high_mois_unit
          if (bio_num == 30) idx_vec_unit <- static_indices_tile$low_mois_unit
          if (bio_num == 32) idx_vec_period <- static_indices_tile$high_mois_period %||% tile_results$mois_periods[, "max_idx"]
          if (bio_num == 33) idx_vec_period <- static_indices_tile$low_mois_period %||% tile_results$mois_periods[, "min_idx"]
          if (bio_num == 24) idx_vec_period <- static_indices_tile$wettest_period %||% tile_results$preci_periods[, "max_idx"]
          if (bio_num == 25) idx_vec_period <- static_indices_tile$driest_period %||% tile_results$preci_periods[, "min_idx"]
          if (bio_num == 26) idx_vec_period <- static_indices_tile$warmest_period %||% tile_results$tempr_periods[, "max_idx"]
          if (bio_num == 27) idx_vec_period <- static_indices_tile$coldest_period %||% tile_results$tempr_periods[, "min_idx"]
          if (bio_num == 34) idx_vec_period <- static_indices_tile$warmest_period %||% tile_results$tempr_periods[, "max_idx"]
          if (bio_num == 35) idx_vec_period <- static_indices_tile$coldest_period %||% tile_results$tempr_periods[, "min_idx"]

          result_var <- NULL
          bio_output_name <- paste0("bio", sprintf("%02d", bio_num))
          # Call BIOS_fun 
          if (bio_num == 1) {
            result_var <- bio01_fun(tavg = tile_results$temperature_avg, 
                                    cell = cell_ids)
          } else if (bio_num == 2) {
            result_var <- bio02_fun(tmin = tile_results$tmin, 
                                    tmax = tile_results$tmax, 
                                    cell = cell_ids)
          } else if (bio_num == 3) {
            result_var <- bio03_fun(bio02V = calculated_bios_in_tile$bio02[, 1, drop = TRUE],
                                    bio07V = calculated_bios_in_tile$bio07[, 1, drop = TRUE],
                                    cell = cell_ids)  
          } else if (bio_num == 4) {
            result_var <- bio04_fun(tavg = tile_results$temperature_avg, 
                                    cell = cell_ids)
          } else if (bio_num == 5) {
            result_var <- bio05_fun(tmax = tile_results$tmax, 
                                    cell = cell_ids, 
                                    index_vector = idx_vec_unit)
          } else if (bio_num == 6) {
            result_var <- bio06_fun(tmin = tile_results$tmin, 
                                    cell = cell_ids, 
                                    index_vector = idx_vec_unit)
          } else if (bio_num == 7) {
            result_var <- bio07_fun(bio05V = calculated_bios_in_tile$bio05[, 1, drop = TRUE], 
                                    bio06V = calculated_bios_in_tile$bio06[, 1, drop = TRUE],
                                    cell = cell_ids)
          } else if (bio_num == 8) {
            result_var <- bio08_fun(tperiod = tile_results$tempr_periods,
                                    pperiod_max_idx = idx_vec_period,
                                    period_length = period_length,
                                    cell = cell_ids)
          } else if (bio_num == 9) {
            result_var <- bio09_fun(tperiod = tile_results$tempr_periods,
                                    pperiod_min_idx = idx_vec_period,
                                    period_length = period_length,
                                    cell = cell_ids)
          } else if (bio_num == 10) {
            result_var <- bio10_fun(tperiod = tile_results$tempr_periods,
                                    tperiod_max_idx = idx_vec_period,
                                    period_length = period_length,
                                    cell = cell_ids)
          } else if (bio_num == 11) {
            result_var <- bio11_fun(tperiod = tile_results$tempr_periods,
                                    tperiod_min_idx = idx_vec_period,
                                    period_length = period_length,
                                    cell = cell_ids)
          } else if (bio_num == 12) {
            result_var <- bio12_fun(precp = tile_results$prec_vals, 
                                    cell = cell_ids)
          } else if (bio_num == 13) {
            result_var <- bio13_fun(precp = tile_results$prec_vals, 
                                    cell = cell_ids, 
                                    index_vector = idx_vec_unit)
          } else if (bio_num == 14) {
            result_var <- bio14_fun(precp = tile_results$prec_vals, 
                                    cell = cell_ids, 
                                    index_vector = idx_vec_unit)
          } else if (bio_num == 15) {
            result_var <- bio15_fun(precp = tile_results$prec_vals,
                                    bio12V = calculated_bios_in_tile$bio12[, 1, drop = TRUE],
                                    n_units = n_units, 
                                    cell = cell_ids)
          
          } else if (bio_num == 16) {
            result_var <- bio16_fun(pperiod = tile_results$preci_periods,
                                    pperiod_max_idx = idx_vec_period,
                                    cell = cell_ids)
          } else if (bio_num == 17) {
            result_var <- bio17_fun(pperiod = tile_results$preci_periods,
                                    pperiod_min_idx = idx_vec_period,
                                    cell = cell_ids)
          } else if (bio_num == 18) {
            result_var <- bio18_fun(pperiod = tile_results$preci_periods,
                                    tperiod_max_idx = idx_vec_period,
                                    cell = cell_ids)
          } else if (bio_num == 19) {
            result_var <- bio19_fun(pperiod = tile_results$preci_periods,
                                    tperiod_min_idx = idx_vec_period,
                                    cell = cell_ids)
          } else if (bio_num == 20) {
            result_var <- bio20_fun(srad = tile_results$srad_vals, 
                                    cell = cell_ids)
          } else if (bio_num == 21) {
            result_var <- bio21_fun(srad = tile_results$srad_vals, 
                                    cell = cell_ids, 
                                    index_vector = idx_vec_unit)
          } else if (bio_num == 22) {
            result_var <- bio22_fun(srad = tile_results$srad_vals, 
                                    cell = cell_ids, 
                                    index_vector = idx_vec_unit)
          } else if (bio_num == 23) {
            result_var <- bio23_fun(srad = tile_results$srad_vals, 
                                    cell = cell_ids)
          } else if (bio_num == 24) {
            result_var <- bio24_fun(speriod = tile_results$srad_periods,
                                    pperiod_max_idx = idx_vec_period,
                                    cell = cell_ids)
          } else if (bio_num == 25) {
            result_var <- bio25_fun(speriod = tile_results$srad_periods,
                                    pperiod_min_idx = idx_vec_period,
                                    cell = cell_ids)
          } else if (bio_num == 26) {
            result_var <- bio26_fun(speriod = tile_results$srad_periods,
                                    tperiod_max_idx = idx_vec_period,
                                    cell = cell_ids)
          } else if (bio_num == 27) {
            result_var <- bio27_fun(speriod = tile_results$srad_periods,
                                    tperiod_min_idx = idx_vec_period,
                                    cell = cell_ids)
          } else if (bio_num == 28) {
            result_var <- bio28_fun(mois = tile_results$mois_vals, 
                                    cell = cell_ids)
          } else if (bio_num == 29) {
            result_var <- bio29_fun(mois = tile_results$mois_vals, 
                                    cell = cell_ids, 
                                    index_vector = idx_vec_unit)
          } else if (bio_num == 30) {
            result_var <- bio30_fun(mois = tile_results$mois_vals, 
                                    cell = cell_ids, 
                                    index_vector = idx_vec_unit)
          } else if (bio_num == 31) {
            result_var <- bio31_fun(mois = tile_results$mois_vals, 
                                    cell = cell_ids)
          } else if (bio_num == 32) {
            result_var <- bio32_fun(speriod = tile_results$mois_periods,
                                    speriod_max_idx = idx_vec_period,
                                    cell = cell_ids)
          } else if (bio_num == 33) {
            result_var <- bio33_fun(speriod = tile_results$mois_periods,
                                    speriod_min_idx = idx_vec_period,
                                    cell = cell_ids)
          } else if (bio_num == 34) {
            result_var <- bio34_fun(speriod = tile_results$mois_periods,
                                    tperiod_max_idx = idx_vec_period,
                                    cell = cell_ids)
          } else if (bio_num == 35) {
            result_var <- bio35_fun(speriod = tile_results$mois_periods,
                                    tperiod_min_idx = idx_vec_period,
                                    cell = cell_ids)
          }
          if (!is.null(result_var)) { 
            calculated_bios_in_tile[[bio_output_name]] <- result_var
            if (bio_num %in% bios) { 
              rio::export(result_var, 
                bios_qs_paths[[bio_output_name]][x])
            }
          }
          }
        
        # --- Cleanup ---
        rm(tile_results, 
          calculated_bios_in_tile, 
          static_indices_tile, 
          climate_matrix, 
          evars_stack_tile, 
          rt0_geom, 
          extracted_data)
        gc()
        return(NULL)
      }, 
      future.seed = TRUE, 
      future.globals = export_vars, 
      future.packages = c("sf", "terra", "exactextractr", "Rfast", "rio", "purrr"))
      # })
      
      message("Parallel computation finished.")
      return(bioclima_dir)
}