#' Export Bioclimatic Variables to GeoTIFF
#'
#' Assembles intermediate `.qs` files into a full vector in memory, performs
#' cell ID mapping if necessary, and writes final GeoTIFF rasters. Derives geometry from
#' metadata saved within the input directory.
#'
#' @param biovardir Character string. The path to the directory containing the
#'        intermediate `.qs` files AND the `template_info.qs` file.
#' @param save_dir Character string. Path for the final GeoTIFF files. Default: "bioclimatic".
#' @param clean_temporary_files Logical. Delete `biovardir` after writing? Default: FALSE.
#'
#' @return None. Writes GeoTIFF files to `save_dir`.
#' @author Luis Osorio-Olvera, Gonzalo E. Pinilla-Buitrago
#'
#' @export
#' @import terra purrr stringr qs rio

write_layers <- function(biovardir, save_dir = "bioclimatic",
                         clean_temporary_files = FALSE){

  # --- 1. Input Validation and Load Template Info ---
  if (!dir.exists(biovardir)) stop("Dir not found: ", biovardir)
  template_info_file <- file.path(biovardir, "template_info.qs")
  if (!file.exists(template_info_file)) stop("template_info.qs not found.")
  message("Loading template geometry information...")
  template_info <- tryCatch({ 
    qs::qread(template_info_file) 
  }, error = function(e) { 
    stop("Failed read template: ",e$message) 
  })
  if (!is.list(template_info) || !all(c("original_geom", "target_geom") %in% names(template_info))) stop("template_info.qs error.")
  validate_geom <- function(geom_list, name) { 
    if (!is.list(geom_list) || !all(c("extent","dimensions","crs") %in% names(geom_list))) stop(name," fmt err.")
    if (length(geom_list$extent) != 4 || length(geom_list$dimensions) != 2 || !is.character(geom_list$crs)) stop(name," fmt err.") 
  }
  validate_geom(template_info$original_geom, "original_geom")
  validate_geom(template_info$target_geom, "target_geom")
  # --- 2. Find and Organize Intermediate BIO Files ---
  message("Finding and organizing intermediate BIO files...")
  bio_paths <- list.files(biovardir, pattern = "bio\\d{1,2}_.*\\.qs$", full.names = TRUE)
  if (length(bio_paths) == 0) stop("No intermediate bioclim '.qs' files found.")
  bioqs_names <- basename(bio_paths)
  bios_names <- stringr::str_extract(pattern = "bio\\d{1,2}", bioqs_names)
  if (all(is.na(bios_names))) stop("Could not extract BIO names.")
  bio_paths_df <- data.frame(bio_paths, bios_names)
  bio_paths_df <- bio_paths_df[!is.na(bio_paths_df$bios_names), ]
  if (nrow(bio_paths_df) == 0) stop("No valid intermediate files found.")
  bio_paths_dfL <- bio_paths_df |> split(bio_paths_df$bios_names)
  writing_order <- stringr::str_extract(pattern = "\\d{1,2}", names(bio_paths_dfL)) |> as.numeric() |> order()
  bio_paths_dfL <- bio_paths_dfL[writing_order]
  rm(bio_paths, bioqs_names, bios_names, bio_paths_df, writing_order)
  gc()
  message("Found ", length(bio_paths_dfL), " BIO variables.")
  # --- 3. Reconstruct Template Rasters ---
  message("Reconstructing template raster geometries...")
  original_template <- tryCatch({ 
    g <- template_info$original_geom
    terra::rast(xmin = g$extent[1], xmax = g$extent[2], ymin = g$extent[3], ymax = g$extent[4], 
                nrows = g$dimensions[1], ncols = g$dimensions[2], crs = g$crs)
  }, error = function(e) {
    stop("Failed create original template: ", e$message)
  })
  target_template <- tryCatch({
    g <- template_info$target_geom
    terra::rast(xmin = g$extent[1], xmax = g$extent[2], ymin = g$extent[3], ymax = g$extent[4], 
                nrows = g$dimensions[1], ncols = g$dimensions[2], crs = g$crs)
  }, error = function(e) {
    stop("Failed create target template: ", e$message)
  })
  rm(template_info)
  gc()
  if (terra::ncell(target_template) == 0) stop("Target template empty.")
  if (terra::ncell(original_template) == 0) stop("Original template empty.")
  crop_occurred <- !all(terra::ext(original_template) == terra::ext(target_template)) || 
                     !all(terra::dim(original_template)[1:2] == terra::dim(target_template)[1:2])
  n_target_cells <- terra::ncell(target_template)
  message(sprintf("Initializing value vector for %d target cells.", n_target_cells))
  rvals0 <- tryCatch({ 
    rep(NA_real_, n_target_cells) 
  }, error = function(e) { 
    stop(sprintf("Failed allocate NA vector: %s", e$message)) 
  })
  rm(rvals0)
  gc()
  
  # --- 4. Process Each Bioclim Variable ---
  seq_along(bio_paths_dfL) |> purrr::walk(function(x) {
    
    bio_paths_df_base <- bio_paths_dfL[[x]]
    biovar_paths <- bio_paths_df_base$bio_paths
    current_bio_name <- bio_paths_df_base$bios_names[1]

    ms <- paste("\nAssembling", current_bio_name, "in memory...")
    message(ms)
    pb <- utils::txtProgressBar(min = 0, max = length(biovar_paths), style = 3, width = 50, char = "=")
    rvals <- rep(NA_real_, n_target_cells)

    if (crop_occurred) {
      message("\nMapping cell IDs...")
      for (i in seq_along(biovar_paths)) {
        bioval <- NULL
        tryCatch({ 
          bioval <- rio::import(biovar_paths[i]) 
        }, error = function(e) {})
        if (is.null(bioval) || !("cell"%in%colnames(bioval)) || nrow(bioval)==0) {
          utils::setTxtProgressBar(pb, i)
          next
        }
        cellID <- bioval$cell
        xy_pix <- terra::xyFromCell(original_template, cellID)
        cellID <- na.omit(terra::cellFromXY(target_template, xy_pix))
        id_na <- attr(cellID,"na.action")
        if(!is.null(id_na)){
          rvals[cellID] <- bioval[[1]][-id_na]
        } else{
          rvals[cellID] <- bioval[[1]]
          
        }
        utils::setTxtProgressBar(pb, i)
      }
    } else {
      message("\nAssigning cell IDs directly...")
      for (i in seq_along(biovar_paths)) {
        bioval <- NULL; tryCatch({ bioval <- rio::import(biovar_paths[i]) }, error=function(e){})
        if (is.null(bioval) || !("cell"%in%colnames(bioval)) || nrow(bioval)==0) {
          utils::setTxtProgressBar(pb, i)
          next
        }
        cellID <- bioval$cell
        value_data <- bioval[[1]]
        valid_idx <- which(cellID > 0 & cellID <= terra::ncell(target_template))
        if (length(valid_idx) > 0) {
          rvals[cellID[valid_idx]] <- value_data[valid_idx]
        }
        utils::setTxtProgressBar(pb, i)
      }
    }
    close(pb)
    message("Assembly complete.")
    outRast <- terra::rast(target_template, names = current_bio_name)
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    output_file <- file.path(save_dir, paste0(current_bio_name, ".tif"))

    message("Writing GeoTIFF using original block logic: ", output_file)
    write_success <- FALSE
    tryCatch({
        terra::writeStart(outRast, filename = output_file, overwrite = TRUE)

        n_target_cols <- terra::ncol(outRast) # Use outRast (target) dimensions
        n_target_rows <- terra::nrow(outRast)
        n_target_cells_total <- terra::ncell(outRast) # Use total target cells
        step_size_write <- ifelse(nchar(n_target_rows) == 1, 1, 20)
        step_size_write <- min(step_size_write, n_target_rows)
        valores_comienzo <- seq(1, n_target_cells_total + 1,
                                by = n_target_cols * step_size_write)
        if(valores_comienzo[length(valores_comienzo)] > n_target_cells_total) {
           valores_comienzo <- valores_comienzo[-length(valores_comienzo)]
        }

        comienzo <- seq(1, n_target_rows, by = step_size_write)

        message(sprintf("Writing in blocks of %d rows (original logic)...", step_size_write))
        pb_write <- utils::txtProgressBar(min = 0, max = length(comienzo), 
        style = 3, width = 50, char = "=")

        for (i in seq_along(comienzo)) {
            start_row_block <- comienzo[i]
            num_rows_in_block <- min(step_size_write, n_target_rows - start_row_block + 1)
            if (num_rows_in_block <= 0) break
            start_cell_block <- terra::cellFromRowCol(outRast, start_row_block, 1)
            end_cell_block <- terra::cellFromRowCol(outRast, start_row_block + num_rows_in_block - 1, n_target_cols)
            values_for_block <- rvals[start_cell_block:end_cell_block]
            terra::writeValues(x = outRast, v = values_for_block,
                               start = start_row_block, # Start row number
                               nrows = num_rows_in_block) # Number of rows in this block
            rm(values_for_block)
            utils::setTxtProgressBar(pb_write, i)
        }
        close(pb_write)
        write_success <- TRUE
    }, error = function(e){
        warning("Write error ", output_file, ": ", e$message, call. = FALSE)
    }, finally = {
        try(terra::writeStop(outRast), silent = TRUE)
    })
    if (write_success) { 
      message("\n", current_bio_name, " written successfully.")
    } else { 
      if (file.exists(output_file)) try(file.remove(output_file), silent=TRUE)
      warning(current_bio_name, " failed write.") }
    rm(rvals, outRast)
    gc()
    return(NULL)
  })
  
  # --- 5. Final Cleanup ---
  rm(bio_paths_dfL, original_template, target_template, crop_occurred, n_target_cells)
  gc()
  if (clean_temporary_files) { 
    message("Cleaning temp dir: ", biovardir)
    unlink(biovardir, recursive=TRUE, force=TRUE)
    message("Temp cleaned.")
  } else { 
    message("Intermediate files kept: ", biovardir) 
  }
  message("\nAll processed GeoTIFFs written to: ", normalizePath(save_dir))
  return(invisible(NULL))
}