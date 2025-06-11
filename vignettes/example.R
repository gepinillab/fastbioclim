library(terra)
# ACTIVATE PROGRESS BAR
progressr::handlers(global = TRUE)

tmin_path <- list.files(
  # "/Users/gepb/Library/CloudStorage/Dropbox-CityCollege/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/tmin",
  # "/Users/Gonzalo/City College Dropbox/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/tmin",
  # "~/Downloads/tmin/",
  "/Users/Gonzalo/Library/CloudStorage/GoogleDrive-gepinillab@iecologia.unam.mx/My Drive/data/raster/chelsa/1981-2010",
  pattern = "tmin.*\\.tif$",
  full.names = TRUE
)

tmax_path <- list.files(
  # "/Users/gepb/Library/CloudStorage/Dropbox-CityCollege/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/tmax",
  # "/Users/Gonzalo/City College Dropbox/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/tmax",
  # "~/Downloads/tmax",
  "/Users/Gonzalo/Library/CloudStorage/GoogleDrive-gepinillab@iecologia.unam.mx/My Drive/data/raster/chelsa/1981-2010",
  pattern = "tmax.*\\.tif$",
  full.names = TRUE
)

prec_path <- list.files(
  # "/Users/gepb/Library/CloudStorage/Dropbox-CityCollege/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/prcp",
  # "/Users/Gonzalo/City College Dropbox/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/prcp",
  # "~/Downloads/prcp/",
  "/Users/Gonzalo/Library/CloudStorage/GoogleDrive-gepinillab@iecologia.unam.mx/My Drive/data/raster/chelsa/1981-2010",
  pattern = "prcp.*\\.tif$",
  full.names = TRUE
)

srad_path <- list.files(
  # "/Users/gepb/Library/CloudStorage/Dropbox-CityCollege/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/prcp",
  # "/Users/Gonzalo/City College Dropbox/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/prcp",
  # "~/Downloads/prcp/",
  "/Users/Gonzalo/Library/CloudStorage/GoogleDrive-gepinillab@iecologia.unam.mx/My Drive/data/raster/chelsa/1981-2010",
  pattern = "srad.*\\.tif$",
  full.names = TRUE
)

cmi_path <- list.files(
  # "/Users/gepb/Library/CloudStorage/Dropbox-CityCollege/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/prcp",
  # "/Users/Gonzalo/City College Dropbox/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/prcp",
  # "~/Downloads/prcp/",
  "/Users/Gonzalo/Library/CloudStorage/GoogleDrive-gepinillab@iecologia.unam.mx/My Drive/data/raster/chelsa/1981-2010",
  pattern = "cmi.*\\.tif$",
  full.names = TRUE
)

# WORLD - M2 (16 Gb)
# [v10] 497 + 710 = 1207
# [v11] 453 + 716 = 1169
# [v12] 426 + 678 = 1104
# [fastbioclim] 558 + 605 = 1163
# [fastbioclim] 564 + 563 = 1126
world <- geodata::world(resolution = 1, level = 0, path = tempdir(), 
version = "latest")

# MEXICO
# M1 (8Gb)
# [seq] 24.875 + 17.172 = 42.047
# [w4] 19.557 + 16.111 = 35.668
# [w4] 21.519 + 19.469 = 40.988
# M2 (16Gb)
# [seq] 17.216 + 11.844 = 29.06
# [w4] 12.766 + 12.814 = 25.58
# [w4] 12.433 + 12.731 = 25.164
gc()
future::plan("multisession", workers = 4)
# future::plan("sequential") 
mex <- AOI::aoi_get(country = "Mexico")
tictoc::tic("MEX calculation")
bioclim_mex_path <- fastbioclim::bioclim_vars(bios = 1:19,
                                                    n_units = 12,
                                                    tmin_path = tmin_path,
                                                    tmax_path = tmax_path,
                                                    prec_path = prec_path,
                                                    user_region = mex, 
                                                    write_raw_vars = FALSE,
                                                    temp_dir = "/Users/Gonzalo/bioclim_qs")
                                                    # temp_dir = "/Users/gepb/bioclim_qs")
tictoc::toc()
gc()
future::plan("sequential") 
tictoc::tic("Write MEX")
fastbioclim::write_layers(input_dir = bioclim_mex_path ,
  save_dir = "/Users/Gonzalo/bioclim_mex",
  # save_dir = "/Users/gepb/bioclim_mex",
  clean_temporary_files = FALSE)
tictoc::toc()
bio01_mex <- terra::rast("/Users/Gonzalo/bioclim_mex/bio12.tif")
# bio01_mex <- rast("/Users/gepb/bioclim_mex/bio1.tif")
plot(bio01_mex)

# BIOCLIMA
# M1 (8Gb)
# 15.31 + 60.55 = 75.86
# 17.774 + 82.204 = 99.978
# M2 (16Gb) | En memoría?
# 3.366 + 22.63 = 25.996
# 3.478 + 24.974 = 28.452
gc()
tictoc::tic("Crop MEX")
tmin_mex <- terra::crop(terra::rast(tmin_path), mex, mask = TRUE)
tmax_mex <- terra::crop(terra::rast(tmax_path), mex, mask = TRUE)
prcp_mex <- terra::crop(terra::rast(prec_path), mex, mask = TRUE)
tictoc::toc()

tictoc::tic("Bioclima MEX")
bios_mex <- fastbioclim::clima(
  bios = 1:19,
  tmax = tmax_mex,
  tmin = tmin_mex,
  prcp = prcp_mex,
  checkNA = FALSE, 
)
tictoc::toc()


# CHECK
# i <- 19
# # r <- rast(paste0("/Users/gepb/bioclim_mex/bio", sprintf("%02d", i), ".tif"))
# r <- rast(paste0("/Users/Gonzalo/bioclim_mex/bio", sprintf("%02d", i), ".tif"))
# r_chk <- r - bios_mex[[paste0("bio", sprintf("%02d", i))]]
# terra::minmax(r_chk)
# plot(r_chk, main = i)

# COLOMBIA
# M1 (8 Gb)
# [seq] X + X = X
# [w4] X + X = X
# [w4] X + X = X
# M2 (16 Gb)
# [seq] 9.978 + 13.229 = 23.207
# [w4] 9.935 + 12.998 = 22.933
# [w4] 9.795 + 13.052 = 22.847
gc()
# future::plan("multisession", workers = 4)
future::plan("sequential") 
col <- AOI::aoi_get(country = "Colombia")
tictoc::tic("COL calculation")
bioclim_col_path <- fastbioclim::bioclim_vars(bios = 1:19,
                                                    n_units = 12,
                                                    tmin_path = tmin_path,
                                                    tmax_path = tmax_path,
                                                    prec_path = prec_path,
                                                    user_region = col, 
                                                    write_raw_vars = FALSE,
                                                    temp_dir = "/Users/Gonzalo/bioclim_qs")
                                                    # temp_dir = "/Users/gepb/bioclim_qs")
tictoc::toc()
gc()
future::plan("sequential") 
tictoc::tic("Write COL")
fastbioclim::write_layers(input_dir = bioclim_col_path ,
  save_dir = "/Users/Gonzalo/bioclim_col",
  # save_dir = "/Users/gepb/bioclim_col",
  clean_temporary_files = FALSE)
tictoc::toc()
bio01_col <- rast("/Users/Gonzalo/bioclim_col/bio01.tif")
# bio01_col <- rast("/Users/gepb/bioclim_col/bio01.tif")
plot(bio01_col)

# BIOCLIMA
# M1 (8Gb)
# X + X = X
# X + X = X
# M2 (16Gb)
# 1.378 + 5.087 = 6.465
# 1.376 + 4.512 = 5.888
gc()
tictoc::tic("Crop COL")
tmin_col <- terra::crop(terra::rast(tmin_path), col, mask = TRUE)
tmax_col <- terra::crop(terra::rast(tmax_path), col, mask = TRUE)
prcp_col <- terra::crop(terra::rast(prec_path), col, mask = TRUE)
tictoc::toc()

tictoc::tic("Bioclima COL")
bios_col <- fastbioclim::clima(
  bios = 1:19,
  tmax = tmax_col,
  tmin = tmin_col,
  prcp = prcp_col,
  checkNA = FALSE, 
)
tictoc::toc()

# SOUTH AMERICA
# M1 (8 Gb)
# [seq] X + X = X
# [w4] X + X = X
# [w4] X + X = X
# M2 (16 Gb)
# [seq] 106.436 + 37.997 = 144.433
# [w4] 39.069 + 60.98 = 100.049
# [w4] 39.167 + 58.267 = 97.434
gc()
future::plan("multisession", workers = 4)
# future::plan("sequential") 
sa <- AOI::aoi_get(country = "South America") %>%
  st_union() %>%  
  st_as_sf()    
plot(sa)
tictoc::tic("SA calculation")
bioclim_sa_path <- fastbioclim::bioclim_vars(bios = 1:19,
                                                    n_units = 12,
                                                    tmin_path = tmin_path,
                                                    tmax_path = tmax_path,
                                                    prec_path = prec_path,
                                                    user_region = sa, 
                                                    write_raw_vars = FALSE,
                                                    temp_dir = "/Users/Gonzalo/bioclim_qs")
                                                    # temp_dir = "/Users/gepb/bioclim_qs")
tictoc::toc()
gc()
future::plan("sequential") 
tictoc::tic("Write SA")
fastbioclim::write_layers(input_dir = bioclim_sa_path ,
  save_dir = "/Users/Gonzalo/bioclim_sa",
  # save_dir = "/Users/gepb/bioclim_sa",
  clean_temporary_files = FALSE)
tictoc::toc()
bio01_sa <- rast("/Users/Gonzalo/bioclim_sa/bio1.tif")
# bio01_sa <- rast("/Users/gepb/bioclim_sa/bio1.tif")
plot(bio01_sa)

# BIOCLIMA
# M1 (8Gb)
# X + X = X
# X + X = X
# M2 (16Gb)
# 44.555 + 339.847 = 384.402
# 45.312 + 321.626 = 366.938
gc()
tictoc::tic("Crop SA")
tmin_sa <- terra::crop(terra::rast(tmin_path), sa, mask = TRUE)
tmax_sa <- terra::crop(terra::rast(tmax_path), sa, mask = TRUE)
prcp_sa <- terra::crop(terra::rast(prec_path), sa, mask = TRUE)
tictoc::toc()

tictoc::tic("Bioclima SA")
bios_sa <- fastbioclim::clima(
  bios = 1:19,
  tmax = tmax_sa,
  tmin = tmin_sa,
  prcp = prcp_sa,
  checkNA = FALSE, 
)
tictoc::toc()

# CHECK STATIC
# library(terra)
# static_random <- terra::rast(tmin_path[[1]]) 
# static_random[!is.na(static_random)] <- sample(1:12, 
#   global(!is.na(static_random) , "sum") |> as.numeric(), replace = TRUE)
# static_06 <- static_random
# static_06[!is.na(static_06)] <- 6
# writeRaster(static_random, "/Users/Gonzalo/fastbioclim_eg/static_random.tif")
# writeRaster(static_06, "/Users/Gonzalo/fastbioclim_eg/static_06.tif")
random_path <- "/Users/Gonzalo/fastbioclim_eg/static_random.tif"
six_path <- "/Users/Gonzalo/fastbioclim_eg/static_06.tif"
mex <- AOI::aoi_get(country = "Mexico")
bioclim_sta_path <- fastbioclim::bioclim_vars(bios = 1:19,
                                                    n_units = 12,
                                                    tmin_path = tmin_path,
                                                    tmax_path = tmax_path,
                                                    prec_path = prec_path,
                                                    user_region = mex, 
                                                    write_raw_vars = FALSE,
                                                    warmest_unit_path = random_path,
                                                    coldest_unit_path = six_path,
                                                    driest_unit_path = random_path,
                                                    wettest_unit_path = six_path,
                                                    warmest_period_path = random_path,
                                                    coldest_period_path = six_path,
                                                    driest_period_path = random_path,
                                                    wettest_period_path = six_path,
                                                    temp_dir = "/Users/Gonzalo/bioclim_qs")
                                                    # temp_dir = "/Users/gepb/bioclim_qs")
fastbioclim::write_layers(input_dir = bioclim_sta_path,
  save_dir = "/Users/Gonzalo/bioclim_sta",
  clean_temporary_files = FALSE)
bio05_mex <- rast("/Users/Gonzalo/bioclim_mex/bio05.tif")
bio05_sta <- rast("/Users/Gonzalo/bioclim_sta/bio05.tif")
plot(bio05_mex)
plot(bio05_sta)
bio06_mex <- rast("/Users/Gonzalo/bioclim_mex/bio06.tif")
bio06_sta <- rast("/Users/Gonzalo/bioclim_sta/bio06.tif")
plot(bio06_mex)
plot(bio06_sta)
i <- 10
bio_sta <- rast(paste0("/Users/Gonzalo/bioclim_sta/bio", sprintf("%02d", i), ".tif"))
bio_mex <- rast(paste0("/Users/Gonzalo/bioclim_mex/bio", sprintf("%02d", i), ".tif"))
plot(bio_sta - bio_mex, main = i)

# MEXICO - 35 BIOS
# M1 (8Gb)
# [seq] X + X = 2X
# [w4] X + X = 2X
# [w4] X + X = 2X
# M2 (16Gb)
# [seq] X + X = 2X
# [w4] X + X = 2X
# [w4] X + X = 2X
gc()
future::plan("multisession", workers = 4)
# future::plan("sequential") 
mex <- AOI::aoi_get(country = "Mexico")
tictoc::tic("MEX calculation")
i <- 31
bioclim_mex_path <- fastbioclim::bioclim_vars(bios = 1:35,
                                              n_units = 12,
                                              tmin_path = tmin_path,
                                              tmax_path = tmax_path,
                                              prec_path = prec_path,
                                              srad_path = srad_path,
                                              mois_path = cmi_path, 
                                              user_region = mex, 
                                              write_raw_vars = FALSE,
                                              temp_dir = "/Users/Gonzalo/bioclim_qs")
                                              # temp_dir = "/Users/gepb/bioclim_qs")
tictoc::toc()
gc()
future::plan("sequential") 
tictoc::tic("Write MEX")
fastbioclim::write_layers(input_dir = bioclim_mex_path ,
  save_dir = "/Users/Gonzalo/bioclim_mex",
  # save_dir = "/Users/gepb/bioclim_mex",
  clean_temporary_files = FALSE)
tictoc::toc()
biocheck <- rast(paste0("/Users/Gonzalo/bioclim_mex/bio", 
                        sprintf("%02d", 27), ".tif"))
# bio01_mex <- rast("/Users/gepb/bioclim_mex/bio1.tif")
plot(biocheck)


### CUSTOM VARS
wind_path <- list.files(
  # "/Users/gepb/Library/CloudStorage/Dropbox-CityCollege/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/tmin",
  # "/Users/Gonzalo/City College Dropbox/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/tmin",
  # "~/Downloads/tmin/",
  "/Users/Gonzalo/Library/CloudStorage/GoogleDrive-gepinillab@iecologia.unam.mx/My Drive/data/raster/chelsa/1981-2010",
  pattern = "wind.*\\.tif$",
  full.names = TRUE
)
tavg_path <- list.files(
  # "/Users/gepb/Library/CloudStorage/Dropbox-CityCollege/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/tmax",
  # "/Users/Gonzalo/City College Dropbox/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/tmax",
  # "~/Downloads/tmax",
  "/Users/Gonzalo/Library/CloudStorage/GoogleDrive-gepinillab@iecologia.unam.mx/My Drive/data/raster/chelsa/1981-2010",
  pattern = "tavg.*\\.tif$",
  full.names = TRUE
)
prec_path <- list.files(
  # "/Users/gepb/Library/CloudStorage/Dropbox-CityCollege/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/prcp",
  # "/Users/Gonzalo/City College Dropbox/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/prcp",
  # "~/Downloads/prcp/",
  "/Users/Gonzalo/Library/CloudStorage/GoogleDrive-gepinillab@iecologia.unam.mx/My Drive/data/raster/chelsa/1981-2010",
  pattern = "prcp.*\\.tif$",
  full.names = TRUE
)
progressr::handlers(global = TRUE)
future::plan("multisession", workers = 4)
mex <- AOI::aoi_get(country = "Mexico")
tictoc::tic()
wind_vars <- fastbioclim::stats_vars(
  variable_path = wind_path,
  n_units = 12, 
  stats = c("mean", "max", "min", "cv_cli", "max_period", "min_period"),
  inter_variable_path = prec_path,
  inter_stats = c("max_inter", "min_inter"),
  prefix_variable = "wind",
  suffix_inter_max = "wettest",
  suffix_inter_min = "driest",
  user_region = mex,
  temp_dir = "/Users/Gonzalo/bioclim_qs"
)
tictoc::toc() # 24.441
tictoc::tic()
wind_vars_more <- fastbioclim::stats_vars(
  variable_path = wind_path,
  n_units = 12, 
  stats = NULL,
  inter_variable_path = tavg_path,
  inter_stats = c("max_inter", "min_inter"),
  prefix_variable = "wind",
  suffix_inter_max = "warmest",
  suffix_inter_min = "coldest",
  user_region = mex,
  temp_dir = "/Users/Gonzalo/bioclim_qs"
)
tictoc::toc() # 20.068
tictoc::tic()
fastbioclim::write_layers(input_dir = wind_vars,
  save_dir = "/Users/Gonzalo/bioclim_wind",
  file_pattern = "wind",
  clean_temporary_files = FALSE)
tictoc::toc() # 3.779
tictoc::tic()
fastbioclim::write_layers(input_dir = wind_vars_more,
  save_dir = "/Users/Gonzalo/bioclim_wind",
  file_pattern = "wind",
  clean_temporary_files = FALSE)
tictoc::toc() # 1.086
r_wind <- terra::rast(list.files("/Users/Gonzalo/bioclim_wind", pattern = "*.tif", full.names = TRUE))
plot(r_wind)
