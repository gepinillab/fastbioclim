tmin_path <- list.files(
  # "/Users/gepb/Library/CloudStorage/Dropbox-CityCollege/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/tmin",
  "/Users/Gonzalo/City College Dropbox/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/tmin",
  #"~/Downloads/tmin/",
  pattern = "\\.tif$",
  full.names = TRUE
)

tmax_path <- list.files(
  # "/Users/gepb/Library/CloudStorage/Dropbox-CityCollege/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/tmax",
  "/Users/Gonzalo/City College Dropbox/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/tmax",
  #"~/Downloads/tmax",
  pattern = "\\.tif$",
  full.names = TRUE
)

prec_path <- list.files(
  # "/Users/gepb/Library/CloudStorage/Dropbox-CityCollege/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/prcp",
  "/Users/Gonzalo/City College Dropbox/Gonzalo Pinilla Buitrago/data/rasters/climate/chelsa_2.1/1981-2010/prcp",
  #"~/Downloads/prcp/",
  pattern = "\\.tif$",
  full.names = TRUE
)

# WORLD - M2 (16 Gb)
# [v10] 497 + 710 = 1207
# [v11] 453 + 716 = 1169
# [v12] 426 + 678 = 1104
# [fastbioclim] 558 + 605 = 1163
# [fastbioclim] 564 + 563 = 1126
# COLOMBIA - M1 (8 Gb)
# [t1] 9.76 + 10.988
# [t1] 9.67 + 10.887
col <- AOI::aoi_get(country = "Colombia")
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
j <- 18
bioclim_mex_path <- fastbioclim::bioclim_vars(bios = j,
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
fastbioclim::write_layers(biovardir = bioclim_mex_path ,
  save_dir = "/Users/Gonzalo/bioclim_mex",
  # save_dir = "/Users/gepb/bioclim_mex",
  clean_temporary_files = FALSE)
tictoc::toc()
bio01_mex <- rast("/Users/Gonzalo/bioclim_mex/bio1.tif")
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
  bios = j,
  tmax = tmax_mex,
  tmin = tmin_mex,
  prcp = prcp_mex,
  checkNA = FALSE, 
)
tictoc::toc()
# plot(bios_mex)
# plot(bio01_mex - bios_mex[[1]])

# CHECK BIOS: 18:19
i <- j
# r <- rast(paste0("/Users/gepb/bioclim_mex/bio", i, ".tif"))
r <- rast(paste0("/Users/Gonzalo/bioclim_mex/bio", i, ".tif"))
r_chk <- r - bios_mex[[paste0("bio", sprintf("%02d", i))]]
terra::minmax(r_chk)
plot(r_chk, main = i)

