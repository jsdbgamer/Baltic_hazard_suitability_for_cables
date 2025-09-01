# ============================================================================ #
# Model 2 : Predictive Anchor-Drag Map (Baltic AOI)                            #
# Trained on BCS East-West Interlink & C-Lion events                           #
# Author: Shaaz                                                                #
# Date: 10/08/2025                                                             #
# ============================================================================ #

# Function of this code:
# 1) Loads and standardising environmental predictors to a common CRS over 
#    the Baltic AOI.
# 2) Builds presence/background samples from two documented anchor-drag events.
# 3) Fits a MaxEnt model with 4-fold CV for performance assessment, then
#    refits on all data ("full" model) to create the final mapping surface.
# 4) Saves probability-like (cloglog) suitability rasters with PNGs 
# 5) Exports a training table for reproducibility.

# Please note the following for re-use:
# - Requires Java and the 'dismo' MaxEnt interface.
# - All predictors are projected to EPSG:32633 (UTM 33N) and masked to the AOI.
# - The CV surface is for validation.


# STEP 1: Libraries required
suppressPackageStartupMessages({
  library(sf)     # vector data + handling CRS
  library(terra)  # raster data 
  library(raster) # RasterStack for dismo::maxent
  library(dismo)  # MaxEnt wrapper
  library(dplyr)  # preparation
})

# Seed
SEED <- 20250725
set.seed(SEED)

# STEP 2: File paths (EDIT THESE ONLY)
base_dir <- "D:/UCL/Modules/Dissertation"  # outside the repo

baltic_shp     <- file.path(base_dir, "Baltic sea shapefile", "New iho", "Clipped_shaped_baltic.shp")
cable_line_shp <- file.path(base_dir, "Cable data", "Baltic cables", "Baltic_cables.shp")
curr_path      <- file.path(base_dir, "daily Baltic Sea physical oceanography model", "CurrentSpeed_Nov2024.tif")
bathy_path     <- file.path(base_dir, "Bathymetry data", "Clipped_bathymetry", "Clipped_Bathy.tif")
slope_path     <- file.path(base_dir, "Bathymetry data", "slope", "Slope_bathy.tif")
BCS_csv        <- file.path(base_dir, "Vessel data", "Incident_v3", "Lat_Long_Nov_17.csv")
clion_csv      <- file.path(base_dir, "Vessel data", "C_Lion_incident", "C_lion_incident.csv")

# Output folders
out_dir      <- file.path(base_dir, "Outputs", "Model2_Baltic")
mx_path_cv   <- file.path(out_dir, "MaxEnt_CV")
mx_path_full <- file.path(out_dir, "MaxEnt_Full")
dir.create(mx_path_cv,   showWarnings = FALSE, recursive = TRUE)
dir.create(mx_path_full, showWarnings = FALSE, recursive = TRUE)

# STEP 3: Load and project predictors 
target_crs <- "EPSG:32633"  # UTM 33N (metres)

# Load the updated rasters
bathy0 <- rast(bathy_path)
slope0 <- rast(slope_path)
curr0  <- rast(curr_path)

# Create one template and match other to it
template <- project(bathy0, target_crs, method = "bilinear")
bathy    <- template
slope    <- project(slope0, template, method = "bilinear")
curr     <- project(curr0,  template, method = "bilinear")

# load the AOI polygon and convert for masking
baltic_poly <- st_read(baltic_shp, quiet = TRUE) |> st_transform(target_crs)
baltic_vect <- vect(baltic_poly)

# Clip
bathy <- mask(bathy, baltic_vect)
slope <- mask(slope, baltic_vect)
curr  <- mask(curr,  baltic_vect)

# STEP 4: Euclidean distance to nearest cable 
cables_sf   <- st_read(cable_line_shp, quiet = TRUE) |> st_transform(target_crs)
cables_vect <- vect(cables_sf)
cables_rast <- rasterize(cables_vect, bathy, field = 1, background = NA)
dist        <- distance(cables_rast) |> mask(baltic_vect)

# STEP 5: Stack predictors
cov_stack <- c(dist, slope, bathy, curr)
names(cov_stack) <- c("dist", "slope", "bathy", "current")

# STEP 6: Presences & background 
BCS_sf  <- read.csv(BCS_csv)  |> select(lon, lat) |> st_as_sf(coords = c("lon","lat"), crs = 4326) |> st_transform(target_crs)
clio_sf <- read.csv(clion_csv) |> select(lon, lat) |> st_as_sf(coords = c("lon","lat"), crs = 4326) |> st_transform(target_crs)
pres_sf <- rbind(BCS_sf, clio_sf)

# Random background across Baltic Sea
bg_n  <- 20000
bg_sf <- st_sample(baltic_poly, size = bg_n, type = "random") |> st_as_sf(crs = st_crs(baltic_poly))

# Remove NA predictors
pres_ok <- complete.cases(terra::extract(cov_stack, vect(pres_sf), ID = FALSE))
bg_ok   <- complete.cases(terra::extract(cov_stack, vect(bg_sf),   ID = FALSE))
pres_sf <- pres_sf[pres_ok, ]
bg_sf   <- bg_sf[bg_ok, ]

# Save training table
pres_vals <- terra::extract(cov_stack, vect(pres_sf), ID = FALSE) |> cbind(drag = 1)
bg_vals   <- terra::extract(cov_stack, vect(bg_sf),   ID = FALSE) |> cbind(drag = 0)
train_tbl <- na.omit(rbind(pres_vals, bg_vals))
write.csv(train_tbl, file.path(out_dir, "training_values_model2.csv"), row.names = FALSE)

# dismo inputs
p_coords  <- as.matrix(st_coordinates(pres_sf))
a_coords  <- as.matrix(st_coordinates(bg_sf))
env_stack <- raster::stack(cov_stack)

# STEP 7: MaxEnt -> 4-fold cross-validation 
mx_cv <- maxent(
  x = env_stack, 
  p = p_coords, 
  a = a_coords,
  path = mx_path_cv,
  args = c("betamultiplier=1.5", "replicates=4", "replicatetype=crossvalidate",
           "responsecurves=true", "jackknife=true", "outputformat=cloglog")
)
print(mx_cv)  # open CV results in HTML

# CV average (cloglog)
risk_cv <- raster::predict(mx_cv, env_stack, args = c("outputformat=cloglog"), progress = "text")

# STEP 8: MaxEnt –> full-data model 
mx_full <- maxent(
  x = env_stack, 
  p = p_coords, 
  a = a_coords,
  path = mx_path_full,
  args = c("responsecurves=true", "jackknife=true", "outputformat=cloglog")
)
print(mx_full)  # open full-data results in HTML

# Final mapping surface used for figures/mapping
risk_full <- raster::predict(mx_full, env_stack, args = c("outputformat=cloglog"), progress = "text")

# STEP 9: Save the created rasters and show a preview
writeRaster(risk_cv,   file.path(out_dir, "Risk_Model2_CV_cloglog.tif"),   overwrite = TRUE)
writeRaster(risk_full, file.path(out_dir, "Risk_Model2_Full_cloglog.tif"), overwrite = TRUE)

png(file.path(out_dir, "Risk_Model2_CV_preview.png"), width = 1100, height = 800)
plot(risk_cv, main = "Model 2 (CV-avg): Predicted Anchor-Drag Risk (cloglog)")
plot(baltic_poly$geometry, add = TRUE, border = "black")
dev.off()

png(file.path(out_dir, "Risk_Model2_Full_preview.png"), width = 1100, height = 800)
plot(risk_full, main = "Model 2 (Full): Predicted Anchor-Drag Risk (cloglog)")
plot(baltic_poly$geometry, add = TRUE, border = "black")
dev.off()

# STEP 10: Session info 
capture.output(sessionInfo(), file = file.path(out_dir, "R_sessionInfo.txt"))

# Optional, only added for deriving the 90th percentile value
# Load raster
r <- rast("D:/UCL/Modules/Dissertation/Outputs/Model2_Baltic/Risk_Model2_Full_cloglog.tif")

# Calculate 90th percentile
q90 <- quantile(values(r), 0.9, na.rm = TRUE)
q90

# ============================================================================ #
# END OF CODE

# ============================================================================ #
