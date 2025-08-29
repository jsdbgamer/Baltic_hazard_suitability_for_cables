# ============================================================================ #
# Model 1B: C-Lion (18 Nov 2024)                                               #
# Proof-of-concept                                                             #
# # Author: Shaaz                                                              #
# Date: 10/08/2025                                                             #
# ============================================================================ #

# Function of this code:
# 1) Loads corridor-scale predictors and aligns them to a single CRS.
# 2) Restricts analysis to a 5.5 km buffer around the BCS cable corridor.
# 3) Rebuilds Euclidean distance to the cable on the working grid.
# 4) Constructs presence (incident points) and corridor-ring background samples
#    (0.5–1.5 km from the cable).
# 5) Fits MaxEnt with 4-fold CV, predicts CV-average surface, then refits on all
#    data to generate the final “full” surface.
# 6) Saves GeoTIFFs, PNGs, training table, and session info.

# Please note the following for re-use:
# - Requires Java and the 'dismo' MaxEnt interface.
# - Corridor-ring background (0.5–1.5 km) avoids near-zero distances and 
#   far-field areas.


# STEP 1: Libraries required
library(sf)     # vector data + handling CRS
library(terra)  # raster data 
library(raster) # RasterStack for dismo::maxent
library(dismo)  # MaxEnt wrapper
library(dplyr)  # preparation

# STEP 2: Base & output directories 
base_dir <- "D:/UCL/Modules/Dissertation" 
out_dir  <- file.path(base_dir, "Outputs", "Model1B_CLion")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# STEP 3: Load predictor rasters 
bathy0 <- rast(file.path(base_dir, "Bathymetry data", "Clipped_bathymetry", "Clipped_Bathy.tif"))
slope0 <- rast(file.path(base_dir, "Bathymetry data", "slope", "Slope_bathy.tif"))
curr0  <- rast(file.path(base_dir, "daily Baltic Sea physical oceanography model", "CurrentSpeed_Clip.tif"))

# STEP 4: Reproject 
target_crs <- "EPSG:32633"
template <- project(bathy0, target_crs, method = "bilinear")
bathy    <- template
slope    <- project(slope0, template, method = "bilinear")
curr     <- project(curr0,  template, method = "bilinear")

# STEP 5: Clip to 5.5 km corridor buffer
corridor_poly <- st_read(file.path(base_dir, "Buffer_split", "C-Lion_cable_Buffer.shp")) |>
  st_transform(target_crs)
corridor_vect <- vect(corridor_poly)

# Align to route extent
bathy <- mask(crop(bathy, corridor_vect), corridor_vect)
slope <- mask(crop(slope, corridor_vect), corridor_vect)
curr  <- mask(crop(curr,  corridor_vect), corridor_vect)

# STEP 6: Rebuild Euclidean distance
cable_sf   <- st_read(file.path(base_dir, "Cable data", "C-Lion shp", "C-Lion_cable.shp")) |>
  st_transform(target_crs)
cable_vect <- vect(cable_sf)
cable_rast <- rasterize(cable_vect, bathy, field = 1, background = NA)
dist       <- distance(cable_rast)
dist       <- mask(crop(dist, corridor_vect), corridor_vect)

# STEP 7: Stack predictors
cov_stack <- c(dist, slope, bathy, curr)
names(cov_stack) <- c("dist", "slope", "bathy", "current")
env_stack <- stack(cov_stack)   

# STEP 8: Background sample (0.5–1.5 km ring)
pres_df <- read.csv(file.path(base_dir, "Vessel data", "C_Lion_incident", "C_lion_incident.csv"))
pres_sf <- st_as_sf(pres_df, coords = c("lon","lat"), crs = 4326) |>
  st_transform(target_crs)

inner_buf <- st_buffer(cable_sf,  500)
outer_buf <- st_buffer(cable_sf, 1500)
ring_sf   <- st_difference(outer_buf, inner_buf) |> st_intersection(corridor_poly)

# Seed
set.seed(20250725)
bg_sf <- st_sample(ring_sf, size = 10000, type = "random") |> st_as_sf(crs = st_crs(ring_sf))

# STEP 9: Extract values
pres_vals <- terra::extract(cov_stack, vect(pres_sf), ID = FALSE) |> cbind(drag = 1)
bg_vals   <- terra::extract(cov_stack, vect(bg_sf),   ID = FALSE) |> cbind(drag = 0)
all_vals  <- na.omit(rbind(pres_vals, bg_vals))

# Save training data
write.csv(all_vals, file.path(out_dir, "training_values_model1B.csv"), row.names = FALSE)

# Convert to dismo inputs
p_coords <- as.matrix(st_coordinates(pres_sf))
a_coords <- as.matrix(st_coordinates(bg_sf))

# STEP 10: MaxEnt -> 4-fold cross-validation 
mx_path <- file.path(out_dir, "MaxEnt_CV")
dir.create(mx_path, showWarnings = FALSE, recursive = TRUE)

set.seed(20250725)
mx_model <- maxent(
  x = env_stack,
  p = p_coords,
  a = a_coords,
  path = mx_path,
  args = c("betamultiplier=1", "replicates=4", "replicatetype=crossvalidate",
    "responsecurves=true", "jackknife=true", "outputformat=cloglog")
)

print(mx_model) # Open CV results in HTML

# CV average (cloglog)
risk_rast <- raster::predict(mx_model, env_stack, args = c("outputformat=cloglog"), progress = "text")

# GeoTIFF
writeRaster(risk_rast, file.path(out_dir, "Risk_Model1B_CLion_CV_cloglog.tif"), overwrite = TRUE)

# Preview of CV run
png(file.path(out_dir, "Risk_Model1B_CLion_CV_preview.png"), width = 1100, height = 800)
plot(risk_rast, main = "Model 1B (C-Lion, CV avg) – Predicted Anchor-Drag Risk (cloglog)")
plot(corridor_poly$geometry, add = TRUE, border = "black")
dev.off()

# STEP 11: MaxEnt –> full-data model 
mx_path_full <- file.path(out_dir, "MaxEnt_Full")
dir.create(mx_path_full, showWarnings = FALSE, recursive = TRUE)

set.seed(20250725)
mx_full <- maxent(
  x    = env_stack,
  p    = p_coords,
  a    = a_coords,
  path = mx_path_full,
  args = c("betamultiplier=1", "outputformat=cloglog",
           "responsecurves=true", "jackknife=true")
)

print(mx_full) # Open full-data results in HTML

# Save output
risk_full <- raster::predict(mx_full, env_stack, args = c("outputformat=cloglog"), progress = "text")

# GeoTIFF
writeRaster(risk_full, file.path(out_dir, "Risk_Model1B_Full_cloglog.tif"), overwrite = TRUE)

# Preview
png(file.path(out_dir, "Risk_Model1B_Full_preview.png"), width = 1100, height = 800)
plot(risk_full, main = "Model 1B (C-Lion, FULL) - Predicted Suitability (cloglog)")
plot(corridor_poly$geometry, add=TRUE, border="black")
dev.off()

# STEP 12: Session info --------------------------------------------------------
capture.output(sessionInfo(), file = file.path(out_dir, "R_sessionInfo.txt"))

# ============================================================================ #
# END OF CODE
# ============================================================================ #