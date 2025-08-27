# =============================================================================
#  Compute surface current speed 
# =============================================================================

# Required libraries
library(terra)
library(ncdf4)  

# Directory
base_dir <- "D:/UCL/Modules/Dissertation"

# Daily files (17 & 18 Nov 2024)
nc_path1 <- file.path(base_dir, "daily Baltic Sea physical oceanography model",
                      "BAL-NEMO_PHY-DailyMeans-20241117.nc")
nc_path2 <- file.path(base_dir, "daily Baltic Sea physical oceanography model",
                      "BAL-NEMO_PHY-DailyMeans-20241118.nc")

# Monthly file (Nov 2024)
nc_month <- file.path(base_dir, "daily Baltic Sea physical oceanography model",
                      "BAL-NEMO_PHY-MonthlyMeans-202411.nc")

# Output paths
out_event  <- file.path(base_dir, "daily Baltic Sea physical oceanography model",
                        "CurrentSpeed_Clip.tif")
out_month  <- file.path(base_dir, "daily Baltic Sea physical oceanography model",
                        "CurrentSpeed_Nov2024.tif")


# 17–18 Nov 2024
# Load u and v components for each day 
uo1 <- rast(nc_path1, subds = "uo")[[1]]
vo1 <- rast(nc_path1, subds = "vo")[[1]]
uo2 <- rast(nc_path2, subds = "uo")[[1]]
vo2 <- rast(nc_path2, subds = "vo")[[1]]

# Same resolution
uo2 <- resample(uo2, uo1, method = "bilinear")
vo2 <- resample(vo2, uo1, method = "bilinear")
vo1 <- resample(vo1, uo1, method = "bilinear")  

# Compute current speed per day 
speed1 <- sqrt(uo1^2 + vo1^2)
speed2 <- sqrt(uo2^2 + vo2^2)
names(speed1) <- "current"
names(speed2) <- "current"

# Average the two daily speed rasters
event_stack  <- c(speed1, speed2)
event_mean   <- app(event_stack, mean, na.rm = TRUE)
names(event_mean) <- "current"

# Save GeoTIFF 
writeRaster(event_mean, out_event, overwrite = TRUE)


# MONTHLY MEAN (November 2024)
# Load monthly u & v components
uo_m <- rast(nc_month, subds = "uo")[[1]]
vo_m <- rast(nc_month, subds = "vo")[[1]]

# Same resolution
vo_m <- resample(vo_m, uo_m, method = "bilinear")

# Compute monthly mean current speed
nov_mean_speed <- sqrt(uo_m^2 + vo_m^2)
names(nov_mean_speed) <- "current"

# Save GeoTIFF
writeRaster(nov_mean_speed, out_month, overwrite = TRUE)

# Quick checks
plot(event_mean, main = "Event-avg current speed (17–18 Nov 2024)")
plot(nov_mean_speed, main = "Monthly mean current speed (Nov 2024)")
summary(values(nov_mean_speed))

