# --------------------------- #
#
# Description:
# Calculate geographic distances between marine sites
#
# Notes:
# R script mainly uses functions from marmap library
#
# --------------------------- #

# Load libraries
library(marmap)


#--------------#
#
# Prepare coordinates
#
#--------------#

# Vector of latitude and longitude
lat = c(40.36, 44.22, 41.30, 36.78, 45.50, 50.47, 56.76, 58.70, 61.43, 54.85)
lon = c(25.00, 13.78, 8.47 , -8.23, -1.86,  0.99, -6.78, 10.47,  4.61, 19.24)

# Vector of site labels
site = c("Med1","Med2","Med3","Atl1","Atl2","Atl3","Atl4","Atl5","Atl6","Bal1")

# Create a dataframe of coordinates
coords = data.frame(Site = site, Lat = lat, Lon = lon)


#--------------#
#
# Calculate least-cost distances
#
#--------------#

# Download bathymetric data
bathy_data = getNOAA.bathy(lon1 = -20, lon2 = 30, lat1 = 35, lat2 = 65, res = 4, keep = TRUE)
summary(bathy_data)

# Get depth of coordinates
# Must be less than ten metres (<-10) to compute distances (may need to adjust coordinates)
get.depth(bathy_data, x = coords$Lon, y = coords$Lat, locator = FALSE)

# Create transition object [long run time]
# Use a minimum depth of -10 to avoid path crossing land masses
# Use a maximum depth of -200 to limit paths to continental shelf
trans1 = trans.mat(bathy_data, min.depth = -10, max.depth = NULL)
# save(trans1, file = "transition_object.RData")
# load("transition_object.RData")

# Compute least-cost paths [long run time]
lc_paths = lc.dist(trans1, subset(coords, select = c("Lon","Lat")), res = "path")
# save(lc_paths, file = "least_cost_paths.RData")
# load("least_cost_paths.RData")

# Plot paths on a map
# Visually check that no path overlaps land
plot.bathy(bathy_data, image= TRUE, land = TRUE, n = 0,
           bpal = list(c(0, max(bathy_data), "grey"),
                       c(min(bathy_data), 0, "royalblue")))
lapply(lc_paths, lines, col = "orange", lwd = 2, lty = 1)

# Compute least-cost distances (km) matrix
lc_dist = lc.dist(trans1, subset(coords, select = c("Lon","Lat")), res = "dist")

# Convert to matrix, rename columns and rows, and export as csv file
lc_mat = as.matrix(lc_dist)
colnames(lc_mat) = as.vector(coords$Site)
rownames(lc_mat) = as.vector(coords$Site)
lc_mat
write.csv(lc_mat, file="lc_distances_km.csv")
