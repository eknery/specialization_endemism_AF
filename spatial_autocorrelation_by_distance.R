setwd("C:/Users/eduar/Documents/GitHub/specialization_endemism_AF")

# loading packages 
library(raster)
library(geosphere)
library(ape)

#loading species occurrences
bg_points_domains=read.table("1_geographic_classification/bg_points_domains.csv", header =T, sep=",")

# loading raster layers
ras1=raster("0_data/rasters/temperature_diurnal_range")
ras2=raster("0_data/rasters/precipitation_seasonality")
ras3=raster("0_data/rasters/solar_radiation")
ras4=raster("0_data/rasters/soil_pH")

ras_stack = stack(ras1,ras2,ras3,ras4)

######################### calculating Moran's I #######################

### row classification according to domain
n_AF = which(bg_points_domains$domains == "AF")
n_out = which(bg_points_domains$domains != "AF")

### extracting env values
env_values = raster::extract(ras_stack, bg_points_domains[,1:2])
af_env_values = env_values[n_AF,]
out_env_values = env_values[n_out,]

### calculating weights from geographic distance 

# inside AF
AF_points = cbind(bg_points_domains[n_AF,1],bg_points_domains[n_AF,2])
AF_distances = distm(x=AF_points, fun=distGeo)
inv_AF_distances <- 1/AF_distances
diag(inv_AF_distances) <- 0

#outside AF
out_points = cbind(bg_points_domains[n_out,1],bg_points_domains[n_out,2])
out_distances = distm(x=out_points, fun=distGeo)
inv_out_distances = 1/out_distances
diag(inv_out_distances) <- 0

### calculating Moran I
af_moran_i = Moran.I(x=af_env_values[,4], weight=inv_AF_distances)

out_moran_i = Moran.I(x=out_env_values[,4], weight=inv_out_distances)
