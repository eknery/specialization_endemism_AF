setwd("C:/Users/eduar/Documents/GitHub/specialization_endemism_AF")

## loading packages
library(raster)
library(rasterdiv)
library(sf)

#loading species occurrences
spp_points=read.table("0_data/spp_points_7km.csv", header =T, sep=",")

# loading raster layers
ras1=raster("0_data/rasters/temperature_diurnal_range")
ras2=raster("0_data/rasters/precipipatation_seasonality")
ras3=raster("0_data/rasters/soil_pH")
ras4=raster("0_data/rasters/solar_radiation")

altitude = raster("0_data/rasters/altitude")

############################# calculating Q from environment #################################

# scaling
ras_list=list(ras1,as.integer(ras2*10^-1),as.integer(ras3 * 10^-3),ras4)

# calculate rao index
rao_env =paRao(ras_list, dist_m="euclidean", window=3, rasterOut = T,
          method="multidimension",alpha=1, rescale=F,  na.tolerance=0.5, simplify=2)

# taking only raster structure
rao_env=rao_env$window.3$alpha.1

#exporting
writeRaster(rao_env, "2_environmental_heterogeneity/rao_rasters/rao_env_3")

########################## assessing Q from species occurrences ###############

rao_env = raster("2_environmental_heterogeneity/rao_rasters/rao_env_3")
spp_geographic_state = read.table("1_geographic_classification/spp_geographic_state.csv", sep=',', h=T)

# extracting Q values
spp_qvals = extract(rao_env,spp_points[,2:3])
spp_qvals = data.frame(spp_points$species,spp_qvals)
spp_qvals = spp_qvals[-which(is.na(spp_qvals)),]

# median Q per species
spp_rao =aggregate(spp_qvals[,2], list(spp_qvals[,1]), median)
spp_geographic_state = spp_geographic_state[match(spp_rao$Group.1, spp_geographic_state$species),]
spp_rao = data.frame(spp_rao, spp_geographic_state$state)
colnames(spp_rao) = c("species","rao", "state")

# exporting
write.table(spp_rao,"2_environmental_heterogeneity/spp_rao.csv",sep=",", quote=F,row.names=F)

############################# calculating Q from altitude #################################
# calculate rao index
altitude_scale = altitude*10^-2

rao_altitude =paRao(altitude_scale, dist_m="euclidean", window=3, rasterOut = T, method="classic",alpha=1, rescale=F,  na.tolerance=0.5, simplify=2)
