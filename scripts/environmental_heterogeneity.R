
## loading packages
library("raster")
library(rasterdiv)
library("sf")

#loading species occurrences
spp_points=read.table("0_data/spp_points_7km.csv", header =T, sep=",")

# loading species geographic distribution
spp_geographic_distribution = read.table("1_geographic_classification/spp_geographic_distribution.csv", sep=',', h=T)

# loading raster layers
ras1=raster("0_data/rasters/temperature_diurnal_range")
ras2=raster("0_data/rasters/precipitation_seasonality")
ras3=raster("0_data/rasters/solar_radiation")
ras4=raster("0_data/rasters/soil_pH")

ras_list = list(ras1,ras2,ras3,ras4)

############################# calculating Q for all environment variables #################################

# calculate rao index
rao_env = paRao(ras_list, 
                dist_m="euclidean",
                window=3, 
                rasterOut = T,
                method="multidimension",
                alpha=1, 
                rescale=T,  
                na.tolerance=0.5, 
                simplify=2)

# taking only raster structure
rao_env_raster = rao_env$window.3$alpha.1

#exporting
writeRaster(rao_env_raster, "2_environmental_heterogeneity/rao_raster/rao_env_raster")

########################## assessing Q from species occurrences ###############

rao_env_raster = raster("2_environmental_heterogeneity/rao_raster/rao_env_raster")

# extracting Q values
spp_qvals = raster::extract(rao_env_raster,spp_points[,2:3])
spp_qvals = data.frame(spp_points$species,spp_qvals)
spp_qvals = spp_qvals[-which(is.na(spp_qvals$spp_qvals)),]

# median Q per species
spp_rao =aggregate(spp_qvals[,2], list(spp_qvals[,1]), median)
spp_geographic_distribution = spp_geographic_distribution[match(spp_rao$Group.1, spp_geographic_distribution$species),]
spp_rao = data.frame(spp_geographic_distribution, spp_rao$x)
colnames(spp_rao) = c("species","state", "rao")

# exporting
write.table(spp_rao,"2_environmental_heterogeneity/spp_rao.csv",sep=",", quote=F,row.names=F)

############################# calculating differences and effects #################
spp_rao = read.table("2_environmental_heterogeneity/spp_rao.csv", header =T, sep=",",  na.strings = "NA", fill=T)

state_overall_rao = aggregate(spp_rao$rao, list(spp_rao$state), mean)
state_dispersion_rao = aggregate(spp_rao$rao, list(spp_rao$state), sd )

(state_overall_rao$x[1] - state_overall_rao$x[2])/ state_overall_rao$x[2]
(state_overall_rao$x[1] - state_overall_rao$x[3])/ state_overall_rao$x[3]

############################# altitude  values #################################

#loading species occurrences
spp_points=read.table("0_data/spp_points_7km.csv", header =T, sep=",")

# altitude raster
altitude = raster("0_data/rasters/altitude")

# overall q-values
spp_rao = read.table("2_environmental_heterogeneity/spp_rao.csv", header =T, sep=",",  na.strings = "NA", fill=T)

### extracting altitude values
spp_alt_values = raster::extract(altitude, spp_points[,2:3])
spp_alt_values = data.frame(spp_points$species, spp_alt_values)
spp_alt_values = spp_alt_values[-which(is.na(spp_alt_values$spp_alt_values)),]

### median value per species
spp_altitude = aggregate(spp_alt_values[,2], by=list(spp_alt_values[,1]), median)
colnames(spp_altitude) = c("species", "altitude")

# exporting
write.table(spp_altitude,"2_environmental_heterogeneity/spp_altitude.csv",sep=",", quote=F,row.names=F)

############################### SUPPLEMENTARY ##################################

ras_list = list(ras1,ras2,as.integer(ras3*10^-2),as.integer(ras4))
ras_list[[3]]@data@names = ras3@data@names
ras_list[[4]]@data@names = ras4@data@names

### calculating Q for each environmental variables
for (i in 2:length(ras_list)){
  ras_name = ras_list[[i]]@data@names
  rao_one_env = paRao(as.integer(ras_list[[i]]), dist_m="euclidean", window=3, rasterOut = T, method="classic",alpha=1, na.tolerance=0.5, simplify=2)
  rao_one_env_raster = rao_one_env$window.3$alpha.1
  rao_one_env_raster@data@names = paste("rao",ras_name, sep="_")
  writeRaster(rao_one_env_raster, paste("2_environmental_heterogeneity/rao_raster_per_variable/rao_",ras_name, sep=""), overwrite=T )
}

### loading Q rasters for each environmental variable
rao_ras_list = list()
raster_names = list.files(path = "2_environmental_heterogeneity/rao_raster_per_variable", pattern = ".gri")

for (i in 1:length(raster_names)){
  one_raster = raster(paste("2_environmental_heterogeneity/rao_raster_per_variable", raster_names[i],sep="/") )
  rao_ras_list[[i]] = one_raster
}

rao_ras_stack = raster::stack(rao_ras_list)

# extracting Q values
spp_qvals = raster::extract(rao_ras_stack,spp_points[,2:3])
spp_qvals = data.frame(spp_points$species,spp_qvals)

# median Q per species
spp_rao =aggregate(spp_qvals[,-1], list(spp_qvals[,1]), function(x){median(x, na.rm=T)})
spp_rao = data.frame(spp_geographic_distribution$state, spp_rao)
colnames(spp_rao)[1:2] = c("state","species")

# exporting
write.table(spp_rao,"2_environmental_heterogeneity/spp_rao_per_variable.csv",sep=",", quote=F,row.names=F)
