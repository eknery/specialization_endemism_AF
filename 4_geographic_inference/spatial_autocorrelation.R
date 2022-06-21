setwd("C:/Users/eduar/Documents/GitHub/specialization_endemism_AF")

# loading packages 
library(raster)

#loading species occurrences
spp_points=read.table("0_data/spp_points_7km.csv", header =T, sep=",",  na.strings = "NA", fill=T)

### loading  geographic distribution
spp_geographic_distribution = read.table("1_geographic_classification/spp_geographic_distribution.csv", sep=',', h=T)

# loading raster layers
ras1=raster("0_data/rasters/temperature_diurnal_range")
ras2=raster("0_data/rasters/precipitation_seasonality")
ras3=raster("0_data/rasters/solar_radiation")
ras4=raster("0_data/rasters/soil_pH")

ras_list = list(ras1,ras2,ras3,ras4)
ras_stack = stack(ras1,ras2,ras3,ras4)

############# calculating local Moran across environmental rasters ############

### function to calculate Moran index across more than one raster
multi_moran = function(raster_list){
  moran_stack = MoranLocal(raster_list[[1]])
  for (i in 2:length(raster_list)){
    one_moran =  MoranLocal(raster_list[[i]])
    moran_stack = stack(moran_stack, one_moran)
  }
  return(moran_stack)
}

# calculating multiple Moran indexes
multi_moran_ras = multi_moran(ras_list)
plot(multi_moran_ras)

# extracting Moran values 
moran_df = raster::extract(multi_moran_ras, spp_points[,2:3])
spp_moran_df = data.frame(spp_points[,1], moran_df)
colnames(spp_moran_df) = c("species", names(ras_stack))

# median value by species
spp_moran = aggregate(spp_moran_df[,2:5], by=list(spp_moran_df[,1]), function(x){median(x,na.rm=T)})
colnames(spp_moran)[1] = "species"

# summary into single value
moran = apply(spp_moran[,2:5], MARGIN = 1, mean)
names(moran) = spp_moran$species

# by geographic distribution
aggregate(moran, by=list(spp_geographic_distribution$state), median)
aggregate(moran, by=list(spp_geographic_distribution$state), IQR)

# exporting
write.table(spp_moran,"4_geographic_inference/spp_moran.csv",sep=",", quote=F,row.names=F)
