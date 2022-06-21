################################ Rao's Q map ###############################

# wd
setwd("C:/Users/eduar/Documents/GitHub/specialization_endemism_AF")

library(raster)

# loading rao raster
rao_env = raster("2_environmental_heterogeneity/rao_rasters/rao_env_raster")

# getting rid of outliers
raster::cellStats(rao_env, stat=hist)
rao_env[rao_env > 0.75] = 0.75

# graphical parameters
col_func = colorRampPalette(c("steelblue1","green4", "orange1", "red3")) 
legend_colors= col_func(20) 

# plotting
tiff("6_graphs/supp_q_map/rao_q_map.tiff", units="in", width=5.5, height=6, res=600)
plot(rao_env, col=legend_colors , colNA="white")
dev.off()
