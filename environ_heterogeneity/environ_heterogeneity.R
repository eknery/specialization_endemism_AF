
library(raster)
library(sf)
library(rasterdiv)
library(ape)
library(phytools)
library(nlme)


# loading mcc phylogeny
mcc=read.tree("mcc_phylo.nwk")

# loading background sites
bg_points=read.table("bg_points.csv", header =T, sep=",")

#loading species occurrences
spp_points=read.table("spp_points.csv", header =T, sep=",")

# loading raster layers
ras1=raster("temp_diurnal_range")
ras2=raster("precip_seasonality")
ras3=raster("soil_pH")
ras4=raster("solar_radiation")

allras=stack(ras1,ras2,ras3,ras4)

#loading AF-shape boundaries
shape_AF=st_read("AF_shape")

######################### classifying sites into AF ################

# checking sample per cell
cell_count=rasterize(bg_points, ras1, fun="count")
sampled_cells=which(cell_count@data@values==1)

# sampling per geographic domain
bg_sf = st_as_sf(bg_points, coords = c("longitude", "latitude"))
st_crs(bg_sf) = st_crs(shape_AF)
  
bg_intersect = st_intersects(bg_sf,shape_AF)

inside = lengths(bg_intersect) == 1
outside = lengths(bg_intersect) == 0

inside_AF=bg_points[inside,]
outside_AF=bg_points[outside,]

######################## extracting environmental background ###################

bg_env=extract(allras,rbind(inside_AF,outside_AF))
region=c(rep("inside_AF",length(inside_AF[,1])),rep("outside_AF",length(outside_AF[,1])))
bg_env=data.frame(region,bg_env)

# checking NAs
bg_nas = c()
for (i in 1:length(bg_env[1,])){
  bg_nas= c(bg_nas, which(is.na(bg_env[,i])) )
}

# writing df
write.table(bg_env,"environmental_background.csv",sep=",", quote=F,row.names=F)

########################## assessing AF-endemism ####################

spp_sf = st_as_sf(spp_points, coords = c("longitude", "latitude"))
st_crs(spp_sf) = st_crs(shape_AF)

spp_intersect = st_intersects(spp_sf,shape_AF)

bool_inside = lengths(spp_intersect) == 1
spp_point_inside =spp_points[bool_inside,]
total_inside = aggregate(spp_point_inside[,1], by=list(spp_point_inside[,1]), length)

total_point = aggregate(spp_points[,1], by=list(spp_points[,1]), length)

distribution_table = data.frame(matrix(NA, nrow=nrow(total_point), ncol=4))
colnames(distribution_table) = c("species","n", "inside_AF", "ratio_AF")
distribution_table$species = total_point$Group.1

for (i in unique(distribution_table$species)){
  inside = total_inside[total_inside[,1] == i,2]
  total = total_point[total_point[,1] == i,2]
  index = distribution_table$species == i
    if ( length(inside) == 0){ inside = 0 }
  distribution_table[index,2] = total
  distribution_table[index,3] = inside
  distribution_table[index,4] = inside/total
}

distribution_table$ratio_AF = round(distribution_table$ratio_AF,2)

AF_endemics = distribution_table$species[distribution_table$ratio_AF > 0.75]

write.table(distribution_table,"distribution_table.csv",sep=",", quote=F,row.names=F)

######################## extracting species environment ###################

spp_env=extract(allras,spp_points[,2:3])
species= spp_points$species
spp_env= data.frame(species,spp_env)

# checking NAs
spp_nas = c()
for (i in 1:length(spp_env[1,])){
  spp_nas= c(spp_nas, which(is.na(spp_env[,i])) )
}

spp_env = spp_env[-spp_nas,]

# writing df
write.table(spp_env,"species_environment.csv",sep=",", quote=F,row.names=F)

####################### calculating Q and testing differences #################

# scaling
ras_list=list(ras1,as.integer(ras2*10^-1),as.integer(ras3 * 10^-3),ras4)

# calculate rao index
rao=paRao(ras_list, dist_m="euclidean", window=3, rasterOut = T,
          method="multidimension",alpha=1, rescale=F,  na.tolerance=0.5, simplify=2)

rao_ras=rao$window.3$alpha.1

######################### assessing Q from bg sites #############################
bg_qvals=extract(rao_ras,bg_points)
bg_rao=data.frame(bg_points$region,bg_qvals)
colnames(bg_rao)=c("region","rao")

write.table(bg_rao,"bg_rao.csv",sep=",", quote=F,row.names=F)

### testing by permutation
# observed difference
mean_bg_rao=aggregate(bg_rao[,2], list(bg_rao[,1]), mean)
diff=mean_bg_rao[1,2] - mean_bg_rao[2,2]

#simulating random differences
sim_diff=rep(NA, 10000)

for (i in 1:10000){
  sim_rao=data.frame(sample(bg_rao$region),bg_rao$rao)
  sim_means=aggregate(sim_rao[,2], list(sim_rao[,1]), mean)
  sim_diff[i]=sim_means[1,2] - sim_means[2,2]
}

# p-value
1 - ( sum(diff > sim_diff) / 10000)


########################## assessing Q from species occurrences ###############

spp_qvals=extract(rao_ras,spp_points[,2:3])
spp_rao=data.frame(spp_points$species,spp_qvals)
spp_rao = spp_rao[-which(is.na(spp_qvals)),]
colnames(spp_rao)=c("species","rao")

# mean Q per species
mean_spp_rao =aggregate(spp_rao[,2], list(spp_rao[,1]), mean)
colnames(mean_spp_rao) = c("species","rao")

distribution = mean_spp_rao$species
distribution[which(distribution %in% AF_endemics)]="AF-endemic"
distribution[which(distribution != "AF-endemic")]="widespread"
mean_spp_rao = data.frame(distribution, mean_spp_rao)

write.table(mean_spp_rao,"mean_spp_rao.csv",sep=",", quote=F,row.names=F)


### pgls testing
q = mean_spp_rao$rao
names(q) = mean_spp_rao$species

signal_q = phylosig(tree=mcc, x= q, method="lambda")
pagel_cor_q = corPagel(signal_q$lambda, phy = mcc, form = ~1,  fixed = T)

# fit pgls
fit_gls_q = gls(log(rao) ~ distribution, data=mean_spp_rao, correlation=pagel_cor_q,  method = "ML")
summary(fit_gls_q)
plot(fit_gls_q)

#checking residuals
res_gls_q = resid(fit_gls_q)[1:66]
hist(res_gls_q)
shapiro.test(res_gls_q)
