### require packages
library(sf)
library(ape)
library(rangemap)

#loading spp coordinates
spp_points=read.table("0_data/spp_points_7km.csv", header =T, sep=",",  na.strings = "NA", fill=T)

#loading spp geographic classification
spp_geographic_distribution = read.table("1_geographic_classification/spp_geographic_distribution.csv", sep=',', h=T)

############################ measuring convex hull area #############################
all_spp_names = unique(spp_points$species)
all_spp_names =all_spp_names[order(all_spp_names)]

spp_range = c()
for (sp_name in all_spp_names){
  sp_points = spp_points[spp_points$species == sp_name, 2:3]
  multipoint = st_multipoint(as.matrix(sp_points))
  convex = sf::st_convex_hull(multipoint)
  spp_range = c( spp_range, st_area(convex))
}

spp_range = data.frame(all_spp_names,  spp_range)
colnames(spp_range) = c("species","range")

#exporting
write.table(spp_range, "4_range_inference/spp_range.csv", sep=',', quote=F, row.names=F)

### setting regimes
spp_range = read.table("4_range_inference/spp_range.csv", sep=',', h=T)
regimes = data.frame(spp_geographic_distribution, spp_range$range)
colnames(regimes)[3] = "range"

# exporting
write.table(regimes, "4_range_inference/regimes_ranges.csv", sep=',', quote=F, row.names=F)

####################### describing convex hull area per distribution ##################
### load
regimes_ranges = read.table("4_range_inference/regimes_ranges.csv", sep=',', h=T)

### summary
aggregate(regimes_ranges$range, by=list(regimes_ranges$state), median)
aggregate(regimes_ranges$range, by=list(regimes_ranges$state), IQR)


######################### sampling effort ##########################

### species per state
aggregate(spp_geographic_distribution$state, by=list(spp_geographic_distribution$state), length)

### state of each species point
point_state = spp_points$species
all_states = sort(unique(spp_geographic_distribution$state))

for (state_name in all_states){
  group = spp_geographic_distribution$species[spp_geographic_distribution$state == state_name]
  point_state[spp_points$species %in% group] = state_name
}

# number of points in each state
aggregate(spp_points$species, by=list(point_state), length)

### state of each species
spp_n_points = aggregate(spp_points$species, by=list(spp_points$species), length)
spp_state = spp_n_points$Group.1

for (state_name in all_states){
  group = spp_geographic_distribution$species[spp_geographic_distribution$state == state_name]
  spp_state[spp_n_points$Group.1 %in% group] = state_name
}

# mean and sd of points in each state
aggregate(spp_n_points$x, by=list(spp_state), mean)
aggregate(spp_n_points$x, by=list(spp_state), sd)
