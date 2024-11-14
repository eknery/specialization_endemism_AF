# loading packages 
library(raster)
library(sf)

#loading species occurrences
bg_points=read.table("0_data/bg_points_7km.csv", header =T, sep=",")

#loading TEOW boundaries
teow =st_read("1_geographic_classification/teow_modified")
sf::sf_use_s2(FALSE) 

############## classifying points into AF and outside ########################

# cropping the Neotropics
teow_id = teow[1]
neo_teow_id = st_crop(teow_id, y=c(xmin=-109, xmax=-33, ymin=-56, ymax=27))

# polygon index
length(neo_teow_id$geometry)
index = 1:length(neo_teow_id$geometry)
neo_teow_id$index = index

#converting points to sf format
bg_sf = st_as_sf(bg_points, coords = c("longitude", "latitude"))
st_crs(bg_sf) = st_crs(neo_teow_id)

# joining points & shape
join = st_join(bg_sf, neo_teow_id)

# replacing id by domain name
domain_names = c("out", "out", "AF", "AF",
                 "out", "out", "out", "out",
                 "out", "out", "out", "out") 
domains = domain_names[join$index]
domains[is.na(domains)] = "out"

bg_points_domains = data.frame(bg_points, domains)

write.table(bg_points_domains, "1_geographic_classification/bg_points_domains.csv", sep=",", quote = F, row.names= F)
