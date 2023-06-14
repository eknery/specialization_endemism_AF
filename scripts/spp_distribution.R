### setting working directory
setwd("C:/Users/eduar/Documents/GitHub/specialization_endemism_AF")

### loading package
library(sf)
library(plyr)
library(ggplot2)

#loading species occurrences
spp_points=read.table("0_data/spp_points_7km.csv", header =T, sep=",")

#loading TEOW boundaries
teow =st_read("1_geographic_classification/teow_modified")

sf::sf_use_s2(FALSE) 

############################### preparing domain's shape ##########################
### cropping the Neotropics
teow_id = teow[1]
neo_teow_id = st_crop(teow_id, y=c(xmin=-109, xmax=-33, ymin=-56, ymax=27))

# polygon index
length(neo_teow_id$geometry)
index = 1:length(neo_teow_id$geometry)
neo_teow_id$index = index

### labels
centroids = st_centroid(st_geometry(neo_teow_id))
centroid_coords = centroids[[1]][1:2]
for(i in 2:length(centroids)){
  centroid_coords = rbind(centroid_coords, centroids[[i]][1:2])
}

labels = data.frame(centroid_coords, index)
colnames(labels) = c("x", "y", "index")

### colors
id_cols= c("#4BE61D", "#E6E61D", "#1D23E6", "#1D23E6", "#18D4D4", "#E63C1D", "#716C70",
"#E6A61D", "white", "#E6E61D", "#4BE61D", "#E61DC4")


### plotting
tiff("6_graphs/supp_DEC_domains/domains.tiff", units="in", width=5, height=6, res=600)
ggplot() +
  geom_sf(data =neo_teow_id,col="black", aes(fill=as.factor(index))) +
  scale_fill_manual(values=id_cols)+
  #geom_text(data = labels, aes(x, y, label = index), colour = "red", size=3)+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=6),legend.position = "none")
dev.off()

############################ overlapping occurrences and shapes ##############

#converting to sf format
spp_sf = st_as_sf(spp_points, coords = c("longitude", "latitude"))
st_crs(spp_sf) = st_crs(neo_teow_id)

# joining occurrences & shape
join = st_join(spp_sf, neo_teow_id)

# replacing id by domain name
domain_names = c("Ama", "And", "AF", "AF", "Pam", "Cat", "Mes", "Cer", "Pat", "And", "Ama", "Cha") 
domains = domain_names[join$index]
domains[is.na(domains)] = "0"

# points in domains
spp_points_domain = data.frame(spp_points$species, domains)
colnames(spp_points_domain) = c("species", "domain")

# organizing into dataframe
all_spp_name = sort(unique(spp_points$species))
order_domain_names = sort(unique(spp_points_domain$domain))
spp_count_domain = rep(0, (length(order_domain_names)+1))

for(i in 1:length(all_spp_name)){
  sp_name = all_spp_name[i]
  one_sp_domain = spp_points_domain[spp_points_domain$species ==sp_name,2]
  one_count = count(one_sp_domain)
  if(length(one_count$x) <= length(order_domain_names)){
    domain_zero_counts = order_domain_names[!order_domain_names %in% one_count$x]
    for (j in 1:length(domain_zero_counts)){
      one_count = rbind(one_count, c(domain_zero_counts[j], 0) )
    }
  }
  one_count = one_count[order(one_count$x),]
  sp_count = c(sp_name, one_count[,2])
  spp_count_domain = rbind(spp_count_domain, sp_count)
}
spp_count_domain = data.frame(spp_count_domain[-1,])

colnames(spp_count_domain) = c("species", order_domain_names)

write.table(spp_count_domain, "1_geographic_classification/spp_count_domain.csv",sep=",", quote=F, row.names = F)

############################ counting domain combinations #################

library(plyr)

#loading species occurrences
spp_count_domain=read.table("1_geographic_classification/spp_count_domain.csv", header =T, sep=",")

spp_presence_domain = spp_count_domain

for (i in 2:ncol(spp_count_domain)){
  spp_presence_domain[spp_count_domain[,i] != 0, i] = 1
}

domain_combinations = count(spp_presence_domain[, -1])

# exporting
write.table(spp_presence_domain, "1_geographic_classification/spp_presence_domain.csv", sep=',', quote = F, row.names=F)
write.table(domain_combinations, "1_geographic_classification/domain_combinations.csv", sep=',', quote = F, row.names=F)


  