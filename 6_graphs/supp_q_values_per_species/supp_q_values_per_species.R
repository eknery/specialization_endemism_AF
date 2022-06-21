############################## Q value per species ####################################
# wd
setwd("C:/Users/eduar/Documents/GitHub/specialization_endemism_AF")

# packages
library(raster)
library(ggplot2)
library(plyr)

#loading species occurrences
spp_points=read.table("0_data/spp_points_7km.csv", header =T, sep=",")

# loading species geographic state
spp_geographic_distribution = read.table("1_geographic_classification/spp_geographic_distribution.csv", sep=',', h=T)
state = sort(unique(spp_geographic_distribution$state))

# loading rao raster
rao_ras= raster("2_environmental_heterogeneity/rao_rasters/rao_env_3")

### all spp names
all_spp_names = spp_geographic_state$species

### q-vals
qvals = raster::extract(rao_ras,spp_points[,2:3])
qvals[- which(is.na(qvals))]
spp_qvals= data.frame( spp_points$species[- which(is.na(qvals))], qvals[- which(is.na(qvals))] )
colnames(spp_qvals) = c("species", "qval")

# qval range
qval_range = range(spp_qvals$qval)

### color code
mycols = c( "#1E88E5", "#FFC107", "#D81B60")
names(mycols) = state

###### plotting 
for(one_state in state){
  # sorting by geographic state
  sub_state = subset(spp_geographic_distribution, state == one_state)
  state_sp_names = sub_state$species
  # hot to separate species into plots
  graph_breaks = seq(from=0, to=length(state_sp_names), by=4)
  graph_breaks = c(graph_breaks, length(state_sp_names) )
  graph_breaks = graph_breaks[-1]
  # taking qvals by state
  state_spp_qvals = spp_qvals[spp_qvals$species %in% state_sp_names,]
  # getting list ready
  plot_list = list()
  # looping over species
  for (i in  1:length(graph_breaks) ){
    last = graph_breaks[i]
    one_set = state_spp_qvals[state_spp_qvals$species %in% state_sp_names[(last-3):last],]
    one_col = mycols[names(mycols) == one_state]
    median_vals = ddply(one_set, "species", summarise, sp_median=median(qval))
    sp_plot = ggplot(one_set, aes(x=qval)) + 
      geom_histogram(binwidth=50, color="black", fill=one_col ) +
      geom_vline(data=median_vals, aes(xintercept=sp_median), color="black", size=1.5, linetype="dashed" )+
      xlim(qval_range) +
      xlab("Rao's Q value")+
      facet_grid(rows= vars(species)) +
      theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"))
    plot_list[[i]] = sp_plot
  }
  
  for (i in 1:length(plot_list)){
    file_name = paste(one_state,"spp_qval_hist", as.character(i), sep="_")
    output_dir = paste("6_graphs/supp_q_values/",file_name,".tiff", sep="")
    tiff(output_dir, units="in", width=3, height=6, res=600)
    print(plot_list[[i]])
    dev.off()
  }
}
