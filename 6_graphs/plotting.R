################################### plotting ################################
setwd("C:/Users/eduar/Documents/GitHub/specialization_endemism_AF")

library(tidyverse)
library(PupillometryR)
library(ggpubr)
library(readr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)

### loading data
spp_rao=read.table("2_environmental_heterogeneity/spp_rao.csv", header =T, sep=",",  na.strings = "NA", fill=T)
spp_hvolumes=read.table("3_hypervolume_inference/spp_hvolumes.csv", header =T, sep=",",  na.strings = "NA", fill=T)
spp_range=read.table("4_geographic_inference/spp_range.csv", header =T, sep=",",  na.strings = "NA", fill=T)

### organizing dataset to plot
spp_dataset = data.frame(spp_rao, spp_hvolumes$hvolume, spp_range$range)
colnames(spp_dataset)[4:5] = c("hvolume","range")

### my colors
mycols = c( "#1E88E5", "#FFC107", "#D81B60")

############################ plotting current species data #########################

tiff("6_graphs/Qvalues_geostate.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= spp_dataset, aes(x=state, y=rao, fill=state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 2, alpha = 0.65) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.25)+
  geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.25) +
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("geographic state")+ ylab("median Q value")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFothers" = "AF and Other\ndomains", "others" = "Other\ndomains"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=8),legend.position = "none")
dev.off()

tiff("6_graphs/hvolumes_geostate.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= spp_dataset, aes(x=state, y=hvolume, fill=state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 2, alpha = 0.65) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.25)+
  geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.25) +
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("geographic state")+ ylab("hypervolume size")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFothers" = "AF and Other\ndomains", "others" = "Other\ndomains"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=8),legend.position = "none")
dev.off()


ggplot(data= spp_dataset, aes(x=state, y=range, fill=state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 2, alpha = 0.65) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.25)+
  geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.25) +
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("geographic state")+ ylab("range size (km2)")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFothers" = "AF non-endemic", "others" = "Others"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=8),legend.position = "none")

############################### plotting phylogenetic data ############################

# mcc tree
mcc = read.tree("0_data/mcc_phylo.nwk")
n_tips = Ntip(mcc)
n_inner_nodes = mcc$Nnode

# DEC probabilities
mcc_node_probs = read.table("5_comparative_analysis/mcc_node_probs.csv",sep=",", h=T)

# DEC states
mcc_node_states = read.table("5_comparative_analysis/mcc_node_states.csv", sep=",", h=T)

# hypervolume data
hv_spp_volumes=read.table("3_hypervolume_inference/hv_spp_volumes.csv", header =T, sep=",",  na.strings = "NA", fill=T)
hvolumes = hv_spp_volumes$volumes
names(hvolumes) = hv_spp_volumes$species

### visual parameters
# tip states
extent_states = mcc_node_states$node_state[1:n_tips]
names(extent_states) = mcc$tip.label

# tip states probs
tip_states_probs = to.matrix(extent_states[tr$tip.label], c("AF","AFothers", "others"))

# ancestral state probs
inner_node_probs = mcc_node_probs[(1+n_tips):(n_tips+n_nodes),]
inner_node_probs = inner_node_probs[,c(1,3,2)]

# state colors
state_cols=c( "#1E88E5", "#FFC107", "#D81B60")
names(state_cols)=c("AF", "AFothers", "others")

# bar colors
bar_cols=extent_states
bar_cols[extent_states == "AF"] = "#1E88E5"
bar_cols[extent_states == "AFothers"] = "#FFC107"
bar_cols[extent_states == "others"] = "#D81B60"
names(bar_cols)=mcc$tip.label

#### plotting map
tiff("6_graphs/dec_mcc_hvols.tiff", units="in", width=4, height=6, res=600)
plotTree.wBars(ladderize(mcc),hvolumes,type="phylogram",fsize=0.5, col=bar_cols, lmethod="plotTree")
tiplabels(pie=tip_states_probs, piecol=state_cols, cex=0.3)
nodelabels(node=(1+n_tips):(n_tips+n_nodes), pie= inner_node_probs, piecol=state_cols, cex=0.6)
axisPhylo(pos=c(0.5), font=3, cex.axis=0.5)
dev.off()

########################## ploting observed hvolumes & model estimates ###########################

### best-estimates
best_estimates_table = read.table("5_comparative_analysis/best_estimates_table.csv", sep=",", h=T)

# organizing best-fit estimates
af = best_estimates_table[-46,c(1,2,7,8)]
afot = best_estimates_table[-46,c(3,4,9,10)]
ot = best_estimates_table[-46,c(5,6,11,12)]
colnames(af) = colnames(afot) = colnames(ot) = c("theta", "se", "alpha", "sigma_sq")
state = c( rep("AF", nrow(af)), rep("AFothers", nrow(afot)), rep("Others", nrow(ot) ) )
model_estimates = data.frame(state, rbind(af, afot, ot))

### BM simulations
bm_simulations = read.table("5_comparative_analysis/bm_simulations.csv", sep=",", h=T)

# organizing bm simulations
model = rep("bm1",nrow(bm_simulations))
bm_simulations = data.frame(model, bm_simulations)
colnames(bm_simulations)[2] = "hvolume"

## observed hvolumes
regimes_hvolumes = read.table("3_hypervolume_inference/regimes_hvolumes.csv", sep=",", h=T)

# summarizing hvolumes by state
means=aggregate(regimes_hvolumes$hvolume, by=list(regimes_hvolumes$state),mean)
sds=aggregate(regimes_hvolumes$hvolume, by=list(regimes_hvolumes$state),sd)
ns=aggregate(regimes_hvolumes$hvolume, by=list(regimes_hvolumes$state),length)
ses=sds[,2]/sqrt(ns[,2])
summary_hvolume=data.frame(means,sds[,2],ses)
colnames(summary_hvolume)=c("state","hvolume","sd","se")

theta_plot= ggplot(data= model_estimates, aes(x=state, y=theta, fill=state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.10), size = 1, alpha = 0.10) +
  geom_boxplot(width = 0.30, outlier.shape = NA, alpha = 0.25)+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("")+ ylab("theta")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFothers" = "AF non-endemic", "others" = "Others"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=6),legend.position = "none")

sigma_plot= ggplot(data= model_estimates, aes(x=state, y=sqrt(sigma_sq), fill=state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.10), size = 1, alpha = 0.10) +
  geom_boxplot(width = 0.30, outlier.shape = NA, alpha = 0.25)+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("geographic state")+ ylab("sigma")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFothers" = "AF non-endemic", "others" = "Others"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=6),legend.position = "none")

hvolume_plot= ggplot(data=bm_simulations,aes(x=hvolume)) +
  geom_density(alpha=0.5,color="lightgray", fill="lightgray")+
  scale_colour_manual(values=mycols)+
  scale_fill_manual(values=mycols)+
  geom_vline(data=summary_hvolume,aes(xintercept = hvolume),linetype="solid",colour=c(mycols),size=0.75)+
  geom_vline(data=summary_hvolume,aes(xintercept = hvolume-se),linetype="dotted",colour=c(mycols),size=0.45)+
  geom_vline(data=summary_hvolume,aes(xintercept = hvolume+se),linetype="dotted",colour=c(mycols),size=0.45)+
  labs(x="mean \n hypervolume", y="density")+
  xlim(c(0,25))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=12),legend.position = "none")

tiff("6_graphs/evo_model_estimates.tiff", units="in", width=2.5, height=6, res=600)
ggarrange(theta_plot,sigma_plot,hvolume_plot, nrow=3,ncol=1)
dev.off()

############################# plotting ancestral species data ###########################

### my colors
mycols = c( "#1E88E5", "#FFC107", "#D81B60")

### ancestral data
anc_range_data = read.table("5_comparative_analysis/anc_range_data_OUMV.csv", sep=",", h=T)

# time boundaries
old_age = round(min(anc_range_data$anc_node_ages), 0)
new_age = round(max(anc_range_data$anc_node_ages), 0)

# dividing time in intervals
intervals = anc_range_data$anc_node_ages
breaks = seq(new_age, old_age, by= (old_age - new_age)/20)
for (i in 1:length(breaks)){
  intervals[which(anc_range_data$anc_node_ages > breaks[i+1] & anc_range_data$anc_node_ages < breaks[i])] = breaks[i+1] 
}
anc_range_data = data.frame(anc_range_data, intervals)

# cleaning biologically unreasoning entries
summary(anc_range_data$anc_range)
hist(anc_range_data$anc_range)

anc_range_data = anc_range_data[anc_range_data$anc_range > 0 & anc_range_data$anc_range < 2000,] 

tiff("6_graphs/range_size_by_time_OUMV.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= anc_range_data, aes(x=intervals, y=anc_range, fill=state) ) +
  #geom_point(aes(color=state),size = 1, alpha = 0.10) +
  geom_smooth(se= T, aes(color=state, fill=state), alpha=0.25, size=1.5)+
  geom_vline(xintercept = -2.58, linetype="dashed",color = "black", size=0.75) +
  scale_colour_manual(values=mycols)+
  scale_fill_manual(values=mycols)+
  #ylim(c(100,300))+
  xlab("time before present (m.y.a)")+ ylab("range size (km2)")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10), legend.position = "none")
dev.off()

############################### SUPPLEMENTARY MATERIAL ###########################

############################# Topology evaluation plot ##############################

tiff("0_data/topology_dissimilarity.tiff", units="in", width=4, height=3, res=600)
ggplot(data = df_rf,aes(x=Axis.1,y=Axis.2,color=locus)) +
  geom_point(size = 1, alpha = 0.5) +
  scale_colour_manual(values=c("red4","orange2","blue4","cyan4","red2","hotpink1"))+
  labs(x="PCo1", y="PCo2")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=13,face="bold"))
dev.off()

################################ Rao's Q map ###############################

# wd
setwd("C:/Users/eduar/Documents/GitHub/specialization_endemism_AF")

library(raster)

# loading rao raster
rao_env = raster("2_environmental_heterogeneity/rao_rasters/rao_env_3")

# getting rid of outliers
raster::cellStats(rao_env, stat=hist)
rao_ras[rao_ras@data@values > 1000] = 1000

# graphical parameters
col_func = colorRampPalette(c("steelblue1","green4", "orange1", "red3")) 
legend_colors= col_func(20) 

# plotting
tiff("6_graphs/supp_q_map/rao_q_map.tiff", units="in", width=5.5, height=6, res=600)
plot(rao_ras, col=legend_colors , colNA="white")
dev.off()

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
spp_geographic_state = read.table("1_geographic_classification/spp_geographic_state.csv", sep=',', h=T)

# loading rao raster
rao_ras= raster("2_environmental_heterogeneity/rao_rasters/rao_ras_3")

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
names(mycols) = c("AF", "AFothers", "others")

###### plotting 
for(one_state in c("AF", "AFothers", "others")){
  # sorting by geographic state
  sub_state = subset(spp_geographic_state, state == one_state)
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
      one_col = mycols[names(mycols)== one_state]
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
