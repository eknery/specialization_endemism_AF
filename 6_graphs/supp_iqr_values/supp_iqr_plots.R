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

### iqr values 
iqr_env_values = read.table("3_hypervolume_inference/iqr_env_values.csv", sep=",", h=T)

# scalling ph
iqr_env_values$soil_pH = iqr_env_values$soil_pH*10^-1

###loading spp geographic classification
spp_geographic_distribution = read.table("1_geographic_classification/spp_geographic_distribution.csv", sep=',', h=T)

# my colors
mycols = c( "#1E88E5", "#FFC107", "#D81B60")

### looping over variable
plot_list = list()
for (i in 2:5){
  one_env = data.frame(iqr_env_values[,c(1,i)], spp_geographic_distribution$state)
  trait_name = colnames(one_env)[2]
  colnames(one_env)[c(2,3)] = c("trait", "state")
  iqr_plot= ggplot(data= one_env, aes(x=state, y=trait, fill=state)) +
    geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 2, alpha = 0.65) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.25)+
    geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.25) +
    scale_fill_manual(values=mycols)+
    scale_colour_manual(values=mycols)+
    xlab("geographic distribution")+ ylab(trait_name)+
    scale_x_discrete(labels=c("AF" = "AF-endemic", "AFother" = "AF and other\ndomains", "other" = "outside AF"))+
    theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=13,face="bold"),axis.text=element_text(size=6),legend.position = "none")
  plot_list[[i-1]] = iqr_plot
}

tiff("6_graphs/supp_iqr_values/iqr_values.tiff", units="in", width=6, height=6, res=600)
  ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]], plot_list[[4]], nrow=2,ncol=2)
dev.off()