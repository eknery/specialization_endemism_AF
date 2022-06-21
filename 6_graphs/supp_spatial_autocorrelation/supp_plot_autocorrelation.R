
# setting working directory
setwd("C:/Users/eduar/Documents/GitHub/specialization_endemism_AF")

# loading packages
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

# loading data
moran_df = read.table("2_environmental_heterogeneity/moran_df.csv",sep=",",h=T)

####################### plotiing Moran for each variable #######################
# my colors
mycols = c( "#1E88E5", "#D81B60")

plot_list = list()
for (i in 2:5){
  var_name = colnames(moran_df)[i]
  q_plot= ggplot(data= moran_df, aes(x=domain, y=moran_df[,i], fill=domain)) +
    #geom_point(aes(color=distribution),position = position_jitter(width = 0.07), size = 2, alpha = 0.65)+
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.25)+
    geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.25) +
    scale_fill_manual(values=mycols)+
    scale_colour_manual(values=mycols)+
    xlab("geographic\ndistribution")+ ylab(var_name)+
    scale_x_discrete(labels=c("AF" = "AF", "out" = "outside AF"))+
    theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=13,face="bold"),axis.text=element_text(size=10),legend.position = "none")
  plot_list[[i-1]] = q_plot
}

tiff("6_graphs/supp_spatial_autocorrelation/supp_autocorrelation.tiff", units="in", width=6, height=6, res=600)
ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]], plot_list[[4]], nrow=2,ncol=2)
dev.off()