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

### rao values 
spp_rao = read.table("2_environmental_heterogeneity/spp_rao_per_variable.csv",sep=",", h=T)

# my colors
mycols = c( "#1E88E5", "#FFC107", "#D81B60")

### looping over variable
plot_list = list()
for (i in 1:4){
  trait_name = colnames(spp_rao)[2+i]
  one_env = data.frame(spp_rao[,c(1,2)], spp_rao[,2+i])
  colnames(one_env)[3] = "trait"
  q_plot= ggplot(data= one_env, aes(x=state, y=trait, fill=state)) +
    geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 2, alpha = 0.65) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.25)+
    geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.25) +
    scale_fill_manual(values=mycols)+
    scale_colour_manual(values=mycols)+
    xlab("geographic distribution")+ ylab(trait_name)+
    scale_x_discrete(labels=c("AF" = "AF-endemic", "AFother" = "AF and other\ndomains", "other" = "outside AF"))+
    theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=10,face="bold"),axis.text=element_text(size=6),legend.position = "none")
  plot_list[[i]] = q_plot
}

tiff("6_graphs/supp_qvalues_per_variables/supp_qvalues_per_variables.tiff", units="in", width=6, height=6, res=600)
  ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]], plot_list[[4]], nrow=2,ncol=2)
dev.off()