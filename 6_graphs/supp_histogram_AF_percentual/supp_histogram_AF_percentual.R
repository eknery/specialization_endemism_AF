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
library(colorspace)

#loading species occurrences
spp_count_domain=read.table("1_geographic_classification/spp_count_domain.csv", header =T, sep=",")

# total count per species
total_count = apply(spp_count_domain[,-1], MARGIN = 1, sum)

# percentual occurrences in the AF
percent_af = round(spp_count_domain$AF/total_count, 2)
spp_percent_af = data.frame(spp_count_domain$species, percent_af)
colnames(spp_percent_af)[1] = "species"

tiff("6_graphs/supp_histogram_AF_percentual/histogram_AF_percentual.tiff", units="in", width=6, height=6, res=600)
ggplot(data=spp_percent_af, aes(percent_af)) + 
  geom_histogram()+
  xlab("percentual occurrence in the AF")+ ylab("number of species")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=13,face="bold"),axis.text=element_text(size=10),legend.position = "none")
dev.off()

