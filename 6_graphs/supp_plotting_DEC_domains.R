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

library(ape)
############################### plotting phylogenetic data ############################

# mcc tree
mcc = read.tree("0_data/mcc_phylo.nwk")
n_tips = Ntip(mcc)
n_inner_nodes = mcc$Nnode

# DEC node probabilities
DEC_node_domain_probs = read.table("5_comparative_analysis/DEC_node_domain_probs.csv",sep=",", h=T)

### visual parameters
# tip states
extent_states = c()
for( i in 1:n_tips){
  one_state = colnames(DEC_node_domain_probs)[DEC_node_domain_probs[i,] == 1]
  extent_states = c(extent_states, one_state)
}
names(extent_states) = mcc$tip.label

# tip state probs
tip_states_probs = DEC_node_domain_probs[1:n_tips,]

# ancestral state probs
inner_node_probs = DEC_node_domain_probs[(1+n_tips):(n_tips+n_inner_nodes),]

# state colors
all_states = colnames(DEC_node_domain_probs)

# basic domain colors
domain_colors =  c("#1D23E6", "#4BE61D", "#E6E61D", "#E63C1D", "#E6A61D", "#E61DC4", "#716C70", "#18D4D4")
names(domain_colors) = all_states[1:length(domain_colors)]
n_prime_colors = length(domain_colors)
total_colors = length(all_states)

for (one_state in all_states[(n_prime_colors+1): total_colors]){
  chars = strsplit(one_state, split="")[[1]]
  n_chars = length(chars)
  initial_name = chars[1]
  for (i in 2:(n_chars-3)){ initial_name = paste(initial_name, chars[i], sep="")}
  final_name = paste(chars[n_chars-2],chars[n_chars-1], chars[n_chars], sep="")
  x = domain_colors[initial_name]
  name_x = names(domain_colors[initial_name])
  y = domain_colors[final_name]
  name_y = names(domain_colors[final_name])
  col1 = as.vector(col2rgb(x))
  col2 = as.vector(col2rgb(y))
  mix_rgb = mixcolor(0.5, sRGB(col1[1], col1[2], col1[3]), sRGB(col2[1], col2[2], col2[3]))
  mix_hex = rgb(mix_rgb@coords[1], mix_rgb@coords[2], mix_rgb@coords[3], maxColorValue=255)
  mix_name = paste(name_x,name_y, sep="")
  names(mix_hex) = mix_name
  domain_colors = c(domain_colors, mix_hex)
}
domain_colors = domain_colors[all_states]


### plotting DEC
tiff("6_graphs/supp_DEC_domains/DEC_tree.tiff", units="in", width=5, height=8, res=600)
plot(ladderize(mcc),type="phylogram", cex=0.5, label.offset= 0.15)
tiplabels(pie=tip_states_probs,  cex=0.3,  piecol=domain_colors)
nodelabels(node=(1+n_tips):(n_tips+n_inner_nodes), pie= inner_node_probs, piecol=domain_colors, cex=1)
axisPhylo(pos=c(0.1), font=2, cex.axis=0.5)
dev.off()
