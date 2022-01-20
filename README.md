# specialization_endemism_AF
The project encompass three R-scripts:

1) environ_heterogeneity, which: 
a) Quantifies georeferrenced points inside a given geographic shape file;
b) Quantifies environmental heterogeneity based on georeferrenced points projected onto raster layers;
c) Tests differences of environmental heterogeneity based on permutation procedure and PGLS model;

2) hypervolume_inference, which:
a) Sets a function to validate hypervolume models based on k-fold and leave-one-out procedures;
b) Infers hypervolumes and environmental niches of species based on hypervolume models;
c) Displays plots contrasting species' hypervolume;

3) comparative_hypervolume, which:
a) Tests differences of hypervolume between species based on PGLS model;
b) Quantifies phylogenetic signal of species' hypervolume;
c) Sets functions to fit multiple evolutionary models and to choose the best-fit model based on AICc;
d) Infers the best-fit evolutionary model of species'hypervolume based on multiple phylogenetic trees;
