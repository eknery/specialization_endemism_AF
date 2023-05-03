# specialization_endemism_AF

### Description:
This projects aimed to compare ecological specialization ('rarity') inside and outside the Atlantic Forest and explore their causes

### Methods:
The dataset includes georeferenced occurences (~9.000) retrieved from [splink](http://www.splink.org.br) database and raster layers from [WorldClim](https://www.worldclim.org/). Ecological estimates were based on the hypervolume method and evolutionary inferences were carried out using the OUwie R package.
The analytical framework is based on three main R scripts:

1. environmental_heterogeneity.R: 
i. Quantifies environmental heterogeneity based on georeferrenced occurrences projected onto raster layers;
ii. Tests differences of environmental heterogeneity based on a permutation procedure and a generalized least-square (GLS) model;

2) hypervolume_inference.R:
a) Sets a function to validate hypervolume models based on TSS under k-fold and leave-one-out procedures;
b) Infers hypervolumes and environmental niches of species based on hypervolume models;
c) Displays plots contrasting species' hypervolume;

3) comparative_analysis.R:
a) Tests differences of hypervolume size inside and outside the Atlantic Forest based on a GLS model;
b) Sets functions to fit multiple evolutionary models and to choose the best-fit model based on corrected Akaike Information Criterion (AICs);
c) Retrieve parameter estimates from the best-fit evolutionary model;
