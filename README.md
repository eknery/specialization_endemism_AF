# specialization_endemism_AF

### Description:
This projects aimed to compare ecological specialization ('rarity') inside and outside the Atlantic Forest and explore their causes

### Methods:
The dataset includes georeferenced occurences (~9.000) retrieved from [splink](http://www.splink.org.br) database and raster layers from [WorldClim](https://www.worldclim.org/). Ecological estimates were based on the hypervolume method and evolutionary inferences were carried out using the OUwie R package.
The analytical framework is based on three main R scripts:

1. environmental_heterogeneity.R: 
    1. Quantifies environmental heterogeneity based on georeferrenced occurrences projected onto raster layers;
    2. Tests differences of environmental heterogeneity based on a permutation procedure and a generalized least-square (GLS) model;

2. hypervolume_inference.R:
    1. Sets a function to validate hypervolume models based on TSS under k-fold and leave-one-out procedures;
    2. Infers hypervolumes and environmental niches of species based on hypervolume models;
    3. Displays plots contrasting species' hypervolume;

3. comparative_analysis.R:
    1. Tests differences of hypervolume size inside and outside the Atlantic Forest based on a GLS model;
    2. Sets functions to fit multiple evolutionary models and to choose the best-fit model based on the corrected Akaike Information Criterion (AICs);
    3. Retrieve parameter estimates from the best-fit evolutionary model and compare estimates across groups;

If you find any of these scripts useful, please refer to my [article](https://academic.oup.com/aob/advance-article-abstract/doi/10.1093/aob/mcad029/7033295?redirectedFrom=fulltext)
