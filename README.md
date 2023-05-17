# specialization_endemism_AF

### Description:
This is the repositoty for my article [Ecological specialization underlies plant endemism](https://academic.oup.com/aob/advance-article-abstract/doi/10.1093/aob/mcad029/7033295?redirectedFrom=fulltext). This projects aimed to explore the relationship between ecological specialization and plant endemism in the Atlantic Forest. The project utilized the Miconia supersect. Discolor as a study system. 

### Methods:
The dataset includes phylogenetic trees, already-curated georeferenced occurences (~9.000), and raster layers. Estimates of ecological specialiazation were based on the hypervolume method. Evolutionary inferences were carried out using the OUwie R package.

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
