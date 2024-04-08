# specialization_endemism_AF

### Description:
This is the repositoty for my article [Ecological specialization underlies plant endemism](https://academic.oup.com/aob/advance-article-abstract/doi/10.1093/aob/mcad029/7033295?redirectedFrom=fulltext). This projects aimed to test the relationship between ecological specialization and plant endemism in the Atlantic Forest. The project utilized the Miconia supersect. Discolor as a study system. 

### Methods:
The dataset includes phylogenetic trees, already-curated georeferenced occurences (~9.000), and raster layers. Estimates of ecological specialiazation were based on the hypervolume method. Evolutionary inferences were carried out using the OUwie R package.

The analytical framework is based on three main R scripts:

1. environmental_heterogeneity.R: 
    1. Quantifies environmental heterogeneity based on georeferrenced occurrences projected onto raster layers;
    2. Tests differences of environmental heterogeneity among geographic distributions based on a generalized least-square (GLS) model;

2. hypervolume_inference.R:
    1. Sets a function to validate hypervolume models based on TSS under k-fold and leave-one-out procedures;
    2. Infers the environmental niche of species based on their environmental values applied to hypervolume models;
    3. Constrasts species' hypervolume by graphical plots;

3. range_inference.R:
    1. Measures species' geographic range based on the convex hull of their georeferrenced occurrences;
    2.  Constrasts species' geographic range by graphical plots;

4. comparative_analyses.R:
    1. Tests differences of species' hypervolume and geographic range between gepgraphic groups based on GLS models;
    2. Fits multiple models and choose the best-fit model describing the evolution of hypervolume and geopgraphic range based on the corrected Akaike Information Criterion (AIC);
    3. Retrieve parameter estimates from the best-fit evolutionary model and compare estimates across groups;
