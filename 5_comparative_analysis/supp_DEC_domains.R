### setting working directory
setwd("C:/Users/eduar/Documents/GitHub/specialization_endemism_AF")

### packages
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)
library(ape)

# reading range data
geog_fn = ("1_geographic_classification/spp_distribution_domains.data")
moref(geog_fn)

# converting phylip format to tipranges
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geog_fn)
tipranges

# setting maximum number of areas occupied for reconstructions
max_range_size = max(rowSums(dfnums_to_numeric(tipranges@df)))

# Initialize DEC model
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# location of the geography text file
BioGeoBEARS_run_object$geogfn = geog_fn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

# Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$min_branchlength = 0.001    

# set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$include_null_range = FALSE    

# computing options
BioGeoBEARS_run_object$num_cores_to_use = 1

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

### mcc tree
mcc_fn = c("0_data/mcc_phylo.nwk")
mcc = read.tree(mcc_fn)
BioGeoBEARS_run_object$trfn = mcc_fn

# fitting DEC
res_DEC = bears_optim_run(BioGeoBEARS_run_object)

# ploting 
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
plot_BioGeoBEARS_results(res_DEC, plotwhat = "pie", label.offset = 0.5, tipcex = 0.5,
                         statecex = 0.7,  plotsplits = F, cornercoords_loc = scriptdir, 
                         include_null_range = F, tr = mcc, tipranges = tipranges)

# node marginal ML
relprobs_matrix = res_DEC$ML_marginal_prob_each_state_at_branch_top_AT_node

# getting area names
areas = getareas_from_tipranges_object(tipranges)
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=F)

# make the list of ranges
ranges_list = NULL
for (i in 1:length(states_list_0based))
{    
  if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
  {
    tmprange = "_"
  } else {
    tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
  }
  ranges_list = c(ranges_list, tmprange)
}

# Look at the ranges list
ranges_list

# naming marginal ML matrix
colnames(relprobs_matrix) = ranges_list

write.table(relprobs_matrix, "5_comparative_analysis/DEC_node_domain_probs.csv", sep=",", row.names=F, quote=F)

