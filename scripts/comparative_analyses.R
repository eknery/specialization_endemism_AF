
### packages
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)
library(ape)
library(phytools)
library(geiger)
library(OUwie)
library(nlme)

### my own functions
source("./5_comparative_analysis/function_fit_evo_models.R")
source("./5_comparative_analysis/function_choose_best.R")

### number of trees
n_trees = length (list.files("0_data/100_trees"))

### mcc tree
mcc_fn = c("0_data/mcc_phylo.nwk")
mcc = read.tree(mcc_fn)

### spp geographic distribution
spp_geographic_distribution = read.table("1_geographic_classification/spp_geographic_distribution.csv", sep=',', h=T)
# named geographic states
state = spp_geographic_distribution$state
names(state) = spp_geographic_distribution$species
  
### reading rao values
spp_rao = read.table("2_environmental_heterogeneity/spp_rao.csv", header =T, sep=",",  na.strings = "NA", fill=T)
#named rao
rao = spp_rao$rao
names(rao) = spp_rao$species

### reading altitude values
spp_altitude = read.table("2_environmental_heterogeneity/spp_altitude.csv", header =T, sep=",",  na.strings = "NA", fill=T)
altitude = spp_altitude$altitude
names(altitude) = spp_altitude$species

### reading hvolume data
regimes_hvolumes = read.table("3_hypervolume_inference/regimes_hvolumes.csv", sep=",", h=T)
# named hvolumes
hvolume = regimes_hvolumes$hvolume
names(hvolume) = regimes_hvolumes$species

### reading range data
regimes_ranges = read.table("4_range_inference/regimes_ranges.csv", header =T, sep=",",  na.strings = "NA", fill=T)
# named ranges
range = regimes_ranges$range
names(range) = regimes_ranges$species

### reading autocorrelation data
spp_moran = read.table("4_range_inference/spp_moran.csv", header =T, sep=",",  na.strings = "NA", fill=T)
# named autocorrelation
moran = apply(spp_moran[,2:5], MARGIN = 1, mean)
names(moran) = spp_moran$species

######################### phylo GLS  ######################

### rao
# comparing models
fit_bm_rao = fitContinuous(phy= mcc, rao,  model="BM")
fit_ou_rao = fitContinuous(phy= mcc, rao,  model="OU")
# aicc
c(fit_bm_rao$opt$aicc, fit_ou_rao$opt$aicc) 
# correlation structure
cor_rao = corMartins(0.31, phy = mcc, fixed = T)

# pgls rao vs geographcic state
fit_gls_rao = gls(log(rao) ~ state, correlation=cor_rao,  method = "REML")
summary(fit_gls_rao)
plot(fit_gls_rao)
#checking residuals
res_rao = resid(fit_gls_rao)[1:66]
hist(res_rao)
shapiro.test(res_rao)


# pgls rao vs altitude
km_altitude =(altitude*10^-3)
fit_gls_rao = gls(log(rao) ~ km_altitude, correlation=cor_rao,  method = "REML")
summary(fit_gls_rao)
plot(fit_gls_rao)
#checking residuals
res_rao = resid(fit_gls_rao)[1:66]
hist(res_rao)
shapiro.test(res_rao)

ssr = sum( (resid(fit_gls_rao)[1:66])^2 )
sst =  sum((log(rao)- mean(log(rao)))^2)
(r2 = 1 - (ssr/sst))


### hvolumes
# comparing models
fit_bm_hv = fitContinuous(phy= mcc, hvolume,  model="BM")
fit_ou_hv = fitContinuous(phy= mcc, hvolume,  model="OU")
# aicc
c(fit_bm_hv$opt$aicc, fit_ou_hv$opt$aicc) 
# correlation matrix
cor_hv = corMartins(2.71, phy = mcc, fixed = F)

# pgls hvolumes ~ geographic distribution
fit_gls_hv = gls(sqrt(hvolume) ~ state + moran, correlation=cor_hv,  method = "REML")
summary(fit_gls_hv)
plot(fit_gls_hv)
#checking residuals
res_hv = resid(fit_gls_hv)[1:66]
hist(res_hv)
shapiro.test(res_hv)

### range size  
# comparing models
fit_bm_range = fitContinuous(phy= mcc, range,  model="BM")
fit_ou_range = fitContinuous(phy= mcc, range,  model="OU")
c(fit_bm_range$opt$aicc, fit_ou_range$opt$aicc) 
# correlation structure
cor_ran = corMartins(2.71, phy = mcc, fixed = F)

# pgls range size vs hvolume 
cubic_range = range^(1/3)
fit_gls_ran = gls(cubic_range ~ state/hvolume , correlation=cor_ran,  method = "REML")
summary(fit_gls_ran)
plot(fit_gls_ran)
#checking residuals
res_ran = resid(fit_gls_ran)[1:66]
hist(res_ran)
shapiro.test(res_ran)

ssr = sum( (resid(fit_gls_ran)[1:66])^2 )
sst =  sum( (cubic_range - mean(cubic_range) )^2 )
r2 = 1 - (ssr/sst)
r2

########################## preparing DEC  #############################
# reading range data
geog_fn = ("1_geographic_classification/spp_distribution_af.data")
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

#### DEC over trees 
for (i in 1:n_trees){ 
  # phylogeny tree location
  trfn = paste("0_data/100_trees/tree", as.character(i), sep="_")
  tr = read.tree(trfn)
  # n tips and nodes
  n_tips = Ntip(tr)
  n_nodes = tr$Nnode
  # inputting tree into DEC
  BioGeoBEARS_run_object$trfn = trfn
  # fitting DEC
  res_DEC = bears_optim_run(BioGeoBEARS_run_object)
  # node marginal ML
  relprobs_matrix = res_DEC$ML_marginal_prob_each_state_at_branch_top_AT_node
  # node states
  state_labels=c("AF", "other", "AFother")
  node_states = get_ML_states_from_relprobs(relprobs=relprobs_matrix, statenames=state_labels, returnwhat = "states", if_ties = "takefirst")
  # getting only ancestral nodes
  anc_node = (n_tips+1):(n_tips+n_nodes)
  state = node_states[anc_node]
  anc_node_states = data.frame(anc_node, state)
  write.table(anc_node_states , paste("5_comparative_analysis/DEC_ancestral_reconstructions/anc_node_states",as.character(i), sep="_"), sep=",", row.names=F, quote=F)
}

######### fitting evolutionary models to traits over trees #########################

spp_trait_regimes = regimes_ranges

all_best_models = data.frame(matrix(NA, nrow=n_trees, ncol=4))
colnames(all_best_models) = c(c("model","llik","aicc","delta_aicc"))
all_best_estimates = vector("list" , n_trees)
all_models = c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA")

for (i in 1:n_trees){ 
  # phylogenetic tree
  tr_fn = paste("0_data/100_trees/tree", as.character(i), sep="_")
  tr = read.tree(tr_fn)
  # DEC ancestral states
  dec_fn = paste("5_comparative_analysis/DEC_ancestral_reconstructions/anc_node_states",as.character(i), sep="_")
  anc_node_states = read.table(dec_fn, sep=",", h=T)
  tr$node.label = anc_node_states$state
  # fitting evolutionary models
  all_fits = fit_evo_models(tree=tr, regimes=spp_trait_regimes, models_to_fit = all_models)
  best_choice = choose_best(all_fits)
  all_best_models[i,] = best_choice$best_fit
  all_best_estimates[[i]] = best_choice$best_estimates
  print(i)
}

best_estimates =c()
for (i in 1:length(all_best_estimates)){
  best_estimates = rbind(best_estimates,all_best_estimates[[i]])
}

write.table(all_best_models,"5_comparative_analysis/range_best_models.csv",sep=",",quote=F,row.names=F)
write.table(best_estimates,"5_comparative_analysis/range_best_estimates.csv",sep=",",quote=F,row.names=F)

########################### fitting and simulating BM evolution ##################

### which one to estimate?
spp_trait_regimes  = regimes_ranges
  
### fitting BM model
bm_fits=vector(mode = "list", length = n_trees)
for (i in 1:n_trees){
  # phylogenetic tree
  trfn = paste("0_data/100_trees/tree", as.character(i), sep="_")
  tr = read.tree(trfn)
  # DEC ancestral states
  decfn = paste("5_comparative_analysis/DEC_ancestral_reconstructions/anc_node_states",as.character(i), sep="_")
  anc_node_states = read.table(decfn, sep=",", h=T)
  tr$node.label = anc_node_states$state
  # fitting BM1 model
  bm_fits[[i]] = fit_evo_models(tree=tr, regimes=spp_trait_regimes , models_to_fit = "BM1")
  print(i)
}

### taking BM1 estimates
bm_estimates=as.data.frame(matrix(NA,nrow=length(bm_fits),ncol=5))
colnames(bm_estimates)=c("model","theta","se","alpha","sigma_sq")
for (i in 1:length(bm_fits)){
  bm_estimates[i,1]=bm_fits[[1]]$fit_metrics$model
  bm_estimates[i,2:3]=bm_fits[[i]]$model_estimates$BM1$theta[1,] # optima & se
  bm_estimates[i,4:5]=bm_fits[[i]]$model_estimates$BM1$solution[,1] # alpha & sigma_sq
}

write.table(bm_estimates,"5_comparative_analysis/range_bm_estimates.csv", sep=",", quote=F, row.names=F)

### simulating BM1 evolution
bm_simulations=c()
for (i in 1:n_trees){
  # phylogenetic tree
  trfn = paste("0_data/100_trees/tree", as.character(i), sep="_")
  tr = read.tree(trfn)
  for (j in 1:100){
    sim_hvolumes = fastBM(tr, a= sample(spp_trait_regimes[,3], 1), sig2=bm_estimates$sigma_sq[i], bounds=range(spp_trait_regimes[,3])*3, internal=F)
    mean_sim_hvolumes= mean(sim_hvolumes)
    bm_simulations = c(bm_simulations, mean_sim_hvolumes)
  }
}


write.table(bm_simulations,"5_comparative_analysis/range_bm_simulations.csv", sep=",", quote=F, row.names=F)


########################### using best-fit models to get ancestral states ###############

### which one to estimate?
spp_trait_regimes  = regimes_ranges

all_best_models = read.table("5_comparative_analysis/range_best_models.csv",sep=",", h=T)
models_to_use = all_best_models$model

for (i in 1:n_trees){ 
  # phylogenetic tree
  tr_fn = paste("0_data/100_trees/tree", as.character(i), sep="_")
  tr = read.tree(tr_fn)
  # DEC ancestral states
  dec_fn = paste("5_comparative_analysis/DEC_ancestral_reconstructions/anc_node_states",as.character(i), sep="_")
  anc_node_states = read.table(dec_fn, sep=",", h=T)
  tr$node.label = anc_node_states$state
  # n tips and nodes
  n_tips = Ntip(tr)
  n_nodes = tr$Nnode
  # ancestral node numbers
  anc_node = (n_tips+1):(n_tips+n_nodes)
  # fitting best-model
  fit=OUwie(phy=tr, data=spp_trait_regimes, model=models_to_use[i], lb=0, ub=Inf) 
  # getting ancestral traits
  ou_anc = OUwie.anc(fit, knowledge=T)
  anc_trait = ou_anc$NodeRecon
  # into dataframe
  anc_data = data.frame(anc_node, anc_trait)
  # exporting
  write.table(anc_data , paste("5_comparative_analysis/OUwie_ancestral_range/anc_range",as.character(i), sep="_"), sep=",", row.names=F, quote=F)
}

### summarizing ancestral trait into dataframe
anc_data= c(0,0,0,0)
for (i in 1:n_trees){ 
  # phylogenetic tree
  trfn = paste("0_data/100_trees/tree", as.character(i), sep="_")
  tr = read.tree(trfn)
  # n tips and nodes
  n_tips = Ntip(tr)
  n_nodes = tr$Nnode
  # DEC ancestral states
  dec_fn = paste("5_comparative_analysis/DEC_ancestral_reconstructions/anc_node_states",as.character(i), sep="_")
  anc_node_states = read.table(dec_fn, sep=",", h=T)
  # ancestral trait
  data_fn = paste( "5_comparative_analysis/OUwie_ancestral_range/anc_range" ,as.character(i), sep="_")
  anc_trait = read.table(data_fn, sep=",", h=T)
  # node ages
  node_ages = round(node.depth.edgelength(tr),5)
  present = round(max(node_ages), 5)
  node_ages = node_ages - present
  # ancestral node ages
  anc_node_ages = node_ages[(n_tips+1):(n_tips+n_nodes)]
  # into a single df
  one_set = data.frame(anc_node_states,anc_node_ages, anc_trait[,2])
  anc_data= rbind(anc_data, one_set)
}

colnames(anc_data)[4] = "anc_range"
anc_data= anc_data[-1,]
write.table(anc_data, "5_comparative_analysis/anc_range_data_OUMV.csv", sep=",", row.names=F, quote=F)

############################## DEC over mcc only #################################

### mcc tree
mcc_fn = c("0_data/mcc_phylo.nwk")
mcc = read.tree(mcc_fn)
BioGeoBEARS_run_object$trfn = mcc_fn

# fitting DEC
res_DEC = bears_optim_run(BioGeoBEARS_run_object)

# plotting 
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
plot_BioGeoBEARS_results(res_DEC, plotwhat = "pie", label.offset = 0.5, tipcex = 0.5,
                         statecex = 0.7,  plotsplits = F, cornercoords_loc = scriptdir, 
                         include_null_range = F, tr = mcc, tipranges = tipranges)

# node marginal ML
relprobs_matrix = res_DEC$ML_marginal_prob_each_state_at_branch_top_AT_node
node_probs_mcc = relprobs_matrix
colnames(node_probs_mcc) = c("AF", "other", "AFother")

# node states
state_labels=c("AF", "other", "AFother")
node_states_mcc= get_ML_states_from_relprobs(relprobs=relprobs_matrix, statenames=state_labels, returnwhat = "states", if_ties = "takefirst")

write.table(node_probs_mcc, "5_comparative_analysis/DEC_node_probs_mcc.csv", sep=",", quote=F, row.names = F)
write.table(node_states_mcc, "5_comparative_analysis/DEC_node_states_mcc.csv", sep=",", quote=F, row.names = F)

