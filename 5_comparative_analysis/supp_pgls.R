### setting working directory
setwd("C:/Users/eduar/Documents/GitHub/specialization_endemism_AF")

### packages
library(ape)
library(phytools)
library(geiger)
library(nlme)

### spp geographic distribution
spp_geographic_distribution = read.table("1_geographic_classification/spp_geographic_distribution.csv", sep=',', h=T)
# named geographic states
state = spp_geographic_distribution$state
names(state) = spp_geographic_distribution$species

### iqr values 
iqr_env_values = read.table("3_hypervolume_inference/iqr_env_values.csv", sep=",", h=T)

### mcc tree
mcc_fn = c("0_data/mcc_phylo.nwk")
mcc = read.tree(mcc_fn)

###################################### PGLS #########################

colnames(iqr_env_values)

one_env = iqr_env_values[,4]
names(one_env) = iqr_env_values$species

# comparing models
fit_bm_hv = fitContinuous(phy= mcc, one_env,  model="BM")
fit_ou_hv = fitContinuous(phy= mcc, one_env,  model="OU")
# aicc
c(fit_bm_hv$opt$aicc, fit_ou_hv$opt$aicc) 
# correlation matrix
cor_hv = corMartins(fit_ou_hv$opt$alpha, phy = mcc, fixed = F)
# pgls
fit_gls_hv = gls(one_env ~ state, correlation=cor_hv,  method = "REML")
summary(fit_gls_hv)
plot(fit_gls_hv)
#checking residuals
res_hv = resid(fit_gls_hv)[1:66]
hist(res_hv)
shapiro.test(res_hv)

aggregate(one_env, by=list(state), median)
boxplot(one_env ~state)

(3.15 -5.95)/5.95
