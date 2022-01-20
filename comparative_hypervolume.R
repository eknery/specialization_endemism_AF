library(ape)
library(phytools)
library(geiger)
library(OUwie)
library(nlme)

### loading data
spp_hypervol=read.table("spp_hypervol.csv",h=T,sep=",")
spp_env_niche=read.table("spp_env_niche.csv",h=T,sep=",")

trees=read.tree("100_rand_phylos.nwk")
mcc=read.tree("mcc_phylo.nwk")

# defining endemics
AF=c("atlantica","baumgratziana","brunnea","budlejoides","capixaba", 
"castaneiflora","cinerascens","discolor","dura","fasciculata",
"formosa","hyemalis","kollmannii","kriegeriana","lymanii","mellina",
"octopetala","penduliflora","petroniana","polyandra","racemifera",
"robusta","ruschiana","setosociliata","shepherdii","valtheri","willdenowii")

# adding distribution vector
distribution=spp_hypervol$species
distribution[which(distribution %in% AF)]="AF-endemic"
distribution[which(distribution != "AF-endemic")]="widespread"
spp_hypervol = data.frame(distribution, spp_hypervol)

# named vector of volumes
hypervolumes=spp_hypervol$hypervolume
names(hypervolumes)=spp_hypervol$species

# named vector of geographic states
states=spp_hypervol$distribution
names(states)=spp_hypervol$species

# dataframe of sd values
sds=aggregate(spp_env_niche[,-1], by=list(spp_env_niche$species), sd)

sds_dist = sds$Group.1
sds_dist[which(sds_dist %in% AF)]="AF-endemic"
sds_dist[which(sds_dist != "AF-endemic")]="widespread"
sds = data.frame(sds_dist, sds)
colnames(sds)[c(1:2)] = c("distribution", "species")

############################# phylogenetic gls ###################################

## infering correlation matrix
signal_hyper = phylosig(tree=mcc, x=sqrt(hypervolumes), method="lambda")
pagel_cor_hyper = corPagel(signal_hyper$lambda, phy = mcc, form = ~1,  fixed = T)

## pgls
fit_gls_hyper = gls(sqrt(hypervolume) ~ distribution, data=spp_hypervol, pagel_cor_hyper,  method = "ML")
summary(fit_gls_hyper)
plot(fit_gls_hyper)

#checking residuals
res_gls_hyper = resid(fit_gls_hyper)[1:66]
hist(res_gls_hyper)
shapiro.fit_results(res_gls_hyper)

### for sd values
gls_list = vector("list", length(sds[,-c(1:2)]) )
res_gls_list = vector("list", length(sds[,-c(1:2)]) )

for (i in 1:length(sds[,-c(1:2)])){
  sd_log = log(sds[,i+2])
  names(sd_log) = sds$species
  sd_distribution = sds$distribution
  signal_sd = phylosig(tree=mcc, x=sd_log, method="lambda")
  pagel_cor_sd = corPagel(signal_sd$lambda, phy = mcc, form = ~1,  fixed = T)
  fit_gls_sd = gls( sd_log ~ sd_distribution, correlation = pagel_cor_sd,  method = "ML")
  gls_list[[i]] = fit_gls_sd
  res_gls_list[[i]] = resid(fit_gls_sd)[1:66]
}

aggregate(sds[,3:6], by= list(sds$distribution), mean)


############################## phylogenetic signal #############################

signal_matrix=data.frame(rep("NA",length(trees)),rep("NA",length(trees)))
colnames(signal_matrix)=c("lambda","p_value")

for(i in 1:100){
	signal=phylosig(tree=trees[[i]], x=hypervolumes, method="lambda", fit_results=T, nsim=1000)
	signal_matrix[i,1]=signal$lambda
	signal_matrix[i,2]=signal$P
	}

signal_matrix$lambda=as.numeric(signal_matrix$lambda)
signal_matrix$p.value=as.numeric(signal_matrix$p.value)

mean(signal_matrix$lambda)
sd(signal_matrix$lambda)

write.table(signal_matrix,"phylo_signal.csv", quote=F, row.names=F)

########################### reconstruction over mcc ###########################

# stochastic maps
simmaps=make.simmap(mcc, x=states, model="ER", nsim=100)
mean_map=summary(simmaps)

#visual parameters
tip_states<-to.matrix(states[mcc$tip.label], c("AF-endemic","widespread"))

state_cols=c("#0000FF","#FF0000")
names(state_cols)=c("AF-endemic","widespread")

bar_cols=mcc$tip.label
bar_cols[bar_cols %in% AF]= "#0000FF7F"
bar_cols[bar_cols!= "#0000FF7F"]= "#FF0000"
names(bar_cols)=mcc$tip.label

# plotting map
tiff("simmap_mcc_vols.tiff", units="in", width=4, height=6, res=600)
plotTree.wBars(ladderize(mcc),hypervolumes,type="phylogram",fsize=0.5,col=bar_cols,lmethod="plotTree")
tiplabels(pie=tip_states,piecol=state_cols,cex=0.4)
nodelabels(node=1:mcc$Nnode+Ntip(mcc), pie=mean_map$ace,piecol=state_cols,cex=0.6)
dev.off()

###################### fitting evolutionary models ############################

### setting function to fit all models
fit_evo_models = function(tree, regimes, models_to_fit){
	states = regimes[,2]
	names(states) = regimes[,1]
  # reconstructing ancestral states
  simmaps=make.simmap(tree, x=states, model="ER", nsim=100)
	mean_map=summary(simmaps)
	anc_states= rep(NA, length(mean_map$ace[,1]))
	for (i in 1:length(mean_map$ace[,1])){
  	bool = mean_map$ace[i,] == max(mean_map$ace[i,])
  	high_prob_state = names(bool)[bool]
  	anc_states[i] = high_prob_state
	}
  tree$node.label=anc_states
  #setting fitting tables
  model_fit_table = data.frame(matrix(NA, nrow= length(models_to_fit), ncol=3))
  colnames(model_fit_table) = c("model","llik","aicc")
  #setting estimate tables
  model_estimate_list = vector("list", length(models_to_fit))
  for (i in 1:length(models_to_fit)){
    # fitting models
  	fit=OUwie(tree,data=regimes,model=models_to_fit[i])
  	# picking fitting metrics
  	model_fit_table[i,] = c(models_to_fit[i],fit$loglik,fit$AICc)
  	# picking model estimates
  	model_estimate = vector("list", 2)
  	model_estimate[[1]] = fit$theta
  	model_estimate[[2]] = fit$solution
  	names(model_estimate) = c("theta","solution")
  	model_estimate_list[[i]] = model_estimate
  	names(model_estimate_list)[i] = models_to_fit[i]
  }
  fit_results = vector("list", 2)
  fit_results[[1]] = model_fit_table
  fit_results[[2]] = model_estimate_list
  names(fit_results) = c("fit_metrics","model_estimates")
  return(fit_results)
}


### setting function to choose best-fit model and estimates
choose_best = function (fit_results){
  if (min(fit_results$fit_metrics$aicc) < 0){
    # exclude erroneous model
    fit_results$fit_metrics = fit_results$fit_metrics[-which(fit_results$fit_metrics$aicc < 0),]
  }
  # calculate delta_aicc 
  delta_aicc = as.numeric(fit_results$fit_metrics$aicc) - as.numeric(min(fit_results$fit_metrics$aicc))
  fit_results$fit_metrics = data.frame(fit_results$fit_metrics, delta_aicc)
  # find first and second lowest delta aicc
  first_delta = fit_results$fit_metrics[fit_results$fit_metrics$delta_aicc == min(fit_results$fit_metrics$delta_aicc),]
  minus_first = fit_results$fit_metrics[-which(fit_results$fit_metrics$delta_aicc == min(fit_results$fit_metrics$delta_aicc)),]
  second_delta =  minus_first[minus_first$delta_aicc == min(minus_first$delta_aicc),]
  # compare delta aicc and pick best-fit model 
  if (second_delta$delta_aicc - first_delta$delta_aicc > 2){
    best_fit_model = first_delta
  } else {
    best_fit_model = second_delta
  }
  # finding best-model estimates
  best_fit_estimates = fit_results$model_estimates[names(fit_results$model_estimates) == best_fit_model$model]
  best_fit_estimates = best_fit_estimates[[1]]
  # pick best-fit estimates
  theta_se = c()
  for(j in 1:nrow(best_fit_estimates$theta) ){
    thse = c(best_fit_estimates$theta[j,])
    names(thse)[1] = paste("theta", as.character(j), sep="_")
    names(thse)[2] = paste("se", as.character(j), sep="_")
    theta_se = c(theta_se, thse)
  }
  alp_sig = c()
  for(j in 1:ncol(best_fit_estimates$solution) ){
    alsi=best_fit_estimates$solution[,j]
    names(alsi)[1] = paste("alpha", as.character(j), sep="_")
    names(alsi)[2] = paste("sigma_sq", as.character(j), sep="_")
    alp_sig = c(alp_sig, alsi)
  }
  best_estimate = c(theta_se,alp_sig)
  # packing things up
  best_fit_metric = best_fit_model[1,]
  best_list = list(best_fit_metric[-5], best_estimate)
  names(best_list) = c("best_fit", "best_estimates")
  return(best_list)
}



### fitting and choosing models across trees 

# setting regimes & models
regimes = data.frame(spp_hypervol$species, spp_hypervol$distribution, spp_hypervol$hypervolume)
colnames(regimes) = c("species","state", "trait")  
models=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA")


# looping across trees
all_best_models = data.frame(matrix(NA, nrow=length(trees), ncol=4))
colnames(all_best_models) = c(c("model","llik","aicc","delat_aicc"))
all_best_estimates = vector("list" , length(trees))

for (i in 1:length(trees)){
  one_tree = trees[[i]]
  all_fit = fit_evo_models(tree= one_tree, regimes, models_to_fit= models)
  best_choice = choose_best(all_fit)
  all_best_models[i,] = best_choice$best_fit
  all_best_estimates[[i]] = best_choice$best_estimates
}

write.table(all_best_models,"all_best_models.csv",sep=",",quote=F,row.names=F)

best_estimates_table=c()
for (i in 1:length(all_best_estimates)){
  best_estimates_table = rbind(best_estimates_table,all_best_estimates[[i]])
}

best_estimates_table = data.frame(best_estimates_table)

write.table(best_estimates_table,"best_estimates_table.csv",sep=",",quote=F,row.names=F)


########################### describing best-fit #############################

af_estimates = best_estimates_table[,c(1:2,5:6)]
colnames(af_estimates) = c("theta","se","alpha","sigma_q")
af_state = rep("AF-endemic", length(af_estimates[,1]))

ws_estimates = best_estimates_table[,c(3:4,7:8)]
colnames(ws_estimates) = c("theta","se","alpha","sigma_q")
ws_state = rep("widespread", length(ws_estimates[,1]))

distribution_state = c(af_state, ws_state)
estimates = rbind(af_estimates, ws_estimates)
estimates = data.frame(distribution_state, estimates)


########################### fitting and simulating BM evolution ##################

### fitting BM model
bm_fits=vector(mode = "list", length = length(trees))

for (i in 1:length(trees)){
	one_tree=trees[[i]]
	simmaps=make.simmap(one_tree, x=states, model="ER", nsim=100)
	mean_map=summary(simmaps)
	anc_states= rep(NA, length(mean_map$ace[,1]))
	for (j in 1:length(mean_map$ace[,1])){
	  bool = mean_map$ace[j,] == max(mean_map$ace[j,])
	  high_prob_state = names(bool)[bool]
	  anc_states[j] = high_prob_state
	}
	one_tree$node.label=anc_states
	fit=OUwie(one_tree,data=regimes,model="BM1")
	bm_fits[[i]]=fit
}


bm_estimates=as.data.frame(matrix(NA,nrow=length(bm_fits),ncol=5))
colnames(bm_estimates)=c("model","theta","se","alpha","sigma_q")

for (i in 1:length(bm_fits)){
  bm_estimates[i,1]=bm_fits[[i]]$model
	bm_estimates[i,2:3]=bm_fits[[i]]$theta[1,] # optima & se
	bm_estimates[i,4:5]=bm_fits[[i]]$solution[,1] # alpha & sigma_q
	}

write.table(bm_estimates,"bm_estimates.csv", sep=",", quote=F, row.names=F)


### simulating evolution under BM
bm_simulations=c("AF-enfdemic","widespread")

for(i in 1:length(trees)){
	one_tree=trees[[i]]
	for (j in 1:100){
  	sim_hypervolumes = fastBM(one_tree, a=sample(hypervolumes,1), sig2=bm_estimates$sigma_q[i], bounds=c(0.0001,Inf), internal=F)
  	af_mean= mean(sim_hypervolumes[which(names(sim_hypervolumes)%in% AF)])
  	ws_mean= mean(sim_hypervolumes[-which(names(sim_hypervolumes)%in% AF)])
  	bm_simulations = rbind(bm_simulations, c(af_mean, ws_mean))
  }
}

bm_simulations = data.frame(bm_simulations[-1,])
colnames(bm_simulations)<-c("AF-endemic","widespread")

write.table(bm_simulations,"bm_simulations.csv", sep=",", quote=F, row.names=F)


### summarizing hypervolume values
means=aggregate(spp_hypervol$hypervolume,by=list(spp_hypervol$distribution),mean)
sds=aggregate(spp_hypervol$hypervolume,by=list(spp_hypervol$distribution),sd)
ses=sds[,2]/c(sqrt(27),sqrt(39))

summary_hypervolume=data.frame(means,sds[,2],ses)
colnames(summary_hypervolume)=c("distribution","hypervolume","sd","se")

simulated_hypervolumes=data.frame(c(rep("AF-endemcic",10000),rep("widespread",10000)),c(bm_simulations[,1],bm_simulations[,2]))
colnames(simulated_hypervolumes)=c("distribution","hypervolume")
simulated_hypervolumes$hypervolume = as.numeric(simulated_hypervolumes$hypervolume)

sum(summary_hypervolume[1,2] > simulated_hypervolumes[,2])
sum(summary_hypervolume[2,2] > simulated_hypervolumes[,2])

sum(summary_hypervolume[2,2] > bm_simulations[,1])

601/10000
################################### ploting #################################

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

theta_plot=ggplot(data=estimates,aes(x=distribution_state, y=theta, fill=distribution_state)) +
  geom_point(aes(y=theta, color=distribution_state), position = position_jitter(width = 0.07), size = 1, alpha = 0.2) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.3)+
  scale_fill_manual(values=c("blue","red"))+
  scale_colour_manual(values=c("blue","red"))+
  labs(x="", y="theta")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10),legend.position = "none")


sigma_plot=ggplot(data=estimates, aes(x=distribution_state, y=sqrt(sigma_q), fill=distribution_state)) +
  geom_point(aes(y=sqrt(sigma_q), color=distribution_state), position = position_jitter(width = 0.07), size = 1, alpha = 0.2) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.3)+
  scale_fill_manual(values=c("blue","red"))+
  scale_colour_manual(values=c("blue","red"))+
  labs(x="", y="sigma")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10),legend.position = "none")


hypervol_plot=ggplot(data=bm_simulations,aes(x=AF.endemic)) +
	geom_density(alpha=0.5,color="lightgray", fill="lightgray")+
	scale_colour_manual(values=c("blue","red"))+
	scale_fill_manual(values=c("blue1","red1"))+
	geom_vline(data=summary_hypervolume,aes(xintercept = hypervolume),linetype="solid",colour=c("blue3","red3"),size=1)+
	geom_vline(data=summary_hypervolume,aes(xintercept = hypervolume-se),linetype="dotted",colour=c("blue3","red3"),size=0.50)+
	geom_vline(data=summary_hypervolume,aes(xintercept = hypervolume+se),linetype="dotted",colour=c("blue3","red3"),size=0.50)+
  xlim(-0.2, 5)+
	labs(x="mean \n hypervolume", y="density")+
	theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=12),legend.position = "none")


tiff("oumv_estimates_hypervol.tiff", units="in", width=2.5, height=6, res=600)
ggarrange(theta_plot,sigma_plot,hypervol_plot,nrow=3,ncol=1)
dev.off()
