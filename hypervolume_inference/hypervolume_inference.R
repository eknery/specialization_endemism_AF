library(hypervolume)


# loading environmental background. 
bg_env= read.table("environmental_background.csv", header =T, sep=",")

# loading species environment
spp_env= read.table("species_environment.csv", header =T, sep=",")


# calculating descriptive stats and scaling to z-scores
env_means = c()
env_sds = c()
for(i in 1:length(spp_env[1,])){
    env_means = c(env_means, mean(spp_env[,i]))
    env_sds = c(env_sds,sd(spp_env[,i]))
}

spp_scale=data.frame(spp_env[,1],scale(spp_env[,-1]))

str(spp_scale)

######################## setting TSS function for hypervolumes #################

hypervolume_tss= function(data, bg, leave_one=FALSE, k, iterations, threshold){
  if (leave_one){
    k=nrow(data)
    iterations=nrow(data)
    folds=kfold(data, k=k)
    stats=data.frame(matrix(NA, nrow=iterations, ncol=4))
    colnames(stats) = c('accuracy','sensitivity', 'specificity', 'tss')
    for (i in 1:iterations){
      test = data[folds == i,]
      train = data[folds != i,]
      ### estimating species bandwidth and training hypervolume
      data_band=estimate_bandwidth(train,method="silverman")
      data_hv=hypervolume_box(train, samples.per.point = 1000, kde.bandwidth = data_band, tree.chunksize = 10000)
      ### predicting true presences
      data_test=hypervolume_inclusion_test(data_hv, test, reduction.factor = 1, fast.or.accurate ="accurate", accurate.method.threshold=threshold)
      a=sum(data_test, na.rm = TRUE)
      c=length(data_test)-sum(data_test, na.rm = TRUE)
      ### sampling background environment with same effort as test set
      random =round(runif(n=nrow(test), min = 1, max = nrow(bg)))
      bg_rand=bg[random,]
      ### predicting virtual absences
      bg_test=hypervolume_inclusion_test(data_hv, bg_rand, reduction.factor = 1, fast.or.accurate ="accurate", accurate.method.threshold=threshold)
      b=sum(bg_test, na.rm = TRUE)
      d=length(bg_test)-sum(bg_test, na.rm = TRUE)
      ### performance statistics
      n=2*length(test[,1])
      accuracy= (a + d)/n
      sensitivity= a/(a+c)
      specificity= d/(b+d)
      tss= sensitivity + (specificity -1)
      stats[i,]=c(accuracy,sensitivity,specificity,tss)
      }
  }else{
    stats=data.frame(matrix(NA, nrow=iterations, ncol=4))
    colnames(stats) = c('accuracy','sensitivity', 'specificity', 'tss')
    for (i in 1:iterations){
      folds=kfold(data, k=k)
      test = data[folds == 1,]
      train = data[folds != 1,]
      ### estimating species bandwidth and training hypervolume
      data_band=estimate_bandwidth(train,method="silverman")
      data_hv=hypervolume_box(train, samples.per.point = 1000, kde.bandwidth = data_band, tree.chunksize = 10000)
      ### predicting true presences
      data_test=hypervolume_inclusion_test(data_hv, test, reduction.factor = 1, fast.or.accurate ="accurate", accurate.method.threshold=threshold)
      a=sum(data_test, na.rm = TRUE)
      c=length(data_test)-sum(data_test, na.rm = TRUE)
      ### sampling background environment with same effort as test set
      random =round(runif(n=nrow(test), min = 1, max = nrow(bg)))
      bg_rand=bg[random,]
      ### predicting virtual absences
      bg_test=hypervolume_inclusion_test(data_hv, bg_rand, reduction.factor = 1, fast.or.accurate ="accurate", accurate.method.threshold=threshold)
      b=sum(bg_test, na.rm = TRUE)
      d=length(bg_test)-sum(bg_test, na.rm = TRUE)
      ### performance statistics
      n=2*length(test[,1])
      accuracy= (a + d)/n
      sensitivity= a/(a+c)
      specificity= d/(b+d)
      tss= sensitivity + (specificity -1)
      stats[i,]=c(accuracy,sensitivity,specificity,tss)
      }
    }
  return(stats)
}

###################### hypervolume model validation ##############################

### listing species for leave-one-out validation
narrow_spp = c("angelana","capixaba","dura","kollmannii","kriegeriana", "mellina","penduliflora","suberosa")

### validating hypervolume models
spp_list= unique(spp_scale[,1])
validation_tables = vector(mode = "list", length = length(spp_list))

for (i in 1:length(spp_list)){
  sp_data = spp_scale[spp_scale$species== spp_list[i],-1]
  if (spp_list[i] %in% narrow_spp){
    validation_tables[[i]]  = hypervolume_tss(data=sp_data, bg=bg_scale[,-1], leave_one=TRUE, threshold = 0.5)
    names(validation_tables)[i] = spp_list[i]
  } else {
    validation_tables[[i]]  = hypervolume_tss(data=sp_data, bg=bg_scale[,-1], k=3, iterations = 10, threshold = 0.5)
    names(validation_tables)[i] = spp_list[i]
  }
} 


# organizing validation tables into dataframe
hypervolume_validation = data.frame()

for(i in 1:length(spp_list)){
     sp_name = names(validation_tables)[i]
     n_rows = nrow(validation_tables[[i]])
     sp_validation = data.frame(rep(sp_name, n_rows) , validation_tables[[i]])
     hypervolume_validation = rbind(hypervolume_validation, sp_validation)
}

colnames(hypervolume_validation)[1] = "species"

write.table(hypervolume_validation,"hv_validation.csv", sep=",", quote=F, row.names=F)

mean_hv_validation = aggregate(hypervolume_validation[,-1], by = list(hypervolume_validation[,1]), mean)
colnames(mean_hv_validation)[1] = "species"

write.table(mean_hv_validation,"mean_hv_validation.csv", sep=",", quote=F, row.names=F)


####################### inferring hypervolumes and niches ############################

spp_hypervol = data.frame(spp_list,rep(NA, length(spp_list)))
colnames(spp_hypervol) = c("species","hypervolume")

envniche_tables = vector(mode = "list", length = length(spp_list))
names(envniche_tables) = spp_list

for (i in 1:length(spp_list)){
  sp_data = spp_scale[spp_scale[,1] == spp_list[i],-1]
  if (spp_list[i] %in% narrow_spp){
    sp_band = estimate_bandwidth(sp_data,method="silverman")
    sp_hv = hypervolume_box(sp_data, samples.per.point = 1000, kde.bandwidth = sp_band, tree.chunksize = 10000)
    spp_hypervol[i,2] = sp_hv@Volume
    rand_samp=round(runif(n=100, min = 1, max = nrow(sp_hv@RandomPoints)))
    envniche_tables[[i]] = sp_hv@RandomPoints[rand_samp,]
  }else{
    sp_band = estimate_bandwidth(sp_data,method="silverman")
    sp_hv = hypervolume_box(sp_data, samples.per.point = 1000, kde.bandwidth = sp_band, tree.chunksize = 10000)
    ### reducing to 0.5 probability quantile
    hv_th = hypervolume_threshold(sp_hv,quantile.requested = 0.5, quantile.requested.type = "probability", uniform.density = TRUE)
    sp_hv_th = hv_th[[1]]
    # getting hypervolume size
    spp_hypervol[i,2] = sp_hv_th@Volume
    # randomly sampling distributions over environmental variables
    rand_samp=round(runif(n=100, min = 1, max = nrow(sp_hv_th@RandomPoints)))
    envniche_tables[[i]] = sp_hv_th@RandomPoints[rand_samp,]
  }
}


write.table(spp_hypervol, "spp_hypervol.csv", sep=",", quote=F, row.names = F)

# organizing environmental variables into dataframe
spp_scale_niche = data.frame()

for(i in 1:length(spp_list)){
  sp_name = names(envniche_tables)[i]
  n_rows = nrow(envniche_tables[[i]])
  one_scale_niche = data.frame(rep(sp_name, n_rows) , envniche_tables[[i]])
  spp_scale_niche = rbind(spp_scale_niche, one_scale_niche)
}

colnames(spp_scale_niche)[1] = "species"

str(spp_scale_niche)

# back transforming environmental variables
scale_niche = spp_scale_niche[,-1]
env_niche = data.frame(matrix(NA,nrow(scale_niche), ncol(scale_niche)))

for (i in 1:length(scale_niche[1,])){
  env_niche[,i] = ((scale_niche[,i] * env_sds[i]) + env_means[i])
}

spp_env_niche = data.frame(spp_scale_niche[,1], env_niche)
colnames(spp_env_niche) = colnames(spp_scale_niche)

write.table(spp_env_niche,"spp_env_niche.csv", sep=",", quote=F, row.names=F)


######################## graphic analyses ################################
library(tidyverse)
library(PupillometryR)
library(ggpubr)
library(readr)
library(tidyr)

# defining endemics
AF=c("atlantica","baumgratziana","brunnea","budlejoides","capixaba", 
     "castaneiflora","cinerascens","discolor","dura","fasciculata",
     "formosa","hyemalis","kollmannii","kriegeriana","lymanii","mellina",
     "octopetala","penduliflora","petroniana","polyandra","racemifera",
     "robusta","ruschiana","setosociliata","shepherdii","valtheri","willdenowii")

# creating vector of distributions
distribution=spp_hypervol$species
distribution[which(spp_hypervol$species %in% AF)]="AF-endemic"
distribution[-which(spp_hypervol$species %in% AF)]="widespread"

spp_hypervol=data.frame(distribution,spp_hypervol)

aggregate(spp_hypervol$hypervolume, by = list(spp_hypervol$distribution), mean)
aggregate(spp_hypervol$hypervolume, by = list(spp_hypervol$distribution), sd)

### hypervolume x species distribution
tiff("hypervol x spp_distribution.tiff", units="in", width=3.5, height=3.5, res=600)
ggplot(data = spp_hypervol,aes(x=distribution,y=hypervolume,fill=distribution)) +
  geom_point(aes(y = hypervolume, color =distribution), position = position_jitter(width = 0.07), size = 2, alpha = 0.5) +
  geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.3)+
  scale_fill_manual(values=c("blue","red"))+
  scale_colour_manual(values=c("blue","red"))+
  ylim(0,5)+
  labs(x="species' distribution", y="hypervolume size")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=12),legend.position = "none")
dev.off()

# creating vector of distributions
spp_dist=spp_env_niche$species 
spp_dist[which(spp_dist %in% AF)]="AF-endemic"
spp_dist[which(spp_dist != "AF-endemic")]="widespread"

spp_env_niche=data.frame(spp_dist,spp_env_niche)
colnames(spp_env_niche)[1]="distribution"

afs=spp_env_niche[which(spp_env_niche$distribution=="AF-endemic"),]
wds=spp_env_niche[which(spp_env_niche$distribution=="widespread"),]


### niche density distribution x species distribution

# temperature
range(spp_env_niche$temp_diurnal_range)

afs_temp=ggplot(data = afs,aes(x=temp_diurnal_range*10^-1,fill=distribution,col=distribution)) +
  geom_density(aes(group=species), alpha=0.1)+
  scale_color_manual(values=c("blue"))+
  scale_fill_manual(values=c("blue"))+
  labs(x="", y="")+
  xlim(5,15)+ ylim(0,4)+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=11,face="bold"),legend.position = "none")

wds_temp=ggplot(data = wds,aes(x=temp_diurnal_range*10^-1,fill=distribution,col=distribution)) +
  geom_density(aes(group=species), alpha=0.1)+
  scale_color_manual(values=c("red"))+
  scale_fill_manual(values=c("red"))+
  labs(x="temperature range", y="kernel density")+
  xlim(5,15)+ ylim(0,4)+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=11,face="bold"),legend.position = "none")

tiff("temp_niche.tiff", units="in", width=6, height=4.5, res=600)
ggarrange(afs_temp,wds_temp,nrow=2,ncol=1,common.legend = F)
dev.off()


# precipitation
range(spp_env_niche$precip_seasonality)

afs_precip=ggplot(data = afs,aes(x=precip_seasonality,fill=distribution,col=distribution)) +
  geom_density(aes(group=species), alpha=0.1)+
  scale_color_manual(values=c("blue"))+
  scale_fill_manual(values=c("blue"))+
  labs(x="", y="")+
  xlim(0,95)+ylim(0,1)+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=11,face="bold"),legend.position = "none")

wds_precip=ggplot(data = wds,aes(x=precip_seasonality,fill=distribution,col=distribution)) +
  geom_density(aes(group=species), alpha=0.1)+
  scale_color_manual(values=c("red"))+
  scale_fill_manual(values=c("red"))+
  labs(x="precipitation seasonality", y="kernel density")+
  xlim(0,95)+ylim(0,1)+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=11,face="bold"),legend.position = "none")

tiff("precip_niche.tiff", units="in", width=6, height=4.5, res=600)
ggarrange(afs_precip,wds_precip,nrow=2,ncol=1,common.legend = F)
dev.off()


#soil pH
range(spp_env_niche$soil_pH)

afs_soil=ggplot(data = afs,aes(x=soil_pH*10^-1,fill=distribution,col=distribution)) +
  geom_density(aes(group=species), alpha=0.1)+
  scale_color_manual(values=c("blue"))+
  scale_fill_manual(values=c("blue"))+
  labs(x="", y="")+
  xlim(3,7)+ylim(0,10)+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=11,face="bold"),legend.position = "none")


wds_soil=ggplot(data = wds,aes(x=soil_pH*10^-1,fill=distribution,col=distribution)) +
  geom_density(aes(group=species), alpha=0.1)+
  scale_color_manual(values=c("red"))+
  scale_fill_manual(values=c("red"))+
  labs(x="soil pH", y="kernel density")+
  xlim(3,7)+ylim(0,10)+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=11,face="bold"),legend.position = "none")


tiff("soil_niche.tiff", units="in", width=6, height=4.5, res=600)
ggarrange(afs_soil,wds_soil,nrow=2,ncol=1,common.legend = F)
dev.off()


# solar radiation
range(spp_env_niche$solar_radiation)*10^-3

afs_solar=ggplot(data = afs,aes(x=solar_radiation*10^-3,fill=distribution,col=distribution)) +
  geom_density(aes(group=species), alpha=0.1)+
  scale_color_manual(values=c("blue"))+
  scale_fill_manual(values=c("blue"))+
  labs(x="", y="")+
  xlim(12,20)+ylim(0,6)+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=11,face="bold"),legend.position = "none")


wds_solar=ggplot(data = wds,aes(x=solar_radiation*10^-3,fill=distribution,col=distribution)) +
  geom_density(aes(group=species), alpha=0.1)+
  scale_color_manual(values=c("red"))+
  scale_fill_manual(values=c("red"))+
  labs(x="solar radiation", y="kernel density")+
  xlim(12,20)+ylim(0,6)+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=11,face="bold"),legend.position = "none")


tiff("solar_niche.tiff", units="in", width=6, height=4.5, res=600)
ggarrange(afs_solar,wds_solar,nrow=2,ncol=1,common.legend = F)
dev.off()

### niche sds

sds_niche=aggregate(spp_env_niche[,-c(1:2)], by=list(spp_env_niche[,2]), sd)

sds_dist=sds_niche$Group.1
sds_dist[which(sds_dist %in% AF)] ="AF-endemic"
sds_dist[which(sds_dist!="AF-endemic")]="widespread"

sds_niche=data.frame(sds_dist,sds_niche)
colnames(sds_niche)[c(1,2)]=c("distribution","species")

#plots 

sd_temp=ggplot(data = sds_niche,aes(x=distribution,y=temp_diurnal_range,fill=distribution)) +
  geom_point(aes(y = temp_diurnal_range, color =distribution), position = position_jitter(width = 0.07), size = 2, alpha = 0.15) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.3)+
  scale_fill_manual(values=c("blue","red"))+
  scale_colour_manual(values=c("blue","red"))+
  labs(x="", y="sd temperature range")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=13,face="bold"),axis.text.x=element_text(size=11),legend.position = "none")

sd_precip=ggplot(data = sds_niche,aes(x=distribution,y=precip_seasonality,fill=distribution)) +
  geom_point(aes(y = precip_seasonality, color =distribution), position = position_jitter(width = 0.07), size = 2, alpha = 0.15) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.3)+
  scale_fill_manual(values=c("blue","red"))+
  scale_colour_manual(values=c("blue","red"))+
  labs(x="", y="sd precipitation")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=13,face="bold"),axis.text.x=element_text(size=11),legend.position = "none")

sd_soil=ggplot(data = sds_niche,aes(x=distribution,y=soil_pH,fill=distribution)) +
  geom_point(aes(y = soil_pH, color =distribution), position = position_jitter(width = 0.07), size = 2, alpha = 0.15) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.3)+
  scale_fill_manual(values=c("blue","red"))+
  scale_colour_manual(values=c("blue","red"))+
  labs(x="", y="sd soil pH")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=13,face="bold"),axis.text.x=element_text(size=11),legend.position = "none")

sd_solar=ggplot(data = sds_niche,aes(x=distribution,y=solar_radiation,fill=distribution)) +
  geom_point(aes(y = solar_radiation, color =distribution), position = position_jitter(width = 0.07), size = 2, alpha = 0.15) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.3)+
  scale_fill_manual(values=c("blue","red"))+
  scale_colour_manual(values=c("blue","red"))+
  labs(x="species' distribution", y="sd solar radiation")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=13,face="bold"),axis.text.x=element_text(size=11),legend.position = "none")

tiff("sds_niche.tiff", units="in", width=6, height=6.5, res=600)
ggarrange(sd_temp, sd_precip, sd_soil, sd_solar,nrow=2,ncol=2,common.legend = F)
dev.off()
