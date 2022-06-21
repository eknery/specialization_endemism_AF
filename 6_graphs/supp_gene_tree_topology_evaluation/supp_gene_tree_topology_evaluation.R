############################# gree tree topology evaluation plot ##############################

tiff("0_data/topology_dissimilarity.tiff", units="in", width=4, height=3, res=600)
ggplot(data = df_rf,aes(x=Axis.1,y=Axis.2,color=locus)) +
  geom_point(size = 1, alpha = 0.5) +
  scale_colour_manual(values=c("red4","orange2","blue4","cyan4","red2","hotpink1"))+
  labs(x="PCo1", y="PCo2")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=13,face="bold"))
dev.off()