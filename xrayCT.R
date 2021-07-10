###############################
#' Part of the code for the analysis of the soil particle 
#' surface area, volume and weight
###############################

setwd('Documents/git/SoilAggregates/')

weight=read.delim('soilWeightXRAYCT.txt')

map=read.delim('~/Documents/git/SoilAggregates/mapFile.txt', header=T, sep=' ')
map$ID=rownames(map)
rownames(map)=NULL

surfaceArea = read.delim('surfaceArea_full.txt')

weight_short=weight %>%
  mutate(id_new=paste(Site,size, Rep, sep='_')) %>%
  select(id_new, weight)

surface_volume_weight=surfaceArea %>%
  left_join(weight_short, by='id_new')

surfaceAreaPlot=surface_volume_weight %>%
  mutate(id.v2=paste(Site, Size, sep='_')) %>%
  ggplot(aes(x=Size, y=SA.Vol_um.1, color= Site)) +
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  geom_point(size=2) +
  theme_classic() +
  facet_wrap(~Soil) +
  theme(legend.position = 'none') +
  labs(y='area/volume (um-1)')
head(surface_volume_weight)

#Controlling for the soil weight
surfaceAreaPlot2=surface_volume_weight %>%
  mutate(id.v2=paste(Site, Size, sep='_'),
         AreaVolumeWeight=SA.Vol_um.1/weight) %>%
  ggplot(aes(x=Size, y=AreaVolumeWeight, color= Site)) +
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  geom_smooth(method="glm", method.args=list(family=gaussian(link="log")), color='grey', se=F)+
  geom_point(size=3, alpha=.7) +
  stat_summary(fun.data=median_hilow, fun.args=list(conf.int=.5),
               geom='pointrange')+
  theme_classic() +
  facet_wrap(~Soil) +
  theme(legend.position = 'none') +
  labs(y='area/volume/weight ((um x mg)-1)')

#' Area Vs weight
surfaceAreaPlot3=surface_volume_weight %>%
  mutate(id.v2=paste(Site, Size, sep='_'),
         AreaWeight=SA_um/(weight*1000)) %>%
  ggplot(aes(x=Size, y=AreaWeight, color= Site)) +
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  geom_point(size=3, alpha=.7) +
  stat_summary(fun.data=median_hilow, fun.args=list(conf.int=.5),
               geom='pointrange', color='grey')+
  theme_classic() +
  facet_wrap(~Soil) +
  theme(legend.position = 'none') +
  labs(y='area/weight (um2 per g soil)', x='Soil particle size (mm)')

#' Volume Vs weight
surfaceAreaPlot4=surface_volume_weight %>%
  mutate(id.v2=paste(Site, Size, sep='_'),
         VolWeight=Vol_um/weight) %>%
  ggplot(aes(x=Size, y=VolWeight, color= Site)) +
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  geom_point(size=3, alpha=.7) +
  stat_summary(fun.data=median_hilow, fun.args=list(conf.int=.5),
               geom='pointrange', color='grey')+
  theme_classic() +
  facet_wrap(~Soil) +
  theme(legend.position = 'none') +
  labs(y='V/weight (um3mg-1)')


#' Calculating the difference index for each sample as proportion to the smallest
#' number (value/min) 
names(surface_volume_weight)

meanClaculations=surface_volume_weight %>%
  group_by(Soil, ID, Size, Rep.y, id_new) %>%
  summarise(meanS=mean(SA_um), # in um^2
            meanSV=mean(SA.Vol_um.1), # in um-1
            meanW =mean(weight), # in mg
            meanSW=meanS/(meanW*1000),
            meanSVW=meanSV/meanW)

M_indeces=meanClaculations %>%
  filter(Soil == 'MRC') %>%
  mutate(minSW=min(meanClaculations$meanSW[meanClaculations$Soil == 'MRC']),
         minSVW=min(meanClaculations$meanSVW[meanClaculations$Soil == 'MRC']),
         indexA=meanSW/minSW,
         indexB=meanSVW/minSVW)

S_indeces=meanClaculations %>%
  filter(Soil == 'SVERC') %>%
  mutate(minSW=min(meanClaculations$meanSW[meanClaculations$Soil == 'SVERC']),
         minSVW=min(meanClaculations$meanSVW[meanClaculations$Soil == 'SVERC']),
         indexA=meanSW/minSW,
         indexB=meanSVW/minSVW)

IndexDF=rbind(M_indeces, S_indeces)
write.table(IndexDF, 'Surface_volume_weight_indeces.txt', quote = F)

weightFraction=read.delim(file = 'soil_sieving_weightFractions.txt')

mrc_weight=weightFraction[weightFraction$Soil == 'MRC',]
sverc_weight=weightFraction[weightFraction$Soil == 'SVERC',]


m_weight_plot=
  ggplot() +
  geom_bar(data=mrc_weight, 
           aes(x=Soil, y=relative_portion_size, fill=Fractions/2),
           stat = 'identity', color='black') +
  scale_fill_gradient(low = "white", high = "#009e73", 
                      breaks=c(0.05/2,0.5/2,1/2,2/2),
                      labels=c("<0.056um","0.5mm","1mm","2mm"),
                      limits = c(0,1))+
  theme_classic()+
  labs(y='Relative weight contribution of\nsoil particle size fractions',
       fill=NULL) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

s_weight_plot= 
  ggplot() + 
  geom_bar(data=sverc_weight,
           aes(x=Soil, y=relative_portion_size, fill=Fractions/2), 
           stat = "identity", color='black') +
  scale_fill_gradient(low = "white", high = "#0072b2", 
                      breaks=c(0.05/2,0.5/2,1/2,2/2),
                      labels=c("<0.056um","0.5mm","1mm" ,"2mm"),
                      limits = c(0,1)) +
  theme_classic() +
  labs(y='Relative weight contribution of\nsoil particle size fractions') +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) 
 
