###############################
#' Part of the code for the analysis of the soil particle 
#' surface area, volume and weight
###############################
library(tidyverse)

setwd('~/Documents/git/SoilAggregates/')

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

#' Using surface area and volume ratio
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

#' Same as before but controlling for the soil weight used in scanning
surfaceAreaPlot2=surface_volume_weight %>%
  mutate(id.v2=paste(Site, Size, sep='_'),
         AreaVolumeWeight=SA.Vol_um.1/weight) %>%
  ggplot(aes(x=Size, y=AreaVolumeWeight, color= Site)) +
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  geom_smooth(method="glm", method.args=list(family=gaussian(link="log")), color='grey', se=T)+
  geom_point(size=3, alpha=.7) +
  stat_summary(fun.data=median_hilow, fun.args=list(conf.int=.5),
               geom='pointrange')+
  theme_classic() +
  facet_wrap(~Soil) +
  theme(legend.position = 'none') +
  labs(y='area/volume/weight ((um x mg)-1)')

# ANOVA
area.vol.DF=surface_volume_weight %>%
  mutate(AreaVolumeWeight=SA.Vol_um.1/weight)

area.vol.aov <- aov(AreaVolumeWeight ~ Site, data = area.vol.DF)
# Summary of the analysis
summary(area.vol.aov)

# t-test for size comparison
stat.test.area.vol <- area.vol.DF %>%
  group_by(Size) %>%
  t_test(AreaVolumeWeight ~ Site) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test.area.vol

#' Area Vs weight of used soil
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

surface_volume_weight %>%
  mutate(id.v2=paste(Site, Size, sep='_'),
         AreaWeight=SA_um/(weight*1000)) %>%
  filter(Site == 'M') %>%
  ggplot(aes(x=Size, y=AreaWeight)) +
  geom_point(size=3, alpha=.7, col="#009e73") +
  stat_summary(fun.data=median_hilow, fun.args=list(conf.int=.5),
               geom='pointrange', color='grey')+
  geom_smooth(method="lm", formula= (y ~ exp(-x/0.33)), se=T, color='grey')+
  theme_classic() +
  theme(legend.position = 'none') +
  labs(y='area/weight (um2 per g soil)', x='Soil particle size (mm)')


df.test.S=surface_volume_weight %>%
  mutate(id.v2=paste(Site, Size, sep='_'),
         AreaWeight=SA_um/(weight*1000)) %>%
  filter(Site == 'S') %>%
  arrange(Size)

df.test.M=surface_volume_weight %>%
  mutate(id.v2=paste(Site, Size, sep='_'),
         AreaWeight=SA_um/(weight*1000)) %>%
  filter(Site == 'M') %>%
  arrange(Size)

#' Below code was developed by Susan Johnston 
#' (https://sejohnston.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/)
ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

ggplotRegression(lm(log(df.test.M$AreaWeight) ~ df.test.M$Size, data = df.test.M))
ggplotRegression(lm(log(df.test.S$AreaWeight) ~ df.test.S$Size, data = df.test.S))


# ANOVA
area.DF=surface_volume_weight %>%
  mutate(AreaWeight=SA_um/(weight*1000))

area.aov <- aov(AreaWeight ~ Site, data = area.DF)
# Summary of the analysis
summary(area.aov)

# t-test for size comparison
stat.test.area <- area.DF %>%
  group_by(Size) %>%
  t_test(AreaWeight ~ Site) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test.area


#' Volume Vs weight of used soil
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
 
#' What would be the total surface are if considering the measured weight 
#' fractions?

#' Mean surface area values
surface_volume_weight %>%
  mutate(area.weight=SA_um/(weight*1000)) %>%
  group_by(Soil, Size) %>%
  summarise(mean.area.weight=mean(area.weight)) %>%
  mutate(Fractions = Size) %>%
  left_join(weightFraction) %>%
  select(-X, -X.1) %>%
  mutate(SizeByArea= mean.area.weight*relative_portion_size) %>%
  group_by(Soil) %>%
  summarise(sumArea=sum(SizeByArea))

surface_volume_weight %>%
  mutate(area.weight=SA_um/(weight*1000)) %>%
  group_by(Soil, Size) %>%
  summarise(mean.area.weight=mean(area.weight)) %>%
  mutate(Fractions = Size) %>%
  left_join(weightFraction) %>%
  select(-X, -X.1) %>%
  mutate(SizeByArea= mean.area.weight*relative_portion_size) %>%
  filter(Soil == 'SVERC',
         Fractions>.4) %>%
  summarise(sumArea=sum(SizeByArea))
