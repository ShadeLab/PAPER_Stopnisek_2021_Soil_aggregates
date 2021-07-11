################################
#' Code for analysis of AOA and AOB amoA qPCR results
################################

#' Import data
AOAqpcr=read.csv('amoA_qPCR/qPCR_AOA.csv') 
AOBqpcr=read.csv("amoA_qPCR/qPCR_AOB.csv")
map_qpcr=read.csv('amoA_qPCR/qPCR_map.csv')
IndexDF=read.table('Surface_volume_weight_indeces.txt')

head(AOAqpcr)
head(map_qpcr)

qpcrAOA_fig=AOAqpcr %>%
  left_join(map_qpcr, by = c( 'Well' = 'well')) %>%
  filter(sample_type =='sample') %>%
  mutate(conc_per_gram=Starting.Quantity..SQ./.75) %>% #0.75g of dry soil was used to isolate DNA
  ggplot(aes(x=as.factor(Size),y=conc_per_gram, fill=Soil)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#009e73", "#0072b2")) + 
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  labs(x='Soil particle size (mm)', 
       y='AOA amoA gene abundance\n(gene copy/g dry soil)') +
  theme_classic() +
  theme(legend.position = c(.2,.8))

qpcrAOB_fig=AOBqpcr %>%
  left_join(map_qpcr, by = c( 'Well' = 'well')) %>%
  filter(sample_type =='sample') %>%
  mutate(conc_per_gram=Starting.Quantity..SQ./.75) %>% #0.75g of dry soil was used to isolate DNA
  ggplot(aes(x=as.factor(Size),y=Starting.Quantity..SQ., fill=Soil)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#009e73", "#0072b2")) + 
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  labs(x='Soil particle size (mm)', 
       y='AOB amoA gene abundance\n(gene copy/g dry soil)') +
  theme_classic() +
  theme(legend.position = 'none')

ggarrange(qpcrAOA_fig, qpcrAOB_fig,
          labels = c("A", "B"))

#' next we will normalize the gene abundance by particle surface and 
#' surface/volume. We need first to calculate mean of the values per soil and 
#' particle size.

qpcr_index=IndexDF %>%
  group_by(Soil, Size) %>%
  summarise(meanSW=mean(meanSW),
            meanSVW=mean(meanSVW)) %>%
  mutate(Size=if_else(Size == 0.18, 0.16, Size))

qpcrPlot_AOA=AOAqpcr %>%
  left_join(map_qpcr, by = c( 'Well' = 'well')) %>%
  filter(sample_type =='sample') %>%
  left_join(qpcr_index) %>%
  filter(!is.na(meanSVW)) %>%
  mutate(normAbundanceSW=Starting.Quantity..SQ./meanSW,
         normAbundanceSVW=Starting.Quantity..SQ./meanSVW,
         conc_per_gram=normAbundanceSW / 0.75) %>% #0.75g of dry soil was used to isolate DNA
  ggplot(aes(x=as.factor(Size),y=conc_per_gram, fill=Soil)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#009e73", "#0072b2")) + 
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  labs(x='Soil particle size (mm)', y='AOA amoA gene abundance per um2') +
  theme_classic() +
  theme(legend.position = c(.2,.8))

qpcrPlot_AOB=AOBqpcr %>%
  left_join(map_qpcr, by = c( 'Well' = 'well')) %>%
  filter(sample_type =='sample') %>%
  left_join(qpcr_index) %>%
  filter(!is.na(meanSVW)) %>%
  mutate(normAbundanceSW=Starting.Quantity..SQ./meanSW,
         normAbundanceSVW=Starting.Quantity..SQ./meanSVW,
         conc_per_gram=normAbundanceSW / 0.75) %>% #0.75g of dry soil was used to isolate DNA
  ggplot(aes(x=as.factor(Size),y=conc_per_gram, fill=Soil)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#009e73", "#0072b2")) + 
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  labs(x='Soil particle size (mm)', y='AOB amoA gene abundance per um2') +
  theme_classic() +
  theme(legend.position = c(.2,.8))

qcprDF= rbind(AOBqpcr, AOAqpcr) %>% 
  left_join(map_qpcr, by = c( 'Well' = 'well')) %>%
  filter(sample_type =='sample') %>%
  left_join(qpcr_index)

qpcrAbundanceFig=qcprDF %>% 
  filter(!is.na(meanSVW)) %>% mutate(normAbundanceSW=Starting.Quantity..SQ./meanSW,
                                     normAbundanceSVW=Starting.Quantity..SQ./meanSVW,
                                     conc_per_gram=normAbundanceSW / 0.75) %>% #0.75g of dry soil was used to isolate DNA
  ggplot(aes(x=Target,y=conc_per_gram, fill=Soil)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#009e73", "#0072b2")) + 
  labs(x=NULL, y=NULL) +
  facet_grid(~Size)+
  ylim(0, 5e05)+
  theme_classic() +
  theme(legend.position = c(.2,.7), 
        axis.text.y = element_blank())


qpcrRatioFig=qcprDF %>% 
  filter(!is.na(meanSVW)) %>% mutate(normAbundanceSW=Starting.Quantity..SQ./meanSW,
                                     normAbundanceSVW=Starting.Quantity..SQ./meanSVW,
                                     conc_per_gram=normAbundanceSW / 0.75) %>%
  select(Soil, Size, rep,Target, conc_per_gram) %>%
  pivot_wider(values_from=conc_per_gram, names_from=Target) %>%
  group_by(Soil, Size) %>%
  summarise(ratioAOA_AOB=mean(AOA)/mean(AOB)) %>%
  ggplot(aes(x=Soil, y=log2(ratioAOA_AOB), color=Soil)) +
  geom_point(size=3) +
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  theme_classic()+
  facet_wrap(~Size, nrow = 1)+
  theme(legend.position = 'none',
        strip.background = element_blank(), 
        strip.text.x = element_blank(),, 
        axis.text.y = element_blank()) +
  labs(x=NULL, y=NULL) 

ggarrange(qpcrAbundanceFig,qpcrRatioFig, nrow = 2,
          heights = c(2,1))


