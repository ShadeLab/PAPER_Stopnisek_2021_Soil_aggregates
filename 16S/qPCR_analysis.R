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

AOAqpcr %>%
  left_join(map_qpcr, by = c( 'Well' = 'well')) %>%
  filter(sample_type =='sample') %>%
  mutate(conc_per_gram=Starting.Quantity..SQ./.75) %>% #0.75g of dry soil was used to isolate DNA
  ggplot(aes(x=as.factor(Size),y=conc_per_gram, fill=Soil)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#009e73", "#0072b2")) + 
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  labs(x='Soil particle size (mm)', y='AOA amoA gene abundance') +
  theme_classic() +
  theme(legend.position = c(.2,.8))

AOBqpcr %>%
  left_join(map_qpcr, by = c( 'Well' = 'well')) %>%
  filter(sample_type =='sample') %>%
  mutate(conc_per_gram=Starting.Quantity..SQ./.75) %>% #0.75g of dry soil was used to isolate DNA
  ggplot(aes(x=as.factor(Size),y=Starting.Quantity..SQ., fill=Soil)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#009e73", "#0072b2")) + 
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  labs(x='Soil particle size (mm)', y='AOB amoA gene abundance') +
  theme_classic() +
  theme(legend.position = 'none')


#' next we will normalize the gene abundance by particle surface and 
#' surface/volume. We need first to calculate mean of the values per soil and 
#' particle size.

qpcr_index=IndexDF %>%
  group_by(Soil, Size) %>%
  summarise(meanSW=mean(meanSW),
            meanSVW=mean(meanSVW))

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
