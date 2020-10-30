library(RColorBrewer)
library(viridis)
library(microbiome)
library(reshape2)
library(ggpubr)
library(broom)
#library(ggfortify)
library(dplyr)
library(phyloseq)
library(decontam)
library(vegan)
library(tidyverse)
library(indicspecies)
library(lemon)
library(ggalluvial)
library(pheatmap)
library(metagMisc)
BiocManager::install("metagenomeSeq")
library(metagenomeSeq)
install.packages("glue")
library(glue)
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
devtools::install_github("jbisanz/qiime2R")

setwd('~/Documents/git/SoilAggregates/')

otu <- read.table('R1_OTU_table.txt', sep='\t', header=T, row.names=1)
map <- read.csv('map.csv', row.names=1)

OTU=otu_table(as.matrix(otu), taxa_are_rows = T)
MAP=sample_data(map)
otuPhyloseq=phyloseq(MAP, OTU)

sample_data(otuPhyloseq)$is.neg <- sample_data(otuPhyloseq)$Type == "control"

contamdf.prev <- isContaminant(otuPhyloseq, method="prevalence", neg="is.neg", threshold=0.5)
keepOTU <- rownames(contamdf.prev[contamdf.prev$contaminant=='FALSE',])
omit.OTUs <- rownames(contamdf.prev[contamdf.prev$contaminant=='TRUE',])

otuPhyloseq_filt_decont <- prune_taxa(keepOTU,otuPhyloseq)
otuPhyloseq_filt_decont_samp <- subset_samples(otuPhyloseq_filt_decont, Type != 'control')

write_phyloseq(otuPhyloseq_filt_decont_samp, type='METADATA', path=".")
map_final <- read.csv('metadata_table.csv')

otu_final <- phyloseq_to_df(otuPhyloseq_filt_decont_samp, addtax=F, addtot=T)
otu_final <- otu_final[complete.cases(otu_final),]
rownames(otu_final) <- otu_final$OTU
otu_final$Total <- NULL
otu_final$OTU <- NULL

otu_final <- otu_final[,order(colnames(otu_final))]
map_final=map_final[order(map_final$id),]
colnames(otu_final)==map_final$id

set.seed(012)

OTU.rare <- t(rrarefy(t(otu_final), min(colSums(otu_final)))) 

rel.abun.all <- decostand(OTU.rare, method = 'total', MARGIN = 2)

set.seed(014)
map <- map_final
otu.BC <- vegdist(t(OTU.rare), method="bray")
otu.J <- vegdist(t(OTU.rare), method="jaccard")
otu.bc.pcoa <- cmdscale(otu.BC, eig=T)
otu.j.pcoa <- cmdscale(otu.J, eig=T)

map$Axis1.BC <- otu.bc.pcoa$points[,1]
map$Axis2.BC <- otu.bc.pcoa$points[,2]
map$Axis1.J <- otu.j.pcoa$points[,1]
map$Axis2.J <- otu.j.pcoa$points[,2]
ax1.bc.otu <- otu.bc.pcoa$eig[1]/sum(otu.bc.pcoa$eig)
ax2.bc.otu <- otu.bc.pcoa$eig[2]/sum(otu.bc.pcoa$eig)
ax1.j.otu <- otu.j.pcoa$eig[1]/sum(otu.j.pcoa$eig)
ax2.j.otu <- otu.j.pcoa$eig[2]/sum(otu.j.pcoa$eig)


pcoa.BC <- ggplot(map, aes(x=Axis1.BC, y=Axis2.BC)) +
  theme_classic() +
  geom_point(aes(col=Site, alpha=Size_fraction), size=4)+
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(legend.position = c(.5,.5),
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))+
  labs(x=paste('PCoA1 (',100*round(ax1.bc.otu,3),'%)',sep=''),y=paste('PCoA2 (',100*round(ax2.bc.otu,3),'%)', sep=''), 
       alpha='Size (mm)', title='Bray-Curtis')

pcoa.J <- ggplot(map, aes(x=Axis1.J, y=Axis2.J)) +
  theme_classic() +
  geom_point(aes(col=Site, alpha=Size_fraction), size=4)+
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(legend.position = 'none',
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))+
  labs(x=paste('PCoA1 (',100*round(ax1.j.otu,3),'%)',sep=''),y=paste('PCoA2 (',100*round(ax2.j.otu,3),'%)', sep=''), 
       alpha='Size (mm)', title="Jaccard")

ggarrange(ggarrange(pcoa.BC, pcoa.J, 
                    labels = c("A", "B"), ncol = 2)) %>%
  ggexport(filename = "PCoA.pdf", width = 5, height = 3)


s <- specnumber(OTU.rare,MARGIN=2)
h <- vegan::diversity(t(OTU.rare), "shannon")
pielou=h/log(s)

map$Richness <- s
map$Shannon <- h
map$Pielou <- pielou 

ggplot(map, aes(x=as.factor(Size_fraction),y=Richness, fill=Site)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = c("#009e73", "#0072b2")) + 
  labs(x='Size (mm)', y='Richness') +
  theme_classic()

library(gplots)
#For MRF
head(MRC1)
otu_venn <- OTU.rare
MRC_otu <- rel.abun.all[,map$Site=='MRC']
MRC_otu <- MRC_otu[rowSums(MRC_otu)>0,]
SVERC_otu <- otu_venn[,map$Site=='SVERC']
SVERC_otu <- SVERC_otu[rowSums(SVERC_otu)>0,]

MRC1 <- MRC_otu[,"MRC1"]
MRC1_otu <- data.frame(otu=names(MRC1), MRC1) %>%
  filter(MRC1>0)
MRC7 <- MRC_otu[,"MRC7"]
MRC7_otu <- data.frame(otu = names(MRC7), MRC7) %>%
  filter(MRC7>0)
MRC4 <- MRC_otu[,"MRC4"]
MRC4_otu <- data.frame(otu = names(MRC4), MRC4) %>%
  filter(MRC4>0)

MRC1_all <- MRC_otu[as.character(MRC1_otu$otu),]
MRC7_all <- MRC_otu[as.character(MRC7_otu$otu),]
MRC4_all <- MRC_otu[as.character(MRC4_otu$otu),]

data.frame(otu=rownames(MRC_otu), MRC_otu) %>%
  gather(id, abun, -otu) %>%
  left_join(map) %>%
  mutate(abun=if_else(abun>0, 1, 0),
         group=if_else(otu %in% rownames(MRC1_all), "2mm", "other")) %>%
  group_by(id, Size_fraction, group) %>%
  summarise(count_otu=sum(abun)) %>%
  ggplot(aes(x=as.factor(Size_fraction), y=count_otu, alluvium=group)) +
  geom_alluvium(aes(fill=group, col=group))+
  theme_classic()+
  labs(y='#OTUs', x='Size (mm)')

data.frame(otu=rownames(MRC_otu), MRC_otu) %>%
  gather(id, abun, -otu) %>%
  left_join(map) %>%
  mutate(abun=if_else(abun>0, 1, 0),
         group=if_else(otu %in% rownames(MRC7_all), "<0.056mm", "other")) %>%
  group_by(id, Size_fraction, group) %>%
  summarise(count_otu=sum(abun)) %>%
  ggplot(aes(x=as.factor(Size_fraction), y=count_otu, alluvium=group)) +
  geom_alluvium(aes(fill=group, col=group))+
  theme_classic()+
  labs(y='#OTUs', x='Size (mm)')

library(pheatmap)

MRC_abun <- OTU.rare[,map$Site=='MRC']
MRC_otu_v1 <- MRC_otu[rowSums(MRC_otu)>0.001,]
pheatmap(as.matrix(MRC_otu_v1), cluster_rows = T, cluster_cols = F, cutree_rows = 7, scale='row',show_rownames = F)

S.rel.abun.all <- decostand(SVERC_otu, method = 'total', MARGIN = 2)

SVERC_otu_v1 <- S.rel.abun.all[rowSums(S.rel.abun.all)>0.0001,]
pheatmap(as.matrix(SVERC_otu_v1), cluster_rows = T, cluster_cols = F, cutree_rows = 7, scale='row',show_rownames = F)
