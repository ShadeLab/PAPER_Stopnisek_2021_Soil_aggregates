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
#library(metagMisc)
library(metagenomeSeq)
library(qiime2R)

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


#' Using data generated in Oct 2020
#' OTU tables created using QIIME2, R1 only.

otu=read.table("otu_table.txt", header = T, row.names = 1)
tax=read.delim("taxonomy.tsv",row.names = 1)
map=read.csv('metadata_table.csv', row.names = 1)

tax=tax[-2]
tax_df <- colsplit(tax$Taxon, '; ', names =  c("Kingdom", "Phylum", "Class", 
                                     "Order", "Family", "Genus", "Species"))
tax_df[1:7] <- lapply(tax_df[1:7], function(x) gsub(".*__", "", x))
rownames(tax_df) <- rownames(tax)

OTU=otu_table(as.matrix(otu), taxa_are_rows = T)
TAX=tax_table(as.matrix(tax_df))
MAP=sample_data(map)

otuPhyloseq=phyloseq(OTU,TAX,MAP)

#Filtering mitochondria, Chloroplast and Unclassified taxa
otuPhyloseq_filt <- otuPhyloseq %>%
  subset_taxa(Kingdom != "Unassigned")

otuPhyloseq_filt <- otuPhyloseq_filt %>%
  subset_taxa((Family != "Mitochondria") | is.na(Family))

otuPhyloseq_filt <- otuPhyloseq_filt %>%
  subset_taxa((Class != "Chloroplast") | is.na(Class))

filtered_otus <- phyloseq_to_df(otuPhyloseq_filt, addtax=T, addtot=T)

#-------------------------------------------------------------------------
#Using decontam package to remove contamination by prevalence
sample_data(otuPhyloseq_filt)$is.neg <- sample_data(otuPhyloseq_filt)$Type == "neg"

contamdf.prev <- isContaminant(otuPhyloseq_filt, method="prevalence", neg="is.neg", threshold=0.5)
keepOTU <- rownames(contamdf.prev[contamdf.prev$contaminant=='FALSE',])
omit.OTUs <- rownames(contamdf.prev[contamdf.prev$contaminant=='TRUE',])

otuPhyloseq_filt_decont <- prune_taxa(keepOTU,otuPhyloseq_filt)
otuPhyloseq_filt_decont_samp <- subset_samples(otuPhyloseq_filt_decont, Type != 'neg')

otu_filtered <- phyloseq_to_df(otuPhyloseq_filt_decont_samp, addtax=F, addtot=T)
otu_filtered_complete <- otu_filtered[complete.cases(otu_filtered),]

tax_filtered <- phyloseq_to_df(otuPhyloseq_filt_decont_samp, addtax=T, addtot=F)
tax_filtered <- tax_filtered[1:8]

rownames(otu_filtered_complete) <- otu_filtered_complete$OTU
otu_filtered_complete$Total <- NULL
otu_filtered_complete$OTU <- NULL

#ordering samples
otu_filtered_complete <- otu_filtered_complete[,order(colnames(otu_filtered_complete))]
map_filtered=map[map$Type == 'sample',]
map_filtered=map_filtered[order(rownames(map_filtered)),]

colnames(otu_filtered_complete) == rownames(map_filtered)
#-------------------------------------------------------------------------
#Rarefiying data to 40000 reads per sample since the minimum is 42629 reads
set.seed(077)
OTU.rare <- t(rrarefy(t(otu_filtered_complete), 40000)) 
rel.abun.all <- decostand(OTU.rare, method = 'total', MARGIN = 2)


#' Alpha diversity measurements
#' 
s <- specnumber(OTU.rare,MARGIN=2)
h <- vegan::diversity(t(OTU.rare), "shannon")
pielou=h/log(s)

map_filtered$Richness <- s
map_filtered$Shannon <- h
map_filtered$Pielou <- pielou 

Richness <- ggplot(map_filtered, aes(x=as.factor(Size_fraction),y=Richness, col=Site, group=Site)) +
  #geom_boxplot() +
  geom_point(size=1)+
  geom_smooth(method = "lm", se = FALSE, aes(col=Site))+
  scale_fill_manual(values = c("#009e73", "#0072b2")) + 
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  labs(x='Size (mm)', y='Richness') +
  theme_classic()

Shannon <- ggplot(map_filtered, aes(x=as.factor(Size_fraction),y=Shannon, col=Site, group=Site)) +
  #geom_boxplot() +
  geom_point(size=1)+
  geom_smooth(method = "lm", se = FALSE, aes(col=Site))+
  scale_fill_manual(values = c("#009e73", "#0072b2")) + 
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  labs(x='Size (mm)', y='Shannon') +
  theme_classic()

ggarrange(Richness, Shannon, labels = c("A", "B"), nrow = 2) %>%
  ggexport(filename = "alpha_div.pdf", width = 4, height = 4)

#' Beta diversity
#' Combined dataset
otu.BC <- vegdist(t(OTU.rare), method="bray")
otu.J <- vegdist(t(OTU.rare), method="jaccard")
otu.bc.pcoa <- cmdscale(otu.BC, eig=T)
otu.j.pcoa <- cmdscale(otu.J, eig=T)

map_filtered$Axis1.BC <- otu.bc.pcoa$points[,1]
map_filtered$Axis2.BC <- otu.bc.pcoa$points[,2]
map_filtered$Axis1.J <- otu.j.pcoa$points[,1]
map_filtered$Axis2.J <- otu.j.pcoa$points[,2]
ax1.bc.otu <- otu.bc.pcoa$eig[1]/sum(otu.bc.pcoa$eig)
ax2.bc.otu <- otu.bc.pcoa$eig[2]/sum(otu.bc.pcoa$eig)
ax1.j.otu <- otu.j.pcoa$eig[1]/sum(otu.j.pcoa$eig)
ax2.j.otu <- otu.j.pcoa$eig[2]/sum(otu.j.pcoa$eig)

pcoa.BC <- ggplot(map_filtered, aes(x=Axis1.BC, y=Axis2.BC)) +
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

pcoa.J <- ggplot(map_filtered, aes(x=Axis1.J, y=Axis2.J)) +
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

ggarrange(pcoa.BC, pcoa.J, 
                    labels = c("A", "B"), ncol = 2) %>%
  ggexport(filename = "PCoA.pdf", width = 5.5, height = 3)

#' Beta diversity on each soil separately
#' SVERC
otu.S=otu_filtered_complete[,as.character(map_filtered$Site)=='SVERC']
mapS=map_filtered[map_filtered$Site=='SVERC',]
OTU.S.rare <- t(rrarefy(t(otu.S), 40000)) 

otu.S.BC <- vegdist(t(OTU.S.rare), method="bray")
otu.S.J <- vegdist(t(OTU.S.rare), method="jaccard")
otu.S.bc.pcoa <- cmdscale(otu.S.BC, eig=T)
otu.S.j.pcoa <- cmdscale(otu.S.J, eig=T)

mapS$Axis1.BC.S <- otu.S.bc.pcoa$points[,1]
mapS$Axis2.BC.S <- otu.S.bc.pcoa$points[,2]
mapS$Axis1.J.S <- otu.S.j.pcoa$points[,1]
mapS$Axis2.J.S <- otu.S.j.pcoa$points[,2]
ax1.bc.otu.S <- otu.S.bc.pcoa$eig[1]/sum(otu.S.bc.pcoa$eig)
ax2.bc.otu.S <- otu.S.bc.pcoa$eig[2]/sum(otu.S.bc.pcoa$eig)
ax1.j.otu.S <- otu.S.j.pcoa$eig[1]/sum(otu.S.j.pcoa$eig)
ax2.j.otu.S <- otu.S.j.pcoa$eig[2]/sum(otu.S.j.pcoa$eig)

pcoa.BC.S <- ggplot(mapS, aes(x=Axis1.BC.S, y=Axis2.BC.S)) +
  theme_classic() +
  geom_point(aes(col=Site, alpha=Size_fraction), size=4)+
  scale_color_manual(values = c("#0072b2")) + 
  scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(legend.position = 'none')+
  labs(x=paste('PCoA1 (',100*round(ax1.bc.otu.S,3),'%)',sep=''),y=paste('PCoA2 (',100*round(ax2.bc.otu.S,3),'%)', sep=''), 
       alpha='Size (mm)', title='Bray-Curtis')

pcoa.J.S <- ggplot(mapS, aes(x=Axis1.J.S, y=Axis2.J.S)) +
  theme_classic() +
  geom_point(aes(col=Site, alpha=Size_fraction), size=4)+
  scale_color_manual(values = c("#0072b2")) + 
  scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(legend.position = 'none')+
  labs(x=paste('PCoA1 (',100*round(ax1.j.otu.S,3),'%)',sep=''),y=paste('PCoA2 (',100*round(ax2.j.otu.S,3),'%)', sep=''), 
       alpha='Size (mm)', title="Jaccard")

#' MRC
otu.M=otu_filtered_complete[,as.character(map_filtered$Site)!='SVERC']
mapM=map_filtered[map_filtered$Site!='SVERC',]
OTU.M.rare <- t(rrarefy(t(otu.M), 40000)) 

otu.M.BC <- vegdist(t(OTU.M.rare), method="bray")
otu.M.J <- vegdist(t(OTU.M.rare), method="jaccard")
otu.M.bc.pcoa <- cmdscale(otu.M.BC, eig=T)
otu.M.j.pcoa <- cmdscale(otu.M.J, eig=T)

mapM$Axis1.BC.M <- otu.M.bc.pcoa$points[,1]
mapM$Axis2.BC.M <- otu.M.bc.pcoa$points[,2]
mapM$Axis1.J.M<- otu.M.j.pcoa$points[,1]
mapM$Axis2.J.M <- otu.M.j.pcoa$points[,2]
ax1.bc.otu.M <- otu.M.bc.pcoa$eig[1]/sum(otu.M.bc.pcoa$eig)
ax2.bc.otu.M <- otu.M.bc.pcoa$eig[2]/sum(otu.M.bc.pcoa$eig)
ax1.j.otu.M <- otu.M.j.pcoa$eig[1]/sum(otu.M.j.pcoa$eig)
ax2.j.otu.M <- otu.M.j.pcoa$eig[2]/sum(otu.M.j.pcoa$eig)

pcoa.BC.M <- ggplot(mapM, aes(x=Axis1.BC.M, y=Axis2.BC.M)) +
  theme_classic() +
  geom_point(aes(col=Site, alpha=Size_fraction), size=4)+
  scale_color_manual(values = c("#009e73")) + 
  scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(legend.position = 'none')+
  labs(x=paste('PCoA1 (',100*round(ax1.bc.otu.M,3),'%)',sep=''),y=paste('PCoA2 (',100*round(ax2.bc.otu.M,3),'%)', sep=''), 
       alpha='Size (mm)', title='Bray-Curtis')

pcoa.J.M <- ggplot(mapM, aes(x=Axis1.J.M, y=Axis2.J.M)) +
  theme_classic() +
  geom_point(aes(col=Site, alpha=Size_fraction), size=4)+
  scale_color_manual(values = c("#009e73")) + 
  scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(legend.position = 'none')+
  labs(x=paste('PCoA1 (',100*round(ax1.j.otu.M,3),'%)',sep=''),y=paste('PCoA2 (',100*round(ax2.j.otu.M,3),'%)', sep=''), 
       alpha='Size (mm)', title="Jaccard")

ggarrange(pcoa.BC.M,pcoa.BC.S ,pcoa.J.M, pcoa.J.S,
          labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2) %>%
  ggexport(filename = "PCoA_site.pdf", width = 4.5, height = 4)

#' Taxonomic analysis
#' 

tax_filtered %>%
  group_by(Kingdom,Phylum, Class, Family, Genus) %>%
  filter(Kingdom == "Archaea") %>%
  summarize(n_taxa=length(unique(OTU)))

map_filtered$sampleID <- rownames(map_filtered)
names(map_filtered)

(fig_Arch <- data.frame(OTU=rownames(rel.abun.all), rel.abun.all) %>%
  gather(sampleID, abun, -OTU) %>%
  left_join(tax_filtered) %>%
  left_join(map_filtered) %>%
  filter(Phylum == 'Crenarchaeota') %>%
  group_by(Phylum, Site, Size_fraction,Replicate, OM, N) %>%
  summarise(n_abun=sum(abun)/length(unique(sampleID))) %>%
  ggplot(aes(y=n_abun, x=OM/N, group=Site, col=Site)) +
  geom_point(aes(alpha=Size_fraction), size=2.5) +
  scale_alpha_continuous(range = c(0.1, 1)) + 
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  labs(x='C:N', y="Relative abundance", alpha = 'Size (mm)')+
  geom_smooth(method = "lm", se = F, aes(col=Site))+
  theme_classic()
)

names(map_filtered)
CN_correlation <- ggplot(map_filtered, aes(x=N, y=OM, col=Site))+
  geom_point(aes(alpha=Size_fraction), size=2.5) +
  scale_alpha_continuous(range = c(0.1, 1)) + 
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  labs(x='N (%)', y= 'Organic Matter (%)', alpha= 'Size (mm)') +
  geom_smooth(method = "lm", se = F, aes(col=Site))+
  theme_classic()+ 
  theme(legend.position = 'none')

NH4 <- ggplot(map_filtered, aes(x=as.factor(Size_fraction), y=NH4, col=Site))+
  geom_point(size=2.5) +
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  labs(x='Size (mm)', y= 'NH4 (ppm)') +
  theme_classic()+
  theme(legend.position = 'none')

NO3 <- ggplot(map_filtered, aes(x=as.factor(Size_fraction), y=NO3, col=Site))+
  geom_point(size=2.5) +
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  labs(x='Size (mm)', y= 'NO3 (ppm)') +
  theme_classic()+
  theme(legend.position = 'none')

CN <- ggplot(map_filtered, aes(x=as.factor(Size_fraction), y=OM/N, col=Site))+
  geom_point(size=2.5) +
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  labs(x='Size (mm)', y= 'C:N') +
  theme_classic() +
  theme(legend.position = 'none')

ggarrange(NH4,NO3,CN, CN_correlation,
          labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2) %>%
  ggexport(filename = "soil_parameters.pdf", width = 4.5, height = 4)
