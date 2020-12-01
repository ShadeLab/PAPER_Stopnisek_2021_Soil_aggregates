library(RColorBrewer)
library(viridis)
library(microbiome)
library(reshape2)
library(ggpubr)
library(broom)
library(dplyr)
library(phyloseq)
library(decontam)
library(vegan)
library(tidyverse)
library(indicspecies)
library(lemon)
library(ggalluvial)
library(pheatmap)
library(metagenomeSeq)
library(qiime2R)
library(metagMisc)

setwd('~/Documents/git/SoilAggregates/')

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
  geom_point(aes(col=Site, size=Size_fraction))+
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  #scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(legend.position = c(.5,.5),
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))+
  labs(x=paste('PCoA1 (',100*round(ax1.bc.otu,3),'%)',sep=''),y=paste('PCoA2 (',100*round(ax2.bc.otu,3),'%)', sep=''), 
       alpha='Size (mm)', title='Bray-Curtis')

pcoa.J <- ggplot(map_filtered, aes(x=Axis1.J, y=Axis2.J)) +
  theme_classic() +
  geom_point(aes(col=Site, size=Size_fraction))+
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  #scale_alpha_continuous(range = c(0.1, 1)) + 
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
  geom_point(aes(col=Site, size=Size_fraction))+
  scale_color_manual(values = c("#0072b2")) + 
  #scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(legend.position = 'none')+
  labs(x=paste('PCoA1 (',100*round(ax1.bc.otu.S,3),'%)',sep=''),y=paste('PCoA2 (',100*round(ax2.bc.otu.S,3),'%)', sep=''), 
       alpha='Size (mm)', title='Bray-Curtis')

pcoa.J.S <- ggplot(mapS, aes(x=Axis1.J.S, y=Axis2.J.S)) +
  theme_classic() +
  geom_point(aes(col=Site, size=Size_fraction))+
  scale_color_manual(values = c("#0072b2")) + 
  #scale_alpha_continuous(range = c(0.1, 1)) + 
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
  geom_point(aes(col=Site, size=Size_fraction))+
  scale_color_manual(values = c("#009e73")) + 
  scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(legend.position = 'none')+
  labs(x=paste('PCoA1 (',100*round(ax1.bc.otu.M,3),'%)',sep=''),y=paste('PCoA2 (',100*round(ax2.bc.otu.M,3),'%)', sep=''), 
       alpha='Size (mm)', title='Bray-Curtis')

pcoa.J.M <- ggplot(mapM, aes(x=Axis1.J.M, y=Axis2.J.M)) +
  theme_classic() +
  geom_point(aes(col=Site, size=Size_fraction))+
  scale_color_manual(values = c("#009e73")) + 
  #scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(legend.position = 'none')+
  labs(x=paste('PCoA1 (',100*round(ax1.j.otu.M,3),'%)',sep=''),y=paste('PCoA2 (',100*round(ax2.j.otu.M,3),'%)', sep=''), 
       alpha='Size (mm)', title="Jaccard")

ggarrange(pcoa.BC.M,pcoa.BC.S ,pcoa.J.M, pcoa.J.S,
          labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2) %>%
  ggexport(filename = "PCoA_site.pdf", width = 4.5, height = 4)

#' Investigating the dynamics in AOA and OAB communities

map_filtered$sampleID <- rownames(map_filtered)
names(map_filtered)

AOA_AOB_plot<- data.frame(OTU=rownames(rel.abun.all), rel.abun.all) %>%
  gather(sampleID, abun, -OTU) %>%
  left_join(tax_filtered) %>%
  left_join(map_filtered) %>%
  mutate(taxo= paste(Phylum, Class,Order, Family, Genus, Species, sep='.')) %>%
  filter(str_detect(string = taxo, pattern = c('Crenarcha','Nitroso','Nitrososphae','Nitrosotal', 'Nitrosomonas', 'Nitrosococcus', 'Nitrosospira', 'Nitrosovibrio', 
         'Nitrosolobus'))) %>%
  group_by(Phylum, Site, Size_fraction,Replicate, NH4, OM, N) %>%
  summarise(n_abun=sum(abun)/length(unique(sampleID))) %>%
  ggplot(aes(y=n_abun, x=NH4, group=Phylum, col=Phylum)) +
  geom_point(aes(size=as.factor(Size_fraction))) +
  scale_color_manual(values = c("black", "#0072b2"), labels = c("AOA (n=9)", "AOB (n=12)")) + 
  labs(x='NH4 (ppm)', y="Relative abundance", size = 'Size (mm)')+
  geom_smooth(method = "lm", se = F, aes(col=Phylum))+
  facet_wrap(~Site, scales = 'free') +
  theme_classic()
  
data.frame(OTU=rownames(rel.abun.all), rel.abun.all) %>%
  gather(sampleID, abun, -OTU) %>%
  left_join(tax_filtered) %>%
  left_join(map_filtered) %>%
  mutate(taxo= paste(Phylum, Class,Order, Family, Genus, Species, sep='.')) %>%
  filter(str_detect(string = taxo, pattern = c('Crenarcha','Nitroso','Nitrososphae','Nitrosotal', 'Nitrosomonas', 'Nitrosococcus', 'Nitrosospira', 'Nitrosovibrio', 
                                               'Nitrosolobus')),
         abun>0) %>%
  group_by(Site, Phylum, Family) %>%
  summarise(n_taxa=length(unique(OTU)))

AOA_AOB_genus_plot <- data.frame(OTU=rownames(rel.abun.all), rel.abun.all) %>%
  gather(sampleID, abun, -OTU) %>%
  left_join(tax_filtered) %>%
  left_join(map_filtered) %>%
  mutate(taxo= paste(Phylum, Class,Order, Family, Genus, Species, sep='.')) %>%
  filter(str_detect(string = taxo, pattern = c('Crenarcha','Nitroso','Nitrososphae','Nitrosotal', 'Nitrosomonas', 'Nitrosococcus', 'Nitrosospira', 'Nitrosovibrio', 
                                               'Nitrosolobus'))) %>%
  group_by(Phylum, Family, Genus, Site, Size_fraction,Replicate, NH4, OM, N) %>%
  summarise(n_abun=sum(abun)/length(unique(sampleID))) %>%
  ggplot(aes(y=n_abun, x=NH4, group=Genus, col=Family)) +
  geom_point(aes(size=as.factor(Size_fraction))) +
  #scale_color_manual(values = c("black", "#0072b2"), labels = c("AOA (n=9)", "AOB (n=12)")) + 
  labs(x='NH4 (ppm)', y="Relative abundance", size = 'Size (mm)')+
  geom_smooth(method = "lm", se = F, aes(col=Family))+
  facet_wrap(Site~Phylum, scales = 'free') +
  theme_classic()

#' Investigating the soil chemical characteristics
CN_correlation <- ggplot(map_filtered, aes(x=N, y=OM, col=Site))+
  geom_point(aes(size=Size_fraction)) +
  #scale_alpha_continuous(range = c(0.1, 1)) + 
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

#' Species turnover analysis
#' Using betapart package
#' Andres Baselga et al

library(betapart)

# pairwise comparison  
BinaryOTU <- 1*((OTU.rare>0)==1)
#OTU.rare rarefied dataset
#otu_filtered_complete raw dataset

m.OTU <- BinaryOTU[,1:21]
m.OTU <- m.OTU[rowSums(m.OTU)>0,]

s.OTU <- BinaryOTU[,22:42]
s.OTU <- s.OTU[rowSums(s.OTU)>0,]

pair.jac=beta.pair(t(BinaryOTU),index.family = 'jaccard')

m.core=betapart.core(t(m.OTU))
m.multi=beta.multi(m.core)

pair.jac.m=beta.pair(t(m.OTU),index.family = 'jaccard')

pair.jac.s=beta.pair(t(s.OTU),index.family = 'jaccard')


# $ beta.jtu is turnover, beta.jne is nestedness, beta.jac is combined Jaccard
#rar.16s.cohort.ex.j$beta.jtu

distmelt <- function(d){
  d.melt <- d %>% as.matrix %>% as.data.frame %>% tibble::rownames_to_column() %>% melt()
  d.melt$variable <- as.character(d.melt$variable)
  return(d.melt[d.melt$variable > d.melt$rowname, ])
}

jaccard.m.melt <- pair.jac.m$beta.jac %>% distmelt()
turnover.m.melt <- pair.jac.m$beta.jtu %>% distmelt()
nestedness.m.melt <- pair.jac.m$beta.jne %>% distmelt()

jaccard.s.melt <- pair.jac.s$beta.jac %>% distmelt()
turnover.s.melt <- pair.jac.s$beta.jtu %>% distmelt()
nestedness.s.melt <- pair.jac.s$beta.jne %>% distmelt()

jaccard.m.melt$Type <- "Total"
turnover.m.melt$Type <- "Turnover"
nestedness.m.melt$Type <- "Nestedness"

jaccard.s.melt$Type <- "Total"
turnover.s.melt$Type <- "Turnover"
nestedness.s.melt$Type <- "Nestedness"


m.melt <- rbind(jaccard.m.melt, turnover.m.melt, nestedness.m.melt)

s.melt <- rbind(jaccard.s.melt, turnover.s.melt, nestedness.s.melt)

m.melt %>% ggplot(aes(y = variable, x = rowname, fill = value)) +
  geom_raster() + facet_grid(~Type) + theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank()) +
  #scale_fill_viridis(limits=c(0,1)) +
  scale_fill_gradient(low = "black", high = "#009e73", limits=c(0,1))
 
s.melt %>% ggplot(aes(y = variable, x = rowname, fill = value)) +
  geom_raster() + facet_grid(~Type) + theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank()) +
  #scale_fill_viridis(limits=c(0,1), option = "A") +
  scale_fill_gradient(low = "black", high = "#0072b2",limits=c(0,1))

# Add columns listing the categories that the samples are in. Match them using grep.
fix_rows <- function(df){
  
  try(
    # try() to add nice labels, but don't mentioned it if there is no $Type
    df$Type <- factor(df$Type, levels = c("Total", "Nestedness", "Turnover")), silent = T
  )
  
  df$row_cat <- 2.000
  df$row_cat[grepl(".min", df$rowname, fixed = T)] <- 0.050
  df$row_cat[grepl(".005", df$rowname, fixed = T)] <- 0.056
  df$row_cat[grepl(".01", df$rowname, fixed = T)] <- 0.180
  df$row_cat[grepl(".02", df$rowname, fixed = T)] <- 0.250
  df$row_cat[grepl(".05", df$rowname, fixed = T)] <- 0.500
  df$row_cat[grepl(".1", df$rowname, fixed = T)] <- 1.000

  df$var_cat <- 2.000
  df$var_cat[grepl(".min", df$variable, fixed = T)] <- 0.050
  df$var_cat[grepl(".005", df$variable, fixed = T)] <- 0.056
  df$var_cat[grepl(".01", df$variable, fixed = T)] <- 0.180
  df$var_cat[grepl(".02", df$variable, fixed = T)] <- 0.250
  df$var_cat[grepl(".05", df$variable, fixed = T)] <- 0.500
  df$var_cat[grepl(".1", df$variable, fixed = T)] <- 1.000

  df$var_cat <- factor(df$var_cat, levels = c(0.050,0.056,0.180,0.250,0.500,1.000,2.000))
  
  return(df)
}

m.melt.adj <- m.melt %>% fix_rows

s.melt.adj <- s.melt %>% fix_rows()


plot_dm <- function(df, ylab = "Sampling Day"){
  return(
    ggplot(df, aes(y = variable, x = rowname, fill = value)) + geom_raster() +
      theme(strip.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = c(.95,.4)) +
      labs(y = ylab, fill = "Binary \n Jaccard \nDistance") +
      facet_grid(var_cat ~ Type + row_cat, scales = "free", switch = "y")
  )
}

m.melt.adj %>%
  plot_dm(ylab = "Size (mm)") +
  scale_fill_viridis(limits=c(0,1))

s.melt.adj %>%
  plot_dm(ylab = "Size (mm)") +
  scale_fill_viridis(limits=c(0,1), option = "A")

#' Correlates soil particle size better with nestedness or turnover?
#' Using adonis() function to test how much variation of nestedness and turnover can be attributed to soil particle size

m.map=map_filtered[map_filtered$Site == 'MRC',]
s.map=map_filtered[map_filtered$Site != 'MRC',]

day.16s <- rbind(
  adonis(pair.jac.m$beta.jac ~ as.factor(Size_fraction), m.map)$aov.tab %>% data.frame %>% filter(Df == '6') %>% data.frame(Site= 'MRC',Type = "Total"), 
  adonis(pair.jac.m$beta.jtu ~ as.factor(Size_fraction), m.map)$aov.tab %>% data.frame %>% filter(Df == "6") %>% data.frame(Site= 'MRC',Type = "Turnover"), 
  adonis(pair.jac.m$beta.jne ~ as.factor(Size_fraction), m.map)$aov.tab %>% data.frame %>% filter(Df == "6") %>% data.frame(Site= 'MRC',Type = "Nestedness"),
  adonis(pair.jac.s$beta.jac ~ as.factor(Size_fraction), s.map)$aov.tab %>% data.frame %>% filter(Df == '6') %>% data.frame(Site= 'SVERC',Type = "Total"), 
  adonis(pair.jac.s$beta.jtu ~ as.factor(Size_fraction), s.map)$aov.tab %>% data.frame %>% filter(Df == "6") %>% data.frame(Site= 'SVERC',Type = "Turnover"), 
  adonis(pair.jac.s$beta.jne ~ as.factor(Size_fraction), s.map)$aov.tab %>% data.frame %>% filter(Df == "6") %>% data.frame(Site= 'SVERC',Type = "Nestedness")
)



#' How many taxa appear at a timepoint that were not in any previous size fraction 
#' 
OTUdf=otu_table(as.matrix(rel.abun.all), taxa_are_rows = T)
taxPhylo=tax_filtered
rownames(taxPhylo)=taxPhylo$OTU
taxPhylo$OTU=NULL
TAXdf=tax_table(as.matrix(taxPhylo))
MAPdf=sample_data(map_filtered)

PhyloData=phyloseq(OTUdf,TAXdf,MAPdf)

M.rare=subset_samples(PhyloData, Site == 'MRC')
S.rare=subset_samples(PhyloData, Site != 'MRC')

M.cohort.size <- merge_samples(M.rare, group = c("Size_fraction"))

total.taxa.M <- data.frame(Size = factor(c("0.05", "0.056", "0.180", "0.250", "0.500", "1.000", "2.000"),
                                      levels = c("0.05", "0.056", "0.180", "0.250", "0.500", "1.000", "2.000")))
total.taxa.M$Site="MRC"
total.taxa.M$`Total OTUs` <- c(
  subset_samples(M.cohort.size, Size_fraction == .05) %>% filter_taxa(function(x) max(x) > 0, TRUE) %>% ntaxa(),
  subset_samples(M.cohort.size, Size_fraction %in% c(.05, .056)) %>% filter_taxa(function(x) max(x) > 0, TRUE) %>% ntaxa(),
  subset_samples(M.cohort.size, Size_fraction %in% c(.05, .056, .180)) %>% filter_taxa(function(x) max(x) > 0, TRUE) %>% ntaxa(),
  subset_samples(M.cohort.size, Size_fraction %in% c(.05, .056, .180, .250)) %>% filter_taxa(function(x) max(x) > 0, TRUE) %>% ntaxa(),
  subset_samples(M.cohort.size, Size_fraction %in% c(.05, .056, .180, .250, .500)) %>% filter_taxa(function(x) max(x) > 0, TRUE) %>% ntaxa(),
  subset_samples(M.cohort.size, Size_fraction %in% c(.05, .056, .180, .250, .500, 1.000)) %>% filter_taxa(function(x) max(x) > 0, TRUE) %>% ntaxa(),
  subset_samples(M.cohort.size, Size_fraction %in% c(.05, .056, .180, .250, .500, 1.000, 2.000)) %>% filter_taxa(function(x) max(x) > 0, TRUE) %>% ntaxa()
)

total.taxa.M$`New OTUs` <- total.taxa.M$`Total OTUs` - lag(total.taxa.M$`Total OTUs`)

# Replace NA (all OTUs are new at the first timepoint)
total.taxa.M$`New OTUs`[is.na(total.taxa.M$`New OTUs`)] <- total.taxa.M$`Total OTUs`[is.na(total.taxa.M$`New OTUs`)]

# calculate observed turnover
total.taxa.M$`Percent New OTUs` <- total.taxa.M$`New OTUs` / total.taxa.M$`Total OTUs`


S.cohort.size <- merge_samples(S.rare, group = c("Size_fraction"))

total.taxa.S <- data.frame(Size = factor(c("0.05", "0.056", "0.180", "0.250", "0.500", "1.000", "2.000"),
                                         levels = c("0.05", "0.056", "0.180", "0.250", "0.500", "1.000", "2.000")))
total.taxa.S$Site="SVERC"

total.taxa.S$`Total OTUs` <- c(
  subset_samples(S.cohort.size, Size_fraction == .05) %>% filter_taxa(function(x) max(x) > 0, TRUE) %>% ntaxa(),
  subset_samples(S.cohort.size, Size_fraction %in% c(.05, .056)) %>% filter_taxa(function(x) max(x) > 0, TRUE) %>% ntaxa(),
  subset_samples(S.cohort.size, Size_fraction %in% c(.05, .056, .180)) %>% filter_taxa(function(x) max(x) > 0, TRUE) %>% ntaxa(),
  subset_samples(S.cohort.size, Size_fraction %in% c(.05, .056, .180, .250)) %>% filter_taxa(function(x) max(x) > 0, TRUE) %>% ntaxa(),
  subset_samples(S.cohort.size, Size_fraction %in% c(.05, .056, .180, .250, .500)) %>% filter_taxa(function(x) max(x) > 0, TRUE) %>% ntaxa(),
  subset_samples(S.cohort.size, Size_fraction %in% c(.05, .056, .180, .250, .500, 1.000)) %>% filter_taxa(function(x) max(x) > 0, TRUE) %>% ntaxa(),
  subset_samples(S.cohort.size, Size_fraction %in% c(.05, .056, .180, .250, .500, 1.000, 2.000)) %>% filter_taxa(function(x) max(x) > 0, TRUE) %>% ntaxa()
)  

total.taxa.S$`New OTUs` <- total.taxa.S$`Total OTUs` - lag(total.taxa.S$`Total OTUs`)

#' Replace NA (all OTUs are new at the first timepoint)
total.taxa.S$`New OTUs`[is.na(total.taxa.S$`New OTUs`)] <- total.taxa.S$`Total OTUs`[is.na(total.taxa.S$`New OTUs`)]

#' calculate observed turnover
total.taxa.S$`Percent New OTUs` <- total.taxa.S$`New OTUs` / total.taxa.S$`Total OTUs`

#' Join both DFs
total.all <- rbind(total.taxa.M,total.taxa.S)

#' Renaming variables
total.all$`Sum of all observed OTUs` <- total.all$`Total OTUs`
total.all$`Sum of previously unobserved OTUs` <- total.all$`New OTUs`
total.all$`Fraction of community composed of\npreviously unobserved OTUs ` <- total.all$`Percent New OTUs`

total.all$`Total OTUs` <- NULL
total.all$`New OTUs` <- NULL
total.all$`Percent New OTUs` <- NULL

total.all.melt <- melt(total.all)

total.all.melt$variable <- factor(total.all.melt$variable, levels = levels(total.all.melt$variable)[c(3,2,1)])

total.all.melt %>%
  ggplot(aes(x = Size, y = value, color = Site)) +
  geom_point() +
  geom_line(aes(group = Site), size = 1, alpha = .5) +
  facet_wrap(~variable, scales = "free_y") +
  scale_color_manual(values = c("#009e73", "#0072b2")) +
  theme_classic() +
  theme(axis.title = element_blank(),
        #        legend.position = c(.92,.6), # legend in the left panel
        legend.position = c(.5,.62), # legend in the right panel
        legend.title = element_blank()
        #, panel.spacing.x = unit(30, "pt") # add spacing between plots for equations
) 



library(DivNet)
if (!requireNamespace("BiocManager", quietly=TRUE)) 
  install.packages("BiocManager")
BiocManager::install("metagenomeSeq")
library(metagenomeSeq)

