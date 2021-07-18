module load GNU/7.3.0-2.30  OpenMPI/3.1.1-CUDA
module load R/3.6.0-X11-20180604

# Analysis in R
# Read processing
setwd("/mnt/home/stopnise/amoA_MiSeq_sequencing")
library(dada2)
library(phyloseq)

path="~/amoA_MiSeq_sequencing/cutadapt_before_DADA2"

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_trim_R1.fq", full.names = TRUE)) #pattern can be modify
fnRs <- sort(list.files(path, pattern="_trim_R2.fq", full.names = TRUE)) #pattern can be modify
if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not mach")
	
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered_230", paste0(sample.names, "_filt_F.fastq.gz"))
filtRs <- file.path(path, "filtered_230", paste0(sample.names, "_filt_R.fastq.gz"))

#Set truncLen and minLen according to your dataset
out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,230), minLen = 230, truncQ=2, maxN=0, maxEE=c(2,2), compress=TRUE, multithread=TRUE, rm.phix=TRUE)

errF <- dada2::learnErrors(filtFs, multithread=TRUE)
errR <- dada2::learnErrors(filtRs, multithread=TRUE)

dadaFs <- dada2::dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada2::dada(filtRs, err=errR, multithread=TRUE)

mergers <- dada2::mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=TRUE)

seqtab <- dada2::makeSequenceTable(mergers, orderBy = "abundance")
dim(seqtab)

table(nchar(getSequences(seqtab)))

#remove chimera
seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#save ASV as fasta
uniquesToFasta(getUniques(seqtab.nochim), fout="ASV.fa", ids=paste0("ASV", seq(length(getUniques(seqtab.nochim)))))

#rename the ASV from read sequence to ASV#
asvTableAOA=seqtab.nochim
colnames(asvTable)=paste0("ASV", seq(length(getUniques(asvTable))))

#How many chimeric reads are there?
sum(asvTableAOA)/sum(seqtab)

# creating summary table
getN <- function(x) sum(dada2::getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(asvTableAOA))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track <- data.frame(sampleID = row.names(track), track)
track = tibble::remove_rownames(track)
track$sampleID = sapply(strsplit(track$sampleID, "_AOA"), `[`, 1)


#create ASV table for phyloseq
samples.out <- rownames(asvTableAOA)
subject <- sapply(strsplit(samples.out, "_AOA"), `[`, 1)
site <- substr(samples.out,1,1)

samdf <- data.frame(Sample=samples.out, Site=site)
write.table(samdf, 'map_AOA.txt', sep='\t')

mapAOA=read.table('map_AOA.txt', sep='\t', header=T, row.names=1)

psAOA <- phyloseq::phyloseq(otu_table(asvTableAOA, taxa_are_rows=FALSE), 
               sample_data(mapAOA))
			   
#rename ASV from sequences to ASV#
dna <- Biostrings::DNAStringSet(taxa_names(psAOA))
names(dna) <- taxa_names(psAOA)
psAOA <- merge_phyloseq(psAOA, dna)
taxa_names(psAOA) <- paste0("ASV", seq(ntaxa(psAOA)))

#' Joining nucleic acid ASV based on their amino acid sequences (translated using seqkit tool)
#' command for seqkit
#' seqkit translate -T 11 -f 2 ASV.fa --clean > ASV_AA.fa
#' where -T 11 represent frame usage for bacteria
#' -f 2 start from2nd base in the sequence
#' --clean removes the stop codon symbol * into X
library(metagMisc)
library(tidyverse)

AOA_table=phyloseq_to_df(psAOA, addtax=F)
head(AOA_table)

AA.aoa=read.table("~/Documents/git/SoilAggregates/AOA/ASV_AA.tab", header=F)
names(AA.aoa)[1]='OTU'
head(AA.aoa)

nucl.aoa=read.table("~/Documents/git/SoilAggregates/AOA/ASV_NA.tab", header=F)

toFilter=read.delim('~/Documents/git/SoilAggregates/AOA/aoaAA_toRemove.txt', header=F)

AA.aoa_taxa=AA.aoa %>% 
  left_join(AOA_table) %>%
  pivot_longer(!c(V2, OTU, ), names_to = "sampleID", values_to = "count") %>%
  filter(!(OTU %in% toFilter$V1)) %>%
  group_by(V2, sampleID) %>%
  summarise(count=sum(count)) %>%
  dplyr::filter(!grepl("X", V2)) %>%
  pivot_wider(names_from = sampleID,
              values_from = count) 

AA.aoa_taxa$AAasv=paste0("AA.ASV", seq(length(AA.aoa_taxa$V2)))
colnames(AA.aoa_taxa)=str_remove_all(colnames(AA.aoa_taxa),"_trim_filt_F.fastq.gz") 

#Writing out the unique and filtered AA sequences for building phylogenetic tree
write.table(data.frame(x=AA.aoa_taxa$AAasv, y=AA.aoa_taxa$V2),"~/Documents/git/SoilAggregates/AOA/AOAamino.txt", sep='\t')

#' **************************************************
#' Statistical analysis
#' **************************************************

#' Creating an ASV table (ASV based on the AA sequeces)
AOAtable=AA.aoa_taxa[-1] 
samples.out <- colnames(AOAtable)
subject <- sapply(strsplit(samples.out, "_AOA"), `[`, 1)
colnames(AOAtable)=subject

AOAtable=as.data.frame(AOAtable) 
rownames(AOAtable)=AOAtable$AAasv
AOAtable$AAasv=NULL

AOAtable=as.matrix(AOAtable)

library(vegan)
set.seed(335)

ASVaoa.rare <- t(rrarefy(t(AOAtable), 25000)) #M1_a and S_5_c not enough reads
ASVaoa.rare <- ASVaoa.rare[,colSums(ASVaoa.rare)>24999]
length(rownames(ASVaoa.rare))
mapAOA=read.table('~/Documents/git/SoilAggregates/AOB/map_AOB.txt', sep='\t', header=T, row.names=1)
rownames(mapAOA) <- sapply(strsplit(rownames(mapAOA), "_AOB"), `[`, 1)
mapAOA <- mapAOA[rownames(mapAOA)%in%colnames(ASVaoa.rare),]
#write.table(mapAOA[1:7], "~/Documents/git/mapFile.txt")

AOAabun <- decostand(ASVaoa.rare, method = 'total', MARGIN = 2)

rownames(mapAOA) == colnames(ASVaoa.rare)

#' Alpha diversity measurements
#' 
AOAs <- specnumber(ASVaoa.rare,MARGIN=2)
AOAh <- vegan::diversity(t(ASVaoa.rare), "shannon")
AOApielou=AOAh/log(AOAs)

mapAOA$Richness <- AOAs
mapAOA$Shannon <- AOAh
mapAOA$Pielou <- AOApielou 

RichnessAOA <- ggplot(mapAOA, aes(x=as.factor(Size),y=Richness, fill=Site)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#009e73", "#0072b2")) + 
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  labs(x='Size (mm)', y='Richness', title= "AOA") +
  theme_classic()

ShannonAOA <- ggplot(mapAOA, aes(x=as.factor(Size),y=Shannon, fill=Site)) +
  geom_boxplot()+
  scale_fill_manual(values = c("#009e73", "#0072b2")) + 
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  labs(x='Size (mm)', y='Shannon', title= "AOA") +
  theme_classic()

#' Alpha diversity statistics
mapAOA %>%
  group_by(Site, Size) %>%
  summarise(minRich=min(Richness),
            maxRich=max(Richness),
            minShan=min(Shannon),
            maxShan=max(Shannon))
#' ANOVA
alpha.site.aoa.aov <- aov(Shannon ~ Site, #replace Richness with Shannon
                      data = mapAOA) 
alpha.size.aoa.aov <- aov(mapAOA$Shannon ~ Size , #replace Richness with Shannon
                      data = mapAOA) 
alpha.size.M.aoa.aov <- aov(Shannon ~ Size ,  #replace Richness with Shannon
                        data = mapAOA[mapAOA$Site == 'M',])
alpha.size.S.aoa.aov <- aov(Shannon ~ Size , #replace Richness with Shannon
                        data = mapAOA[mapAOA$Site == 'S',])

# Summary of the analysis
summary(alpha.site.aoa.aov)
summary(alpha.size.aoa.aov)
summary(alpha.size.M.aoa.aov)
summary(alpha.size.S.aoa.aov)

#' 
#' Are richness and Shannon index different at same soil size when comparing 
#' two soils?
stat.test.rich.aoa <- mapAOA %>%
  group_by(Size) %>%
  t_test(Richness ~ Site) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test.rich.aoa #Yes in all instances the two soils are sign. different

stat.test.shan.aoa <- mapAOA %>%
  group_by(Size) %>%
  t_test(Shannon ~ Site) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test.shan.aoa # Same is true for the Shannon index

#' Is alpha diversity different with soil particle size?
stat.test.rich.size.aoa <- mapAOA %>%
  group_by(Site) %>%
  t_test(Richness ~ Size) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test.rich.size.aoa # No statistical significant differences

#' Beta diversity
#' Combined dataset
AOAasv.b.BC <- vegdist(t(ASVaoa.rare), method="bray")
AOAasv.b.J <- vegdist(t(ASVaoa.rare), method="jaccard")
AOAasv.b.bc.pcoa <- cmdscale(AOAasv.b.BC, eig=T)
AOAasv.b.j.pcoa <- cmdscale(AOAasv.b.J, eig=T)

mapAOA$Axis1.BC <- AOAasv.b.bc.pcoa$points[,1]
mapAOA$Axis2.BC <- AOAasv.b.bc.pcoa$points[,2]
mapAOA$Axis1.J <- AOAasv.b.j.pcoa$points[,1]
mapAOA$Axis2.J <- AOAasv.b.j.pcoa$points[,2]
AOAax1.bc.asv.b <- AOAasv.b.bc.pcoa$eig[1]/sum(AOAasv.b.bc.pcoa$eig)
AOAax2.bc.asv.b <- AOAasv.b.bc.pcoa$eig[2]/sum(AOAasv.b.bc.pcoa$eig)
AOAax1.j.asv.b <- AOAasv.b.j.pcoa$eig[1]/sum(AOAasv.b.j.pcoa$eig)
AOAax2.j.asv.b <- AOAasv.b.j.pcoa$eig[2]/sum(AOAasv.b.j.pcoa$eig)

pcoa.BC.AOAasv.b <- ggplot(mapAOA, aes(x=Axis1.BC, y=Axis2.BC)) +
  theme_classic() +
  geom_point(aes(col=Site, size=Size))+
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  #scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(#legend.position = c(.5,.5),
        legend.position = 'none',
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))+
  labs(x=paste('PCoA1 (',100*round(AOAax1.bc.asv.b,3),'%)',sep=''),y=paste('PCoA2 (',100*round(AOAax2.bc.asv.b,3),'%)', sep=''), 
       alpha='Size (mm)', title='Bray-Curtis - AOA')

pcoa.J.AOAasv.j=ggplot(mapAOA, aes(x=Axis1.J, y=Axis2.J)) +
  theme_classic() +
  geom_point(aes(col=Site, size=Size))+
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  theme(#legend.position = c(.5,.5),
        legend.position = 'none',
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))+
  labs(x=paste('PCoA1 (',100*round(AOAax1.j.asv.b,3),'%)',sep=''),y=paste('PCoA2 (',100*round(AOAax2.j.asv.b,3),'%)', sep=''), 
       alpha='Size (mm)', title='Jaccard - AOA')

#' Beta diversity on separated dataset
#' MRC
asvAOA.M=ASVaoa.rare[,as.character(mapAOA$Site)=='M']
mapAOA.M=mapAOA[mapAOA$Site=='M',]
asvAOA.M=asvAOA.M[rowSums(asvAOA.M)>0,] # remove ASV not present in the dataset 

asvAOA.M.BC <- vegdist(t(asvAOA.M), method="bray")
asvAOA.M.J <- vegdist(t(asvAOA.M), method="jaccard")
asvAOA.M.bc.pcoa <- cmdscale(asvAOA.M.BC, eig=T)
asvAOA.M.j.pcoa <- cmdscale(asvAOA.M.J, eig=T)

mapAOA.M$Axis1.BC.M <- asvAOA.M.bc.pcoa$points[,1]
mapAOA.M$Axis2.BC.M <- asvAOA.M.bc.pcoa$points[,2]
mapAOA.M$Axis1.J.M <- asvAOA.M.j.pcoa$points[,1]
mapAOA.M$Axis2.J.M <- asvAOA.M.j.pcoa$points[,2]
ax1.bc.asvAOA.M <- asvAOA.M.bc.pcoa$eig[1]/sum(asvAOA.M.bc.pcoa$eig)
ax2.bc.asvAOA.M <- asvAOA.M.bc.pcoa$eig[2]/sum(asvAOA.M.bc.pcoa$eig)
ax1.j.asvAOA.M <- asvAOA.M.j.pcoa$eig[1]/sum(asvAOA.M.j.pcoa$eig)
ax2.j.asvAOA.M <- asvAOA.M.j.pcoa$eig[2]/sum(asvAOA.M.j.pcoa$eig)

pcoa.BC.asvAOA.M <- ggplot(mapAOA.M, aes(x=Axis1.BC.M, y=Axis2.BC.M, size=Size, col=Site)) +
  theme_classic() +
  geom_point()+
  scale_color_manual(values = c("#009e73")) + 
  theme(legend.position =c(.4,.8))+
  labs(x=paste('PCoA1 (',100*round(ax1.bc.asvAOA.M,3),'%)',sep=''),
       y=paste('PCoA2 (',100*round(ax2.bc.asvAOA.M,3),'%)', sep=''), 
       alpha='Size (mm)')

pcoa.J.asvAOA.M <- ggplot(mapAOA.M, aes(x=Axis1.J.M, y=Axis2.J.M, size=Size, col=Site)) +
  theme_classic() +
  geom_point()+
  scale_color_manual(values = c("#009e73")) + 
  theme(#legend.position =c(.4,.8),
        legend.position= 'none')+
  labs(x=paste('PCoA1 (',100*round(ax1.j.asvAOA.M,3),'%)',sep=''),
       y=paste('PCoA2 (',100*round(ax2.j.asvAOA.M,3),'%)', sep=''), 
       alpha='Size (mm)', title='AOA - Jaccard')

#' SVERC
asvAOA.S=ASVaoa.rare[,as.character(mapAOA$Site)=='S']
mapAOA.S=mapAOA[mapAOA$Site=='S',]
asvAOA.S=asvAOA.S[rowSums(asvAOA.S)>0,] # remove ASV not present in the dataset 

asvAOA.S.BC <- vegdist(t(asvAOA.S), method="bray")
asvAOA.S.J <- vegdist(t(asvAOA.S), method="jaccard")
asvAOA.S.bc.pcoa <- cmdscale(asvAOA.S.BC, eig=T)
asvAOA.S.j.pcoa <- cmdscale(asvAOA.S.J, eig=T)

mapAOA.S$Axis1.BC.S <- asvAOA.S.bc.pcoa$points[,1]
mapAOA.S$Axis2.BC.S <- asvAOA.S.bc.pcoa$points[,2]
mapAOA.S$Axis1.J.S <- asvAOA.S.j.pcoa$points[,1]
mapAOA.S$Axis2.J.S <- asvAOA.S.j.pcoa$points[,2]
ax1.bc.asvAOA.S <- asvAOA.S.bc.pcoa$eig[1]/sum(asvAOA.S.bc.pcoa$eig)
ax2.bc.asvAOA.S <- asvAOA.S.bc.pcoa$eig[2]/sum(asvAOA.S.bc.pcoa$eig)
ax1.j.asvAOA.S <- asvAOA.S.j.pcoa$eig[1]/sum(asvAOA.S.j.pcoa$eig)
ax2.j.asvAOA.S <- asvAOA.S.j.pcoa$eig[2]/sum(asvAOA.S.j.pcoa$eig)

pcoa.BC.asvAOA.S <- ggplot(mapAOA.S, aes(x=Axis1.BC.S, y=Axis2.BC.S, size=Size, col=Site)) +
  theme_classic() +
  geom_point()+
  scale_color_manual(values = c("#0072b2")) + 
  theme(legend.position ='none')+
  labs(x=paste('PCoA1 (',100*round(ax1.bc.asvAOA.S,3),'%)',sep=''),
       y=paste('PCoA2 (',100*round(ax2.bc.asvAOA.S,3),'%)', sep=''), 
       alpha='Size (mm)')

pcoa.J.asvAOA.S <- ggplot(mapAOA.S, aes(x=Axis1.J.S, y=Axis2.J.S, size=Size, col=Site)) +
  theme_classic() +
  geom_point()+
  scale_color_manual(values = c("#0072b2")) + 
  theme(#legend.position =c(.4,.8),
        legend.position='none') +
  labs(x=paste('PCoA1 (',100*round(ax1.j.asvAOA.S,3),'%)',sep=''),
       y=paste('PCoA2 (',100*round(ax2.j.asvAOA.S,3),'%)', sep=''), 
       alpha='Size (mm)', title='Jaccard - AOA')

.#' Testing the effect of size, site and chemical parameters on the community 
#' using PERMANOVA
#' Factors: Site, Size, OM, NO3, NH4, N
set.seed(002)
#' Full dataset (M and S) 
adonis(AOAasv.b.BC~mapAOA$Site) # R2=0.838, p=0.001

#' M and S separated 
adonis(asvAOA.S.BC~mapAOA.S$Size) # R2=0.551, p=0.001
adonis(asvAOA.M.BC~mapAOA.M$Size) # R2=0.033, p=0.455

#' Chemical parameters for full SVERC
adonis(asvAOA.S.BC~mapAOA.S$OM)  # R2=0.271, p=0.001
adonis(asvAOA.S.BC~mapAOA.S$N)   # R2=0.123, p=0.075
adonis(asvAOA.S.BC~mapAOA.S$NO3) # R2=0.293, p=0.003
adonis(asvAOA.S.BC~mapAOA.S$NH4) # R2=0.223, p=0.011

#' Chemical analysis for MRC where 2mm soil particles are removed because no 
#' chemistry was done for them
asvAOA.M.BC.2mm <- vegdist(t(asvAOA.M[,-c(15,16,17)]), method="bray") # removing 2 mm samples
mapAOA.M.2=mapAOA.M %>% filter(Size<2)
adonis(asvAOA.M.BC.2mm~mapAOA.M.2$OM) # R2=0.257, p=0.038
adonis(asvAOA.M.BC.2mm~mapAOA.M.2$N)  # R2=0.263, p=0.031
adonis(asvAOA.M.BC.2mm~mapAOA.M.2$NO3)# R2=0.079, p=0.248
adonis(asvAOA.M.BC.2mm~mapAOA.M.2$NH4)# R2=0.074, p=0.272

#############################
#' Calculating beta dispersal

#'MRC 
## Calculate multivariate dispersions
mod.M.aoa <- betadisper(asvAOA.M.BC, mapAOA.M$Size)
mod.M.aoa

## Perform test
anova(mod.M.aoa)

## Permutation test for F
permutest(mod.M.aoa, pairwise = TRUE, permutations = 99)

#' Tukey's Honest Significant Differences
(mod.M.aoa.HSD <- TukeyHSD(mod.M.aoa))
plot(mod.M.aoa.HSD)

#' Plot the groups and distances to centroids on the
## first two PCoA axes with data ellipses instead of hulls
plot(mod.M.aoa, ellipse = TRUE, hull = FALSE, seg.lty = "dashed") # 1 sd data ellipse

#' Draw a boxplot of the distances to centroid for each group
boxplot(mod.M.aoa)

#' SVERC 
## Calculate multivariate dispersions
mod.S.aoa <- betadisper(asvAOA.S.BC, mapAOA.S$Size)
mod.S.aoa

#' Perform test
anova(mod.S.aoa)

#' Permutation test for F
permutest(mod.S.aoa, pairwise = TRUE, permutations = 99)

#' Tukey's Honest Significant Differences
(mod.S.aoa.HSD <- TukeyHSD(mod.S.aoa))
plot(mod.S.aoa.HSD)

#' Plot the groups and distances to centroids on the
#' first two PCoA axes with data ellipses instead of hulls
plot(mod.S.aoa, ellipse = TRUE, hull = FALSE, seg.lty = "dashed") # 1 sd data ellipse

#' Draw a boxplot of the distances to centroid for each group
boxplot(mod.S.aoa)


#######################################################
#' Calculating the phylogenetic diversity (Faith's diversity index (PD))
library(btools)

TREEaoa = read_tree("AOA/AOA_amino.tre")
AOAp=phyloseq(otu_table(as.matrix(ASVaoa.rare), taxa_are_rows = T), sample_data(mapAOA), TREEaoa)

PD.aoa <- estimate_pd(AOAp) %>%
  rownames_to_column('sample_ID')

mapPD.aoa=mapAOA
mapPD.aoa$sample_ID=rownames(mapPD.aoa)

PD.aoa %>%
  left_join(mapPD.aoa) %>%
  ggplot(aes(x=factor(Size), y=PD, color=Site)) +
  geom_boxplot() +
  theme_classic()+
  scale_color_manual(values = c("#009e73", "#0072b2")) +
  labs(x=NULL) 

PD.aoa.df=PD.aoa %>%
  left_join(mapPD.aoa)

#' ANOVA
PD.aoa.site.aov <- aov(PD ~ Site, data = PD.aoa.df) 
PD.aoa.size.aov <- aov(PD ~ Size , data = PD.aoa.df)
PD.aoa.size.M.aov <- aov(PD ~ Size , data = PD.aoa.df[PD.aoa.df$Site == 'M',])
PD.aoa.size.S.aov <- aov(PD ~ Size, data = PD.aoa.df[PD.aoa.df$Site == 'S',])

# Summary of the analysis
summary(PD.aoa.site.aov)
summary(PD.aoa.size.aov)
summary(PD.aoa.size.M.aov)
summary(PD.aoa.size.S.aov)
#' Significant difference between sites but nothing else

#' Does phylogenetic signal/depth differ between soils and different soil 
#' particle sizes?
stat.test.pd.aoa <- PD.aoa.df %>%
  group_by(Size) %>%
  t_test(PD ~ Site) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test.pd.aoa  
# Difference between the MRC and SVERC only at 0.25 mm


