module load GNU/7.3.0-2.30  OpenMPI/3.1.1-CUDA
module load R/3.6.0-X11-20180604

# Analysis in R
setwd("/mnt/home/stopnise/amoA_MiSeq_sequencing")
library(dada2)
library(phyloseq)

path="~/amoA_MiSeq_sequencing/cutadapt_before_DADA2"
#list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_trim_R1.fq", full.names = TRUE)) #pattern can be modify
fnRs <- sort(list.files(path, pattern="_trim_R2.fq", full.names = TRUE)) #pattern can be modify
if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not mach")
	
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)

# Place filtered files in filtered/ subdirectory
#filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt_F.fastq.gz"))
#filtRs <- file.path(path, "filtered", paste0(sample.names, "_filt_R.fastq.gz"))

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
trackAOA=track
head(trackAOA)

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

AA.aoa=read.table("Documents/git/SoilAggregates/AOA/ASV_AA.tab", header=F)
names(AA.aoa)[1]='OTU'
head(AA.aoa)

nucl.aoa=read.table("Documents/git/SoilAggregates/AOA/ASV_NA.tab", header=F)

toFilter=read.delim('Documents/git/SoilAggregates/AOA/aoaAA_toRemove.txt', header=F)

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
write.table(data.frame(x=AA.aoa_taxa$AAasv, y=AA.aoa_taxa$V2),"Documents/git/SoilAggregates/AOA/AOAamino.txt", sep='\t')

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

mapAOA=read.table('Documents/git/SoilAggregates/AOB/map_AOB.txt', sep='\t', header=T, row.names=1)
rownames(mapAOA) <- sapply(strsplit(rownames(mapAOA), "_AOB"), `[`, 1)
mapAOA <- mapAOA[rownames(mapAOA)%in%colnames(ASVaoa.rare),]

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

RichnessAOA <- ggplot(mapAOA, aes(x=as.factor(Size),y=Richness, col=Site, group=Site)) +
  #geom_boxplot() +
  geom_point(size=1)+
  geom_smooth(method = "lm", se = FALSE, aes(col=Site))+
  scale_fill_manual(values = c("#009e73", "#0072b2")) + 
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  labs(x='Size (mm)', y='Richness') +
  theme_classic()

Shannon <- ggplot(mapAOA, aes(x=as.factor(Size),y=Shannon, col=Site, group=Site)) +
  #geom_boxplot() +
  geom_point(size=1)+
  geom_smooth(method = "lm", se = FALSE, aes(col=Site))+
  scale_fill_manual(values = c("#009e73", "#0072b2")) + 
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  labs(x='Size (mm)', y='Shannon') +
  theme_classic()


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
  geom_point(aes(col=Site, size=Size, shape=Rep))+
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  #scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(legend.position = c(.5,.5),
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))+
  labs(x=paste('PCoA1 (',100*round(AOAax1.bc.asv.b,3),'%)',sep=''),y=paste('PCoA2 (',100*round(AOAax2.bc.asv.b,3),'%)', sep=''), 
       alpha='Size (mm)', title='Bray-Curtis')


TREE = read_tree("Documents/git/SoilAggregates/AOA/AOA_amino.tre")
AOAphyloseq=phyloseq(otu_table(as.matrix(AOAtable), taxa_are_rows = T), sample_data(mapAOA), TREE)
AOAphyloseq_norm=phyloseq_transform_css(AOAphyloseq)

plot_tree(AOAphyloseq_norm, ladderize="left", nodelabf=nodeplotboot(), color="Site", size="Size") 

AOA.taxa.cor <- taxa.env.correlation(AOAphyloseq, grouping_column="Site", method="pearson", pvalue.threshold=0.05,
                                     padjust.method="BH", adjustment=5, num.taxa=195, select.variables=NULL)
p <- plot_taxa_env(AOA.taxa.cor)
print(p)

