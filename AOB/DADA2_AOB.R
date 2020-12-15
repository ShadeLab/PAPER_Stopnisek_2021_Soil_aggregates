module load GNU/7.3.0-2.30  OpenMPI/3.1.1-CUDA
module load R/3.6.0-X11-20180604

# Analysis in R
setwd("/mnt/home/stopnise/amoA_MiSeq_sequencing/AOBpart")
library(dada2)
library(phyloseq)

path="~/amoA_MiSeq_sequencing/AOBpart/cutadapt_before_DADA2"
outpath="/mnt/research/ShadeLab/WorkingSpace/Stopnisek/Soil_aggregates"

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

filtFs <- file.path(outpath, "filtered_240", paste0(sample.names, "_filt_F.fastq.gz"))
filtRs <- file.path(outpath, "filtered_240", paste0(sample.names, "_filt_R.fastq.gz"))

#Set truncLen and minLen according to your dataset
#AOB amoA assembly: truncLen=c(229,229), minLen = 229 
out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,240), 
                            minLen = 240, truncQ=2, maxN=0, maxEE=c(2,2), compress=TRUE, 
                            rm.phix=TRUE,  multithread=TRUE)
head(out)

errF <- dada2::learnErrors(filtFs, multithread=TRUE)
errR <- dada2::learnErrors(filtRs, multithread=TRUE)

dadaFs <- dada2::dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada2::dada(filtRs, err=errR, multithread=TRUE)

mergers <- dada2::mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- dada2::makeSequenceTable(mergers, orderBy = "abundance")
dim(seqtab)

table(nchar(dada2::getSequences(seqtab)))

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 449:453] 

#remove chimera
seqtab.nochim <- dada2::removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#save ASV as fasta
dada2::uniquesToFasta(dada2::getUniques(seqtab.nochim), fout="ASV.fa", 
                      ids=paste0("ASV", seq(length(dada2::getUniques(seqtab.nochim)))))

#rename the ASV from read sequence to ASV#
asvTable=seqtab.nochim
colnames(asvTable)=paste0("ASV", seq(length(dada2::getUniques(asvTable))))
write.table(t(asvTable), "Documents/git/SoilAggregates/AOB/ASVna_table_AOB.txt", sep='\t') #ASV table output for AOB

#How many chimeric reads are there?
sum(asvTable)/sum(seqtab2)

# creating summary table
getN <- function(x) sum(dada2::getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(asvTable))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
head(track)

write.delim(x = track, 'Documents/git/SoilAggregates/AOB_dada2_readStat.txt')

library(gt)
library(dplyr)

track = as.data.frame(track)
track <- data.frame(sampleID = row.names(track), track)
track = tibble::remove_rownames(track)
track$sampleID = sapply(strsplit(track$sampleID, "_AOB"), `[`, 1)

track %>%
  gt(rowname_col = "sampleID") %>%
  tab_stubhead(label = "Sample ID") %>%
  tab_header(
    title = md("**Read processing Statistics**"),
    subtitle = "Using DADA2 package in R and amoA AOB amplicons") %>%
  tab_row_group(
    group = "Montcalm",
    rows = 1:21) %>%
  tab_row_group(
    group = "Saginaw Valley",
    rows = 22:42) %>%
  cols_label(
    input = "Input",
    filtered = "Quality\nfiltered",
    denoisedF = "Denoised forward",
    denoisedR = "Denoised reverse",
    merged = "Merged",
    nonchim = "Without chimera")

#create ASV table for phyloseq

samples.out <- rownames(asvTable)
subject <- sapply(strsplit(samples.out, "_AOB"), `[`, 1)
site <- substr(samples.out,1,1)

samdf <- data.frame(Sample=samples.out, Site=site)
write.table(samdf, 'map_AOB.txt', sep='\t') #add outside R soil size and replicate

mapAOB=read.table('Documents/git/SoilAggregates/AOB/map_AOB.txt', sep='\t', header=T, row.names=1)

ps <- phyloseq::phyloseq(otu_table(asvTable, taxa_are_rows=FALSE), 
                         sample_data(mapAOB))

#rename ASV from sequences to ASV#
library(Biostrings)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#' Joining nucleic acid ASV based on their amino acid sequences (translated using seqkit tool)
#' command for seqkit
#' seqkit translate -T 11 -f 2 ASV.fa --clean > ASV_AA.fa
#' where -T 11 represent frame usage for bacteria
#' -f 2 start from2nd base in the sequence
#' --clean removes the stop codon symbol * into X
#' Also convert .fa to .tab using seqkit fx2tab comand
#' 

AOB_table=phyloseq_to_df(ps, addtax=F)
head(AOB_table)

AAaob=read.table("Documents/git/SoilAggregates/ASV_AA_AOB.tab", header=F)
names(AAaob)[1]='OTU'

AA_aob_taxa=AAaob %>% 
  left_join(AOB_table) %>%
  pivot_longer(!c(V2, OTU), names_to = "sampleID", values_to = "count") %>%
  filter(!(OTU %in% c("ASV924","ASV959",'ASV1021','ASV1059'))) %>%
  group_by(V2, sampleID) %>%
  summarize(count=sum(count)) %>%
  dplyr::filter(!grepl("X", V2)) %>%
  pivot_wider(names_from = sampleID,
              values_from = count) 

AA_aob_taxa$AAasv=paste0("AA.ASV", seq(length(AA_aob_taxa$V2)))
colnames(AA_aob_taxa)=str_remove_all(colnames(AA_aob_taxa),"_trim_filt_F.fastq.gz") 

#Writing out amino acid sequences:
write.table(data.frame(x=AA_aob_taxa$AAasv, y=AA_aob_taxa$V2),"Documents/git/SoilAggregates/AOB/AOBamino.txt", sep='\t')
#' Writing out the amino acid ASV table
write_delim(AA_aob_taxa, 'Documents/git/SoilAggregates/AOB/AOB_AA_table.txt') 

#' **************************************************
#' Statistical analysis
#' **************************************************

#' Creating an ASV table (ASV based on the AA sequeces)
AA_aob_taxa=read.delim('Documents/git/SoilAggregates/AOB/AOB_AA_table.txt',sep = ' ')
mapAOB_raw=read.table('Documents/git/SoilAggregates/AOB/map_AOB.txt', sep='\t', header=T, row.names=1)

AOBtable=AA_aob_taxa[-1] #remove the sequences

#' need to rename the colnames and remove anything following "_AOB"
samples.out <- colnames(AOBtable)
subject <- sapply(strsplit(samples.out, "_AOB"), `[`, 1)
colnames(AOBtable)=subject
rownames(AOBtable)=AOBtable$AAasv
AOBtable$AAasv=NULL
AOBtable=as.matrix(AOBtable)

#' same for the map file
rownames(mapAOB) <- sapply(strsplit(rownames(mapAOB), "_AOB"), `[`, 1)

library(vegan)
set.seed(533)

ASVaob.rare <- t(rrarefy(t(AOBtable), 40000)) 
AOBabun <- decostand(ASVaob.rare, method = 'total', MARGIN = 2)

rownames(mapAOB) == colnames(AOBtable)

#' Alpha diversity measurements
#' 
AOBs <- specnumber(ASVaob.rare,MARGIN=2)
AOBh <- vegan::diversity(t(ASVaob.rare), "shannon")
AOBpielou=AOBh/log(AOBs)

mapAOB$Richness <- AOBs
mapAOB$Shannon <- AOBh
mapAOB$Pielou <- AOBpielou 

#' Beta diversity
#' Combined dataset
asv.b.BC <- vegdist(t(ASVaob.rare), method="bray")
asv.b.J <- vegdist(t(ASVaob.rare), method="jaccard")
asv.b.bc.pcoa <- cmdscale(asv.b.BC, eig=T)
asv.b.j.pcoa <- cmdscale(asv.b.J, eig=T)

mapAOB$Axis1.BC <- asv.b.bc.pcoa$points[,1]
mapAOB$Axis2.BC <- asv.b.bc.pcoa$points[,2]
mapAOB$Axis1.J <- asv.b.j.pcoa$points[,1]
mapAOB$Axis2.J <- asv.b.j.pcoa$points[,2]
ax1.bc.asv.b <- asv.b.bc.pcoa$eig[1]/sum(asv.b.bc.pcoa$eig)
ax2.bc.asv.b <- asv.b.bc.pcoa$eig[2]/sum(asv.b.bc.pcoa$eig)
ax1.j.asv.b <- asv.b.j.pcoa$eig[1]/sum(asv.b.j.pcoa$eig)
ax2.j.asv.b <- asv.b.j.pcoa$eig[2]/sum(asv.b.j.pcoa$eig)

pcoa.BC.asv.b <- ggplot(mapAOB, aes(x=Axis1.BC, y=Axis2.BC)) +
  theme_classic() +
  geom_point(aes(col=Site, size=Size))+
  scale_color_manual(values = c("#009e73", "#0072b2")) + 
  #scale_alpha_continuous(range = c(0.1, 1)) + 
  theme(legend.position = c(.5,.5),
        legend.background = element_rect(fill="transparent"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))+
  labs(x=paste('PCoA1 (',100*round(ax1.bc.asv.b,3),'%)',sep=''),y=paste('PCoA2 (',100*round(ax2.bc.asv.b,3),'%)', sep=''), 
       alpha='Size (mm)', title='Bray-Curtis')

#' Using phyloseq for the analysis

TREE = read_tree("Documents/git/SoilAggregates/AOB/AOB_amino.tre")
AOBp=phyloseq(otu_table(as.matrix(AOBtable), taxa_are_rows = T), sample_data(mapAOB), TREE)
AOBphyloseq_norm=phyloseq_transform_css(AOBp)

plot_tree(AOBphyloseq_norm, ladderize="left", nodelabf=nodeplotboot(), color="Site", size="Size") 



AOB.taxa.cor <- taxa.env.correlation(AOBphyloseq_norm, grouping_column="Site", method="pearson", pvalue.threshold=0.05,
                                     padjust.method="BH", adjustment=5, num.taxa=195, select.variables=NULL)
p <- plot_taxa_env(AOB.taxa.cor)
print(p)

