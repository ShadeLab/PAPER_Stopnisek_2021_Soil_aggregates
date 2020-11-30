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
#AOB amoA assembly: truncLen=c(229,229), minLen = 229 
#AOA amoA gap: truncLen=c(200,200), minLen = 200 
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
asvTable=seqtab.nochim
colnames(asvTable)=paste0("ASV", seq(length(getUniques(asvTable))))

#How many chimeric reads are there?
sum(asvTable)/sum(seqtab)

# creating summary table
getN <- function(x) sum(dada2::getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(asvTable))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#create ASV table for phyloseq

samples.out <- rownames(asvTable)
subject <- sapply(strsplit(samples.out, "_AOA"), `[`, 1)
site <- substr(samples.out,1,1)

samdf <- data.frame(Sample=samples.out, Site=site)
write.table(samdf, 'map_AOA.txt', sep='\t')
mapAOA=read.table('map_AOA.txt', sep='\t', header=T, row.names=1)


ps <- phyloseq::phyloseq(otu_table(asvTable, taxa_are_rows=FALSE), 
               sample_data(mapAOA))
			   
#rename ASV from sequences to ASV#
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

pdf('AOA_richness.pdf')
plot_richness(ps, x="Size", measures=c("Observed","Shannon", "Simpson"), color="Site")
dev.off()


#' Joining nucleic acid ASV based on their amino acid sequences (translated using seqkit tool)
#' command for seqkit
#' seqkit translate -T 11 -f 2 ASV.fa --clean > ASV_AA.fa
#' where -T 11 represent frame usage for bacteria
#' -f 2 start from2nd base in the sequence
#' --clean removes the stop codon symbol * into X

AOA_table=phyloseq_to_df(ps, addtax=F)
head(AOA_table)

AA=read.table("ASV_AA.tab", header=F)
names(AA)[1]='OTU'
head(AA)

nucl=read.table("ASV_NA.tab", header=F)

AA_taxa=AA %>% 
  left_join(AOA_table) %>%
  pivot_longer(!c(V2, OTU, ), names_to = "sampleID", values_to = "count") %>%
  group_by(V2, sampleID) %>%
  summarise(count=sum(count)) %>%
  dplyr::filter(!grepl("X", V2)) %>%
  pivot_wider(names_from = sampleID,
              values_from = count) 

AA_taxa$AAasv=paste0("AA.ASV", seq(length(AA_taxa$V2)))
colnames(AA_taxa)=str_remove_all(colnames(AA_taxa),"_trim_filt_F.fastq.gz") 


