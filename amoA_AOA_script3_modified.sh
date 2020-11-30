#!/bin/bash
# 2017-09 Axel Aigle

#
# Script_3_gap
# 
# set "list-of-your-samples-names": if your samples are named, list them as "A B C" etc., if your samples are numbered just type {0..n} instead of "list-of-your-samples-names"

# optimized for our system on HPCC
sample=M_056_a_AOA_S74
startdir=cutadapt_before_DADA2/filtered
enddir=processing_AOA
#for sample in list-of-your-samples-names ; do 
for sample in M_056_a_AOA_S74; do 
  gunzip -c ${startdir}/${sample}_trim_filt_F.fastq.gz > ${enddir}/${sample}_R1.fastq
  gunzip -c ${startdir}/${sample}_trim_filt_R.fastq.gz > ${enddir}/${sample}_R2.fastq
  seqkit fq2fa ${enddir}/${sample}_R1.fastq -o ${enddir}/${sample}_R1_1.fasta
  seqkit fq2fa ${enddir}/${sample}_R2.fastq -o ${enddir}/${sample}_R2_1.fasta
  seqkit seq -r -p ${enddir}/${sample}_R2_1.fasta -o ${enddir}/${sample}_R2_2.fasta #creating reverse complement
  seqkit fx2tab ${enddir}/${sample}_R1_1.fasta > ${enddir}/${sample}_R1_1.tab
  seqkit fx2tab ${enddir}/${sample}_R2_2.fasta > ${enddir}/${sample}_R2_2.tab
  
  awk -F "\t" '{print $1}' ${enddir}/${sample}_R1_1.tab > ${enddir}/${sample}_R1_2.tab #takes the names
  awk 'BEGIN {FS="/t" } NR { $1=NR} { print}' ${enddir}/${sample}_R1_2.tab > ${enddir}/${sample}_R1_3.tab #numbering the lines
  sed 's/^/AOA-'${sample}'-/' ${enddir}/${sample}_R1_3.tab > ${enddir}/${sample}_R1_4.tab #rename sequence names, adding AOA in front and number to the back
  awk -F "\t" '{print $2}' ${enddir}/${sample}_R1_1.tab > ${enddir}/${sample}_R1_5.tab #only forward sequences
  awk -F "\t" '{print $2}' ${enddir}/${sample}_R2_2.tab > ${enddir}/${sample}_R2_3.tab #only reverse complement sequences
  paste ${enddir}/${sample}_R1_5.tab ${enddir}/${sample}_R2_3.tab > ${enddir}/${sample}_1.tab #combine reads 
  sed 's/\t//g' ${enddir}/${sample}_1.tab > ${enddir}/${sample}_2.tab  #remove tab between the reads
  paste ${enddir}/${sample}_R1_4.tab ${enddir}/${sample}_2.tab > ${enddir}/${sample}_reads.tab #merging names and reads
  
  seqkit tab2fx ${enddir}/${sample}_reads.tab > ${enddir}/${sample}_reads.fasta #creating fasta from the previous file
  
  /mnt/home/stopnise/amoA_MiSeq_sequencing/usearch64 -cluster_fast ${enddir}/${sample}_reads.fasta -id 1 -centroids ${enddir}/${sample}_reads_id_1.fasta -sizeout #dereplicate reads (only unique reads)
  
  transeq -sequence ${enddir}/${sample}_reads_id_1.fasta -outseq ${enddir}/${sample}_reads_id_1_1.fasta -frame 3 -table 11 -clean #translate sequences 
  
  grep -v X ${enddir}/${sample}_reads_id_1_1.fasta > ${enddir}/${sample}_reads_id_1_2.fasta # taking out the names AGAIN
  
  /mnt/home/stopnise/amoA_MiSeq_sequencing/usearch64 -fastx_truncate ${enddir}/${sample}_reads_id_1_1.fasta -trunclen 133 -fastaout ${enddir}/${sample}_reads_id_1_3.fasta # truncate AA reads to 133AA
  
  seqkit fx2tab ${enddir}/${sample}_reads_id_1_3.fasta > ${enddir}/${sample}_reads_id_1_4.tab #from .fa to .tab
  awk -F '\t' '{print $1}' ${enddir}/${sample}_reads_id_1_4.tab > ${enddir}/${sample}_reads_id_1_5.tab # taking out the names AGAIN :o
  
  sed 's/_3//' ${enddir}/${sample}_reads_id_1_5.tab > ${enddir}/${sample}_reads_id_1_6.tab #remove _3 from the names
  
  seqkit fx2tab ${enddir}/${sample}_reads_id_1.fasta > ${enddir}/${sample}_reads_id_1.tab # .fa to .tab of the 
  
  awk -F '\t' '{print $2, $1}' OFS='\t' ${enddir}/${sample}_reads_id_1.tab > ${enddir}/${sample}_reads_id_1_2.tab # reverse order in tab(first sequence than seq ID)
  
  awk 'NR==FNR {h[$2] = $1; next} {print $1, h[$1]}' OFS='\t' ${enddir}/${sample}_reads_id_1_2.tab ${enddir}/${sample}_reads_id_1_6.tab > ${enddir}/${sample}_trans_reads_id_1.tab # comparing the files (NA seq and AA seqIDs) and output the consensus reads.
  
  seqkit tab2fx ${enddir}/${sample}_trans_reads_id_1.tab > ${enddir}/${sample}_trans_reads_id_1.fasta #convert tab to fa
  
  /mnt/home/stopnise/amoA_MiSeq_sequencing/usearch64 -unoise3 ${enddir}/${sample}_trans_reads_id_1.fasta -zotus ${enddir}/${sample}_zotus.fasta -tabbedout ${enddir}/${sample}_unoise3.txt -minsize 2 # instead of renaming seq IDs into Zotus it can be kept with the original seqID
  
  seqkit fx2tab head ${enddir}/${sample}_zotus.fasta > ${enddir}/${sample}_zotus.tab #convert .fa to .tab
  
  sed 's/Zotu/amp/' ${enddir}/${sample}_zotus.tab > ${enddir}/${sample}_zotus_1.tab # rename Zotu# to amp#
  
  awk -F '\t' '{print $1, $3}' OFS='\t' ${enddir}/${sample}_unoise3.txt > ${enddir}/${sample}_unoise3_1.tab # remove 2nd column (denoise)
  awk 'NR==FNR {h[$2] = $1; next} {print $1, h[$1]}' OFS='\t' ${enddir}/${sample}_unoise3_1.tab ${enddir}/${sample}_zotus_1.tab > ${enddir}/${sample}_list_non_singletons.tab # compare unoise3 mapping names with the zotu headers (what for don't know?!)
  
  awk -F '\t' '{print $2}' ${enddir}/${sample}_list_non_singletons.tab > ${enddir}/${sample}_list_non_singletons_1.tab #remove amp column == same output as if this would be done from the zotu.fa file
  
  seqkit fx2tab ${enddir}/${sample}_trans_reads_id_1.fasta > ${enddir}/${sample}_valids_reads_id_1.tab
  
  awk -F '\t' '{print $2, $1}' OFS='\t' ${enddir}/${sample}_valids_reads_id_1.tab > ${enddir}/${sample}_valids_reads_id_1_1.tab
  awk 'NR==FNR {h[$2] = $1; next} {print $1, h[$1]}' OFS='\t' ${enddir}/${sample}_valids_reads_id_1_1.tab ${enddir}/${sample}_list_non_singletons_1.tab > ${enddir}/${sample}_Reads_non_singletons.tab # joining information of sequence and ID from 
  
  seqkit tab2fx ${enddir}/${sample}_Reads_non_singletons.tab > ${enddir}/${sample}_Reads_non_singletons_gap.fasta
  
  rm ${enddir}/${sample}_R1_*.*
  rm ${enddir}/${sample}_R2_*.*
  rm ${enddir}/${sample}_reads.fasta  
  mkdir ${enddir}/Intermediate_reads
  mkdir ${enddir}/DADA2_filtered
  mv ${enddir}/${sample}_*.fastq.gz DADA2_filtered
  rm ${enddir}/${sample}_R*.fastq
  rm ${enddir}/${sample}_reads_id_1_*.fasta
  rm ${enddir}/${sample}_zotus.fasta
  mv ${enddir}/${sample}_reads_id_1.fasta Intermediate_reads
  mv ${enddir}/${sample}_trans_reads_id_1.fasta Intermediate_reads
  rm *.tab
  rm *.txt
done
