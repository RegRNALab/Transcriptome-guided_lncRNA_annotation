# Strict Method - Identification of lncRNAs in *Arabidopsis thaliana*

List of used programs:

* [gffread](https://github.com/gpertea/gffread)
* [EMBOSS/v6.6.0](http://emboss.sourceforge.net/apps/cvs/emboss/apps/infoseq.html)    - infoseq function
* [kentUtils/v302](https://github.com/ENCODE-DCC/kentUtils)                           - faSomeRecords, gtfToGenePred and genePredToBed functions
* [TransDecoder/v5.3.0](https://github.com/TransDecoder/TransDecoder/releases)        - TransDecoder.LongOrfs and TransDecoder.Predict functions
* [ncbi-blast+/v2.6.0](https://github.com/ncbi/blast_plus_docs)                       - blastp, blastx and blastdbcmd functions
* [hmmer/v3.1b2](https://www.howtoinstall.me/ubuntu/18-04/hmmer/)                     - hmmscan function
* [SignalP/v4.1](https://vcru.wisc.edu/simonlab/bioinformatics/programs/install/signalp.htm) - signalp function
* [infernal/v1.1.2](http://eddylab.org/infernal/)                                     - cmsearch function
* [get_intron_gtf_from_bed.pl](https://github.com/riverlee/IntronGTF)                 - the whole code

List of input data:

* At_stringtie_merged_team.gtf                : Our annotation in GTF format on the Araport11 genome.
* araport_biotypes.tsv                        : Transcript list with gene model type based on the Araport11 genome release.
* Araport11_gene_type-ed.txt                  : Gene list with gene model type based on the Araport11 genome release. (from https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FAraport11_genome_release)
* Araport11_GFF3_genes_transposons.201606.gtf : Araport11 annotation in GTF format with Col-0 Mitochondrial annotations. (from https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FAraport11_genome_release)
* GreeNC_Araport11_lncRNAs.bed                : File generated from GreeNC IDs for *Arabidopsis thaliana* in BED format.


## 1) Transcripts over 200 nt

### Convert GTF to FASTA

```
gffread At_stringtie_merged_team.gtf -g At_allgenomes.fa -w At_stringm_seq.fasta --keep-genes
```

### Select transcripts over 200 nt

```
infoseq At_stringm_seq.fasta -only -name -type -length | sort -k 3n > At_transcripts_length
awk '$3 >= 200 {print}' At_transcripts_length > At_transcripts_more_than_200
awk '{print $1}' At_transcripts_more_than_200 | tail -n +2 | sort -V > At_ID_list_more_than_200
```

### Obtain the FASTA file 

```
faSomeRecords At_stringm_seq.fasta At_ID_list_more_than_200 At_filter_200nt_st.fasta
```

## 2) Eliminate all sequences with homology to proteins

### Predicting the ORF and amino acid sequence

```
TransDecoder.LongOrfs -m 100 -S -t At_filter_200nt_st.fasta
TransDecoder.Predict -t At_filter_200nt_st.fasta
```

### Obtain all translated ORF or nucleotide sequence had homology to proteins in the Uniprot database (evalue <= 1e10-6)

```
blastx -query At_filter_200nt_st.fasta -db uniprot-sprot -num_threads 8 -max_target_seqs 100 -strand "plus" -outfmt "6 qseqid sseqid qlen slen qstart sstart qend send score evalue bitscore length positive qcovs" -evalue 1e-6 > At_blastx.outfmt6

awk 'BEGIN{FS="\t"}{if($10<=(1e-6)){print $AF}}' At_blastx.outfmt6 > At_blastx_evalua10_6 
awk '{print $1}' At_blastx_evalua10_6 | sort -V | uniq > At_blastx_strand_IDs

blastp -query At_filter_200nt_st.fasta.transdecoder.pep -db uniprot-sprot -num_threads 8 -max_target_seqs 100 -outfmt "6 qseqid sseqid qlen slen qstart sstart qend send score evalue bitscore length positive qcovs" -evalue 1e-6 > At_blastp.outfmt6 

awk 'BEGIN{FS="\t"}{if($10<=(1e-6)){print $AF}}' At_blastp.outfmt6 > At_blastp_evalua10_6

# Identify by strand
awk 'BEGIN { OFS = "\t" } { $13 = $5 - $7 } 1' At_blastp_evalua10_6 | awk 'BEGIN { OFS = "\t" } { $14 = $6 - $8 } 1' | cut -f1,2,13,14 > At_blastp.outfmt6_strand 
awk '$3<(-1) && $4<(-1)' At_blastp.outfmt6_strand > At_blastp.outfmt6_strand_negative
awk '{print $1}' At_blastp.outfmt6_strand_negative | sed 's/.$//' | sed 's/.$//' | sed 's/.$//' | sort -V | uniq > At_blastp_IDs

diff At_blastp_IDs At_blastx_strand_IDs --left-column | grep \> | cut -f2 -d' ' > At_blastxyp_IDs
```

### Obtain all sequences related to protein domains

```
hmmscan -o At_pfam.log_new --domtblout At_TrinotatePFAM.out  --incE 0.0001 --incdomE 0.0001 --cpu 8 Pfam-A.hmm At_filter_200nt_st.fasta.transdecoder.pep
grep "(+)," At_pfam.log_new > At_pfam.log_strand_plus
awk '{print $7}' At_pfam.log_strand_plus | sed 's/:.*//' > At_hmmer_IDs_strand_plus

grep -Ff At_hmmer_IDs_strand_plus At_TrinotatePFAM.out | awk '{print $1,$3,$4,$6,$7}' > At_hmmer_strand
awk 'BEGIN{FS="\t"}{if($5<=(1e-6)){print $AF}}' At_hmmer_strand > At_hmmer_strand_evalua10_6
awk '{print $3}' At_hmmer_strand_evalua10_6 | sed 's/.$//' | sed 's/.$//' | sed 's/.$//' | sort -V| uniq > list_pfam_200nt
```

### Obtain all sequences related to signal peptide 

```
signalp -f short -n signalp.out At_filter_200nt_st.fasta.transdecoder.pep
awk '$6 >= 0.45 {print}' signalp.out  | awk '{print $1}' | sed '1d' | sed 's/.$//'| sed 's/.$//' | sed 's/.$//' | sed 's/##sequence-n//'| sort -V | uniq > At_signalp_IDs
```

### Remove all sequences related to proteins

```
cat At_blastxyp_IDs list_pfam_200nt At_signalp_IDs > IDs_filtros
sort -V IDs_filtros | uniq > IDs_uniq_filtros

diff IDs_uniq_filtros At_ID_list_more_than_200 --suppress-common-lines| grep \> | cut -f2 -d' ' > At_filter5_IDs
```

### Obtain FASTA file 

``` 
faSomeRecords At_filter_200nt_st.fasta At_filter5_IDs At_filter5_200nt.fasta
```

## 3) Remove sequence with an ORF greater than 100 nt

```
getorf -minsize 300 -find 0 At_filter5_200nt.fasta At_ORFs_300_stop-stop.fasta
infoseq At_ORFs_300_stop-stop.fasta -only -name -type -length | sort -k 3n > At_ORF_300_nt_length
grep -v '[[:space:]]100[[:space:]]' At_ORF_300_nt_length > At_ORF_morethan_300_nt_length
awk '{print $1}' At_ORF_morethan_300_nt_length | tail -n +2 | sort -V > sort_At_ORF_300_nt_length
sed 's/.$//' sort_At_ORF_300_nt_length |sed 's/.$//' | sort -V | uniq > At_ORF100aa_IDs_ed
diff At_ORF100aa_IDs_ed At_filter5_IDs --suppress-common-lines| grep \> | cut -f2 -d' ' > ORF100_At_filter6_IDs
```

### Obtain FASTA file 

```
faSomeRecords At_filter5_200nt.fasta ORF100_At_filter6_IDs ORF100_At_filter6.fasta
```

## 4) Eliminate sequences with homology to non-redundant protein (nr)

```
blastx -query ORF100_At_filter6.fasta -db /data/secuencias/NR/nr -num_threads 8 -max_target_seqs 100 -strand "plus" -outfmt "6 qseqid sseqid qlen slen qstart sstart qend send score evalue bitscore length positive qcovs" -evalue 1e-6 -out At_blastx_filtro7_nr_ORF100aa.outmt6

awk 'BEGIN{FS="\t"}{if($10<=(1e-6)){print $AF}}' At_blastx_filtro7_nr_ORF100aa.outmt6 > At_blastxNR_evalua10_6
awk '{print $1}' At_blastxNR_evalua10_6 | sort -V | uniq > At_blastxNR_100aa_IDs
```

## 5) Hypothetical protein search

```
for i in $(cat At_blastxNR_100aa_IDs); do grep $i At_blastx_filtro7_nr_ORF100aa.outmt6 | sort -nrk9,9 | head -1 | cut -d ' ' -f3 | awk '{print $1,$2,$9}' >> At_BlastxNR_IDS_maxvalue; done

awk '{print $2}' At_BlastxNR_IDS_maxvalue | sed 's/ref|//' |  sed 's/gb|//' |  sed 's/dbj|//' |  sed 's/emb|//' | sed 's/prf|//' | sed 's/pir|//' | sed 's/sp|//' | sed 's/|//g' | sort -V | uniq > At_blastx_input
blastdbcmd -db /data/secuencias/NR/nr -entry_batch At_blastx_input  -target_only -out At_BlastxNR_NCBI_IDs
grep ">" At_BlastxNR_NCBI_IDs > At_BlastxNR_NCBI_info_only
```

### Get all categories related to hypothetical proteins according to NCBI

```
grep "hypothetical protein" At_BlastxNR_NCBI_info_only > At_BlastxNR_NCBI_hypoteticalp_info
grep "unnamed protein product" At_BlastxNR_NCBI_info_only > At_BlastxNR_NCBI_unnamedp_info
grep "unknown" At_BlastxNR_NCBI_info_only > At_BlastxNR_NCBI_unknown_info
grep "putative protein" At_BlastxNR_NCBI_info_only > At_BlastxNR_NCBI_putativep_info
grep "predicted protein" At_BlastxNR_NCBI_info_only > At_BlastxNR_NCBI_predictedp_info
grep "similar to" At_BlastxNR_NCBI_info_only > At_BlastxNR_NCBI_similartop_info
```

Also we try with "uncharacterized protein" , "novel" and "unnamed gene", but we did not find anything.

```
# Obtain IDs
cat At_BlastxNR_NCBI_hypoteticalp_info At_BlastxNR_NCBI_unnamedp_info At_BlastxNR_NCBI_unknown_info At_BlastxNR_NCBI_putativep_info At_BlastxNR_NCBI_predictedp_info At_BlastxNR_NCBI_similartop_info> At_hypothetical_protein_info_blastxNR
awk '{print $1}' At_hypothetical_protein_info_blastxNR | sed 's/>//' | sort -V | uniq > At_hypothetical_protein_info_blastxNR_IDs

# Select by maxvalue
grep -Ff At_hypothetical_protein_info_blastxNR_IDs At_BlastxNR_IDS_maxvalue > At_hypothetical_protein_info_maxvalue
awk '{print $1}' At_hypothetical_protein_info_maxvalue > At_hypothetical_protein_BlastxNR_transcript_IDs
```

### Retain hypothetical proteins related to lncRNAs

```
grep -vf At_hypothetical_protein_BlastxNR_transcript_IDs At_blastxNR_100aa_IDs | grep -v "AT3G09922.1" > At_blastxNR_100aa_IDs_ed
diff At_blastxNR_100aa_IDs_ed ../ORF100_At_filter6_IDs --suppress-common-lines| grep \> | cut -f2 -d' ' > ORF100_At_filter7_IDs
```

Note: At_hypothetical_protein_BlastxNR_transcript_IDs contain 972 IDs related to hypothetical proteins.

### Obtain FASTA file 

```
faSomeRecords ORF100_At_filter6.fasta ORF100_At_filter7_IDs ORF100_At_filter7.fasta
```

## 6) Delete tRNA, rRNA and mitochondria sequences

### Remove tRNA and rRNA sequences

```
# Prepare database
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.5/Rfam.cm.gz
mv Rfam.cm.gz Rfam_14.5.cm.gz
gzip -d Rfam_14.5.cm.gz

cmfetch Rfam_14.5.cm tRNA > tRNA_Rfam.cm
cmfetch Rfam_14.5.cm SSU_rRNA_eukarya >SSU_rRNA_eukarya_Rfam.cm
cmfetch Rfam_14.5.cm LSU_rRNA_eukarya > LSU_rRNA_eukarya_Rfam.cm
cmfetch Rfam_14.5.cm 5S_rRNA > 5S_rRNA_Rfam.cm
cmfetch Rfam_14.5.cm 5_8S_rRNA > 5_8S_rRNA_Rfam.cm

cat *_Rfam.cm > housekeeping_rfam_14.5.cm

# Detection using Infernal
cmsearch --cpu 32 --cut_ga --tblout Inferna_housekeeping_rfam.tsv housekeeping_rfam_14.5.cm  ORF100_At_filter7.fasta
awk '{print $1}' At_Inferna_housekeeping_rfam.tsv | sed 1,2d | head -n-10 > At_tRNA_ORF100aa_IDs

# Detection with BLAST
grep -w "rRNA" araport_biotypes.tsv | awk '{print $1}' > Araport_rRNA_IDs
grep -w "tRNA" araport_biotypes.tsv| awk '{print $1}' > Araport_tRNA_IDs

faSomeRecords At_stringm_seq.fasta Araport_rRNA_IDs At_rRNA.fasta
faSomeRecords At_stringm_seq.fasta Araport_tRNA_IDs At_tRNA.fasta

makeblastdb -in At_rRNA.fasta -dbtype nucl
blastn -query ORF100_At_filter7.fasta -db At_rRNA.fasta -num_threads 8 -max_target_seqs 1 -outfmt "6 qseqid sseqid qlen slen qstart sstart qend send pident score evalue length positive qcovs" -out At_team_rRNA_ORF100.outfmt6

blastn -query ORF100_At_filter7.fasta -db At_tRNA.fasta -num_threads 8 -max_target_seqs 1 -outfmt "6 qseqid sseqid qlen slen qstart sstart qend send pident score evalue length positive qcovs" -out At_team_tRNA_ORF100.outfmt6

# Obtain IDs
awk '{print $1}' At_team_rRNA_ORF100.outfmt6 > At_rRNA_ORF100aa_IDs_blastn
awk '{print $1}' At_team_tRNA_ORF100.outfmt6 > At_tRNA_ORF100aa_IDs_blastn
cat At_rRNA_ORF100aa_IDs_blastn At_tRNA_ORF100aa_IDs_blastn | sort -V | uniq > At_tRNA_rRNA_ORF100aa_IDs_blastn

grep -vf At_tRNA_ORF100aa_IDs BlastXNR_At_ORF100a_new/ORF100_At_filter7_IDs | grep -vf At_tRNA_rRNA_ORF100aa_IDs_blastn > ORF100_At_filter8_IDs 

```

### Remove sequences from mitochondria

```
grep -v "ATM" ORF100_At_filter8_IDs > ORF100_At_filter8_clean_IDs

```

### Obtain FASTA file 

```
faSomeRecords ORF100_At_filter7.fasta ORF100_At_filter8_clean_IDs ORF100_At_filter8.fasta
```

### Obtain BED and GTF files

```
# Convert to BED file
for file in At_stringtie_merged_team.gtf; do clean=$(echo $file | sed 's/\.gtf//'); gtfToGenePred $file ${clean}.genePred; done
for file in At_stringtie_merged_team.genePred;  do clean=$(echo $file | sed 's/\.genePred//'); genePredToBed $file ${clean}.bed; done

# Obtain GTF file
grep -Ff ORF100_At_filter8_clean_IDs ../../At_stringtie_merged_team.gtf | sort -k 1,1 -k2,2n > At_team_lncRNAs.gtf


```

## 7) Delete sequences with introns over 6000 nt

```
# Obtain exons
grep '[[:space:]]exon[[:space:]]' At_team_lncRNAs.gtf > At_exonlncRNAs.gtf

# Convert to BED file
for file in At_exonlncRNAs.gtf; do clean=$(echo $file | sed 's/\.gtf//'); ./gtfToGenePred $file ${clean}.genePred; done
for file in At_exonlncRNAs.genePred;  do clean=$(echo $file | sed 's/\.genePred//'); ./genePredToBed $file ${clean}.bed; done
sort -k 1,1 -k2,2n At_exonlncRNAs.bed > At_exon_lncRNAs_sorted.bed 

# Obtain introns
perl ./get_intron_gtf_from_bed.pl At_exon_lncRNAs_sorted.bed At_lncRNAs_intrones.gtf
awk '{print $1,$4,$5,$10,$6,$7}' At_lncRNAs_intrones.gtf | sed 's/;//' | sed 's/transcript://' | perl -pi -e 's/[""]//g' > At_introns_lncRNAs.bed

# Obtain length of introns
awk 'BEGIN { OFS = "\t" } { $7 = $3 - $2 } 1' At_introns_lncRNAs.bed | cut -f4,7 > At_introns_lncRNAs_length.txt
sort -k2 -n At_introns_lncRNAs_length.txt > At_introns_lncRNAs_length_sorted.txt

# Obtain IDs and delete from the file
awk 'BEGIN {OFS=FS="\t"} $2 > 6000' At_introns_lncRNAs_length_sorted.txt | awk '{print $1}' | grep "MSTRG" | sort -V | uniq  > At_team_introns_lncRNAs_more_than_6000_IDs
grep -vf At_team_introns_lncRNAs_more_than_6000_IDs At_team_lncRNAs.gtf > At_team_lncRNAs_clean.gtf

# Obtain BED file
for file in At_team_lncRNAs_clean.gtf; do clean=$(echo $file | sed 's/\.gtf//'); ./gtfToGenePred $file ${clean}.genePred; done
for file in At_team_lncRNAs_clean.genePred;  do clean=$(echo $file | sed 's/\.genePred//'); ./genePredToBed $file ${clean}.bed; done

# Obtain IDs
awk '{print $4}' At_team_lncRNAs_clean.bed | sort -V | uniq > At_team_lncRNAs_IDs
```

## 8) Manual screening to eliminate annotated sequences and transmembrane proteins

### Analyze the isoforms of the genes

```
clusterBed -s -i ../At_team_lncRNAs_clean.bed > cluster_At_team_ORF100_lncRNAs.bed
awk 'BEGIN{FS="\t"; OFS=FS} {print $1,$2,$3,$4,$13}' cluster_At_team_ORF100_lncRNAs.bed > cluster_At_team_lncRNAs_transcript_shorted.bed

# ---- First revision and corrections to the file (In R)------
library("dplyr")
library("tidyverse")

clusters <- read.table("cluster_At_team_lncRNAs_transcript_filtered.bed", colClasses =c("factor", "numeric", "numeric", "factor", "factor"), col.names = c("chrom", "start", "end", "gene", "clust")) 
lncRNAs_new_genes <- read.table("At_lncRNas_faltantes_IDs", header =F, colClasses ="factor")
numb <- nrow(lncRNAs_new_genes)
numb_plus <- sum(numb,5508) -1 # ultimo numero anotado 5507
lncRNAs_new_genes_df<- data.frame(lncRNAs_new_genes, 5508:numb_plus)
colnames(lncRNAs_new_genes_df) <- c("gene", "clust")
lncRNAs_new_genes_df$clust <- as.factor(lncRNAs_new_genes_df$clust)

clusters_genes <-select(clusters, "gene", "clust")
clusters_df = rbind(clusters_genes, lncRNAs_new_genes_df) 
clust_group <- factor(clusters_df$clust)

# check isoforms
lncRNAs_by_isoforms <- table(clust_group)[table(clust_group) >1]
lncRNAs_by_genes <- table(clust_group)[table(clust_group) ==1]
lncRNAs_by_isoforms_clusters <- c(x=names(lncRNAs_by_isoforms))

# Obtain gene isoforms
lncRNAs_by_isoforms_out<-NULL
lncRNAs_by_isoforms_outdb<-NULL

# cluster by groups
for (i in 1:length(lncRNAs_by_isoforms_clusters)) {
  lncRNAs_by_isoforms_out <- clusters_df %>% group_by(clust) %>% filter(clust == lncRNAs_by_isoforms_clusters[i])
  lncRNAs_by_isoforms_outdb <- rbind(lncRNAs_by_isoforms_outdb,lncRNAs_by_isoforms_out)
}

lncRNAs_by_isoforms_db <- as.data.frame(lncRNAs_by_isoforms_outdb)

# change the name in rows
rownames(lncRNAs_by_isoforms_db) <- lncRNAs_by_isoforms_db$gene
head(lncRNAs_by_isoforms_db)

write.table(lncRNAs_by_isoforms_db, "./At_team_cluster_isoforms.txt", sep="\t", row.names = FALSE)
#-----

# Note: Saves the contents of the lncRNAs_by_isoforms_db file to another file without characters, At_team_cluster_isoforms_IDs.
perl -pi -e 's/[""]//g' At_team_cluster_isoforms_IDs
head -n-1 At_team_cluster_isoforms_IDs > At_team_cluster_isoforms_ed_IDs
awk '{print $1}' At_team_cluster_isoforms_ed_IDs | sort -V | uniq > At_team_cluster_isoforms_uniq_IDs
grep -wf At_team_cluster_isoforms_uniq_IDs cluster_At_team_lncRNAs_transcript_filtered.bed > At_team_cluster_isoforms.bed

# Note: Transcriptome coordinates were analyzed to confirm and discard the isoforms of each gene.

# ----Second revision and corrections to the file (In R)------
library("dplyr")
library("tidyverse")

# Load file
clusters_clean_data <- read.table("./At_team_cluster_isoforms_ed.bed", colClasses =c("factor", "numeric", "numeric", "factor", "factor"), col.names = c("chrom", "start", "end", "gene", "clust")) 
clusters_clean_data_select <-select(clusters_clean_data, "gene", "clust")
clust_group_clean <- factor(clusters_clean_data_select$clust)

lncRNAs_by_isoforms_clean <- table(clust_group_clean)[table(clust_group_clean) >1]
lncRNAs_by_genes_clean <- table(clust_group_clean)[table(clust_group_clean) ==1]

lncRNAs_by_isoforms_clusters_clean <- c(x=names(lncRNAs_by_isoforms_clean))

# Obtain gene isoforms
lncRNAs_by_isoforms_out_clean<-NULL
lncRNAs_by_isoforms_outdb_clean<-NULL

for (i in 1:length(lncRNAs_by_isoforms_clusters_clean)) {
  lncRNAs_by_isoforms_out_clean <- clusters_clean_data_select %>% group_by(clust) %>% filter(clust == lncRNAs_by_isoforms_clusters_clean[i])
  lncRNAs_by_isoforms_outdb_clean <- rbind(lncRNAs_by_isoforms_outdb_clean,lncRNAs_by_isoforms_out_clean)
}

lncRNAs_by_isoforms_db_clean <- as.data.frame(lncRNAs_by_isoforms_outdb_clean)
rownames(lncRNAs_by_isoforms_db_clean) <- lncRNAs_by_isoforms_db_clean$gene

write.table(lncRNAs_by_isoforms_db_clean, "./At_lncRNAs_by_isoforms_clusters_clean.txt", row.names = F, col.names = F)
# --------------

# Note: Saves the contents of the At_lncRNAs_by_isoforms_clusters_clean.txt file to another file without characters, At_team_cluster_isoforms_class_ed.

```

### Manual screening From tx2gene

```
awk '{print $1}' tx2gene_lncRNAs_team.txt | grep -v "transcriptID" | grep -v "MSTRG"| sort -V | uniq > At_team_lncRNas_transcriptIDs
grep -f At_team_lncRNas_transcriptIDs araport_biotypes.tsv > At_team_lncRNas_transcriptIDs_Araport_results
grep "mRNA" At_team_lncRNas_transcriptIDs_Araport_results > At_team_lncRNas_transcriptIDs_Araport_results_mRNA

# Detection of mRNA in the data
awk '{print $1}' At_team_lncRNas_transcriptIDs_Araport_results_mRNA > At_team_lncRNas_transcriptIDs_Araport_results_mRNA_IDs

# mRNA related to hypothetical proteins
grep -f At_team_lncRNas_transcriptIDs_Araport_results_mRNA_IDs At_hipotetical_protein_BlastxNR_transcript_IDs > At_team_lncRNas_transcriptIDs_Araport_results_mRNA_HP_IDs

# hypothetical in our data
grep -vf At_team_lncRNas_transcriptIDs_Araport_results_mRNA_HP_IDs At_team_lncRNas_transcriptIDs_Araport_results_mRNA_IDs > At_team_lncRNas_transcriptIDs_Araport_results_mRNA_ANALIZAR_IDs

grep -f At_team_lncRNas_transcriptIDs_Araport_results_mRNA_ANALIZAR_IDs At_team_lncRNAs_clean.bed > At_team_lncRNAs_ANALIZAR.bed

# lncRNAs previously reported in GreeNC

awk '{print $4}' GreeNC_Araport11_lncRNAs.bed | sort -V | uniq > GreeNC_lncRNAs_IDs
sed 's/.$//' At_team_lncRNas_transcriptIDs_Araport_results_mRNA_HP_IDs | sed 's/.$//' | sort -V | uniq > At_team_lncRNas_transcriptIDs_Araport_results_mRNA_HP_genes_IDs
grep -Ff GreeNC_lncRNAs_IDs At_team_lncRNas_transcriptIDs_Araport_results_mRNA_HP_genes_IDs > GreeNC_lncRNAs_HP_IDs
grep -Ff GreeNC_lncRNAs_HP_IDs At_team_lncRNAs_clean.bed > At_team_lncRNAs_ANALIZAR_HP.bed

# Others IDs
grep -vf GreeNC_lncRNAs_HP_IDs At_team_lncRNas_transcriptIDs_Araport_results_mRNA_HP_IDs > At_team_lncRNas_transcriptIDs_Araport_results_mRNA_HP_others_IDs
grep -Ff At_team_lncRNas_transcriptIDs_Araport_results_mRNA_HP_others_IDs At_team_lncRNAs_clean.bed > At_team_lncRNAs_ANALIZAR_HP_others.bed

```

With this manual review we obtained the file *At_team_transcriptIDs_deleted_sorted*, this file contains 97 IDs.

```
grep -vf At_team_transcriptIDs_deleted_sorted At_team_lncRNAs_clean.bed > At_team_lncRNAs_clean_v2.bed
```