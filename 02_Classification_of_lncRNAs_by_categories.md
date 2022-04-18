# Classification of lncRNAs by categories

List of used programs:

* [kentUtils/v302](https://github.com/ENCODE-DCC/kentUtils) - gtfToGenePred and genePredToBed functions
* [get_intron_gtf_from_bed.pl](https://github.com/riverlee/IntronGTF) - the whole code
* [BEDtools](https://bedtools.readthedocs.io/en/latest/index.html) - intersectBed function

List of input data:

* Araport11_gene_type-ed.txt                  : Gene list with gene model type based on the Araport11 genome release. (from https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FAraport11_genome_release)
* Araport11_GFF3_genes_transposons.201606.gtf : Araport11 annotation in GTF format with Col-0 Mitochondrial annotations. (from https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FAraport11_genome_release)
* At_team_lncRNAs_clean_v2.bed                : File generated using the Strict Method (SM) for the identification of lncRNAs in *Arabidopsis thaliana* in BED format.


## 1) Obtain the positions of exons and introns from protein-coding genes

```
# Obtain protein exons 
grep "protein_coding" Araport11_gene_type-ed.txt > Araport11_mRNAs_gene_IDs
grep "AT[1-5]" Araport11_mRNAs_gene_IDs | sort -V | uniq > Araport_mRNAs_gene_IDs_nucleo
awk '{print $1}' Araport_mRNAs_gene_IDs_nucleo > Araport_mRNAs_gene_nucleo_IDs
awk '{print $4}' sort-Coding_Araport11_josian | sed 's/.$//' | sed 's/.$//' | sort -V | uniq | wc -l

cat  Araport11_GFF3_genes_transposons.201606.gtf | grep '[[:space:]]exon[[:space:]]' | awk '{print $1,$4,$5,$10,$6,$7}'| sed 's/;//' | perl -pi -e 's/[""]//g' | sort -k 1,1 -k2,2n > At_exon_pcoding_sorted.bed

grep -Ff Araport_mRNAs_gene_nucleo_IDs At_exon_pcoding_sorted.bed > At_coding_transcript_exon_IDs_Araport11.bed

# Obtain protein introns
grep "exon" Araport11_GFF3_genes_transposons.201606.gtf > Araport11_exones.gtf
for file in Araport11_exones.gtf; do clean=$(echo $file | sed 's/\.gtf//'); gtfToGenePred $file ${clean}.genePred; done
for file in Araport11_exones.genePred;  do clean=$(echo $file | sed 's/\.genePred//'); genePredToBed $file ${clean}.bed12; done
perl get_intron_gtf_from_bed.pl Araport11_exones.bed12 Araport11_intrones.gtf
cat  Araport11_intrones.gtf | grep '[[:space:]]exon[[:space:]]' | awk '{print $1,$4,$5,$10,$6,$7}'| sed 's/;//' | perl -pi -e 's/[""]//g' > At_intron_pcoding.bed 

# Delete hypothetical proteins
grep -vf At_hipotetical_protein_BlastxNR_transcript_IDs At_coding_transcript_exon_IDs_Araport11.bed > At_coding_transcript_exon_IDs_Araport11_sinhp.bed
grep -vf At_hipotetical_protein_BlastxNR_transcript_IDs At_coding_transcript_intron_IDs_Araport11.bed > At_coding_transcript_intron_IDs_Araport11_sinhp.bed

# Change spaces to tabs 
perl -p -i -e 's/ /\t/g' At_coding_transcript_exon_IDs_Araport11_sinhp.bed
perl -p -i -e 's/ /\t/g' At_coding_transcript_intron_IDs_Araport11_sinhp.bed

```

## 2) Classification of lncRNAs

```
# > sense-exonic
intersectBed -wo -f 0.1 -s -a At_team_lncRNAs_clean_v2.bed -b At_coding_transcript_exon_IDs_Araport11_sinhp.bed > At_team_lncRNAs_sense.bed
awk '{print $1,$2,$3,$4,$6,$13,$14,$15,$16,$18,$19}' At_team_lncRNAs_sense.bed > At_team_lncRNAs_sense_summary.bed
awk '{print $4}' At_team_lncRNAs_sense_summary.bed | sort -V | uniq > At_team_lncRNAs_sentido_IDs

# > NAT
intersectBed -wo -f 0.1 -S -a At_team_lncRNAs_clean_v2.bed -b At_coding_transcript_exon_IDs_Araport11_sinhp.bed > At_team_lncRNAs_NAT.bed
awk '{print $1,$2,$3,$4,$6,$13,$14,$15,$16,$18,$19}' At_team_lncRNAs_NAT.bed > At_team_lncRNAs_NAT_summary.bed
awk '{print $4}' At_team_lncRNAs_NAT_summary.bed | sort -V | uniq > At_team_lncRNAs_antisentido_IDs

# > intronic
intersectBed -wo -f 1 -a At_team_lncRNAs_clean_v2.bed  -b At_coding_transcript_intron_IDs_Araport11_sinhp.bed > At_team_lncRNAs_intronic.bed
awk '{print $1,$2,$3,$4,$6,$13,$14,$15,$16,$18,$19}' At_team_lncRNAs_intronic.bed > At_team_lncRNAs_intronic_summary.bed
awk '{print $4}' At_team_lncRNAs_intronic_summary.bed | sort -V | uniq > At_team_lncRNAs_intronicos_IDs

# > intronic nonsense in the strand
grep -w "[.]" At_team_lncRNAs_intronicos_summary.bed > At_team_lncRNAs_intronicos_SIN_sentido_summary.bed
awk '{print $4}' At_team_lncRNAs_intronicos_SIN_sentido_summary.bed | sort -V | uniq > At_team_lncRNAs_intronicos_SIN_sentido_IDs

# > intergenic
intersectBed -wo -v -a At_team_lncRNAs_clean_v2.bed -b At_coding_transcript_exon_IDs_Araport11_sinhp.bed At_coding_transcript_intron_IDs_Araport11_sinhp.bed > At_team_lncRNAs_intergenicos.bed
awk '{print $1,$2,$3,$4,$6}' At_team_lncRNAs_intergenicos.bed > At_team_lncRNAs_intergenicos_summary.bed
awk '{print $4}' At_team_lncRNAs_intergenicos_summary.bed | sort -V | uniq > At_team_lncRNAs_intergenicos_IDs

```

### Analyze the categorized lncRNAs

```
#--- Load data in R------
library(VennDiagram)

# > intergenic
At_team_lncRNAs_intergenicos_IDs          <- read.table("At_team_lncRNAs_intergenicos_IDs", header= F)

# > intronic
At_team_lncRNAs_intronicos_IDs             <- read.table("At_team_lncRNAs_intronicos_IDs", header= F)
At_team_lncRNAs_intronicos_SIN_sentido_IDs <- read.table("At_team_lncRNAs_intronicos_SIN_sentido_IDs", header= F)

# > sense-exonic
At_team_lncRNAs_sentido_IDs                <- read.table("At_team_lncRNAs_sentido_IDs", header= F)

# > NAT
At_team_lncRNAs_antisentido_IDs            <- read.table("At_team_lncRNAs_antisentido_IDs", header= F)

ORF100_At_filter8_clean_IDs                <- read.table("./cluster_bed_genes/At_team_lncRNAs_clean_v2_IDs", header= F)

At_team_lncRNAs_clasificados  <- c(At_team_lncRNAs_antisentido_IDs$V1,At_team_lncRNAs_sentido_IDs$V1,At_team_lncRNAs_intronicos_IDs$V1, At_team_lncRNAs_intergenicos_IDs$V1)

venn.diagram(x = list(At_team_lncRNAs_clasificados, ORF100_At_filter8_clean_IDs$V1),
  "At_team_lncRNA_clasificados_noclasificados.png" , 
  fill = c("blue", 'red'),
  category.names = c("Clasificados", "TODOS"))

# Obtain the IDs that are not classified
At_team_noanotados_overlap <- calculate.overlap(x= list("TODOS" = ORF100_At_filter8_clean_IDs$V1, "Clasificados" = At_team_lncRNAs_clasificados))
At_team_NOanotados         <- subset(ORF100_At_filter8_clean_IDs,!(At_team_noanotados_overlap$a1 %in% At_team_noanotados_overlap$a3)) 
At_team_compartidos        <- At_team_noanotados_overlap$a3

write.table(At_team_NOanotados$V1, "./At_team_NOanotados_noclasificados.txt", row.names = F, col.names = F)
#------

# Create a new file without characters
# Note :  At_team_NOanotados_IDs contain the same information of At_team_NOanotados_noclasificados.txt.
perl -pi -e 's/[""]//g' At_team_NOanotados_IDs 
head -n-1 At_team_NOanotados_IDs > At_team_NO_anotados_IDs
grep -Ff At_team_NO_anotados_IDs At_team_lncRNAs_clean_v2.bed > At_team_lncRNAs_NO_anotados.bed

intersectBed -wo -S -a At_team_lncRNAs_NO_anotados.bed -b At_coding_transcript_intron_IDs_Araport11_sinhp.bed > At_team_lncRNA_NO_anotados_intronico.bed
intersectBed -wo -S -a At_team_lncRNAs_NO_anotados.bed -b At_coding_transcript_exon_IDs_Araport11_sinhp.bed  > At_team_lncRNA_NO_anotados_exones.bed

awk '{print $4}' At_team_lncRNA_NO_anotados_exones.bed | sort -V | uniq > At_team_lncRNAs_NO_anotados_ex_IDs
awk '{print $4}' At_team_lncRNA_NO_anotados_intronico.bed | sort -V | uniq > At_team_lncRNAs_NO_anotados_int_IDs

cat At_team_lncRNAs_NO_anotados_ex_IDs At_team_lncRNAs_NO_anotados_int_IDs | sort -V| uniq > At_team_lncRNAs_NO_anotados_ex_inIDS
diff At_team_lncRNAs_NO_anotados_ex_inIDS At_team_NO_anotados_IDs --suppress-common-lines | grep \> | cut -f2 -d' ' > At_team_lncRNAs_44_NOanno_IDs

# > NAT
grep -v "AT5G00460.1" At_team_lncRNAs_NO_anotados_ex_inIDS > At_team_lncRNAs_NO_anotados_ex_inIDS_ed # Correct annotation of intergenic to NAT
cat At_team_lncRNAs_antisentido_IDs At_team_lncRNAs_NO_anotados_ex_inIDS_ed | sort -V | uniq > At_team_NAT_IDs

# > intronic
cat At_team_lncRNAs_intronicos_antisentido_IDs At_team_lncRNAs_intronicos_sentido_IDs | sort -V | uniq > At_team_lncRNAs_intronic_IDs

# > sense-exonic
cat At_team_lncRNAs_sentido_IDs At_team_lncRNAs_44_NOanno_IDs | sort -V | uniq > At_team_sense_exonic_IDs

#--- Analyze the distribution of classifications (in R)----
At_team_sense_exonic_IDs <- read.table("At_team_sense_exonic_IDs", header= F)
At_team_NAT_IDs <- read.table("At_team_NAT_IDs", header= F)
At_team_lncRNAs_intergenicos_IDs <- read.table("At_team_lncRNAs_intergenicos_IDs", header= F)
At_team_lncRNAs_intronic_IDs  <- read.table("At_team_lncRNAs_intronic_IDs", header= F)

# Graph
venn.diagram(x = list(At_team_lncRNAs_intronic_IDs$V1, 
  At_team_sense_exonic_IDs$V1, At_team_NAT_IDs$V1, At_team_lncRNAs_intergenicos_IDs$V1),
  "At_team_Clasificacion_comparacionIDS_Araport11.png" , 
  fill = c("blue", 'red', 'green', 'orange'),
  category.names = c("Intronic", "Sense-exonic", "NAT", "Int"))

At_team_clasificados_overlap <- calculate.overlap(x= list("Intronic" = At_team_lncRNAs_intronic_IDs$V1, "sense_exonic" = At_team_sense_exonic_IDs$V1, "NAT" = At_team_NAT_IDs$V1, "Intergenic" = At_team_lncRNAs_intergenicos_IDs$V1))

At_team_13compartidos_NAT_Intronic <- At_team_clasificados_overlap$a4
At_team_4compartidos_NAT_senseexonic <- At_team_clasificados_overlap$a13
write.table(At_team_13compartidos_NAT_Intronic, "./Comparaciones_database_a13.txt", row.names = F, col.names = F)
#------

# Create a new file without characters
# Note :  At_team_20_repetidos_IDs contain the same information of Comparaciones_database_a13.txt.
perl -pi -e 's/[""]//g' At_team_20_repetidos_IDs
head -n-1 At_team_20_repetidos_IDs > At_team_20_repetidos_IDs_ed 

# Remove shared IDs between NAT and intronic from the intronic category
grep -vf At_team_20_repetidos_IDs_ed At_team_lncRNAs_intronic_IDs > At_team_lncRNAs_intronic_depurado_IDs
grep -v "AT3G09922.1" At_team_sense_exonic_IDs > At_team_sense_exonic_IDs_ed
cat gen_AT5G00460 gen_AT3G09922 At_team_lncRNAs_intergenicos_IDs > At_team_lncRNAs_intergenicos_IDs_ed

At_team_lncRNAs_intronic_depurado_IDs  <- read.table("At_team_lncRNAs_intronic_depurado_IDs", header= F)

#--- Analyze the distribution of classifications (in R)----
At_team_lncRNAs_intronic_depurado_IDs  <- read.table("At_team_lncRNAs_intronic_depurado_IDs", header= F)
venn.diagram(x = list(At_team_lncRNAs_intronic_depurado_IDs$V1, 
  At_team_sense_exonic_IDs$V1, At_team_NAT_IDs$V1, At_team_lncRNAs_intergenicos_IDs$V1),
  "At_team_Clasificacion_comparacionIDS_Araport11_depurados.png" , 
  fill = c("blue", 'red', 'green', 'orange'),
  category.names = c("Intronic", "Sense-exonic", "NAT", "Int"))
#------

```

### Classify by genes

```
grep -wf At_team_lncRNAs_clean_v2_IDs At_team_lncRNas_class_table.txt > At_team_lncRNAs_class_table_complete.txt

# > intronic
grep -Ff At_team_lncRNAs_intronic_depurado_IDs ./cluster_bed_genes/At_team_lncRNAs_class_table_complete.txt > ./cluster_bed_genes/At_team_intronic_stranded_lncRNAs_class_table.txt
grep -Ff At_team_lncRNAs_intronicos_SIN_sentido_IDs ./cluster_bed_genes/At_team_lncRNAs_class_table_complete.txt > ./cluster_bed_genes/At_team_intronic_wstranded_lncRNAs_class_table.txt
cat ./cluster_bed_genes/At_team_intronic_stranded_lncRNAs_class_table.txt ./cluster_bed_genes/At_team_intronic_wstranded_lncRNAs_class_table.txt > At_team_intronic_lncRNAs_table.txt
awk '{print $2}' At_team_intronic_lncRNAs_table.txt | sort -V | uniq > At_team_lncRNAs_intronic_depurado_GENES_IDs

# > sense-exonic
grep -Ff At_team_sense_exonic_IDs_ed ./cluster_bed_genes/At_team_lncRNAs_class_table_complete.txt > At_team_sense_exonic_lncRNAs_table.txt
awk '{print $2}' At_team_sense_exonic_lncRNAs_table.txt | sort -V | uniq > At_team_sense_exonic_GENES_IDs

# > NAT
grep -Ff At_team_NAT_IDs ./cluster_bed_genes/At_team_lncRNAs_class_table_complete.txt > At_team_NAT_lncRNAs_table.txt
awk '{print $2}' At_team_NAT_lncRNAs_table.txt | sort -V | uniq > At_team_NAT_GENES_IDs

# > intergenic
grep -Ff At_team_lncRNAs_intergenicos_IDs_ed ./cluster_bed_genes/At_team_lncRNAs_class_table_complete.txt > At_team_intergenic_lncRNAs_table.txt
awk '{print $2}' At_team_intergenic_lncRNAs_table.txt | sort -V | uniq > At_team_lncRNAs_intergenicos_GENES_IDs


```

### Obtain BED file

```
# > intronic
grep -Ff At_team_lncRNAs_intronic_depurado_IDs At_team_lncRNAs.bed > At_team_intronic_lncRNAs.bed
grep -Ff At_team_lncRNAs_intronicos_SIN_sentido_IDs At_team_lncRNAs.bed > At_team_wsense_intronic_lncRNAs.bed
awk '$(NF+1) = "intronic_lncRNA"' At_team_intronic_lncRNAs.bed > At_team_intronic_lncRNAs_ed.bed
awk '$(NF+1) = "intronic_lncRNA"' At_team_wsense_intronic_lncRNAs.bed > At_team_wsense_intronic_lncRNAs_ed2.bed

# > sense-exonic
grep -Ff At_team_sense_exonic_IDs_ed At_team_lncRNAs.bed > At_team_sense_exonic_lncRNAs.bed
awk '$(NF+1) = "sense_exonic_lncRNA"' At_team_sense_exonic_lncRNAs.bed > At_team_sense_exonic_lncRNAs_ed.bed

# > NAT
grep -Ff At_team_NAT_IDs At_team_lncRNAs.bed > At_team_NAT_lncRNAs.bed
awk '$(NF+1) = "NAT_lncRNA"' At_team_NAT_lncRNAs.bed > At_team_NAT_lncRNAs_ed.bed


# > intergenic
grep -Ff At_team_lncRNAs_intergenicos_IDs_ed At_team_lncRNAs_clean_v2.bed > At_team_intergenic_lncRNAs.bed
awk '$(NF+1) = "intergenic_lncRNA"' At_team_intergenic_lncRNAs.bed > At_team_intergenic_lncRNAs_ed.bed

cat At_team_NAT_lncRNAs_ed.bed At_team_intronic_lncRNAs_ed.bed At_team_wsense_intronic_lncRNAs_ed2.bed At_team_sense_exonic_lncRNAs_ed.bed At_team_intergenic_lncRNAs_ed.bed | sort -k1,1 -k2,2n > At_team_lncRNA_clasification.bed

```
