
### Description ###

This repository contains a compilation of programs and scripts for do transcriptome-guided annotation and functional classification of long non-coding RNAs in *Arabidopsis thaliana*

### Citation ###

**Jose Antonio Corona-Gomez, Evelia Lorena Coss-Navarrete, Irving Jair Garcia-Lopez , Jaime Alejandro Pérez-Patiño and Selene L. Fernandez-Valverde**
*Transcriptome-guided annotation and functional classification of long non-coding RNAs in Arabidopsis thaliana*
Scientific Reports 12:14063, 2022

[https://doi.org/10.1038/s41598-022-18254-0](https://doi.org/10.1038/s41598-022-18254-0)


#### SOFTWARE REQUIREMENTS ###

Summary of set up:

The programs are powered by R and need next libraries: 
 
>(WGCNA, DESeq2)

External programs:

> [FastQC v0.11.2](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
> [MultiQC v1.0](https://multiqc.info/)
> [Trimmomatic v0.32](http://www.usadellab.org/cms/?page=trimmomatic)
> [Kallisto v0.44.0](https://pachterlab.github.io/kallisto/)
> [STAR v2.7.2.b](https://github.com/alexdobin/STAR)
> [EMBOSS v6.6.0](https://www.ebi.ac.uk/Tools/emboss/)
> [HMMER v3.1b2](http://hmmer.org/)
> [signalP v4.1](https://services.healthtech.dtu.dk/service.php?SignalP-4.1)
> [blastx](https://www.ncbi.nlm.nih.gov/books/NBK569861/)
> [infernal v1.1.2](http://eddylab.org/infernal/)
> [BLASTn](https://www.ncbi.nlm.nih.gov/books/NBK569861/)
> [CPC v0.9.r2](http://cpc.gao-lab.org/)
> [BedTools](https://bedtools.readthedocs.io/en/latest/)
---

---

### Transcriptomes ###

It is necessary to have transcripts of good quality and depth if you want to look for lncRNAs.
And have reads of transcriptome in fastq format example: *Arabidopsis thaliana* 


---
#### Filtering, assembly, and quantification of transcripts across all transcriptomes

We assessed the quality of all transcriptomes using FastQC v0.11.2 and MultiQC. Low-quality reads and adapters were removed using Trimmomatic v0.32 . All quality filter reads were aligned to the A. thaliana TAIR10 genome using STAR v2.7.2.b (--alignMatesGapMax 120000). The resulting alignments were assembled using StringTie v1.3.4 (-f 0.3 -m 50 -a 10 -j 15 -c 2.5), using the Araport11 annotation as a reference. The resulting transcripts were joined using the merge function (-c 2.5 -f 0.3) of the StringTie v1.3.4 program. 

---
#### lncRNA identification: Strict Method - Identification of lncRNAs in *Arabidopsis thaliana*

To identify the lncRNAs, we first generated the amino acid sequence for all transcripts using TransDecoder v5.3.0. We then applied nine sequential filters based on previous studies. We refer to this process as the Strict Method (SM). 

>First, (1) we selected all autosomal transcripts ≥ 200 nt using the infoseq program of EMBOSS v6.6.0. We eliminated sequences whose translated ORF or nucleotide sequence had homology to proteins in the Uniprot database as measured by the (2) blastp (e-value ≤ 1e-6) or (3) blastx (e-value ≤ 1e-6, strand=”plus”) program, respectively. We subsequently removed sequences with (4) identifiable protein domains found in the base of Pfam (v33.0) using the HMMER v3.1b2 program (e-value ≤ 1e-6) or (5) with identifiable signal peptides using signalP v4.1 (D-cutoff: 0.45). For any reminder sequences, (6) we removed those that had an ORF > 100 aa using the program getorf of EMBOSS v6.6.0. We did an additional filtering step of all sequences with homology to non-redundant proteins (nr) annotated in the NCBI database using BLASTx (evalue ≤ 1e-6, strand = "plus"). For each remaining transcript, we identified the best blast hit against the ‘nr’ database with a percentage of identity above 70% (pident ≥ 70.000). For each best hit, we used the blastdbcmd function to obtain the information related to the ID. The transcripts annotated in NCBI as: "hypothetical protein" (in Refseq), "similar to" (NCBI's annotation pipeline), "putative protein", "unknown (unknown protein, unknown, partial, unknown)", "predicted protein” and “unnamed protein product” (https://www.ncbi.nlm.nih.gov/books/NBK3840/) were retained. tRNAs and rRNAs were identified using infernal v1.1.2 and the covariance models in the Rfam database . We additionally compared sequences with tRNAs and rRNAs reported in A. thaliana using BLASTn (evalue ≤ 1e-6, strand = "plus"). All sequences identified as tRNAs or rRNAs were discarded. Finally, we eliminated transcripts with introns > 6,000 bp.


Detailed information and code in file: 

> `01_Strict_Method.md`


---

#### Classification of lncRNAs by categories

LncRNAs are generally classified by their positional relationship to other genes. We used the following non-overlapping categories, based on the GENCODE annotation:

+ a) Intergenic lncRNAs (lincRNAs): lncRNAs found in intergenic regions.
+ b) Natural antisense lncRNA (NAT): lncRNAs that totally or partially overlap an exon of another gene in the complementary chain.
+ c) Sense-exonic lncRNAs: lncRNAs that totally or partially overlap the exon of another gene with the same direction of transcription (transcribed from the same DNA strand). 
+ d) Intronic lncRNAs: lncRNAs found within the intron of another gene without overlapping any of its exons, including those on the same chain or complementary to the superimposed gene.


It is worth mentioning that all the isoforms of the overlapping gene are considered for all these categories. To know with which genes our lncRNAs overlap, we used the annotation of Araport11 and BedTools intersectBed (sense_exonic lncRNAs -wo -f 0.1 -s, NAT -wo -f 0.1 -S, intronic -wo -f 1 and lincRNAs -wo -v). Finally, all final annotations were inspected by visualizing them in the UCSC Genome Browser. 

Detailed information and code in file: 

> `02_Classification_of_lncRNAs_by_categories.md`

---

#### Quantification of lncRNAs by tissue and stage

To identify lncRNAs specific to a tissue or stage of development, we calculated the value of SPM (Specificity Measure) based on the work methodology described in (Julca et al. 2020). The SPM value of a gene is obtained by dividing its mean TPM value in a category with the sum of the TPM values in all categories, resulting in a value from 0 to 1 where genes that are tissue or stage-specific have values close to 1. Only the top 5% values close to 1 in all transcriptomes were considered tissue-specific or developmental stage-specific.

See the example of how it is done in the following R markdown:

> `03_SPM_git.Rmd`

---

#### Generation of coding and non-coding gene co-expression networks

Gene co-expression networks were generated as described in the script: 

> `04_WGCNA_annotation_At.R`

---

### Contact ###

[https://regrnalab.github.io/](https://regrnalab.github.io/)
