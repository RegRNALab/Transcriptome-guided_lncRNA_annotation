
#### Final script lncRNAs Annotation paper 13/04/2022

### WGCNA


#### Install packages 
#BiocManager::install("dynamicTreeCut")
#BiocManager::install("fastcluster")
#BiocManager::install("WGCNA")
#BiocManager::install("foreach")
#BiocManager::install("iterators")
#BiocManager::install("doParallel")
#BiocManager::install("Hmisc")
#BiocManager::install("Formula")
#BiocManager::install("acepack")
#BiocManager::install("base64enc")
#BiocManager::install("latticeExtra")
#BiocManager::install("jpeg")
#BiocManager::install("htmlTable")
#BiocManager::install("checkmate")
#BiocManager::install("htmlwidgets")
#BiocManager::install("htmltools")
#BiocManager::install("knitr")
#BiocManager::install("impute")
#BiocManager::install("preprocessCore")
#BiocManager::install("janitor")
#BiocManager::install("lubridate")

#### Directory
setwd("C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene")
#load ("D:/R-studio/script_limpio_WGCNA/")
### Cargar paquetes necesarios 
library (GenomicRanges)
library (DESeq2)
library (WGCNA)
library (janitor)
library (dplyr)
library (tidyr)

#### Updating table 224 transcriptomes 

### Load count table

countData <- read.csv(file = 'KallistoExpr_gene.csv', header = TRUE)
head (countData) [1:10,1:10]
dim (countData)

### Removing NAs 
names(countData)[1] <- "Geneid"
countData <- countData[!(is.na(countData$Geneid) | countData$Geneid==""), ]
### Rounding
countData[,-1] <-round (countData[,-1],0)
countData [1:10,1:5]

### Columns in rownames 
countData2 <- countData[,-1]
rownames(countData2) <- countData[,1]
countData <- countData2
rm (countData2)
### Mostrar primeras columnas y filas 
countData [1:5,1:5]
dim (countData)
tail (countData) [1:5,1:5]


### Convert SRR in data frame to comparing with trancriptomes table
colnams_df <- data.frame (Run = colnames(countData))
dim (colnams_df)

table_at_wgcna <- read.csv (file = "C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/Transcriptomas_A_thaliana_wgcna_usados.csv")
nrow(table_at_wgcna)
#head (table_at_wgcna)
#table_at_wgcna [grepl("*_seedling", table_at_wgcna$code),]

### Generating file coldata for Deseq2 
### Add factors to experiments   
coldata <- merge(colnams_df, table_at_wgcna ,by="Run")
colnames(coldata)
nrow (coldata)

col_order <- coldata$Run
### Order columns 
countData <- countData[,col_order]
## Add to each column the name of 
### Agregar a cada columna el nombre del corresponding treatment
colnames(countData) <- coldata$code

### Columns in factors
coldata  <- mutate_at(coldata , vars(time_hrs, tissue, treatment), as.factor)

### Object Deseq2 for vst transformation

dds_wgcna <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = coldata,
  design = ~ tissue)
dds_wgcna

### Remove genes with few counts
ds_wgcna <- dds_wgcna[rowSums(counts(dds_wgcna)) > 10, ]

### transforming (Variance Stabilizing transfomrmation)
vsd_wgcna <- varianceStabilizingTransformation(dds_wgcna, blind=FALSE)
#class (vsd_wgcna)
### Manipulate object 
vsd_wgcna <- assay(vsd_wgcna)

### transpose matrix 
datExpr <- t(vsd_wgcna)    
head (datExpr[,1:10])
dim(datExpr)
rm(dds_wgcna)
rm(vsd_wgcna)

### Save matrix as csv 
write.csv(datExpr, file = "C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/datExpr.csv",
          row.names = TRUE)
### Load matrix 
datExpr <- read.csv("C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/datExpr.csv", row.names = 1)
head (datExpr[,1:10])
dim (datExpr)


###Next we cluster the samples (in contrast to clustering genes that will come later) 
#to see if there are any obvious outliers.
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)



sampleTrees <- hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(17,13)
pdf(file = "C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/sampleTrees.pdf",
    width = 20 , height = 10);
par(cex = 1);
par(mar = c(0,4,2,0))
plot(sampleTrees, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2,
     cex.axis = 2, cex.main = 3)
dev.off();

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTrees, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 120, col = "red");
# Determine cluster under the line
clust <- cutreeStatic(sampleTrees, cutHeight = 120, minSize = 10)
table(clust)

######################################## WGCNA_automatic ###########################
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "coexpression_ara.RData");
#The variable lnames contains the names of loaded variables.
lnames


# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 4,  networkType = "signed", corFnc = bicor)
sft
# Plot the results:

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
pdf(file = "C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/scale_and_mean_connectivity.pdf",
    width = 12 , height = 10)
par(mfrow = c(1,2));
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


############################ Order samples by tissue and stage #################################
Trans_tm_tj <- read.csv("C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/Transcriptomas_ordenados_tiempo_tejido.csv", row.names = 1)

#### Removing sample pollen
Trans_tm_tj[grep("^pollen", Trans_tm_tj$tissue),]
Trans_tm_tj <- Trans_tm_tj[!(Trans_tm_tj$tissue=="pollen"),]


datExpr[1:9,1:9]
rownames(datExpr)
#### Rownames as column
datExpr_order <- data.frame(names = rownames(datExpr), datExpr)
datExpr_order[1:9,1:9]
datExpr_order <- datExpr_order[ order(match(datExpr_order$names, Trans_tm_tj$code)), ]
datExpr_order$names <- NULL

#### without sample pollen
write.csv(datExpr_order, file = "C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/datExpr_stage.csv",
          row.names = TRUE)


#### Network construction 
getwd()
datExpr_stage <- read.csv("C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/datExpr_stage.csv", row.names = 1)
datExpr <- datExpr_stage

net_unmer <- blockwiseModules(datExpr, power = 12,
                              TOMType = "signed", minModuleSize = 50,
                              networkType = "signed",
                              maxBlockSize = 8000,
                              mergeCutHeight = FALSE, deepSplit = 2,
                              corType = "bicor",
                              numericLabels = TRUE, pamRespectsDendro = FALSE,
                              saveTOMs = TRUE,
                              saveTOMFileBase = "TOM_wgcna_unmerged", 
                              verbose = 4)


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting

unmergedColors <- labels2colors(net_unmer$colors)
unmergedColors

unique(mergedColors)
# Plot the dendrogram and the module colors underneath
pdf(file = "C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/clusterDendogram.pdf",
    width = 20 , height = 10);
plotDendroAndColors(net_unmer$dendrograms[[1]], unmergedColors[net_unmer$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()



############################################# Gene enrinchment topGO
####### load packages for annotation 
library(igraph)
library(DBI)
library(sqldf)
library(digest)
library(AnnotationDbi)
library(GO.db)
library(topGO)
library(ALL)
library (ath1121501.db)
library(Rcpp)


########################## Load gene universe from datExpr matrix 
universe <- colnames(datExpr)
unique (moduleLabels_unmer)
unique (net_unmer$colors)
unique(bwModuleColors_unmer)
moduleLabels_unmer <- net_unmer$colors
#bwLabels_unmer

mode(universe)
str(universe)
mode (moduleLabels_unmer)
dim (moduleLabels_unmer)
str(moduleLabels_unmer)

moduleLabels_unmer[moduleLabels_unmer = 1]

############# top node 100 

GOresults<-data.frame()
for(module in unique(moduleLabels_unmer))
{
  genes<-universe[net_unmer$colors==module]
  #genes<-universe[net_unmer$colors=="1"]
  
  geneList <- factor(as.integer(universe %in% genes))
  names(geneList) <- universe
  
  #pdf(file=paste("C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/topgo/cME_",module,".pdf", sep=""))
  
  # topGO analysis
  if(exists("enrich")) {remove(enrich)}
  for(on in c("MF","BP","CC"))
  {
    print(on)
    # Make topGO object
    GOdata_wgcna <- new ("topGOdata", ontology = on, allGenes = geneList , nodeSize = 10, annot = annFUN.org, mapping = "org.At.tair.db")
    #GOdata <- new("topGOdata", ontology = on, allGenes = geneList, nodeSize = 5, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    # fisher test
    result <- runTest(GOdata_wgcna, algorithm = "classic", statistic = "fisher")
    results.table <- GenTable(GOdata_wgcna, result, topNodes = 100)
    colnames(results.table) <- c("GO.ID", "Term", "Annotated", "Significant", "Expected", "result1")
    results.table$result1 <- gsub('<','', results.table$result1)
    
    #results.table <- data.frame(sapply(results.table, pattern = "<", replacement = ""))
    results.table$result1 <- as.numeric(results.table$result1)
    # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, FDR <= 5%
    results.table$qval.bh<-p.adjust(results.table[,"result1"],method="BH")
    # label ontology type
    results.table$ontology<-on
    # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, consider FDR <= 5% in future
    keep <- results.table[as.numeric(results.table[,"qval.bh"])<0.01,]
    if(exists("enrich")) enrich<- rbind(enrich, keep)
    if(!exists("enrich")) enrich<- keep
    
    # draw figure for GO terms pval<=0.05 before FDR correction
    #if(is.na(sigNo<-length(keep$ontology))){next}
    #showSigOfNodes(GOdata_wgcna, score(result), firstSigNodes = sigNo, useInfo = "all")
    #mtext(on, line=-1)
  }
  #dev.off()
  if(dim(enrich)[1]>0)
  {
    enrichME<-enrich
    enrichME$MEs_mer=module
    GOresults<-rbind(GOresults,enrichME)   }
}
write.table(GOresults,
            file="C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/topgo/consensus_GOresults_unmer_tnode100.txt", sep="\t", row.names=FALSE)


#################### lncRNAs in each module 

#### New annotation 7270 lncRNAs
lnc_bed_at <- read.table("C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/At_team_lncRNA_clasification_revmanual.bed",
                         header=FALSE, sep = "")

colnames(lnc_bed_at) <- c("chrom", "chrSt", "chrE", "names", "score", "strand", "tickS",
                          "thickE", "itemRgb", "Exons", "blockSizes", "blockStarts",
                          "class")

### Numbers modules
colors_dat <- as.data.frame(unmergedColors)
numbers_dat <- as.data.frame(moduleLabels_unmer)
all_modules <- cbind (numbers_dat, colors_dat)


colnames(all_modules) <- c("module_number", "module_color")

#write.table (write.table(all_modules,
 #                        file="C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/genes_each_module.txt",
  #                       sep="\t", row.names=TRUE))
all_modules <- read.table(file="C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/genes_each_module.txt", header = T)

#### Classify the genes in each module
all_modules$names <- rownames(all_modules)
rownames(all_modules) <- NULL
head (all_modules)
uni_modules <-  unique(all_modules)
### Order modules
uni_modules <- uni_modules[order(uni_modules$module_number),]
head (uni_modules)


module_gene_freq <- as.data.frame(table(all_modules$module_number))
module_gene_freq_color <- as.data.frame(table(all_modules$module_color))

module_gene_freq_color <- module_gene_freq_color[order(module_gene_freq_color$Freq, decreasing = T),]
### Gene frequency per module
table_freq_modules <- cbind (module_gene_freq, module_gene_freq_color)[,c(1,3,4)]
colnames(table_freq_modules) <- c("mod_number", "mod_color", "freq_genes")
head (table_freq_modules)


#### Merge the names of modules with names of lncRNAs classified

lnc_bed_at_class <- lnc_bed_at[,c(4,13)]
head (lnc_bed_at_class)
library (stringr)
#### Remove decimals in IDs of genes
lnc_bed_at_class$names <- str_sub(lnc_bed_at_class$names, end=-3)
lnc_modules <- merge(all_modules, lnc_bed_at_class,
                     by= "names",
                     all = FALSE)

### Delete repeated genes (isoforms)
lnc_modules <- unique (lnc_modules)

### Annotation updated
write.table(lnc_modules , file="C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/lnc_modules_gene.txt",
            sep="\t", col.names = TRUE, row.names=FALSE)

lnc_modules <- read.table("C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/lnc_modules_gene.txt", header = T)

#### Frequency of lncRNAs

module_lnc_freq <- as.data.frame(table(lnc_modules$moduleLabels_unmer))
line_43 <-data.frame(as.factor(43), as.integer(0))
names(line_43)<-c("Var1","Freq")

module_lnc_freq_2 <- rbind(module_lnc_freq, line_43)
module_lnc_freq_2 <- module_lnc_freq_2[c(1:43,46,44,45),]

table_freqlnc_modules <- cbind(table_freq_modules, module_lnc_freq_2)[,c(1:3,5)]
colnames(table_freqlnc_modules) <- c("mod_number", "mod_color", "freq_genes", "freq_lnc")

write.table (write.table(table_freqlnc_modules,
                         file="C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/table_freqlnc_modules.txt",
                         sep="\t", row.names=TRUE))


### Calculate eigengenes for modules named with numbers

MEs_unmer_number <- moduleEigengenes(datExpr,moduleLabels_unmer)$eigengenes;



#write.table(MEs_unmer_number, file="C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/MEs_unmer_number.txt",
 #           sep="\t", col.names = TRUE, row.names= TRUE)

MEs_unmer_number <- read.table("C:/Users/Hells/Desktop/wgcna_anolnc/kal_gene/MEs_unmer_number.txt", header = T)


