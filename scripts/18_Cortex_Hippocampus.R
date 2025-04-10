### Goals
#1. Translating kalisto transcript counts to gene counts (STEP1-3)
#2. Cortex input into DESeq2, normalize counts, diff exp, fsgea (STEP4)
#3. Hippocampus input into DESeq2, normalize counts, diff exp, fsgea (STEP5)
#4. Positive selection overlap (STEP6)


###Step1: First set up all your libraries
library(stringr)
library(readr)
library(assertr)
library(tximport)
library(GenomicFeatures)
library(DESeq2)
library(ggplot2)
library(regionReport)
library(rhdf5)
library(edgeR)
library(tibble)
library( "genefilter" )
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(UpSetR)
library(TCseq)
library(cluster)
library(EnhancedVolcano)
library(fgsea)
library(biomaRt)

### STEP2: Create a mechanism for getting transcript to gene
TOGAsortxdb <- makeTxDbFromGFF("~/Sorex_Genome2/data/0_refs/GCF_027595985.1_mSorAra2.pri_genomic.gtf",format = "auto")
TOGAsor_k <- keys(TOGAsortxdb, keytype="TXNAME")       
TOGAsortx2gene <- AnnotationDbi::select(TOGAsortxdb, TOGAsor_k, "GENEID", "TXNAME")
TOGAsortx2gene <- TOGAsortx2gene[!duplicated(TOGAsortx2gene[,1]),]
TOGAsortx2gene <- na.omit(TOGAsortx2gene)

### STEP3: Import kallisto quantifications and write out
#Hippocampus
TOGA_subset_h <- read.table("~/Sorex_Genome2/data/00_transcriptomics/ids/ids_hippocampus.txt", header = T)
TOGA_subset_files_h <-file.path("~/Sorex_Genome2/analysis/5_expression/Hippocampus/TranscriptAbundances", TOGA_subset_h$Sample_name, "abundance.tsv")
names(TOGA_subset_files_h) <- paste0("sample_", TOGA_subset_h$Sample_name)
#check all files exist
all(file.exists(TOGA_subset_files_h))
TOGA_h.count.tsv <- tximport(TOGA_subset_files_h, type = "kallisto", tx2gene = TOGAsortx2gene, ignoreAfterBar=TRUE)
TOGA_h.tpm.tsv <- tximport(TOGA_subset_files_h, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = TOGAsortx2gene, ignoreAfterBar=TRUE)
#write transript abundances out
write.table(TOGA_h.tpm.tsv$abundance, "~/Sorex_Genome2/analysis/5_expression/Hippocampus/TranscriptAbundances/hippoc.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(TOGA_h.count.tsv$abundance, "~/Sorex_Genome2/analysis/5_expression/Hippocampus/TranscriptAbundances/hippoc.count.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

#Cortex
TOGA_subset_c <- read.table("~/Sorex_Genome2/data/00_transcriptomics/ids/ids_cortex.txt", header = T)
TOGA_subset_files_c <-file.path("~/Sorex_Genome2/analysis/5_expression/Cortex/TranscriptAbundances", TOGA_subset_c$Sample_name, "abundance.tsv")
names(TOGA_subset_files_c) <- paste0("sample_", TOGA_subset_c$Sample_name)
all(file.exists(TOGA_subset_files_c))
TOGA_c.count.tsv <- tximport(TOGA_subset_files_c, type = "kallisto", tx2gene = TOGAsortx2gene, ignoreAfterBar=TRUE)
TOGA_c.tpm.tsv <- tximport(TOGA_subset_files_c, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = TOGAsortx2gene, ignoreAfterBar=TRUE)
write.table(TOGA_c.tpm.tsv$abundance, "~/Sorex_Genome2/analysis/5_expression/Cortex/TranscriptAbundances/cortex.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(TOGA_c.count.tsv$abundance, "~/Sorex_Genome2/analysis/5_expression/Cortex/TranscriptAbundances/cortex.count.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)


## STEP4: Cortex, Normalization, PCA, DESEQ2, fsgea
#Begin by setting up design matrix
colnames(TOGA_c.count.tsv$counts) <- c("Stg1_1","Stg1_2","Stg1_3","Stg1_4","Stg1_5","Stg2_1","Stg2_2","Stg2_3","Stg2_4","Stg3_1","Stg3_2","Stg3_3","Stg3_4","Stg3_5","Stg4_1","Stg4_2","Stg4_3","Stg4_4","Stg4_5","Stg5_1","Stg5_2","Stg5_3","Stg5_4","Stg5_5")
c_stages <- factor(c(TOGA_subset_c$Run))
c_organs <- factor(c(TOGA_subset_c$Condition))
c_full1 <- factor(c(TOGA_subset_c$Sample_name))
c_sex <- factor(c(TOGA_subset_c$Sex))
cor_stages_organ_frame <-cbind(as.data.frame(c_stages),as.data.frame(c_organs),as.data.frame(c_full1),as.data.frame(c_sex))
#Now to normalize our reads
dds_cor_all <- DESeqDataSetFromMatrix(round(TOGA_c.count.tsv$counts), DataFrame(cor_stages_organ_frame), ~c_sex + c_stages)
mcols(dds_cor_all) <- cbind(mcols(dds_cor_all), row.names(TOGA_c.count.tsv$counts))
rownames(dds_cor_all) <- row.names(TOGA_c.count.tsv$counts)
dds_cor_all <- DESeq(dds_cor_all)
#And look at pcas and heatmaps of counts
vst_dds_cor_all <- vst(dds_cor_all)
pcaData_cor_all<- plotPCA(vst_dds_cor_all,intgroup=c("c_stages","c_organs"), ntop=1000, returnData=TRUE)
ggplot(pcaData_cor_all, aes(x = PC1, y = PC2, color = factor(c_stages))) +
  geom_point(size=3)+
  theme_bw()
###cortexHeatmap###
corsampleDists <- dist(t(assay(vst_dds_cor_all)))
corsampleDistMatrix <- as.matrix(corsampleDists)
colnames(corsampleDistMatrix) <- NULL
#make the heatmap
pheatmap(corsampleDistMatrix, clustering_distance_rows=corsampleDists,
         clustering_distance_cols = corsampleDists, color = colorRampPalette(rev(brewer.pal(n = 9, name ="Reds")))(255))
#make a heatmap of genes now
cortopVarGenes <- head( order( rowVars( assay(vst_dds_cor_all) ), decreasing=TRUE ), 15 )
heatmap.2( assay(vst_dds_cor_all)[ cortopVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
)
#DESEQ2
#CortexStage4vs2
cor24res <- results(dds_cor_all, contrast = c("c_stages","Stage4","Stage2"))
cor24resSig <- subset(cor24res,cor24res$padj<.05)
cor24up <- subset(cor24resSig,(cor24resSig$log2FoldChange)>=0)
cor24down <- subset(cor24resSig,(cor24resSig$log2FoldChange)<=0)
write.table(cor24res, "~/Sorex_Genome2/analysis/5_expression/Cortex/DESEQ2/Stage4vsStage2_DESeq.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(cor24resSig, "~/Sorex_Genome2/analysis/5_expression/Cortex/DESEQ2/Stage4vsStage2_Sig.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(cor24up, "~/Sorex_Genome2/analysis/5_expression/Cortex/DESEQ2/Stage4vsStage2Upsig.txt", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(cor24down, "~/Sorex_Genome2/analysis/5_expression/Cortex/DESEQ2/Stage4vsStage2Downsig.txt", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
####CREATE VOLCANO PLOT FOR FIGURE
cor24resX <-  cor24res
for (i in 1:length(cor24resX$padj)) {
  if  (cor24resX$padj[i]<1e-5 & !is.na (cor24resX$padj[i])) {
    cor24resX$padj[i] <- 1e-5
  }
  if (cor24resX$log2FoldChange[i]>4 & !is.na (cor24resX$log2FoldChange[i])) {
    cor24resX$log2FoldChange[i] <- 4
  }
  if (cor24resX$log2FoldChange[i]< -4 & !is.na (cor24resX$log2FoldChange[i])) {
    cor24resX$log2FoldChange[i] <- -4
  }
}

# create custom key-value pairs for different cell-types
# this can be achieved with nested ifelse statements
keyvals.shape <- ifelse(
  abs(cor24resX$log2FoldChange) == 4, 17,
  ifelse(cor24resX$padj==1e-5, 17,
         16))
keyvals.shape[is.na(keyvals.shape)] <- 1
names(keyvals.shape)[keyvals.shape == 16] <- '<log1e-15'
names(keyvals.shape)[keyvals.shape == 17] <- 'Off-Graph'
###
keyvals <- ifelse(
  cor24resX$padj > 0.05, 'grey',
  ifelse(cor24resX$log2FoldChange <= -1.58, 'red',
         ifelse(cor24resX$log2FoldChange >= 1.58, 'blue',
                ifelse(cor24resX$log2FoldChange >= 0, 'blue',
                       'red'))))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'grey'] <- 'NotSig'
names(keyvals)[keyvals == 'red'] <- 'DownRegulated'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated High Effect'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated'
names(keyvals)[keyvals == 'red'] <- 'Downregulated High Effect'
####
corPSG_over <- Reduce(intersect, list(shrewtotalPSG$p_adjusted,rownames(cor24resSig)))
EnhancedVolcano(cor24resX,
                lab = rownames(cor24resX),
                xlim=c(-4 ,4),
                ylim=c(0,5.5),
                x = 'log2FoldChange',
                y = 'padj',
                shapeCustom = keyvals.shape,
                selectLab = corPSG_over,
                pCutoff = .05,
                FCcutoff = 10,
                colCustom = keyvals,
                pointSize = 4.0,
                legendPosition = 'none',
                drawConnectors = TRUE,
                gridlines.major  = FALSE,
                widthConnectors = 0.5)

#fsgea
cor_Stage4vsStage2_DESeq <- read_delim("~/Sorex_Genome2/analysis/5_expression/Cortex/DESEQ2/Stage4vsStage2_DESeq.tsv", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)
COR_hres <- cor_Stage4vsStage2_DESeq %>% 
  dplyr::select(Gene, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(Gene) %>% 
  summarize(stat=mean(stat))
library(tidyverse)
COR_hranks <- deframe(COR_hres)
fgsea_cor_Stage4vsStage2_DESeq <- fgsea(pathways=gmtPathways("~/Sorex_Genome2/data/000_miscFiles/c2.cp.kegg.v2023.1.Hs.symbols.gmt.txt"), COR_hranks) %>% 
  as_tibble() %>% 
  arrange(padj)
fgsea_cor_Stage4vsStage2_DESeqTidy <- fgsea_cor_Stage4vsStage2_DESeq %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table:
fgsea_cor_Stage4vsStage2_DESeqTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()
corexp_Stage4vs2_fsgea <-ggplot(subset(fgsea_cor_Stage4vsStage2_DESeqTidy,padj<0.05 ), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES>0)) +
  coord_flip() +
  labs(x="KEGGPathwayLiver", y="Normalized Enrichment Score") + 
  scale_fill_manual(values =  c("red", "blue"))+
  theme_bw()
corexp_Stage4vs2_fsgea
ggsave("~/Dehnels_Seasonal_RNAseq3/data/Liver/DifferentialExp/Stage3vs1/fsgea_SumvsWinter.png", livexp_Stage3vs1_fsgea,width = 7.1, height = 6, dpi =300,)
fgsea_cor_Stage4vsStage2_DESeqTidy2 <- apply(fgsea_cor_Stage4vsStage2_DESeqTidy,2,as.character)
write.table(fgsea_cor_Stage4vsStage2_DESeqTidy2, file='~/Dehnels_Seasonal_RNAseq3/data/Liver/DifferentialExp/Stage3vs1/fsgea_liverStage3vs1Tidy.tsv', quote=FALSE, sep='\t')


#Hippocampus
## STEP5: Hippocampus, Normalization, PCA, DESEQ2, fsgea
#Begin by setting up design matrix
colnames(TOGA_h.count.tsv$counts) <- c("Stg1_1","Stg1_2","Stg1_3","Stg1_4","Stg1_5","Stg2_1","Stg2_2","Stg2_3","Stg2_4","Stg3_1","Stg3_2","Stg3_3","Stg3_4","Stg3_5","Stg4_1","Stg4_2","Stg4_3","Stg4_4","Stg4_5","Stg5_1","Stg5_2","Stg5_3","Stg5_4","Stg5_5")
h_stages <- factor(c(TOGA_subset_h$Run))
h_organs <- factor(c(TOGA_subset_h$Condition))
h_full1 <- factor(c(TOGA_subset_h$Sample_name))
h_sex <- factor(c(TOGA_subset_h$Sex))
hip_stages_organ_frame <-cbind(as.data.frame(h_stages),as.data.frame(h_organs),as.data.frame(h_full1),as.data.frame(h_sex))
#Now to normalize our reads
dds_hip_all <- DESeqDataSetFromMatrix(round(TOGA_h.count.tsv$counts), DataFrame(hip_stages_organ_frame), ~h_sex + h_stages)
mcols(dds_hip_all) <- cbind(mcols(dds_hip_all), row.names(TOGA_h.count.tsv$counts))
rownames(dds_hip_all) <- row.names(TOGA_h.count.tsv$counts)
dds_hip_all <- DESeq(dds_hip_all)
#And look at pcas and heatmaps of counts
vst_dds_hip_all <- vst(dds_hip_all)
pcaData_hip_all<- plotPCA(vst_dds_hip_all,intgroup=c("h_stages","h_organs"), ntop=1000, returnData=TRUE)
ggplot(pcaData_hip_all, aes(x = PC1, y = PC2, color = factor(h_stages))) +
  geom_point(size=3)+
  theme_bw()
###hippocHeatmap###
hipsampleDists <- dist(t(assay(vst_dds_hip_all)))
hipsampleDistMatrix <- as.matrix(hipsampleDists)
colnames(hipsampleDistMatrix) <- NULL
#make the heatmap
pheatmap(hipsampleDistMatrix, clustering_distance_rows=hipsampleDists,
         clustering_distance_cols = hipsampleDists, color = colorRampPalette(rev(brewer.pal(n = 9, name ="Reds")))(255))
#make a heatmap of genes now
hiptopVarGenes <- head( order( rowVars( assay(vst_dds_hip_all) ), decreasing=TRUE ), 15 )
heatmap.2( assay(vst_dds_hip_all)[ hiptopVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
)

#HippocampusStage4vs2
hip24res <- results(dds_hip_all, contrast = c("h_stages","Stage4","Stage2"))
hip24resSig <- subset(hip24res,hip24res$padj<.05)
hip24up <- subset(hip24resSig,(hip24resSig$log2FoldChange)>=0)
hip24down <- subset(hip24resSig,(hip24resSig$log2FoldChange)<=0)
write.table(hip24res, "~/Sorex_Genome2/analysis/5_expression/Hippocampus/DESEQ2/Stage4vsStage2_DESeq.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(hip24resSig, "~/Sorex_Genome2/analysis/5_expression/Hippocampus/DESEQ2/Stage4vsStage2_Sig.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(hip24up, "~/Sorex_Genome2/analysis/5_expression/Hippocampus/DESEQ2/Stage4vsStage2Upsig.txt", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(hip24down, "~/Sorex_Genome2/analysis/5_expression/Hippocampus/DESEQ2/Stage4vsStage2Downsig.txt", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
#heatmap
####CREATE VOLCANO PLOT FOR FIGURE
hip24resX <-  hip24res
for (i in 1:length(hip24resX$padj)) {
  if  (hip24resX$padj[i]<1e-5 & !is.na (hip24resX$padj[i])) {
    hip24resX$padj[i] <- 1e-5
  }
  if (hip24resX$log2FoldChange[i]>4 & !is.na (hip24resX$log2FoldChange[i])) {
    hip24resX$log2FoldChange[i] <- 4
  }
  if (hip24resX$log2FoldChange[i]< -4 & !is.na (hip24resX$log2FoldChange[i])) {
    hip24resX$log2FoldChange[i] <- -4
  }
}

# create custom key-value pairs for different cell-types
# this can be achieved with nested ifelse statements
keyvals.shape <- ifelse(
  abs(hip24resX$log2FoldChange) == 4, 17,
  ifelse(hip24resX$padj==1e-5, 17,
         16))
keyvals.shape[is.na(keyvals.shape)] <- 1
names(keyvals.shape)[keyvals.shape == 16] <- '<log1e-15'
names(keyvals.shape)[keyvals.shape == 17] <- 'Off-Graph'
###
keyvals <- ifelse(
  hip24resX$padj > 0.05, 'grey',
  ifelse(hip24resX$log2FoldChange <= -1.58, 'red',
         ifelse(hip24resX$log2FoldChange >= 1.58, 'blue',
                ifelse(hip24resX$log2FoldChange >= 0, 'blue',
                       'red'))))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'grey'] <- 'NotSig'
names(keyvals)[keyvals == 'red'] <- 'DownRegulated'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated High Effect'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated'
names(keyvals)[keyvals == 'red'] <- 'Downregulated High Effect'
####
hipPSG_over <- Reduce(intersect, list(shrewtotalPSG$p_adjusted,rownames(hip24resSig)))
EnhancedVolcano(hip24resX,
                lab = rownames(hip24resX),
                xlim=c(-4 ,4),
                ylim=c(0,5.5),
                x = 'log2FoldChange',
                y = 'padj',
                shapeCustom = keyvals.shape,
                selectLab = hipPSG_over,
                pCutoff = .05,
                FCcutoff = 10,
                colCustom = keyvals,
                pointSize = 4.0,
                legendPosition = 'none',
                drawConnectors = TRUE,
                gridlines.major  = FALSE,
                widthConnectors = 0.5)

#fsgea
hip_Stage4vsStage2_DESeq <- read_delim("~/Sorex_Genome2/analysis/5_expression/Hippocampus/DESEQ2/Stage4vsStage2_DESeq.tsv", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)
hip_hres <- hip_Stage4vsStage2_DESeq %>% 
  dplyr::select(Gene, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(Gene) %>% 
  summarize(stat=mean(stat))
library(tidyverse)
hip_hranks <- deframe(hip_hres)
fgsea_hip_Stage4vsStage2_DESeq <- fgsea(pathways=gmtPathways("~/Sorex_Genome2/data/000_miscFiles/c2.cp.kegg.v2023.1.Hs.symbols.gmt.txt"), hip_hranks) %>% 
  as_tibble() %>% 
  arrange(padj)
fgsea_hip_Stage4vsStage2_DESeqTidy <- fgsea_hip_Stage4vsStage2_DESeq %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table:
fgsea_hip_Stage4vsStage2_DESeqTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()
hipexp_Stage4vs2_fsgea <-ggplot(subset(fgsea_hip_Stage4vsStage2_DESeqTidy,padj<0.05 ), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES>0)) +
  coord_flip() +
  labs(x="KEGGPathwayLiver", y="Normalized Enrichment Score") + 
  scale_fill_manual(values =  c("red", "blue"))+
  theme_bw()
hipexp_Stage4vs2_fsgea
#plotfsgea for both
which(fgsea_hip_Stage4vsStage2_DESeqTidy$pathway == "KEGG_NOTCH_SIGNALING_PATHWAY")
fgsea_hip_Stage4vsStage2_DESeqTidy$leadingEdge
plotCounts(dds_hip_all, gene="GLUT4", intgroup="h_stages")
#
fgsea_hip_Stage4vsStage2_DESeqTidy <- fgsea_hip_Stage4vsStage2_DESeqTidy %>%
  mutate(source = "hip")

fgsea_cor_Stage4vsStage2_DESeqTidy <- fgsea_cor_Stage4vsStage2_DESeqTidy %>%
  mutate(source = "cor")

# Combine the two tibbles
combined_tibble <- bind_rows(fgsea_hip_Stage4vsStage2_DESeqTidy, fgsea_cor_Stage4vsStage2_DESeqTidy)
ggplot(subset(combined_tibble, padj < 0.05), aes(reorder(pathway, NES), NES)) +
  geom_point(aes(shape = source, fill = NES > 0), size = 7, color = "black") +
  coord_flip() +
  labs(x = "KEGG Pathway Liver", y = "Normalized Enrichment Score") +
  scale_fill_manual(values = c("red", "blue"), labels = c("Negative NES", "Positive NES")) +
  scale_shape_manual(values = c(21, 22), labels = c("hip", "cor")) +
  theme_bw() +
  theme(legend.position = "right")
ggplot(subset(combined_tibble, padj < 0.05), aes(reorder(pathway, NES), NES)) +
  geom_point(aes(shape = source, fill = NES > 0), size = 10, color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("red", "blue"), labels = c("Negative NES", "Positive NES")) +
  scale_shape_manual(values = c(21, 22), labels = c("hip", "cor")) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  )

#STEP 6: Comparing data sets with PSGs
#import data sets
shrewtotalPSG <- read_table("~/ShrewProjects/Sorex_Genome2/analysis/3_hyphy/shrewtotal.txt")
shrewspecificPSG <-read_table("~/ShrewProjects/Sorex_Genome2/analysis/3_hyphy/shrewspecific.txt")
hypothal_evo <- read_table("~/Sorex_Genome2/data/000_miscFiles/Hypothal_data_2024/shrew_sigShrewBS", 
                           col_names = FALSE)
hyp24up <- read_table("~/Sorex_Genome2/data/000_miscFiles/Hypothal_data_2024/hyp24up")
hyp24down <- read_table("~/Sorex_Genome2/data/000_miscFiles/Hypothal_data_2024/hyp24down")

#Total Shrews
##Comparison 1a: Total Shrew + Hypothalamus Evolution
Reduce(intersect, list(shrewtotalPSG$p_adjusted,hypothal_evo$X3)) #MLST8, VEGFA, SPHK2, GRIP2, FANCI"
#Comparison 1b: Total Shrew + Cycling Hypothalamus
Reduce(intersect, list(shrewtotalPSG$p_adjusted,hyp24down$Gene)) #KCNK12
Reduce(intersect, list(shrewtotalPSG$p_adjusted,hyp24up$Gene)) #"KNDC1" "PARP4"
#Comparison 1c: Total Shrew + Cycling Cortex
Reduce(intersect, list(shrewtotalPSG$p_adjusted,rownames(cor24resSig))) #SOX9, ASPHD1
#Comparison 1d: Total Shrew + Cycling Hippocampus
Reduce(intersect, list(shrewtotalPSG$p_adjusted,rownames(hip24resSig))) #SOX9


#Shrew Specific
##Comparison 1a: Total Shrew + Hypothalamus Evolution
Reduce(intersect, list(shrewspecificPSG$p_adjusted,hypothal_evo$X3)) #MLST8, VEGFA, SPHK2, GRIP2, FANCI"
#Comparison 2b: Total Shrew + Cycling Hypothalamus
Reduce(intersect, list(shrewspecificPSG$p_adjusted,hyp24down$Gene)) #KCNK12
Reduce(intersect, list(shrewspecificPSG$p_adjusted,hyp24up$Gene)) #"KNDC1" "PARP4"
#Comparison 2c: Total Shrew + Cycling Cortex
Reduce(intersect, list(shrewspecificPSG$p_adjusted,rownames(cor24resSig))) #SOX9, ASPHD1
#Comparison 2d: Total Shrew + Cycling Hippocampus
Reduce(intersect, list(shrewspecificPSG$p_adjusted,rownames(hip24resSig))) #SOX9

###Create upset

gene_listPSGupset <- list(PSG = shrewtotalPSG$p_adjusted,SasPSG = shrewspecificPSG$p_adjusted,SeasonalCortex = rownames(cor24resSig),SeasonalHippocampus = rownames(hip24resSig),SeasonalHypothalamus =rownames(hyp24resSig), ComparativeHypothalamus =hypothal_evo$X3)
upset_gene_listPSGupset <- fromList(gene_listPSGupset)
upset(upset_gene_listPSGupset,nsets = 6,order.by = "freq",)
svg(filename = "~/Sorex_Genome2/analysis/3_hyphy/plot_listPSGupset.svg", width = 10, height = 8)


#overlap between diff exp
#all
Reduce(intersect, list(rownames(cor24down),rownames(hip24down),hyp24down$Gene)) #"DDIT4"    "NFKBIA"   "TP53INP1"
Reduce(intersect, list(rownames(cor24up),rownames(hip24up),hyp24up$Gene))

#corhip
Reduce(intersect, list(rownames(cor24resSig),rownames(hip24resSig))) #21 overlap total
Reduce(intersect, list(rownames(cor24down),rownames(hip24down))) #all in same direction
Reduce(intersect, list(rownames(cor24up),rownames(hip24up))) #all in same direction
#corhyp
Reduce(intersect, list(rownames(cor24down),hyp24down$Gene))# "DDIT4"    "LIMS2"    "NFKBIA"   "STARD8"   "TP53INP1" "ZFP36"
Reduce(intersect, list(rownames(cor24up),hyp24up$Gene)) #"ALOX5AP" "CASQ2"   "GAS8"    "GPR68"   "PRDM12"  "ZFHX3"
#hiphyp
Reduce(intersect, list(rownames(hip24down),hyp24down$Gene))# "CLDN5"    "DDIT4"    "NFKBIA"   "TP53INP1"
Reduce(intersect, list(rownames(hip24up),hyp24up$Gene))




#STEP7 Plotting Notch and Sox9
interesting_cor <- DESeq2::counts(dds_cor_all, normalized=TRUE)
interesting_hip <- DESeq2::counts(dds_hip_all, normalized=TRUE)
#NOTCH1
plotCounts(dds_hip_all, gene="NOTCH1", intgroup="h_stages")
plotCounts(dds_cor_all, gene="NOTCH1", intgroup="c_stages")
NOTCH1_cor <-as.data.frame(interesting_cor[rownames(interesting_cor) %in% "NOTCH1", ])
NOTCH1_hip <-as.data.frame(interesting_hip[rownames(interesting_hip) %in% "NOTCH1", ])
NOTCH1_cor$Stage <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
NOTCH1_hip$Stage <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
colnames(NOTCH1_cor) <- c("Expression", "Stage")
colnames(NOTCH1_hip) <- c("Expression", "Stage")
#plot
summary_NOTCH1_cor <- NOTCH1_cor %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_NOTCH1_cor$Organ <- c(rep("Cortex",5))
summary_NOTCH1_hip <- NOTCH1_hip %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_NOTCH1_hip$Organ <- c(rep("Hippocampus",5))
combined_NOTCH1 <- bind_rows(summary_NOTCH1_hip, summary_NOTCH1_cor)
# Step 3: Plot the data
NOTCH1_plot <- ggplot(combined_NOTCH1, aes(x = Stage, y = mean_expression, shape = Organ, group = Organ)) +
  geom_line(size = 2, color = "blue") +
  geom_point(size = 8, fill = "blue", color = "blue") +
  scale_shape_manual(values = c("Cortex" = 16, "Hippocampus" = 15)) + # 21 = circle, 22 = square
  theme_classic() +
  labs(x = "Stage", y = "Mean Expression", shape = "Organ") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
NOTCH1_plot
ggsave("~/Sorex_Genome2/analysis/5_expression/Hippocampus/DESEQ2/NOTCH1_plot.png", NOTCH1_plot,width = 6, height = 8.5, dpi =300,)

#SOX9
plotCounts(dds_hip_all, gene="NOTCH2", intgroup="h_stages")
plotCounts(dds_cor_all, gene="NOTCH2", intgroup="c_stages")
SOX9_cor <-as.data.frame(interesting_cor[rownames(interesting_cor) %in% "SOX9", ])
SOX9_hip <-as.data.frame(interesting_hip[rownames(interesting_hip) %in% "SOX9", ])
SOX9_cor$Stage <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
SOX9_hip$Stage <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
colnames(SOX9_cor) <- c("Expression", "Stage")
colnames(SOX9_hip) <- c("Expression", "Stage")
summary_SOX9_cor <- SOX9_cor %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_SOX9_cor$Organ <- c(rep("Cortex",5))
summary_SOX9_hip <- SOX9_hip %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_SOX9_hip$Organ <- c(rep("Hippocampus",5))
combined_SOX9 <- bind_rows(summary_SOX9_hip, summary_SOX9_cor)
SOX9_plot <- ggplot(combined_SOX9, aes(x = Stage, y = mean_expression, shape = Organ, group = Organ)) +
  geom_line(size = 1.5, color = "blue") +
  geom_point(size = 4, fill = "blue", color = "blue") +
  scale_shape_manual(values = c("Cortex" = 16, "Hippocampus" = 15)) + # 21 = circle, 22 = square
  theme_classic() +
  labs(x = "Stage", y = "Mean Expression", shape = "Organ") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
SOX9_plot
ggsave("~/Sorex_Genome2/analysis/5_expression/Hippocampus/DESEQ2/SOX9_plot.png", SOX9_plot,width = 6, height = 8.5, dpi =300,)
subset(hyp24res, rownames(hyp24res) == "SOX9")

###SOX9 AA visualization
#amino acid boces
SOX9boxes_4R <- read_delim("~/Sorex_Genome2/analysis/5_expression/1_for_visualizations/SOX9/SOX9boxes_4R_toga.txt", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

SOX9snp_data<-read_table("~/Sorex_Genome2/analysis/5_expression/1_for_visualizations/SOX9/SOX9_rfigure_inputTOGA.txt", 
                          col_names = FALSE)

SOX9boxes_4R$size <- ifelse(SOX9boxes_4R$domain == "NO", 15, 20)  # "NO" segments thinner (10), others thicker (20)
####
SOX9boxes_4R
# Plot with mutations as lollipops (corrected y positioning)
ggplot(SOX9boxes_4R, aes(x = start_position, xend = end_position, y = 1, yend = 1, color = domain, size = size)) +
  # Segment bar
  geom_segment() +
  # Mutation lollipops with adjusted position
  geom_segment(
    data = SOX9snp_data,
    aes(x = X2, xend = X2, y = 1 + 10 / 2, yend = 1 + 10 / 2 + 0.5),  # Stick starts at top of the bars
    color = "black",  # Stick color
    inherit.aes = FALSE
  ) +
  geom_point(
    data = SOX9snp_data,
    aes(x = X2, y = 1 + 10 / 2 + 0.5),  # Head just above the stick
    color = "red",  # Head color
    size = 3,
    inherit.aes = FALSE
  ) +
  # Custom colors for segments
  scale_color_manual(values = c("NO" = "grey", "HMG" = "#FFC20A")) +
  scale_size_identity() +  # Use the size column directly
  theme_minimal() +
  labs(
    title = "Segment Bar with Mutations for Sorex araneus",
    x = "Position",
    y = NULL,  # Remove y-axis label
    color = "Segment Label"
  ) +
  theme(
    axis.text.y = element_blank(),   # Remove y-axis text
    axis.ticks.y = element_blank(), # Remove y-axis ticks
    panel.grid = element_blank(),    # Remove grid lines
    legend.position = "none"
  )


###STEP 8 box and tick plots for FANCI, PALB2, FAAP100
####FANCI Boxes
#amino acid boces
FANCIboxes_4R <- read_delim("~/Sorex_Genome2/analysis/5_expression/1_for_visualizations/FANCI/FANCIboxes_4R_.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

FANCIsnp_data<-read_table("~/Sorex_Genome2/analysis/5_expression/1_for_visualizations/FANCI/FANCI_rfigure_input.txt", 
                         col_names = TRUE)

FANCIboxes_4R$size <- ifelse(FANCIboxes_4R$domain %in% c("NO", "I", "D"), 15, 20)####
FANCIboxes_4R$end_position <- FANCIboxes_4R$end_position +1
# Plot with mutations as lollipops (corrected y positioning)
FANCIsnp_data
ggplot(FANCIboxes_4R, aes(x = start_position, xend = end_position, y = 1, yend = 1, color = domain, size = size)) +
  # Segment bar
  geom_segment() +
  # Mutation lollipops with adjusted position
  geom_segment(
    data = FANCIsnp_data,
    aes(x = boxes, xend = boxes, y = 1 + 10 / 2, yend = 1 + 10 / 2 + 0.5),  # Stick starts at top of the bars
    color = "black",  # Stick color
    inherit.aes = FALSE
  ) +
  geom_point(
    data = FANCIsnp_data,
    aes(x = boxes, y = 1 + 10 / 2 + 0.5),  # Head just above the stick
    color = "red",  # Head color
    size = 3,
    inherit.aes = FALSE
  ) +
  # Custom colors for segments
  scale_color_manual(values = c("NO" = "grey", "S" = "#88CCEE", "SC" ="#332288", "HD"="#E69F00", "INSD"="#117733", "I"="#117733", "DELD"="#882255", "D"="#882255")) +
  scale_size_identity() +  # Use the size column directly
  theme_minimal() +
  labs(
    title = "Segment Bar with Mutations for Sorex araneus",
    x = "Position",
    y = NULL,  # Remove y-axis label
    color = "Segment Label"
  ) +
  theme(
    axis.text.y = element_blank(),   # Remove y-axis text
    axis.ticks.y = element_blank(), # Remove y-axis ticks
    panel.grid = element_blank(),    # Remove grid lines
    legend.position = "none"
  )





####FAAP100 Boxes
#amino acid boces
FAAP100boxes_4R <- read_delim("~/Sorex_Genome2/analysis/5_expression/1_for_visualizations/FAAP100/FAAP100boxes_4R_gaps.txt", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)

FAAP100snp_data<-read_table("~/Sorex_Genome2/analysis/5_expression/hyphy/meme2/FAAP100/FAAP100_rfigure_input.txt", 
                          col_names = TRUE)

FAAP100boxes_4R$size <- ifelse(FAAP100boxes_4R$domain %in% c("NO", "I", "D"), 15, 20)####
FAAP100boxes_4R$end_position <- FAAP100boxes_4R$end_position +1
# Plot with mutations as lollipops (corrected y positioning)
FAAP100snp_data
ggplot(FAAP100boxes_4R, aes(x = start_position, xend = end_position, y = 1, yend = 1, color = domain, size = size)) +
  # Segment bar
  geom_segment() +
  # Mutation lollipops with adjusted position
  geom_segment(
    data = FAAP100snp_data,
    aes(x = boxes, xend = boxes, y = 1 + 10 / 2, yend = 1 + 10 / 2 + 0.5),  # Stick starts at top of the bars
    color = "black",  # Stick color
    inherit.aes = FALSE
  ) +
  geom_point(
    data = FAAP100snp_data,
    aes(x = boxes, y = 1 + 10 / 2 + 0.5),  # Head just above the stick
    color = "red",  # Head color
    size = 3,
    inherit.aes = FALSE
  ) +
  # Custom colors for segments
  scale_color_manual(values = c("NO" = "grey", "S" = "#88CCEE", "FANCAA"="#C38FFF", "INSD"="#117733", "I"="#117733", "DELD"="#882255", "D"="#882255")) +
  scale_size_identity() +  # Use the size column directly
  theme_minimal() +
  labs(
    title = "Segment Bar with Mutations for Sorex araneus",
    x = "Position",
    y = NULL,  # Remove y-axis label
    color = "Segment Label"
  ) +
  theme(
    axis.text.y = element_blank(),   # Remove y-axis text
    axis.ticks.y = element_blank(), # Remove y-axis ticks
    panel.grid = element_blank(),    # Remove grid lines
    legend.position = "none"
  )





####PALB2 Boxes
#amino acid boces
PALB2boxes_4R <- read_delim("~/Sorex_Genome2/analysis/5_expression/1_for_visualizations/PALB2/PALB2boxes_4R_gap.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

PALB2snp_data<-read_table("~/Sorex_Genome2/analysis/5_expression/1_for_visualizations/PALB2/PALB2_rfigure_input.txt", 
                            col_names = TRUE)

PALB2boxes_4R$size <- ifelse(PALB2boxes_4R$domain %in% c("NO", "I", "D"), 15, 20)####
PALB2boxes_4R$end_position <- PALB2boxes_4R$end_position +1
# Plot with mutations as lollipops (corrected y positioning)
PALB2snp_data
ggplot(PALB2boxes_4R, aes(x = start_position, xend = end_position, y = 1, yend = 1, color = domain, size = size)) +
  # Segment bar
  geom_segment() +
  # Mutation lollipops with adjusted position
  geom_segment(
    data = PALB2snp_data,
    aes(x = boxes, xend = boxes, y = 1 + 10 / 2, yend = 1 + 10 / 2 + 0.5),  # Stick starts at top of the bars
    color = "black",  # Stick color
    inherit.aes = FALSE
  ) +
  geom_point(
    data = PALB2snp_data,
    aes(x = boxes, y = 1 + 10 / 2 + 0.5),  # Head just above the stick
    color = "red",  # Head color
    size = 3,
    inherit.aes = FALSE
  ) +
  # Custom colors for segments
  scale_color_manual(values = c("NO" = "grey", "WD40"="#63BEBC", "INSD"="#117733", "I"="#117733", "DELD"="#882255", "D"="#882255")) +
  scale_size_identity() +  # Use the size column directly
  theme_minimal() +
  labs(
    title = "Segment Bar with Mutations for Sorex araneus",
    x = "Position",
    y = NULL,  # Remove y-axis label
    color = "Segment Label"
  ) +
  theme(
    axis.text.y = element_blank(),   # Remove y-axis text
    axis.ticks.y = element_blank(), # Remove y-axis ticks
    panel.grid = element_blank(),    # Remove grid lines
    legend.position = "none"
  )
