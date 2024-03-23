


library("BiocManager")
# You can define a list of samples to query and download providing relative TCGA barcodes.
listSamples <- c(
  "TCGA-E9-A1NG-11A-52R-A14M-07","TCGA-BH-A1FC-11A-32R-A13Q-07",
  "TCGA-A7-A13G-11A-51R-A13Q-07","TCGA-BH-A0DK-11A-13R-A089-07",
  "TCGA-E9-A1RH-11A-34R-A169-07","TCGA-BH-A0AU-01A-11R-A12P-07",
  "TCGA-C8-A1HJ-01A-11R-A13Q-07","TCGA-A7-A13D-01A-13R-A12P-07",
  "TCGA-A2-A0CV-01A-31R-A115-07","TCGA-AQ-A0Y5-01A-11R-A14M-07"
)

# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  barcode = listSamples )

# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
BRCA.Rnaseq.SE <- GDCprepare(query)

BRCAMatrix <- assay(BRCA.Rnaseq.SE,"unstranded")
# For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
BRCA.RNAseq_CorOutliers <- TCGAanalyze_Preprocessing(BRCA.Rnaseq.SE)



library(TCGAbiolinks)

if (!require("EDASeq", quietly = TRUE))
  install.packages("EDASeq")
install.packages("EDASeq")
library("EDASeq")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EDASeq")


# normalization of genes
dataNorm <- TCGAanalyze_Normalization(
  tabDF = BRCA.RNAseq_CorOutliers, 
  geneInfo =  geneInfoHT
)

# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile", 
  qnt.cut =  0.25
)

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(
  barcode = colnames(dataFilt),
  typesample = c("NT")
)

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(
  barcode = colnames(dataFilt), 
  typesample = c("TP")
)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")


# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[,samplesNT],
  mat2 = dataFilt[,samplesTP],
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT"
)

# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(
  FC_FDR_table_mRNA = dataDEGs,
  typeCond1 = "Tumor",
  typeCond2 = "Normal",
  TableCond1 = dataFilt[,samplesTP],
  TableCond2 = dataFilt[,samplesNT]
)
head(dataDEGsFiltLevel)

# LogFC değeri sıfırdan büyük olan satırları seçme
logFC_greater_than_zero <- dataDEGsFiltLevel[dataDEGsFiltLevel$logFC > 0, ]

# Sonuçları görüntüleme
head(logFC_greater_than_zero)

library(ggplot2)

install.packages("ggpubr")
library(ggpubr)

# ggpubr ve ggplot2 paketlerini yükle
library(ggpubr)
library(ggplot2)

# ggpubr ve ggplot2 paketlerini yükle
library(ggpubr)
library(ggplot2)

# Veri görselleştirme kütüphanesini yükle
library(ggplot2)

# 'logFC' sütunundaki değerleri al
logFC_values <- dataDEGsFiltLevel$logFC

# Histogram oluştur ve 0'dan büyük ve küçük değerleri farklı renklerde göster
ggplot(dataDEGsFiltLevel, aes(x = mRNA)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  labs(title = "logFC Değerlerinin Histogramı",
       x = "logFC",
       y = "Frekans") +
  theme_minimal()

#


head(dataDEGsFiltLevel)

class(dataDEGsFiltLevel)

library(ggplot2)



gene_names <- dataDEGsFiltLevel$mRNA[1:150]

logfc_values<- dataDEGsFiltLevel$logFC[1:150]

data_frame <- data.frame(gene_names, logfc_values)

data_frame
