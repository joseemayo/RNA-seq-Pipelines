# 11/08/21 - José Mayorga
#
# This script is designed to analyze an RNA-seq dataset specific to my project.
#
# ** Requirements **
# 
# Create the following folders in the working directory of the script:
#
# Output_Files
# Plots

library(sva)
library(plyr)
library(edgeR)
library(ggsci)
library(gtools)
library(ggpubr)
library(plotrix)
library(kohonen)
library(tidyverse)
library(VennDiagram)
library(scatterplot3d)


#### Accounting for Batch to Batch Variation with ComBat-seq ####

# Import Gene Count Matrix (change name as needed)
geneCountMatrix <- read.csv('JM_YB_RNA_seq_gene_count_matrix.csv')

# Remove microRNA
geneCountMatrix <- geneCountMatrix[substr(geneCountMatrix$gene_id, 1, 1) != 'a',]

# Change this sample info dataframe as needed with your sample names.
sampleInfo <- data.frame(sample_info = c('sorted_1.yb2016.t2', 'sorted_2.yb2016.t2', 'sorted_3.yb2016.t2',
                           'sorted_4.yb2016.t2', 'sorted_5.yb2016.t2', 'sorted_6.yb2016.t2',
                           'sorted_7.yb2016.t2', 'sorted_8.yb2016.t2', 'sorted_9.yb2016.t2',
                           'sorted_10.yb2016.t2', 'sorted_11.yb2016.t2', 'sorted_12.yb2016.t2',
                           'sorted_13.yb2016.t2', 'sorted_14.yb2016.t2', 'sorted_15.yb2016.t2',
                           'sorted_16.yb2016.t2', 'sorted_17.yb2016.t2', 'sorted_18.yb2016.t2',
                           'sorted_19.yb2016.t2', 'sorted_20.yb2016.t2', 'sorted_21.yb2016.t2',
                           'sorted_22.yb2016.t2', 'sorted_23.yb2016.t2', 'sorted_24.yb2016.t2',
                           'sorted_25.yb2016.t2', 'sorted_26.yb2016.t2', 'sorted_27.yb2016.t2',
                           'sorted_28.yb2016.t2', 'sorted_29.yb2016.t2', 'sorted_30.yb2016.t2',
                           'sorted_31.yb2016.t2', 'sorted_32.yb2016.t2', 'sorted_33.yb2016.t2',
                           'sorted_34.yb2016.t2', 'sorted_35.yb2016.t2', 'sorted_36.yb2016.t2',
                           'sorted_37.yb2016.t2', 'sorted_38.yb2016.t2', 'sorted_39.yb2016.t2'),
           sample_name = c('Naive_0H_1', 'Naive_0H_2', 'Naive_0H_3',
                           'Mock_1H_1', 'Mock_1H_2', 'Mock_1H_3', 'Mock_6H_1', 'Mock_6H_2', 'Mock_6H_3',
                           'Mock_24H_1', 'Mock_24H_2', 'Mock_24H_3', 'Mock_48H_1', 'Mock_48H_2', 'Mock_48H_3',
                           'Vir_1H_1', 'Vir_1H_2', 'Vir_1H_3', 'Vir_6H_1', 'Vir_6H_2', 'Vir_6H_3',
                           'Vir_24H_1', 'Vir_24H_2', 'Vir_24H_3', 'Vir_48H_1', 'Vir_48H_2', 'Vir_48H_3',
                           'Avr_1H_1', 'Avr_1H_2', 'Avr_1H_3', 'Avr_6H_1', 'Avr_6H_2', 'Avr_6H_3',
                           'Avr_24H_1', 'Avr_24H_2', 'Avr_24H_3', 'Avr_48H_1', 'Avr_48H_2', 'Avr_48H_3'))

# Adding Sample Names
rownames(geneCountMatrix) <- geneCountMatrix$gene_id
geneCountMatrix <- geneCountMatrix %>% select(-c(gene_id))
all(colnames(geneCountMatrix) == sampleInfo$sample_info)
colnames(geneCountMatrix) <- sampleInfo$sample_name

# Removing Zero Expression Genes
geneCountMatrix$sum <- rowSums(geneCountMatrix)
geneCountMatrix <- geneCountMatrix[geneCountMatrix$sum != 0,]
geneCountMatrix <- geneCountMatrix %>% select(-c(sum))

# Save Raw Data with fixed names
write.csv(geneCountMatrix %>% rownames_to_column('gene'), 'Output_Files/RNA_seq_Raw_Gene_Count_Matrix.csv', row.names = F)

# Generating batch information. Rows should match samples.
batchInfo <- rep(c(1,2,3),13)

# Calling ComBat-seq with no Covariates
adjustedGeneCountMatrix <- data.frame(ComBat_seq(as.matrix(geneCountMatrix), 
                                                 batch = batchInfo))

# Sort by gene
adjustedGeneCountMatrix <- rownames_to_column(adjustedGeneCountMatrix, 'gene')
adjustedGeneCountMatrix <- adjustedGeneCountMatrix[order(adjustedGeneCountMatrix$gene),]

write.csv(adjustedGeneCountMatrix, 'Output_Files/RNA_seq_CS_Normalized.csv', row.names = F)

rm(list = ls())

#


#### TMM Normalizing Data using edgeR ####

# Import Gene Count Matrix, Raw or Normalized
rawData <- read.csv('Output_Files/RNA_seq_Raw_Gene_Count_Matrix.csv', row.names = 1)
csData <- read.csv('Output_Files/RNA_seq_CS_Normalized.csv', row.names = 1)

# sample metadata
sampleInfo <- data.frame(
  sample_name = c('Naive_0H_1', 'Naive_0H_2', 'Naive_0H_3',
                  'Mock_1H_1', 'Mock_1H_2', 'Mock_1H_3', 'Mock_6H_1', 'Mock_6H_2', 'Mock_6H_3',
                  'Mock_24H_1', 'Mock_24H_2', 'Mock_24H_3', 'Mock_48H_1', 'Mock_48H_2', 'Mock_48H_3',
                  'Vir_1H_1', 'Vir_1H_2', 'Vir_1H_3', 'Vir_6H_1', 'Vir_6H_2', 'Vir_6H_3',
                  'Vir_24H_1', 'Vir_24H_2', 'Vir_24H_3', 'Vir_48H_1', 'Vir_48H_2', 'Vir_48H_3',
                  'Avr_1H_1', 'Avr_1H_2', 'Avr_1H_3', 'Avr_6H_1', 'Avr_6H_2', 'Avr_6H_3',
                  'Avr_24H_1', 'Avr_24H_2', 'Avr_24H_3', 'Avr_48H_1', 'Avr_48H_2', 'Avr_48H_3'),
  condition = c('Naive_0H', 'Naive_0H', 'Naive_0H',
                'Mock_1H', 'Mock_1H', 'Mock_1H', 'Mock_6H', 'Mock_6H', 'Mock_6H',
                'Mock_24H', 'Mock_24H', 'Mock_24H', 'Mock_48H', 'Mock_48H', 'Mock_48H',
                'Vir_1H', 'Vir_1H', 'Vir_1H', 'Vir_6H', 'Vir_6H', 'Vir_6H',
                'Vir_24H', 'Vir_24H', 'Vir_24H', 'Vir_48H', 'Vir_48H', 'Vir_48H',
                'Avr_1H', 'Avr_1H', 'Avr_1H', 'Avr_6H', 'Avr_6H', 'Avr_6H',
                'Avr_24H', 'Avr_24H', 'Avr_24H', 'Avr_48H', 'Avr_48H', 'Avr_48H')
  )

# Operates on multiple files for simplicity
normalize_Data <- function(geneCountMatrix, sampleInfo, name){
  
  # Generate DGE Object
  DGEObject <- DGEList(counts = geneCountMatrix, 
                       group = sampleInfo$condition)
  
  # TMM normalization
  DGEObject <- calcNormFactors(DGEObject, method = 'TMM')
  tmmNorm <- as.data.frame(cpm(DGEObject))
  
  # Sort by gene
  tmmNorm <- rownames_to_column(tmmNorm, 'gene')
  tmmNorm <- tmmNorm[order(tmmNorm$gene),]
  
  write.csv(tmmNorm, paste0('Output_Files/', name, '.csv'), row.names = F)
}

normalize_Data(rawData, sampleInfo, 'RNA_seq_Raw_TMM_Normalized')
normalize_Data(csData, sampleInfo, 'RNA_seq_CS_TMM_Normalized')

rm(list = ls())
#


#### Visualizing PR1 Expression ####

# Raw Data
rawCountsData <- read.csv('Output_Files/RNA_seq_Raw_Gene_Count_Matrix.csv')
comBatSeqDataRaw <- read.csv('Output_Files/RNA_seq_CS_Normalized.csv')

# Normalized Data
tmmCountsData <- read.csv('Output_Files/RNA_seq_Raw_TMM_Normalized.csv')
comBatSeqDatatmm <- read.csv('Output_Files/RNA_seq_CS_TMM_Normalized.csv')

# Selecting for PR1 (AT2G14610)
rawCountsData <- rawCountsData[rawCountsData$gene == 'AT2G14610',] %>%
  mutate(Mock_0H_1 = Naive_0H_1, Vir_0H_1 = Naive_0H_1, Avr_0H_1 = Naive_0H_1, 
         Mock_0H_2 = Naive_0H_2, Vir_0H_2 = Naive_0H_2, Avr_0H_2 = Naive_0H_2,
         Mock_0H_3 = Naive_0H_3, Vir_0H_3 = Naive_0H_3, Avr_0H_3 = Naive_0H_3) %>%
  select(-c(Naive_0H_1, Naive_0H_2, Naive_0H_3))

comBatSeqDataRaw <- comBatSeqDataRaw[comBatSeqDataRaw$gene == 'AT2G14610',] %>%
  mutate(Mock_0H_1 = Naive_0H_1, Vir_0H_1 = Naive_0H_1, Avr_0H_1 = Naive_0H_1, 
         Mock_0H_2 = Naive_0H_2, Vir_0H_2 = Naive_0H_2, Avr_0H_2 = Naive_0H_2,
         Mock_0H_3 = Naive_0H_3, Vir_0H_3 = Naive_0H_3, Avr_0H_3 = Naive_0H_3) %>%
  select(-c(Naive_0H_1, Naive_0H_2, Naive_0H_3))

tmmCountsData <- tmmCountsData[tmmCountsData$gene == 'AT2G14610',] %>%
  mutate(Mock_0H_1 = Naive_0H_1, Vir_0H_1 = Naive_0H_1, Avr_0H_1 = Naive_0H_1, 
         Mock_0H_2 = Naive_0H_2, Vir_0H_2 = Naive_0H_2, Avr_0H_2 = Naive_0H_2,
         Mock_0H_3 = Naive_0H_3, Vir_0H_3 = Naive_0H_3, Avr_0H_3 = Naive_0H_3) %>%
  select(-c(Naive_0H_1, Naive_0H_2, Naive_0H_3))

comBatSeqDatatmm <- comBatSeqDatatmm[comBatSeqDatatmm$gene == 'AT2G14610',] %>%
  mutate(Mock_0H_1 = Naive_0H_1, Vir_0H_1 = Naive_0H_1, Avr_0H_1 = Naive_0H_1, 
         Mock_0H_2 = Naive_0H_2, Vir_0H_2 = Naive_0H_2, Avr_0H_2 = Naive_0H_2,
         Mock_0H_3 = Naive_0H_3, Vir_0H_3 = Naive_0H_3, Avr_0H_3 = Naive_0H_3) %>%
  select(-c(Naive_0H_1, Naive_0H_2, Naive_0H_3))

# Preparing Data for Plotting
rawCountsData$norm_method <- 'Raw Counts'
comBatSeqDataRaw$norm_method <- 'ComBat-seq Raw'

tmmCountsData$norm_method <- 'TMM Counts'
comBatSeqDatatmm$norm_method <- 'ComBat-seq -> TMM'

# Gathering Data
revRawCountsData <- rawCountsData %>%
  gather('treatment','value',-norm_method,-gene) %>%
  separate('treatment', into = c('Treatment','Time','rep'))

revComBatSeqDataRaw <- comBatSeqDataRaw %>%
  gather('treatment','value',-norm_method,-gene) %>%
  separate('treatment', into = c('Treatment','Time','rep'))

revTmmCountsData <- tmmCountsData %>%
  gather('treatment','value',-norm_method,-gene) %>%
  separate('treatment', into = c('Treatment','Time','rep'))

revComBatSeqDatatmm <- comBatSeqDatatmm %>%
  gather('treatment','value',-norm_method,-gene) %>%
  separate('treatment', into = c('Treatment','Time','rep'))

# Combining data and merging replicates by mean value
allData <- rbind(revRawCountsData, revComBatSeqDataRaw, revTmmCountsData, revComBatSeqDatatmm) %>%
  group_by(gene, norm_method, Treatment, Time) %>%
  summarise(mean = mean(value))

# Setting factors for plotting order
allData$Treatment <- factor(allData$Treatment, levels = c('Mock','Vir','Avr'))
allData$Time <- factor(allData$Time, levels = c('0H','1H','6H','24H','48H'))
allData$norm_method <- factor(allData$norm_method, 
                              levels = c('Raw Counts', 'ComBat-seq Raw', 'TMM Counts', 'ComBat-seq -> TMM'))

# Plotting and saving plot
pr1Plot <- ggplot(allData, aes(x = Time, y = log2(mean), group = norm_method, colour = norm_method)) +
  geom_line(size = 1.5) +
  scale_colour_aaas() +
  theme_linedraw() + 
  theme(legend.text = element_text(size = 12),
        strip.background = element_rect(fill = 'gray'),
        strip.text = element_text(colour = 'black', size = 12),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid(. ~ Treatment) +
  ggtitle('PR1 Gene Expression') +
  labs(y = expression(log[2]('Mean Expression Value')), x = '', colour = '')

ggsave(paste0('Plots/PR1_Expression_Plot.png'), pr1Plot, 
       width = 200, height = 100, units = 'mm')

rm(list = ls())
#

#### Visualizing PR2 Expression ####

# Raw Data
rawCountsData <- read.csv('Output_Files/RNA_seq_Raw_Gene_Count_Matrix.csv')
comBatSeqDataRaw <- read.csv('Output_Files/RNA_seq_CS_Normalized.csv')

# Normalized Data
tmmCountsData <- read.csv('Output_Files/RNA_seq_Raw_TMM_Normalized.csv')
comBatSeqDatatmm <- read.csv('Output_Files/RNA_seq_CS_TMM_Normalized.csv')

# Selecting for PR2 (AT3G57260)
rawCountsData <- rawCountsData[rawCountsData$gene == 'AT3G57260',] %>%
  mutate(Mock_0H_1 = Naive_0H_1, Vir_0H_1 = Naive_0H_1, Avr_0H_1 = Naive_0H_1, 
         Mock_0H_2 = Naive_0H_2, Vir_0H_2 = Naive_0H_2, Avr_0H_2 = Naive_0H_2,
         Mock_0H_3 = Naive_0H_3, Vir_0H_3 = Naive_0H_3, Avr_0H_3 = Naive_0H_3) %>%
  select(-c(Naive_0H_1, Naive_0H_2, Naive_0H_3))

comBatSeqDataRaw <- comBatSeqDataRaw[comBatSeqDataRaw$gene == 'AT3G57260',] %>%
  mutate(Mock_0H_1 = Naive_0H_1, Vir_0H_1 = Naive_0H_1, Avr_0H_1 = Naive_0H_1, 
         Mock_0H_2 = Naive_0H_2, Vir_0H_2 = Naive_0H_2, Avr_0H_2 = Naive_0H_2,
         Mock_0H_3 = Naive_0H_3, Vir_0H_3 = Naive_0H_3, Avr_0H_3 = Naive_0H_3) %>%
  select(-c(Naive_0H_1, Naive_0H_2, Naive_0H_3))

tmmCountsData <- tmmCountsData[tmmCountsData$gene == 'AT3G57260',] %>%
  mutate(Mock_0H_1 = Naive_0H_1, Vir_0H_1 = Naive_0H_1, Avr_0H_1 = Naive_0H_1, 
         Mock_0H_2 = Naive_0H_2, Vir_0H_2 = Naive_0H_2, Avr_0H_2 = Naive_0H_2,
         Mock_0H_3 = Naive_0H_3, Vir_0H_3 = Naive_0H_3, Avr_0H_3 = Naive_0H_3) %>%
  select(-c(Naive_0H_1, Naive_0H_2, Naive_0H_3))

comBatSeqDatatmm <- comBatSeqDatatmm[comBatSeqDatatmm$gene == 'AT3G57260',] %>%
  mutate(Mock_0H_1 = Naive_0H_1, Vir_0H_1 = Naive_0H_1, Avr_0H_1 = Naive_0H_1, 
         Mock_0H_2 = Naive_0H_2, Vir_0H_2 = Naive_0H_2, Avr_0H_2 = Naive_0H_2,
         Mock_0H_3 = Naive_0H_3, Vir_0H_3 = Naive_0H_3, Avr_0H_3 = Naive_0H_3) %>%
  select(-c(Naive_0H_1, Naive_0H_2, Naive_0H_3))

# Preparing Data for Plottig
rawCountsData$norm_method <- 'Raw Counts'
comBatSeqDataRaw$norm_method <- 'ComBat-seq Raw'

tmmCountsData$norm_method <- 'TMM Counts'
comBatSeqDatatmm$norm_method <- 'ComBat-seq -> TMM'

revRawCountsData <- rawCountsData %>%
  gather('treatment','value',-norm_method,-gene) %>%
  separate('treatment', into = c('Treatment','Time','rep'))

revComBatSeqDataRaw <- comBatSeqDataRaw %>%
  gather('treatment','value',-norm_method,-gene) %>%
  separate('treatment', into = c('Treatment','Time','rep'))

revTmmCountsData <- tmmCountsData %>%
  gather('treatment','value',-norm_method,-gene) %>%
  separate('treatment', into = c('Treatment','Time','rep'))

revComBatSeqDatatmm <- comBatSeqDatatmm %>%
  gather('treatment','value',-norm_method,-gene) %>%
  separate('treatment', into = c('Treatment','Time','rep'))

allData <- rbind(revRawCountsData, revComBatSeqDataRaw, revTmmCountsData, revComBatSeqDatatmm) %>%
  group_by(gene, norm_method, Treatment, Time) %>%
  summarise(mean = mean(value))

allData$Treatment <- factor(allData$Treatment, levels = c('Mock','Vir','Avr'))
allData$Time <- factor(allData$Time, levels = c('0H','1H','6H','24H','48H'))
allData$norm_method <- factor(allData$norm_method, 
                              levels = c('Raw Counts', 'ComBat-seq Raw', 'TMM Counts', 'ComBat-seq -> TMM'))

pr2Plot <- ggplot(allData, aes(x = Time, y = log2(mean), group = norm_method, colour = norm_method)) +
  geom_line(size = 1.5) +
  scale_colour_aaas() +
  theme_linedraw() + 
  theme(legend.text = element_text(size = 12),
        strip.background = element_rect(fill = 'gray'),
        strip.text = element_text(colour = 'black', size = 12),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid(. ~ Treatment) +
  ggtitle('PR2 Gene Expression') +
  labs(y = expression(log[2]('Mean Expression Value')), x = '', colour = '')

ggsave(paste0('Plots/PR2_Expression_Plot.png'), pr2Plot, 
       width = 200, height = 100, units = 'mm')

rm(list = ls())
#

#### Visualizing PR5 Expression ####

# Raw Data
rawCountsData <- read.csv('Output_Files/RNA_seq_Raw_Gene_Count_Matrix.csv')
comBatSeqDataRaw <- read.csv('Output_Files/RNA_seq_CS_Normalized.csv')

# Normalized Data
tmmCountsData <- read.csv('Output_Files/RNA_seq_Raw_TMM_Normalized.csv')
comBatSeqDatatmm <- read.csv('Output_Files/RNA_seq_CS_TMM_Normalized.csv')

# Selecting for PR5 (AT1G75040)
rawCountsData <- rawCountsData[rawCountsData$gene == 'AT1G75040',] %>%
  mutate(Mock_0H_1 = Naive_0H_1, Vir_0H_1 = Naive_0H_1, Avr_0H_1 = Naive_0H_1, 
         Mock_0H_2 = Naive_0H_2, Vir_0H_2 = Naive_0H_2, Avr_0H_2 = Naive_0H_2,
         Mock_0H_3 = Naive_0H_3, Vir_0H_3 = Naive_0H_3, Avr_0H_3 = Naive_0H_3) %>%
  select(-c(Naive_0H_1, Naive_0H_2, Naive_0H_3))

comBatSeqDataRaw <- comBatSeqDataRaw[comBatSeqDataRaw$gene == 'AT1G75040',] %>%
  mutate(Mock_0H_1 = Naive_0H_1, Vir_0H_1 = Naive_0H_1, Avr_0H_1 = Naive_0H_1, 
         Mock_0H_2 = Naive_0H_2, Vir_0H_2 = Naive_0H_2, Avr_0H_2 = Naive_0H_2,
         Mock_0H_3 = Naive_0H_3, Vir_0H_3 = Naive_0H_3, Avr_0H_3 = Naive_0H_3) %>%
  select(-c(Naive_0H_1, Naive_0H_2, Naive_0H_3))

tmmCountsData <- tmmCountsData[tmmCountsData$gene == 'AT1G75040',] %>%
  mutate(Mock_0H_1 = Naive_0H_1, Vir_0H_1 = Naive_0H_1, Avr_0H_1 = Naive_0H_1, 
         Mock_0H_2 = Naive_0H_2, Vir_0H_2 = Naive_0H_2, Avr_0H_2 = Naive_0H_2,
         Mock_0H_3 = Naive_0H_3, Vir_0H_3 = Naive_0H_3, Avr_0H_3 = Naive_0H_3) %>%
  select(-c(Naive_0H_1, Naive_0H_2, Naive_0H_3))

comBatSeqDatatmm <- comBatSeqDatatmm[comBatSeqDatatmm$gene == 'AT1G75040',] %>%
  mutate(Mock_0H_1 = Naive_0H_1, Vir_0H_1 = Naive_0H_1, Avr_0H_1 = Naive_0H_1, 
         Mock_0H_2 = Naive_0H_2, Vir_0H_2 = Naive_0H_2, Avr_0H_2 = Naive_0H_2,
         Mock_0H_3 = Naive_0H_3, Vir_0H_3 = Naive_0H_3, Avr_0H_3 = Naive_0H_3) %>%
  select(-c(Naive_0H_1, Naive_0H_2, Naive_0H_3))

# Preparing Data for Plottig
rawCountsData$norm_method <- 'Raw Counts'
comBatSeqDataRaw$norm_method <- 'ComBat-seq Raw'

tmmCountsData$norm_method <- 'TMM Counts'
comBatSeqDatatmm$norm_method <- 'ComBat-seq -> TMM'

revRawCountsData <- rawCountsData %>%
  gather('treatment','value',-norm_method,-gene) %>%
  separate('treatment', into = c('Treatment','Time','rep'))

revComBatSeqDataRaw <- comBatSeqDataRaw %>%
  gather('treatment','value',-norm_method,-gene) %>%
  separate('treatment', into = c('Treatment','Time','rep'))

revTmmCountsData <- tmmCountsData %>%
  gather('treatment','value',-norm_method,-gene) %>%
  separate('treatment', into = c('Treatment','Time','rep'))

revComBatSeqDatatmm <- comBatSeqDatatmm %>%
  gather('treatment','value',-norm_method,-gene) %>%
  separate('treatment', into = c('Treatment','Time','rep'))

allData <- rbind(revRawCountsData, revComBatSeqDataRaw, revTmmCountsData, revComBatSeqDatatmm) %>%
  group_by(gene, norm_method, Treatment, Time) %>%
  summarise(mean = mean(value))

allData$Treatment <- factor(allData$Treatment, levels = c('Mock','Vir','Avr'))
allData$Time <- factor(allData$Time, levels = c('0H','1H','6H','24H','48H'))
allData$norm_method <- factor(allData$norm_method, 
                              levels = c('Raw Counts', 'ComBat-seq Raw', 'TMM Counts', 'ComBat-seq -> TMM'))

pr5Plot <- ggplot(allData, aes(x = Time, y = log2(mean), group = norm_method, colour = norm_method)) +
  geom_line(size = 1.5) +
  scale_colour_aaas() +
  theme_linedraw() + 
  theme(legend.text = element_text(size = 12),
        strip.background = element_rect(fill = 'gray'),
        strip.text = element_text(colour = 'black', size = 12),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid(. ~ Treatment) +
  ggtitle('PR5 Gene Expression') +
  labs(y = expression(log[2]('Mean Expression Value')), x = '', colour = '')

ggsave(paste0('Plots/PR5_Expression_Plot.png'), pr5Plot, 
       width = 200, height = 100, units = 'mm')

rm(list = ls())
#

#### Visualizing Read Distribution ####

# Raw Data
rawCountsData <- read.csv('Output_Files/RNA_seq_Raw_Gene_Count_Matrix.csv')
comBatSeqDataRaw <- read.csv('Output_Files/RNA_seq_CS_Normalized.csv')

# Normalized Data
tmmCountsData <- read.csv('Output_Files/RNA_seq_Raw_TMM_Normalized.csv')
comBatSeqDatatmm <- read.csv('Output_Files/RNA_seq_CS_TMM_Normalized.csv')

sampleInfo <- data.frame(
  sample_name = c('Naive_0H_1', 'Naive_0H_2', 'Naive_0H_3',
                  'Mock_1H_1', 'Mock_1H_2', 'Mock_1H_3', 'Mock_6H_1', 'Mock_6H_2', 'Mock_6H_3',
                  'Mock_24H_1', 'Mock_24H_2', 'Mock_24H_3', 'Mock_48H_1', 'Mock_48H_2', 'Mock_48H_3',
                  'Vir_1H_1', 'Vir_1H_2', 'Vir_1H_3', 'Vir_6H_1', 'Vir_6H_2', 'Vir_6H_3',
                  'Vir_24H_1', 'Vir_24H_2', 'Vir_24H_3', 'Vir_48H_1', 'Vir_48H_2', 'Vir_48H_3',
                  'Avr_1H_1', 'Avr_1H_2', 'Avr_1H_3', 'Avr_6H_1', 'Avr_6H_2', 'Avr_6H_3',
                  'Avr_24H_1', 'Avr_24H_2', 'Avr_24H_3', 'Avr_48H_1', 'Avr_48H_2', 'Avr_48H_3'))

#barplot(1:4, col = c('grey40','grey80','#ED0000FF','#00468BFF'))

plotColors <- c(rep('grey80', 3), rep('grey40', 12), rep('#ED0000FF', 12), rep('#00468BFF', 12))

# Distribution Raw Counts Data
png('Plots/Raw_Data_Read_Distribution.png', width = 1000)
boxplot(log10(rawCountsData[2:40] + 1), las = 2, col = plotColors, xaxt = 'n', ylab = 'log10()')
text(x = (1:39) - 1, y = -0.75, sampleInfo$sample_name, xpd = TRUE, srt = 45)
dev.off()

# Distribution of ComBat-seq Raw Data
png('Plots/CS_Raw_Data_Read_Distribution.png', width = 1000)
boxplot(log10(comBatSeqDataRaw[2:40] + 1), las = 2, col = plotColors, xaxt = 'n', ylab = 'log10()')
text(x = (1:39) - 1, y = -0.75, sampleInfo$sample_name, xpd = TRUE, srt = 45)
dev.off()

# Distribution TMM Counts Data
png('Plots/TMM_Data_Read_Distribution.png', width = 1000)
boxplot(log10(tmmCountsData[2:40] + 1), las = 2, col = plotColors, xaxt = 'n', ylab = 'log10()')
text(x = (1:39) - 1, y = -0.65, sampleInfo$sample_name, xpd = TRUE, srt = 45)
dev.off()

# Distribution of ComBat-seq TMM Data
png('Plots/CS_TMM_Data_Read_Distribution.png', width = 1000)
boxplot(log10(comBatSeqDatatmm[2:40] + 1), las = 2, col = plotColors, xaxt = 'n', ylab = 'log10()')
text(x = (1:39) - 1, y = -0.65, sampleInfo$sample_name, xpd = TRUE, srt = 45)
dev.off()

rm(list = ls())

#


#### 2D and 3D PCA Analysis ####

# Raw Data
rawCountsData <- read.csv('Output_Files/RNA_seq_Raw_Gene_Count_Matrix.csv', row.names = 1)
comBatSeqDataRaw <- read.csv('Output_Files/RNA_seq_CS_Normalized.csv', row.names = 1)

# Normalized Data
tmmCountsData <- read.csv('Output_Files/RNA_seq_Raw_TMM_Normalized.csv', row.names = 1)
comBatSeqDatatmm <- read.csv('Output_Files/RNA_seq_CS_TMM_Normalized.csv', row.names = 1)

# Checking if all data labels appear the same.
all(rownames(rawCountsData) == rownames(comBatSeqDataRaw) && rownames(rawCountsData) == rownames(tmmCountsData) 
    && rownames(tmmCountsData) == rownames(comBatSeqDatatmm))

# Automatically Calculated and Plots PCA Data
plotPCA <- function(dataSet, name){
  
  # Remove genes with exhibiting no variance
  dataSet <- dataSet[, which(apply(dataSet, 2, var) != 0)]
  
  # Calculate PCA
  myPca <- prcomp(dataSet, scale. = TRUE)
  
  # Calculate Contributions contribution
  variance <- (myPca$sdev)^2
  varPercent <- variance/sum(variance) * 100
  
  # Visualize Contributions
  #barplot(varPercent, xlab='PC', ylab='Percent Variance',
  #          names.arg=1:length(varPercent), las=1, col='gray') 
  #abline(h = 1 / 39 * 100, col='red')
  
  # Extract PC Scores 
  pcaData <- myPca$x
  pcaData <- rownames_to_column(as.data.frame(pcaData), 'category')
  pcaData <- pcaData %>%
    separate(category, into = c('treatment','time','rep'))
  
  # Prepare data for plotting
  pc1Cont <- round(varPercent[1], 2) # PC1 Contribution
  pc2Cont <- round(varPercent[2], 2) # PC2 Contribution
  pc3Cont <- round(varPercent[3], 2) # PC2 Contribution
  pcaData$treatment <- factor(pcaData$treatment, levels = c('Naive','Mock','Vir','Avr'))
  pcaData$time <- factor(pcaData$time, levels = c('0H','1H','6H','24H','48H'))
  
  colors <- c('grey80','black','#ED0000FF','#00468BFF')
  
  ## 2D Plot ##
  twoDPca <- ggplot(pcaData, aes(x = PC1, y = PC2, colour = treatment, shape = time)) +
    geom_point(size = 3) +
    geom_text(aes(label = rep), hjust = 0.5, vjust = -1, colour = 'black') +
    scale_colour_manual(values = colors) +
    theme_classic() +
    labs(y = paste0('PC2 (',pc2Cont,'%)'), 
         x = paste0('PC1 (',pc1Cont,'%)'),
         colour = 'Treatment',
         shape = 'Time')
  
  ggsave(paste0('Plots/2D_PCA_',name,'.png'), twoDPca, 
         width = 200, height = 150, units = 'mm')
  
  ## 3D Plot ##
  # Scaling Colors
  pallete <- colors
  colors <- colors[as.numeric(pcaData$treatment)]
  
  # Scaling Shapes
  shape <- c(21:25)
  shape <- shape[as.numeric(pcaData$time)]
  
  png(paste0('Plots/3D_PCA_',name,'.png'), width = 600)
  with(pcaData, {
    s3d <- scatterplot3d(PC1, PC2, PC3, 
                         xlab = paste0('PC1 (', pc1Cont, '%)'),
                         ylab = '',
                         zlab = paste0('PC3 (', pc3Cont, '%)'),
                         pch = shape,
                         color = colors)
    
    s3d.coords <- s3d$xyz.convert(PC1, PC2, PC3) # convert 3D coords to 2D projection
    
    legend('topright', legend = levels(pcaData$treatment),
           col =  pallete, pch = 19)
    
    legend('topleft', legend = levels(pcaData$time), pch = c(21:25))
    
    
    text(s3d.coords$x, s3d.coords$y, labels = rep, cex = 1, pos = 3)           
    
    # Fixing y-lab
    dims <- par('usr')
    x <- dims[1]+ 0.90 * diff(dims[1:2])
    y <- dims[3]+ 0.2 * diff(dims[3:4])
    text(x, y, paste0('PC2 (', pc2Cont, '%)'), srt = 50)
    
  })
  dev.off()
}

# Transpose Data and call function
plotPCA(t(rawCountsData), 'Raw_Data')
plotPCA(t(comBatSeqDataRaw), 'CS_Raw_Data')
plotPCA(t(tmmCountsData), 'TMM_Data')
plotPCA(t(comBatSeqDatatmm), 'CS_TMM')

rm(list = ls())

#


#### Identifying DEG using edgeR ####

# Import Data
rawCountsData <- read.csv('Output_Files/RNA_seq_Raw_Gene_Count_Matrix.csv', row.names = 1)
comBatSeqDataRaw <- read.csv('Output_Files/RNA_seq_CS_Normalized.csv', row.names = 1)

# Sample metadata
sampleInfo <- data.frame(
  sample_name = c('Naive_0H_1', 'Naive_0H_2', 'Naive_0H_3',
                  'Mock_1H_1', 'Mock_1H_2', 'Mock_1H_3', 'Mock_6H_1', 'Mock_6H_2', 'Mock_6H_3',
                  'Mock_24H_1', 'Mock_24H_2', 'Mock_24H_3', 'Mock_48H_1', 'Mock_48H_2', 'Mock_48H_3',
                  'Vir_1H_1', 'Vir_1H_2', 'Vir_1H_3', 'Vir_6H_1', 'Vir_6H_2', 'Vir_6H_3',
                  'Vir_24H_1', 'Vir_24H_2', 'Vir_24H_3', 'Vir_48H_1', 'Vir_48H_2', 'Vir_48H_3',
                  'Avr_1H_1', 'Avr_1H_2', 'Avr_1H_3', 'Avr_6H_1', 'Avr_6H_2', 'Avr_6H_3',
                  'Avr_24H_1', 'Avr_24H_2', 'Avr_24H_3', 'Avr_48H_1', 'Avr_48H_2', 'Avr_48H_3'),
  condition = c('Naive_0H', 'Naive_0H', 'Naive_0H',
                'Mock_1H', 'Mock_1H', 'Mock_1H', 'Mock_6H', 'Mock_6H', 'Mock_6H',
                'Mock_24H', 'Mock_24H', 'Mock_24H', 'Mock_48H', 'Mock_48H', 'Mock_48H',
                'Vir_1H', 'Vir_1H', 'Vir_1H', 'Vir_6H', 'Vir_6H', 'Vir_6H',
                'Vir_24H', 'Vir_24H', 'Vir_24H', 'Vir_48H', 'Vir_48H', 'Vir_48H',
                'Avr_1H', 'Avr_1H', 'Avr_1H', 'Avr_6H', 'Avr_6H', 'Avr_6H',
                'Avr_24H', 'Avr_24H', 'Avr_24H', 'Avr_48H', 'Avr_48H', 'Avr_48H')
)

determineDEGs <- function(data, name, FDR, criteria){
  
  # Generate DGE Object
  DGEObject <- DGEList(counts = data, group = sampleInfo$condition)
  
  # TMM normalization
  y <- calcNormFactors(DGEObject, method = 'TMM')
  
  # Make grouping vector
  group <- factor(sampleInfo$condition, 
                  levels = c('Naive_0H',
                             'Mock_1H','Mock_6H','Mock_24H','Mock_48H',
                             'Vir_1H','Vir_6H','Vir_24H','Vir_48H',
                             'Avr_1H','Avr_6H','Avr_24H','Avr_48H'))
  
  # Make design matrix
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  rownames(design) <- colnames(y$counts)
  
  # Estimate dispersion across tags
  y <- estimateDisp(y, design)
  
  # Fit model to counts
  fit <- glmQLFit(y, design)
  
  # Test matrix *** If you are only interested in a small subset of comparisons you can make 
  # the desired modifications ***
  myContrasts <- makeContrasts(
    Mock_1H = Mock_1H - Naive_0H,
    Mock_6H = Mock_6H - Naive_0H,
    Mock_24H = Mock_24H - Naive_0H,
    Mock_48H = Mock_48H - Naive_0H,
    
    Vir_1H = Vir_1H - Naive_0H,
    Vir_6H = Vir_6H - Naive_0H,
    Vir_24H = Vir_24H - Naive_0H,
    Vir_48H = Vir_48H - Naive_0H,
    
    Avr_1H = Avr_1H - Naive_0H,
    Avr_6H = Avr_6H - Naive_0H,
    Avr_24H = Avr_24H - Naive_0H,
    Avr_48H = Avr_48H - Naive_0H,
    
    Vir_vs_Mock_1H = (Vir_1H - Naive_0H) - (Mock_1H - Naive_0H),
    Vir_vs_Mock_6H = (Vir_6H - Naive_0H) - (Mock_6H - Naive_0H),
    Vir_vs_Mock_24H = (Vir_24H - Naive_0H) - (Mock_24H - Naive_0H),
    Vir_vs_Mock_48H = (Vir_48H - Naive_0H) - (Mock_48H - Naive_0H),
    
    Avr_vs_Mock_1H = (Avr_1H - Naive_0H) - (Mock_1H - Naive_0H),
    Avr_vs_Mock_6H = (Avr_6H - Naive_0H) - (Mock_6H - Naive_0H),
    Avr_vs_Mock_24H = (Avr_24H - Naive_0H) - (Mock_24H - Naive_0H),
    Avr_vs_Mock_48H = (Avr_48H - Naive_0H) - (Mock_48H - Naive_0H),
    
    Avr_vs_Vir_1H = (Avr_1H - Naive_0H) - (Vir_1H - Naive_0H),
    Avr_vs_Vir_6H = (Avr_6H - Naive_0H) - (Vir_6H - Naive_0H),
    Avr_vs_Vir_24H = (Avr_24H - Naive_0H) - (Vir_24H - Naive_0H),
    Avr_vs_Vir_48H = (Avr_48H - Naive_0H) - (Vir_48H - Naive_0H),
    levels = design)
  
  # Test for DEGs using differential timepoints
  V1TopTags <- topTags(glmTreat(fit, contrast = myContrasts[,'Vir_vs_Mock_1H']), n = Inf)
  V6TopTags <- topTags(glmTreat(fit, contrast = myContrasts[,'Vir_vs_Mock_6H']), n = Inf)
  V24TopTags <- topTags(glmTreat(fit, contrast = myContrasts[,'Vir_vs_Mock_24H']), n = Inf)
  V48TopTags <- topTags(glmTreat(fit, contrast = myContrasts[,'Vir_vs_Mock_48H']), n = Inf)
  
  A1TopTags <- topTags(glmTreat(fit, contrast = myContrasts[,'Avr_vs_Mock_1H']), n = Inf)
  A6TopTags <- topTags(glmTreat(fit, contrast = myContrasts[,'Avr_vs_Mock_6H']), n = Inf)
  A24TopTags <- topTags(glmTreat(fit, contrast = myContrasts[,'Avr_vs_Mock_24H']), n = Inf)
  A48TopTags <- topTags(glmTreat(fit, contrast = myContrasts[,'Avr_vs_Mock_48H']), n = Inf)
  
  A1_V1TopTags <- topTags(glmTreat(fit, contrast = myContrasts[,'Avr_vs_Vir_1H']), n = Inf)
  A6_V6TopTags <- topTags(glmTreat(fit, contrast = myContrasts[,'Avr_vs_Vir_6H']), n = Inf)
  A24_V24TopTags <- topTags(glmTreat(fit, contrast = myContrasts[,'Avr_vs_Vir_24H']), n = Inf)
  A48_V48TopTags <- topTags(glmTreat(fit, contrast = myContrasts[,'Avr_vs_Vir_48H']), n = Inf)
  
  # Selecting based on FDR <= criteria
  V1TopTags <- V1TopTags[V1TopTags$table$FDR <= FDR,]
  V6TopTags <- V6TopTags[V6TopTags$table$FDR <= FDR,]
  V24TopTags <- V24TopTags[V24TopTags$table$FDR <= FDR,]
  V48TopTags <- V48TopTags[V48TopTags$table$FDR <= FDR,]
  
  A1TopTags <- A1TopTags[A1TopTags$table$FDR <= FDR,]
  A6TopTags <- A6TopTags[A6TopTags$table$FDR <= FDR,]
  A24TopTags <- A24TopTags[A24TopTags$table$FDR <= FDR,]
  A48TopTags <- A48TopTags[A48TopTags$table$FDR <= FDR,]
  
  A1_V1TopTags <- A1_V1TopTags[A1_V1TopTags$table$FDR <= FDR,]
  A6_V6TopTags <- A6_V6TopTags[A6_V6TopTags$table$FDR <= FDR,]
  A24_V24TopTags <- A24_V24TopTags[A24_V24TopTags$table$FDR <= FDR,]
  A48_V48TopTags <- A48_V48TopTags[A48_V48TopTags$table$FDR <= FDR,]
  
  # Selecting for  genes with FC  > | Value |
  Vir1_Mock <- rownames(V1TopTags$table[V1TopTags$table$logFC > criteria | V1TopTags$table$logFC < -criteria,])
  Vir6_Mock <- rownames(V6TopTags$table[V6TopTags$table$logFC > criteria | V6TopTags$table$logFC < -criteria,])
  Vir24_Mock <- rownames(V24TopTags$table[V24TopTags$table$logFC > criteria | V24TopTags$table$logFC < -criteria,])
  Vir48_Mock <- rownames(V48TopTags$table[V48TopTags$table$logFC > criteria | V48TopTags$table$logFC < -criteria,])
  
  Avr1_Mock <- rownames(A1TopTags$table[A1TopTags$table$logFC > criteria | A1TopTags$table$logFC < -criteria,])
  Avr6_Mock <- rownames(A6TopTags$table[A6TopTags$table$logFC > criteria | A6TopTags$table$logFC < -criteria,])
  Avr24_Mock <- rownames(A24TopTags$table[A24TopTags$table$logFC > criteria | A24TopTags$table$logFC < -criteria,])
  Avr48_Mock <- rownames(A48TopTags$table[A48TopTags$table$logFC > criteria | A48TopTags$table$logFC < -criteria,])
  
  Avr1_Vir <- rownames(A1_V1TopTags$table[A1_V1TopTags$table$logFC > criteria | A1_V1TopTags$table$logFC < -criteria,])
  Avr6_Vir <- rownames(A6_V6TopTags$table[A6_V6TopTags$table$logFC > criteria | A6_V6TopTags$table$logFC < -criteria,])
  Avr24_Vir <- rownames(A24_V24TopTags$table[A24_V24TopTags$table$logFC > criteria | A24_V24TopTags$table$logFC < -criteria,])
  Avr48_Vir <- rownames(A48_V48TopTags$table[A48_V48TopTags$table$logFC > criteria | A48_V48TopTags$table$logFC < -criteria,])
  
  # If you are interested in retaining the data, simply remove the 'select' function
  #Avr1VirDEGs <- data.frame(data.frame(A1_V1TopTags[row.names(A1_V1TopTags) %in% Avr1_Vir,]) %>% select(-c('unshrunk.logFC', 'logCPM', 'PValue')),
  #                          comparison = 'Avr1_Vir')
  Avr6VirDEGs <- data.frame(data.frame(A6_V6TopTags[row.names(A6_V6TopTags) %in% Avr6_Vir,]) %>% select(-c('unshrunk.logFC', 'logCPM', 'PValue')),
                            comparison = 'Avr6_Vir')
  
  # Selecting for all data within the Avr-Vir 6 hpi pair
  #Avr6VirDEGs <- data.frame(data.frame(A6_V6TopTags) %>% select(-c('unshrunk.logFC', 'logCPM')))
  #write.csv(Avr6VirDEGs, 'Output_Files/all_DEG_data_Avr6_Vir.csv')
  
  Avr24VirDEGs <- data.frame(data.frame(A24_V24TopTags[row.names(A24_V24TopTags) %in% Avr24_Vir,]) %>% select(-c('unshrunk.logFC', 'logCPM', 'PValue')),
                             comparison = 'Avr24_Vir')
  Avr48VirDEGs <- data.frame(data.frame(A48_V48TopTags[row.names(A48_V48TopTags) %in% Avr48_Vir,]) %>% select(-c('unshrunk.logFC', 'logCPM', 'PValue')), 
                             comparison = 'Avr48_Vir')
  
  # Avr1VirDEGs is empty, so it is omitted
  allDEGsAvrVir <- rbind(Avr6VirDEGs, Avr24VirDEGs, Avr48VirDEGs)
  
  write.csv(allDEGsAvrVir, 'Output_Files/YB_RNA_Seq_DEG_CS_TMM_0.01_1_Avr_Vir.csv')
  
  # Saving Avr6_Vir DEGs
  #write.csv(data.frame(gene = Avr6_Vir), 'Output_Files/DEG_CS_TMM_0.01_1_Avr6_Vir.csv', row.names = F)
  
  # Count DEGs 
  DEGenes <- data.frame(
    name = c('Vir vs Mock_1H','Vir vs Mock_6H','Vir vs Mock_24H','Vir vs Mock_48H',
             'Avr vs Mock_1H','Avr vs Mock_6H','Avr vs Mock_24H','Avr vs Mock_48H',
             'Avr vs Vir_1H','Avr vs Vir_6H','Avr vs Vir_24H','Avr vs Vir_48H'),
    
    count = c(length(Vir1_Mock), length(Vir6_Mock), length(Vir24_Mock), length(Vir48_Mock),
              length(Avr1_Mock), length(Avr6_Mock), length(Avr24_Mock), length(Avr48_Mock),
              length(Avr1_Vir), length(Avr6_Vir), length(Avr24_Vir), length(Avr48_Vir)))
  
  DEGenes <- DEGenes %>% separate(name, into = c('condition','time'), sep = '_')
  DEGenes$condition <- factor(DEGenes$condition, levels = c('Vir vs Mock','Avr vs Mock','Avr vs Vir'))
  DEGenes$time <- factor(DEGenes$time, levels = c('1H','6H','24H','48H'))
  
  # All DEGs 
  allDEGs <- data.frame(
    gene = unique(c(Vir1_Mock, Vir6_Mock, Vir24_Mock, Vir48_Mock,
                    Avr1_Mock, Avr6_Mock, Avr24_Mock, Avr48_Mock,
                    Avr1_Vir, Avr6_Vir, Avr24_Vir, Avr48_Vir)))
  
  # All DEGs 
  avrVirDEGs <- data.frame(gene = unique(c(Avr1_Vir, Avr6_Vir, Avr24_Vir, Avr48_Vir)))
  
  # Plotting DEGs per group and Time
  DEGPlot <- ggplot(DEGenes, aes(x = condition, y = count, fill = condition)) + 
    geom_bar(stat = 'identity') +
    theme_linedraw() +
    theme(strip.background = element_rect(fill = 'gray'),
          strip.text = element_text(colour = 'black'),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 12)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)), limits = c(0,8000)) +
    scale_fill_manual(values = c('grey40','#00468BFF','#ED0000FF')) +
    ggtitle(label = 'DEGs', 
            subtitle = paste0('FDR <= ', FDR, '; log2() fold change >= ', criteria, '\n', 
                              round(nrow(avrVirDEGs) / nrow(allDEGs), 3) * 100)) +
    labs(x = '', y = '# of DEGs') +
    facet_wrap(. ~ time, nrow = 1)
  
  ggsave(paste0('Plots/DEG_', name, '_' , FDR, '_', criteria ,'_Plot.png'), DEGPlot, 
         width = 150, height = 100, units = 'mm')
  
  # Adding Data Label and saving DEGs
  avrVirDEGs$dynamics <- 'Avr_vs_Vir'
  
  write.csv(avrVirDEGs, paste0('Output_Files/DEG_', name, '_' , FDR, '_', criteria ,'.csv'), row.names = F)
}

# Calling the function with various FDR and criteria (log2FC selection)
#determineDEGs(rawCountsData, 'RAW_TMM', 0.05, 1)
#determineDEGs(comBatSeqDataRaw, 'CS_TMM', 0.05, 1)
determineDEGs(rawCountsData, 'RAW_TMM', 0.01, 1)
determineDEGs(comBatSeqDataRaw, 'CS_TMM', 0.01, 1)

#determineDEGs(rawCountsData, 'RAW_TMM', 0.05, 0.5)
#determineDEGs(comBatSeqDataRaw, 'CS_TMM', 0.05, 0.5)
#determineDEGs(rawCountsData, 'RAW_TMM', 0.01, 0.5)
#determineDEGs(comBatSeqDataRaw, 'CS_TMM', 0.01, 0.5)

rm(list = ls())

#

