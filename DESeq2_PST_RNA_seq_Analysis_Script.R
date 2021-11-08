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

library(ggsci)
library(gtools)
library(DESeq2)
library(tidyverse)
library(VennDiagram)
library(scatterplot3d)

#
#### Compiling Reads to Create Raw Counts Matrix ####

s1 <- read.table('HT-seq_Quantification_Files/sorted_1.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s2 <- read.table('HT-seq_Quantification_Files/sorted_2.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s3 <- read.table('HT-seq_Quantification_Files/sorted_3.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s4 <- read.table('HT-seq_Quantification_Files/sorted_4.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s5 <- read.table('HT-seq_Quantification_Files/sorted_5.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s6 <- read.table('HT-seq_Quantification_Files/sorted_6.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s7 <- read.table('HT-seq_Quantification_Files/sorted_7.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s8 <- read.table('HT-seq_Quantification_Files/sorted_8.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s9 <- read.table('HT-seq_Quantification_Files/sorted_9.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s10 <- read.table('HT-seq_Quantification_Files/sorted_10.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s11 <- read.table('HT-seq_Quantification_Files/sorted_11.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s12 <- read.table('HT-seq_Quantification_Files/sorted_12.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s13 <- read.table('HT-seq_Quantification_Files/sorted_13.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s14 <- read.table('HT-seq_Quantification_Files/sorted_14.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s15 <- read.table('HT-seq_Quantification_Files/sorted_15.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s16 <- read.table('HT-seq_Quantification_Files/sorted_16.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s17 <- read.table('HT-seq_Quantification_Files/sorted_17.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s18 <- read.table('HT-seq_Quantification_Files/sorted_18.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s19 <- read.table('HT-seq_Quantification_Files/sorted_19.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s20 <- read.table('HT-seq_Quantification_Files/sorted_20.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s21 <- read.table('HT-seq_Quantification_Files/sorted_21.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s22 <- read.table('HT-seq_Quantification_Files/sorted_22.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s23 <- read.table('HT-seq_Quantification_Files/sorted_23.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s24 <- read.table('HT-seq_Quantification_Files/sorted_24.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s25 <- read.table('HT-seq_Quantification_Files/sorted_25.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s26 <- read.table('HT-seq_Quantification_Files/sorted_26.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s27 <- read.table('HT-seq_Quantification_Files/sorted_27.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s28 <- read.table('HT-seq_Quantification_Files/sorted_28.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s29 <- read.table('HT-seq_Quantification_Files/sorted_29.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s30 <- read.table('HT-seq_Quantification_Files/sorted_30.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s31 <- read.table('HT-seq_Quantification_Files/sorted_31.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s32 <- read.table('HT-seq_Quantification_Files/sorted_32.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s33 <- read.table('HT-seq_Quantification_Files/sorted_33.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s34 <- read.table('HT-seq_Quantification_Files/sorted_34.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s35 <- read.table('HT-seq_Quantification_Files/sorted_35.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s36 <- read.table('HT-seq_Quantification_Files/sorted_36.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s37 <- read.table('HT-seq_Quantification_Files/sorted_37.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s38 <- read.table('HT-seq_Quantification_Files/sorted_38.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')
s39 <- read.table('HT-seq_Quantification_Files/sorted_39.yb2016.t2.count') %>% subset(substr(V1,1,1) == 'A')

countsData <- cbind(s1, s2$V2, s3$V2, s4$V2, s5$V2, s6$V2, s7$V2, s8$V2, s9$V2, s10$V2, 
      s11$V2, s12$V2, s13$V2, s14$V2, s15$V2, s16$V2, s17$V2, s18$V2, s19$V2, s20$V2, 
      s21$V2, s22$V2, s23$V2, s24$V2, s25$V2, s26$V2, s27$V2, s28$V2, s29$V2, s30$V2,
      s31$V2, s32$V2, s33$V2, s34$V2, s35$V2, s36$V2, s37$V2, s38$V2, s39$V2)

colnames(countsData) <- c('gene','Naive_0H_1','Naive_0H_2','Naive_0H_3',
                          'Mock_1H_1','Mock_1H_2','Mock_1H_3','Mock_6H_1','Mock_6H_2','Mock_6H_3',
                          'Mock_24H_1','Mock_24H_2','Mock_24H_3','Mock_48H_1','Mock_48H_2','Mock_48H_3',
                          'Vir_1H_1','Vir_1H_2','Vir_1H_3','Vir_6H_1','Vir_6H_2','Vir_6H_3',
                          'Vir_24H_1','Vir_24H_2','Vir_24H_3','Vir_48H_1','Vir_48H_2','Vir_48H_3',
                          'Avr_1H_1','Avr_1H_2','Avr_1H_3','Avr_6H_1','Avr_6H_2','Avr_6H_3',
                          'Avr_24H_1','Avr_24H_2','Avr_24H_3','Avr_48H_1','Avr_48H_2','Avr_48H_3')

write.csv(countsData, 'Output_Files/YB_RNA_seq_Raw_Counts_Matrix.csv', row.names = F)

rm(list = ls())
#


#### Normalizing data with DESeq2 ####

# Example for DEGs
countsData <- read.csv('Output_Files/YB_RNA_seq_Raw_Counts_Matrix.csv', row.names = 1)

# sample metadata
sampleInfo <- data.frame(
  
  sample_name = c('Naive_0H_1', 'Naive_0H_2', 'Naive_0H_3',
                  'Mock_1H_1', 'Mock_1H_2', 'Mock_1H_3', 'Mock_6H_1', 'Mock_6H_2', 'Mock_6H_3',
                  'Mock_24H_1', 'Mock_24H_2', 'Mock_24H_3', 'Mock_48H_1', 'Mock_48H_2', 'Mock_48H_3',
                  'Vir_1H_1', 'Vir_1H_2', 'Vir_1H_3', 'Vir_6H_1', 'Vir_6H_2', 'Vir_6H_3',
                  'Vir_24H_1', 'Vir_24H_2', 'Vir_24H_3', 'Vir_48H_1', 'Vir_48H_2', 'Vir_48H_3',
                  'Avr_1H_1', 'Avr_1H_2', 'Avr_1H_3', 'Avr_6H_1', 'Avr_6H_2', 'Avr_6H_3',
                  'Avr_24H_1', 'Avr_24H_2', 'Avr_24H_3', 'Avr_48H_1', 'Avr_48H_2', 'Avr_48H_3'),
  
  treatment_time = factor(c('Naive_0H', 'Naive_0H', 'Naive_0H',
                            'Mock_1H', 'Mock_1H', 'Mock_1H', 'Mock_6H', 'Mock_6H', 'Mock_6H',
                            'Mock_24H', 'Mock_24H', 'Mock_24H', 'Mock_48H', 'Mock_48H', 'Mock_48H',
                            'Vir_1H', 'Vir_1H', 'Vir_1H', 'Vir_6H', 'Vir_6H', 'Vir_6H',
                            'Vir_24H', 'Vir_24H', 'Vir_24H', 'Vir_48H', 'Vir_48H', 'Vir_48H',
                            'Avr_1H', 'Avr_1H', 'Avr_1H', 'Avr_6H', 'Avr_6H', 'Avr_6H',
                            'Avr_24H', 'Avr_24H', 'Avr_24H', 'Avr_48H', 'Avr_48H', 'Avr_48H'),
                          levels = c('Naive_0H',
                                     'Mock_1H','Mock_6H','Mock_24H','Mock_48H',
                                     'Vir_1H','Vir_6H','Vir_24H','Vir_48H',
                                     'Avr_1H','Avr_6H','Avr_24H','Avr_48H')),
  
  replicate = factor(rep(c(1,2,3),13), levels = c(1 ,2 ,3))
)

# Checking if sample info matches
all(colnames(countsData) == sampleInfo$sample_name)
all(sampleInfo$sample_name == paste0(sampleInfo$treatment_time, '_', sampleInfo$replicate))

# Design is used for GLM: We are focused on the interaction
# between Treatment/Time (treatment_time) and we want to account for 
# batch to batch variation (replicate).
# The last position in the formula is what is tested form while those before
# are controlled for. 
dds <- DESeqDataSetFromMatrix(countsData, sampleInfo, design = ~ replicate + treatment_time)

dds <- estimateSizeFactors(dds)
normalizedData <- counts(dds, normalized = T)

write.csv(normalizedData, 'Output_Files/DESeq2_Normalized_YB_RNA_seq.csv')

rm(list = ls())
#


#### Visualizing PR1 Expression ####

# Import Data
rawCountsData <- read.csv('Output_Files/YB_RNA_seq_Raw_Counts_Matrix.csv')
normData <- read.csv('Output_Files/DESeq2_Normalized_YB_RNA_seq.csv')

# Selecting for PR1 (AT2G14610)
rawCountsData <- rawCountsData[rawCountsData$gene == 'AT2G14610',] %>%
  mutate(Mock_0H_1 = Naive_0H_1, Vir_0H_1 = Naive_0H_1, Avr_0H_1 = Naive_0H_1, 
         Mock_0H_2 = Naive_0H_2, Vir_0H_2 = Naive_0H_2, Avr_0H_2 = Naive_0H_2,
         Mock_0H_3 = Naive_0H_3, Vir_0H_3 = Naive_0H_3, Avr_0H_3 = Naive_0H_3) %>%
  select(-c(Naive_0H_1, Naive_0H_2, Naive_0H_3))

normData <- normData[normData$gene == 'AT2G14610',] %>%
  mutate(Mock_0H_1 = Naive_0H_1, Vir_0H_1 = Naive_0H_1, Avr_0H_1 = Naive_0H_1, 
         Mock_0H_2 = Naive_0H_2, Vir_0H_2 = Naive_0H_2, Avr_0H_2 = Naive_0H_2,
         Mock_0H_3 = Naive_0H_3, Vir_0H_3 = Naive_0H_3, Avr_0H_3 = Naive_0H_3) %>%
  select(-c(Naive_0H_1, Naive_0H_2, Naive_0H_3))


# Preparing Data for Plottig
rawCountsData$norm_method <- 'Raw Counts'
normData$norm_method <- 'DESeq2 Norm'

revRawCountsData <- rawCountsData %>%
  gather('treatment','value',-norm_method,-gene) %>%
  separate('treatment', into = c('Treatment','Time','rep'))

revnormData <- normData %>%
  gather('treatment','value',-norm_method,-gene) %>%
  separate('treatment', into = c('Treatment','Time','rep'))


allData <- rbind(revRawCountsData, revnormData) %>%
  group_by(gene, norm_method, Treatment, Time) %>%
  summarise(mean = mean(value))

allData$Treatment <- factor(allData$Treatment, levels = c('Mock','Vir','Avr'))
allData$Time <- factor(allData$Time, levels = c('0H','1H','6H','24H','48H'))
allData$norm_method <- factor(allData$norm_method, 
                              levels = c('Raw Counts', 'DESeq2 Norm'))

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

# Import Data
rawCountsData <- read.csv('Output_Files/YB_RNA_seq_Raw_Counts_Matrix.csv')
normData <- read.csv('Output_Files/DESeq2_Normalized_YB_RNA_seq.csv')

# Selecting for PR2 (AT3G57260)
rawCountsData <- rawCountsData[rawCountsData$gene == 'AT3G57260',] %>%
  mutate(Mock_0H_1 = Naive_0H_1, Vir_0H_1 = Naive_0H_1, Avr_0H_1 = Naive_0H_1, 
         Mock_0H_2 = Naive_0H_2, Vir_0H_2 = Naive_0H_2, Avr_0H_2 = Naive_0H_2,
         Mock_0H_3 = Naive_0H_3, Vir_0H_3 = Naive_0H_3, Avr_0H_3 = Naive_0H_3) %>%
  select(-c(Naive_0H_1, Naive_0H_2, Naive_0H_3))

normData <- normData[normData$gene == 'AT3G57260',] %>%
  mutate(Mock_0H_1 = Naive_0H_1, Vir_0H_1 = Naive_0H_1, Avr_0H_1 = Naive_0H_1, 
         Mock_0H_2 = Naive_0H_2, Vir_0H_2 = Naive_0H_2, Avr_0H_2 = Naive_0H_2,
         Mock_0H_3 = Naive_0H_3, Vir_0H_3 = Naive_0H_3, Avr_0H_3 = Naive_0H_3) %>%
  select(-c(Naive_0H_1, Naive_0H_2, Naive_0H_3))

# Preparing Data for Plotting
rawCountsData$norm_method <- 'Raw Counts'
normData$norm_method <- 'DESeq2 Norm'

revRawCountsData <- rawCountsData %>%
  gather('treatment','value',-norm_method,-gene) %>%
  separate('treatment', into = c('Treatment','Time','rep'))

revnormData <- normData %>%
  gather('treatment','value',-norm_method,-gene) %>%
  separate('treatment', into = c('Treatment','Time','rep'))

allData <- rbind(revRawCountsData, revnormData) %>%
  group_by(gene, norm_method, Treatment, Time) %>%
  summarise(mean = mean(value))

allData$Treatment <- factor(allData$Treatment, levels = c('Mock','Vir','Avr'))
allData$Time <- factor(allData$Time, levels = c('0H','1H','6H','24H','48H'))
allData$norm_method <- factor(allData$norm_method, 
                              levels = c('Raw Counts', 'DESeq2 Norm'))

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

# Import Data
rawCountsData <- read.csv('Output_Files/YB_RNA_seq_Raw_Counts_Matrix.csv')
normData <- read.csv('Output_Files/DESeq2_Normalized_YB_RNA_seq.csv')

# Selecting for PR5 (AT1G75040)
rawCountsData <- rawCountsData[rawCountsData$gene == 'AT1G75040',] %>%
  mutate(Mock_0H_1 = Naive_0H_1, Vir_0H_1 = Naive_0H_1, Avr_0H_1 = Naive_0H_1, 
         Mock_0H_2 = Naive_0H_2, Vir_0H_2 = Naive_0H_2, Avr_0H_2 = Naive_0H_2,
         Mock_0H_3 = Naive_0H_3, Vir_0H_3 = Naive_0H_3, Avr_0H_3 = Naive_0H_3) %>%
  select(-c(Naive_0H_1, Naive_0H_2, Naive_0H_3))

normData <- normData[normData$gene == 'AT1G75040',] %>%
  mutate(Mock_0H_1 = Naive_0H_1, Vir_0H_1 = Naive_0H_1, Avr_0H_1 = Naive_0H_1, 
         Mock_0H_2 = Naive_0H_2, Vir_0H_2 = Naive_0H_2, Avr_0H_2 = Naive_0H_2,
         Mock_0H_3 = Naive_0H_3, Vir_0H_3 = Naive_0H_3, Avr_0H_3 = Naive_0H_3) %>%
  select(-c(Naive_0H_1, Naive_0H_2, Naive_0H_3))

# Preparing Data for Plottig
rawCountsData$norm_method <- 'Raw Counts'
normData$norm_method <- 'DESeq2 Norm'

revRawCountsData <- rawCountsData %>%
  gather('treatment','value',-norm_method,-gene) %>%
  separate('treatment', into = c('Treatment','Time','rep'))

revnormData <- normData %>%
  gather('treatment','value',-norm_method,-gene) %>%
  separate('treatment', into = c('Treatment','Time','rep'))

allData <- rbind(revRawCountsData, revnormData) %>%
  group_by(gene, norm_method, Treatment, Time) %>%
  summarise(mean = mean(value))

allData$Treatment <- factor(allData$Treatment, levels = c('Mock','Vir','Avr'))
allData$Time <- factor(allData$Time, levels = c('0H','1H','6H','24H','48H'))
allData$norm_method <- factor(allData$norm_method, 
                              levels = c('Raw Counts', 'DESeq2 Norm'))

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
rawCountsData <- read.csv('Output_Files/YB_RNA_seq_Raw_Counts_Matrix.csv')
normData <- read.csv('Output_Files/DESeq2_Normalized_YB_RNA_seq.csv')

sampleInfo <- read.csv('RNA_seq_sample_info.csv')

#barplot(1:4, col = c('grey40','grey80','#ED0000FF','#00468BFF'))

plotColors <- c(rep('grey80', 3), rep('grey40', 12), rep('#ED0000FF', 12), rep('#00468BFF', 12))

# Distribution Raw Counts Data
png('Plots/Raw_Data_Read_Distribution.png', width = 1000)
boxplot(log10(rawCountsData[2:40] + 1), las = 2, col = plotColors, xaxt = 'n', ylab = 'log10()')
text(x = (1:39) - 1, y = -0.75, sampleInfo$sample_name, xpd = TRUE, srt = 45)
dev.off()

# Distribution of ComBat-seq Raw Data
png('Plots/Norm_Data_Read_Distribution.png', width = 1000)
boxplot(log10(normData[2:40] + 1), las = 2, col = plotColors, xaxt = 'n', ylab = 'log10()')
text(x = (1:39) - 1, y = -0.75, sampleInfo$sample_name, xpd = TRUE, srt = 45)
dev.off()

rm(list = ls())
#


#### 2D and 3D PCA Analysis ####

# Set Genes to Rownames
rawCountsData <- read.csv('Output_Files/YB_RNA_seq_Raw_Counts_Matrix.csv', row.names = 1)
normData <- read.csv('Output_Files/DESeq2_Normalized_YB_RNA_seq.csv', row.names = 1)


# Checking if all data labels appear the same.
all(rownames(rawCountsData) == rownames(normData))

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
plotPCA(t(normData), 'DESeq2_Norm')

rm(list = ls())
#


#### Performing DEG Analysis with DESeq2 ####

# Example for DEGs
countsData <- read.csv('Output_Files/YB_RNA_seq_Raw_Counts_Matrix.csv', row.names = 1)

# sample metadata
sampleInfo <- data.frame(
  
  sample_name = c('Naive_0H_1', 'Naive_0H_2', 'Naive_0H_3',
                  'Mock_1H_1', 'Mock_1H_2', 'Mock_1H_3', 'Mock_6H_1', 'Mock_6H_2', 'Mock_6H_3',
                  'Mock_24H_1', 'Mock_24H_2', 'Mock_24H_3', 'Mock_48H_1', 'Mock_48H_2', 'Mock_48H_3',
                  'Vir_1H_1', 'Vir_1H_2', 'Vir_1H_3', 'Vir_6H_1', 'Vir_6H_2', 'Vir_6H_3',
                  'Vir_24H_1', 'Vir_24H_2', 'Vir_24H_3', 'Vir_48H_1', 'Vir_48H_2', 'Vir_48H_3',
                  'Avr_1H_1', 'Avr_1H_2', 'Avr_1H_3', 'Avr_6H_1', 'Avr_6H_2', 'Avr_6H_3',
                  'Avr_24H_1', 'Avr_24H_2', 'Avr_24H_3', 'Avr_48H_1', 'Avr_48H_2', 'Avr_48H_3'),
  
  treatment_time = factor(c('Naive_0H', 'Naive_0H', 'Naive_0H',
                            'Mock_1H', 'Mock_1H', 'Mock_1H', 'Mock_6H', 'Mock_6H', 'Mock_6H',
                            'Mock_24H', 'Mock_24H', 'Mock_24H', 'Mock_48H', 'Mock_48H', 'Mock_48H',
                            'Vir_1H', 'Vir_1H', 'Vir_1H', 'Vir_6H', 'Vir_6H', 'Vir_6H',
                            'Vir_24H', 'Vir_24H', 'Vir_24H', 'Vir_48H', 'Vir_48H', 'Vir_48H',
                            'Avr_1H', 'Avr_1H', 'Avr_1H', 'Avr_6H', 'Avr_6H', 'Avr_6H',
                            'Avr_24H', 'Avr_24H', 'Avr_24H', 'Avr_48H', 'Avr_48H', 'Avr_48H'),
                          levels = c('Naive_0H',
                                     'Mock_1H','Mock_6H','Mock_24H','Mock_48H',
                                     'Vir_1H','Vir_6H','Vir_24H','Vir_48H',
                                     'Avr_1H','Avr_6H','Avr_24H','Avr_48H')),
  
  replicate = factor(rep(c(1,2,3),13), levels = c(1 ,2 ,3))
)

#ddscabu = DESeqDataSetFromMatrix(rccabu, cdcabu, design= ~ condicabu)
# Design is used for GLM: We are focused on the interaction
# between Treatment/Time (treatment_time) and we want to account for 
# batch to batch variation (replicate).
# The last position in the formula is what is tested form, while those before
# are controlled for. 
dds <- DESeqDataSetFromMatrix(countsData, sampleInfo, design = ~ replicate + treatment_time)

# Creating DESeq Object
dds <- DESeq(dds)

# Collecting DEGs 1 Hour
DEGs1H = results(dds, contrast = c('treatment_time','Avr_1H','Vir_1H')) %>%
  data.frame() %>%
  subset(., padj <= 0.01 & abs(log2FoldChange) > 1) %>%
  mutate(kinetics = ifelse(log2FoldChange > 1, 'Upregulated','Downregulated'),
         time = '1H') %>%
  rownames_to_column('gene')

# Collecting DEGs 6 Hour
DEGs6H = results(dds, contrast = c('treatment_time','Avr_6H','Vir_6H')) %>%
  data.frame() %>%
  subset(., padj <= 0.01 & abs(log2FoldChange) > 1) %>%
  mutate(kinetics = ifelse(log2FoldChange > 1, 'Upregulated','Downregulated'),
         time = '6H') %>%
  rownames_to_column('gene')

# Collecting DEGs 24 Hour
DEGs24H = results(dds, contrast = c('treatment_time','Avr_24H','Vir_24H')) %>%
  data.frame() %>%
  subset(., padj <= 0.01 & abs(log2FoldChange) > 1) %>%
  mutate(kinetics = ifelse(log2FoldChange > 1, 'Upregulated','Downregulated'),
         time = '24H') %>%
  rownames_to_column('gene')

# Collecting DEGs 48 Hour
DEGs48H = results(dds, contrast = c('treatment_time','Avr_48H','Vir_48H')) %>%
  data.frame() %>%
  subset(., padj <= 0.01 & abs(log2FoldChange) > 1) %>%
  mutate(kinetics = ifelse(log2FoldChange > 1, 'Upregulated','Downregulated'),
         time = '48H') %>%
  rownames_to_column('gene')

# Generating a Venn Diagram indicating overlap between Avr and Vir 6, 24, 48
venn.diagram(
  x = list(DEGs6H$gene, DEGs24H$gene, DEGs48H$gene),
  category = c('A-V 6H','A-V 24H','A-V 48H'),
  filename = 'Plots/Avr_Vir_Overlap_DESeq2.png',
  output = TRUE,

# Circles
  lwd = 2,
  lty = 'blank',
  fill = c('black','blue','forestgreen'),
# Numbers
  cex = 2,
  fontface = 'bold',
  fontfamily = 'sans',

# Set names
  cat.cex = 2,
  cat.fontface = 'bold',
  cat.default.pos = 'outer',
  cat.fontfamily = 'sans'
)

allDEGs <- rbind(DEGs1H, DEGs6H, DEGs24H, DEGs48H)

write.csv(allDEGs, 'Output_Files/YB_RNA_Seq_DESeq2_Avr_Vir_DEGs.csv', row.names = F)

## Calculating DEGs for all treatments and times ##
allDEGs$comp <- 'Avr vs Vir'

# Avirulent #
# Collecting DEGs 1 Hour
DEGs1H = results(dds, contrast = c('treatment_time','Avr_1H','Mock_1H')) %>%
  data.frame() %>%
  subset(., padj <= 0.01 & abs(log2FoldChange) > 1) %>%
  mutate(kinetics = ifelse(log2FoldChange > 1, 'Upregulated','Downregulated'),
         time = '1H', comp = 'Avr vs Mock') %>%
  rownames_to_column('gene')

# Collecting DEGs 6 Hour
DEGs6H = results(dds, contrast = c('treatment_time','Avr_6H','Mock_6H')) %>%
  data.frame() %>%
  subset(., padj <= 0.01 & abs(log2FoldChange) > 1) %>%
  mutate(kinetics = ifelse(log2FoldChange > 1, 'Upregulated','Downregulated'),
         time = '6H', comp = 'Avr vs Mock') %>%
  rownames_to_column('gene')

# Collecting DEGs 24 Hour
DEGs24H = results(dds, contrast = c('treatment_time','Avr_24H','Mock_24H')) %>%
  data.frame() %>%
  subset(., padj <= 0.01 & abs(log2FoldChange) > 1) %>%
  mutate(kinetics = ifelse(log2FoldChange > 1, 'Upregulated','Downregulated'),
         time = '24H', comp = 'Avr vs Mock') %>%
  rownames_to_column('gene')

# Collecting DEGs 48 Hour
DEGs48H = results(dds, contrast = c('treatment_time','Avr_48H','Mock_48H')) %>%
  data.frame() %>%
  subset(., padj <= 0.01 & abs(log2FoldChange) > 1) %>%
  mutate(kinetics = ifelse(log2FoldChange > 1, 'Upregulated','Downregulated'),
         time = '48H', comp = 'Avr vs Mock') %>%
  rownames_to_column('gene')

# Appending Avr vs Mock DEGs
allDEGs <- rbind(allDEGs, DEGs1H, DEGs6H, DEGs24H, DEGs48H)

# Virulent #
# Collecting DEGs 1 Hour
DEGs1H = results(dds, contrast = c('treatment_time','Vir_1H','Mock_1H')) %>%
  data.frame() %>%
  subset(., padj <= 0.01 & abs(log2FoldChange) > 1) %>%
  mutate(kinetics = ifelse(log2FoldChange > 1, 'Upregulated','Downregulated'),
         time = '1H', comp = 'Vir vs Mock') %>%
  rownames_to_column('gene')

# Collecting DEGs 6 Hour
DEGs6H = results(dds, contrast = c('treatment_time','Vir_6H','Mock_6H')) %>%
  data.frame() %>%
  subset(., padj <= 0.01 & abs(log2FoldChange) > 1) %>%
  mutate(kinetics = ifelse(log2FoldChange > 1, 'Upregulated','Downregulated'),
         time = '6H', comp = 'Vir vs Mock') %>%
  rownames_to_column('gene')

# Collecting DEGs 24 Hour
DEGs24H = results(dds, contrast = c('treatment_time','Vir_24H','Mock_24H')) %>%
  data.frame() %>%
  subset(., padj <= 0.01 & abs(log2FoldChange) > 1) %>%
  mutate(kinetics = ifelse(log2FoldChange > 1, 'Upregulated','Downregulated'),
         time = '24H', comp = 'Vir vs Mock') %>%
  rownames_to_column('gene')

# Collecting DEGs 48 Hour
DEGs48H = results(dds, contrast = c('treatment_time','Vir_48H','Mock_48H')) %>%
  data.frame() %>%
  subset(., padj <= 0.01 & abs(log2FoldChange) > 1) %>%
  mutate(kinetics = ifelse(log2FoldChange > 1, 'Upregulated','Downregulated'),
         time = '48H', comp = 'Vir vs Mock') %>%
  rownames_to_column('gene')

# Appending Vir vs Mock DEGs
allDEGs <- rbind(allDEGs, DEGs1H, DEGs6H, DEGs24H, DEGs48H)

# Plotting DEGs per Time
plotData <- allDEGs %>%
  select(c('gene','time','comp')) %>%
  group_by(time, comp) %>%
  summarise(count = n())

plotData$comp <- factor(plotData$comp, levels = c('Vir vs Mock', 'Avr vs Mock', 'Avr vs Vir'))
plotData$time <- factor(plotData$time, levels = c('1H','6H','24H','48H'))

# Plotting DEGs per group and Time
DEGPlot <- ggplot(plotData, aes(x = comp, y = count, fill = comp)) + 
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
          subtitle = 'DESeq2; p-value <= 0.01; log2FC > 1 or < -1') +
  labs(x = '', y = '# of DEGs') +
  facet_wrap(~time, nrow = 1)

ggsave('Plots/DEGs_Comp_Plot.png', DEGPlot, width = 150, height = 100, units = 'mm')

rm(list = ls())
#

