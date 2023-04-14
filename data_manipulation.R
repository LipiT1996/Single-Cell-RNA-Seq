#script to manipulate gene expression data
# setwd("~/Desktop/Core Bioinformatics/scripts")

#load the libraries
library(dplyr)
library(tidyverse)
library(GEOquery)

#reading the data

data <- read.csv(file = "../Data/GSE183947_fpkm.csv")
dim(data)

#get metadata
gse <- getGEO(GEO = 'GSE183947', GSEMatrix = TRUE)
gse

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

#modify the data

metadata_subset <- select(metadata, c(1 ,10,11,17)) #novice way to extract specific info.

#lets just tidyverse and dplyr

metadata_modified<- metadata %>%
  select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1, metastasis = characteristics_ch1.1)%>%
  mutate(tissue = gsub("tissue: ", "", tissue))%>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))
  
head(data)  

#reshaping wide format to long format

data_long <- data %>%
  rename(gene = X) %>%
  gather(key = 'samples',value = 'FPKM', -gene) 
  
#join data frames: data_long + metadata_modified

data_long <- data_long %>%
  left_join(., metadata_modified, by = c('samples' = 'description'))

#explore and filter specific information

brca_data <- data_long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  group_by(gene, tissue) %>%
  summarise(mean_FPKM = mean(FPKM), 
            median_FPKM = median(FPKM)) %>%
  arrange(-mean_FPKM) 