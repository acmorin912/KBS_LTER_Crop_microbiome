---
  title: "PLB 847 Group Project - LTER Soybean Roots and their Fungal Community"
author: "Moriah Young"
date: "2022-11-29"
output: pdf_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
rm(list = ls()) # clear working environment
# set working directory - don't need this, working directory is git repository
#setwd("/Users/moriahyoung/Desktop/PLB 847/PLB 847 Group Project/")
# load packages
library(tidyverse)
library(stringr)
library(plotrix) #std.error
library(phyloseq) 
library(tibble)
library(plyr)
library(decontam)
library(Biostrings)
# isContamination
#library(microViz)
# Three tables are needed for the phyloseq object
# OTU
# Taxonomy
# Samples
# Sequences (as far as I can tell, this isn't necessary???)
otu_mat <- read.delim("otu_table_ITS_UPARSE_R1.txt") # OTUs - fyi, mat = matrix
tax_mat <- read.delim("constax_taxonomy.txt") # taxonomy
samples_df <- read.delim("root_fun_map_soybean_2021_1.txt") # metadata
sequences <- readDNAStringSet(otus_R1.fasta, format= "fasta", seek.first.rec = TRUE)
# change column name in "otu_mat" dataframe to match the "taxa" dataframe column
colnames(otu_mat)[colnames(otu_mat) == "X.OTU.ID"] <- "OTU_ID" 
# make it so the first column becomes the row names
# this is setting it up for the phyloseq object
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("OTU_ID") 
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("OTU_ID") 
# remove the _1 that comes after every taxa name
tax_mat[] <- lapply(tax_mat, function(x) sub("_1", "", x, fixed = TRUE))
## ignore!
# not sure if we want to do this BUT 
## define function - this will make all the empty cells in the taxa data frame now "unknown"
# or do we make these NAs which is done later on in the script?
#empty_as_unknown <- function(x){
#    if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
#    ifelse(as.character(x)!="", x, "Unknown")
#}
#tax_mat <- tax_mat %>% mutate_each(list(empty_as_unknown)) # call the function
colnames(samples_df)[colnames(samples_df) == "X.SampleID"] <- "SampleID" # get rid of that X. that shows up
samples_df <- samples_df %>% 
  tibble::column_to_rownames("SampleID") 
#samples_df2 <- samples_df[-c("NC1root", "PC1root","NC1root"),]
# turn data frames into matrices
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

physeq_object_roots <- phyloseq(OTU, TAX, samples, sequences)
physeq_object_roots

sample_data(physeq_object_roots)$is.neg <- sample_data(physeq_object_roots)$Sample_or_Control == "Control Sample"
contamdf.prev_roots <- isContaminant(physeq_object_roots, method="prevalence", neg="is.neg")
table(contamdf.prev_roots$contaminant)

ps.noncontam_roots <- prune_taxa(!contamdf.prev_roots$contaminant, physeq_object_roots)
# with contaminants removed
otu_table(ps.noncontam_roots)

ps.noncontam_roots <- subset_samples(ps.noncontam_roots, experiment%in%c("obj_1","obj_2"))
otu_table(ps.noncontam_roots) <- otu_table(ps.noncontam_roots)[which(rowSums(otu_table(ps.noncontam_roots)) >= 1),]
ps.noncontam_roots
otu_rare_roots <- as.data.frame(otu_table(ps.noncontam_roots))
metadata_roots_rare <- as.data.frame(sample_data(ps.noncontam_roots))



#roots
otu_fungi_roots <- as.data.frame(otu_table(ps.noncontam_roots))
taxa_fungi_roots <- as.data.frame(as.matrix(tax_table(ps.noncontam_roots)))
metadata_fungi_roots <- as.data.frame(as.matrix(sample_data(ps.noncontam_roots)))


#below permanova
options(scipen = 999) 
library(vegan)
library(RVAideMemoire)
model.matrix(~ Description * Management, data=metadata_fungi_roots)
model.matrix(~ Description + Management + Description : Management, data=metadata_fungi_roots)

adonis(t(otu_fungi_roots) ~ Description * Management, data=metadata_fungi_roots, permutations=9999) # by = "margin"
adonis(t(otu_fungi_roots) ~ Description + Management + Description : Management, data=metadata_fungi_roots, permutations=9999) 
adonis_fungi_roots <- adonis(t(otu_fungi_roots) ~ Description + Management + Description : Management, data=metadata_fungi_roots, permutations=9999)


vegan::vegdist(t(otu_fungi_roots), method="bray") -> dist_otu_fungi_roots
permdisp_otu_fungi_roots_M <- betadisper(dist_otu_fungi_roots, metadata_fungi_roots$Management)
permdisp_otu_fungi_roots_GS<- betadisper(dist_otu_fungi_roots, metadata_fungi_roots$Description)
anova(permdisp_otu_fungi_roots_M, permutations = 9999)
permutest(permdisp_otu_fungi_roots_M, permutations = 9999, pairwise = T)
plot(permdisp_otu_fungi_roots_M)
plot(TukeyHSD(permdisp_otu_fungi_roots_M), las=1)
boxplot(permdisp_otu_fungi_roots_M)
anova(permdisp_otu_fungi_roots_GS, permutations = 9999)
permutest(permdisp_otu_fungi_roots_GS, permutations = 9999, pairwise = T)
plot(permdisp_otu_fungi_roots_GS)
plot(TukeyHSD(permdisp_otu_fungi_roots_GS), las=1)
boxplot(permdisp_otu_fungi_roots_GS)


