        #Ashlynn Morin
    #THE IMPACTS OF MANAGEMENT REGIME AND CROP ROTATION ON MICROBIAL COMMUNITIES
#In the USA, a growing movement towards a three-crop rotation involving maize (Zea mays), wheat (Triticum aestivum), and soybeans (Glycine max)
#are in a broader push towards sustainable agriculture practices, such as reducing chemical inputs and preservation of soils.
#This study, utilizing the Kellogg Biological Station (KBS) Long-Term Ecological Research (LTER) site, aims to understand how differing
#agricultural management practices affect the microbial assembly longitudinally. Annually, we collected leaf, stem, root, and soil samples
#from conventional, no-tillage, and organically managed plots at different growth stages. Through 16S and Internal Transcribed Spacer (ITS)
#region rRNA amplicon data for bacteria and fungi respectively, we analyzed a full crop rotation of microbial data to visualize the
#community dynamics and stability of microbial communities across varying agricultural management systems. Previously, Dr. Reid Longley using
#the 2018 soybean microbiome data had come to the significant conclusion that management regime did affect the microbiome,
#an insight expected to continue across the rest of the years. Continuation and expansion of this previous study aims to provide
#a more prolonged approach and understanding of ‘sustainable’ agriculture as it continues to gain traction amongst farmers.
#A study such as this is not only adding to the knowledge of what differences a long-term management systems does to microbial communities,
#but also laying groundwork for future and more comprehensive studies that address the more multifaceted nature of agroecosystems and their long-term impacts.

      #Started 07-15-24

#load in Packages
library(dplyr)
library(tidyverse)
library(phyloseq)
library(vegan)
library(S4Vectors)
library(decontam)
library(devtools)
library(ggpubr)
library(tidyr)
library(tibble)
library("indicspecies")
library("ComplexHeatmap")
library("circlize")
library(Biostrings)
library(yaml)
library(colorspace)
library(ggplot2)
library(indicspecies)
library(randomForest)

#set working directory
setwd("~/Desktop/Paper/R_code")

###Creation of phyloseq object
  #ASV
ASV_fungi<-read.delim("Fungi_ASVtable_180bp.txt", row.names = 1)
head(ASV_fungi)
ps_ASV<-otu_table(ASV_fungi, taxa_are_rows = TRUE)
  #metadata
Meta_fungi<-read.delim("fungi_metadata.txt", row.names = 1)
head(Meta_fungi)
ps_meta<-sample_data(Meta_fungi)
  #taxonomy
Tax_fungi<- read.delim("fungi_constax_taxonomy_180.txt", header = TRUE, row.names = 1)
head(Tax_fungi)
ps_tax<-tax_table(as.matrix(Tax_fungi))
  #fasta
fasta_fungi<- readDNAStringSet("asv_180bp.fasta", format = "fasta", seek.first.rec = TRUE, use.names = TRUE)
head(fasta_fungi)
    ##THE OBJECT##
ps_object<- phyloseq(ps_ASV,
                     ps_meta,
                     ps_tax,
                     fasta_fungi)
ps_object
tax_table(ps_object)
sample_data(ps_object)

#format the taxonomy table
tax_table(ps_object)
tax_table(ps_object)[, "Kingdom"] <- gsub("_1", "", tax_table(ps_object)[, "Kingdom"])
tax_table(ps_object)[, "Phylum"] <- gsub("_1", "", tax_table(ps_object)[, "Phylum"])
tax_table(ps_object)[, "Class"] <- gsub("_1", "", tax_table(ps_object)[, "Class"])
tax_table(ps_object)[, "Order"] <- gsub("_1", "", tax_table(ps_object)[, "Order"])
tax_table(ps_object)[, "Family"] <- gsub("_1", "", tax_table(ps_object)[, "Family"])
tax_table(ps_object)[, "Genus"] <- gsub("_1", "", tax_table(ps_object)[, "Genus"])
tax_table(ps_object)[, "Species"] <- gsub("_1", "", tax_table(ps_object)[, "Species"])
tax_table(ps_object)
ps_object <- subset_taxa(ps_object, Phylum!="Chloroplast")
ps_object <- subset_taxa(ps_object, Class!="Chloroplast")
ps_object <- subset_taxa(ps_object, Order!="Chloroplast")
ps_object <- subset_taxa(ps_object, Family!="Chloroplast")
ps_object <- subset_taxa(ps_object, Genus!="Chloroplast")
tax_table(ps_object)
ps_object <- subset_taxa(ps_object, Phylum!="Mitochondria")
ps_object <- subset_taxa(ps_object, Class!="Mitochondria")
ps_object <- subset_taxa(ps_object, Order!="Mitochondria")
ps_object <- subset_taxa(ps_object, Family!="Mitochondria")
ps_object <- subset_taxa(ps_object, Genus!="Mitochondria")
tax_table(ps_object)
ps_object <- subset_taxa(ps_object, Kingdom!="Anthophyta")
ps_object <- subset_taxa(ps_object, Kingdom!="Alveolata")
ps_object <- subset_taxa(ps_object, Kingdom!="Ichthyosporia")
ps_object <- subset_taxa(ps_object, Kingdom!="Protista")
ps_object <- subset_taxa(ps_object, Kingdom!="Metazoa")
ps_object <- subset_taxa(ps_object, Kingdom!="Rhizaria")
ps_object <- subset_taxa(ps_object, Kingdom!="Viridiplantae")
sort(unique(as.data.frame(tax_table(ps_object))$Kingdom))
ps_object

#remove obj_1 (previous fungicide experiment that was included)
ps_obj1<-subset_samples(ps_object, Experiment%in%c("obj_1"))
sample_data(ps_obj1)

#subset by year (for ease)
ps_2018<-subset_samples(ps_obj1, Year%in%c("2018"))
  sample_data(ps_2018)
ps_2019<-subset_samples(ps_obj1, Year%in%c("2019"))
  sample_data(ps_2019)
ps_2020<-subset_samples(ps_obj1, Year%in%c("2020"))
  sample_data(ps_2020)
ps_2021<-subset_samples(ps_obj1, Year%in%c("2021"))
  sample_data(ps_2021)
  
  ###Check library distribution
#2018
df_ps_2018 <- as.data.frame(sample_data(ps_2018))
df_ps_2018$LibrarySize_leaf <- sample_sums(ps_2018)
df_ps_2018 <- df_ps_2018[order(df_ps_2018$LibrarySize_leaf),]
df_ps_2018$Index <- seq(nrow(df_ps_2018))
ggplot(data=df_ps_2018, aes(x=Index, y=LibrarySize_leaf, color=Sample_or_Control)) + geom_point()
#2019
df_ps_2018 <- as.data.frame(sample_data(ps_2018))
df_ps_2018$LibrarySize_leaf <- sample_sums(ps_2018)
df_ps_2018 <- df_ps_2018[order(df_ps_2018$LibrarySize_leaf),]
df_ps_2018$Index <- seq(nrow(df_ps_2018))
ggplot(data=df_ps_2018, aes(x=Index, y=LibrarySize_leaf, color=Sample_or_Control)) + geom_point()
#2020
df_ps_2020 <- as.data.frame(sample_data(ps_2020))
df_ps_2020$LibrarySize_leaf <- sample_sums(ps_2020)
df_ps_2020 <- df_ps_2020[order(df_ps_2020$LibrarySize_leaf),]
df_ps_2020$Index <- seq(nrow(df_ps_2020))
ggplot(data=df_ps_2020, aes(x=Index, y=LibrarySize_leaf, color=Sample_or_Control)) + geom_point()
#2021
df_ps_2021 <- as.data.frame(sample_data(ps_2021))
df_ps_2021$LibrarySize_leaf <- sample_sums(ps_2021)
df_ps_2021 <- df_ps_2021[order(df_ps_2021$LibrarySize_leaf),]
df_ps_2021$Index <- seq(nrow(df_ps_2021))
ggplot(data=df_ps_2021, aes(x=Index, y=LibrarySize_leaf, color=Sample_or_Control)) + geom_point()

#filter by prevalence
  #2018
sample_data(ps_2018)$is.neg <- sample_data(ps_2018)$Sample_or_Control == "Control_Sample"
contamdf.prev_2018 <- isContaminant(ps_2018, method="prevalence", neg="is.neg")
table(contamdf.prev_2018$contaminant)
  #2019
sample_data(ps_2019)$is.neg <- sample_data(ps_2019)$Sample_or_Control == "Control_Sample"
contamdf.prev_2019 <- isContaminant(ps_2019, method="prevalence", neg="is.neg")
table(contamdf.prev_2019$contaminant)
  #2020
sample_data(ps_2020)$is.neg <- sample_data(ps_2020)$Sample_or_Control == "Control_Sample"
contamdf.prev_2020 <- isContaminant(ps_2020, method="prevalence", neg="is.neg")
table(contamdf.prev_2020$contaminant)
  #2021
sample_data(ps_2021)$is.neg <- sample_data(ps_2021)$Sample_or_Control == "Control_Sample"
contamdf.prev_2021 <- isContaminant(ps_2021, method="prevalence", neg="is.neg")
table(contamdf.prev_2021$contaminant)

###Remove contaminants and negative controls
  #2018
ps.noncontam_2018 <- prune_taxa(!contamdf.prev_2018$contaminant, ps_2018) #remove contaminants
ps.noncontam_2018<-subset_samples(ps.noncontam_2018,Sample_or_Control%in%c("True_Sample")) #remove negative control
ps.noncontam_2018
otu_table(ps.noncontam_2018) <- otu_table(ps.noncontam_2018)[which(rowSums(otu_table(ps.noncontam_2018)) >= 1),]
ps.noncontam_2018
otu_table(ps.noncontam_2018)
sample_data(ps.noncontam_2018)
ps.noncontam_2018
  #2019
ps.noncontam_2019 <- prune_taxa(!contamdf.prev_2019$contaminant, ps_2019) #remove contaminants
ps.noncontam_2019<-subset_samples(ps.noncontam_2019,Sample_or_Control%in%c("True_Sample")) #remove negative control
ps.noncontam_2019
otu_table(ps.noncontam_2019) <- otu_table(ps.noncontam_2019)[which(rowSums(otu_table(ps.noncontam_2019)) >= 1),]
ps.noncontam_2019
otu_table(ps.noncontam_2019)
sample_data(ps.noncontam_2019)
ps.noncontam_2019
  #2020
ps.noncontam_2020 <- prune_taxa(!contamdf.prev_2020$contaminant, ps_2020) #remove contaminants
ps.noncontam_2020<-subset_samples(ps.noncontam_2020,Sample_or_Control%in%c("True_Sample")) #remove negative control
ps.noncontam_2020
otu_table(ps.noncontam_2020) <- otu_table(ps.noncontam_2020)[which(rowSums(otu_table(ps.noncontam_2020)) >= 1),]
ps.noncontam_2020
otu_table(ps.noncontam_2020)
sample_data(ps.noncontam_2020)
ps.noncontam_2020
  #2021
ps.noncontam_2021 <- prune_taxa(!contamdf.prev_2021$contaminant, ps_2021) #remove contaminants
ps.noncontam_2021<-subset_samples(ps.noncontam_2021,Sample_or_Control%in%c("True_Sample")) #remove negative control
ps.noncontam_2021
otu_table(ps.noncontam_2021) <- otu_table(ps.noncontam_2021)[which(rowSums(otu_table(ps.noncontam_2021)) >= 1),]
ps.noncontam_2021
otu_table(ps.noncontam_2021)
sample_data(ps.noncontam_2021)
ps.noncontam_2021

##Export otu table to check samples
# Following removal of contaminants identified by the decontam package, removing any samples which had 
# less than 1000 reads to avoid biasing beta diversity analyses
# removing samples with less than 1000 reads, including samples that had less than 1000 read following sampling
# using csv files to check
  #2018
write.csv(otu_table(ps.noncontam_2018), file = "filtering/filtering_low_2018.csv")
sum_reads_2018 <- rowSums(otu_table(ps.noncontam_2018))
samples_to_remove_2018<- names(which(sum_reads_2018 <1000))
ps.filtered_2018 <- prune_samples(!sample_names(ps.noncontam_2018) %in% samples_to_remove_2018, ps.noncontam_2018)
write.csv(otu_table(ps.filtered_2018), file = "filtering/supposed_filtered_low_2018.csv")
otu_table(ps.filtered_2018) <- subset(otu_table(ps.filtered_2018),
                                      select = -c(sample1695,sample1801,sample1814,sample1828,sample1829,sample1895,sample1917,
                                                  sample1947,sample1948,sample1994,sample2035,sample2038,sample2061,sample2096,
                                                  sample2150,sample2154,sample2156,sample2159,sample2162,sample2163,sample2164,
                                                  sample2165,sample2171,sample2172,sample2176,sample2177,sample2179,sample2181,
                                                  sample2186,sample2187,sample2188,sample2191,sample2192,sample2198,sample2198,
                                                  sample2214,sample2218,sample2223,sample2225,sample2229,sample2234,sample2437,
                                                  sample899,sample2207))
write.csv(otu_table(ps.filtered_2018), file = "filtering/doublechecked_filtered_low_2018.csv")
  #2019
write.csv(otu_table(ps.noncontam_2019), file = "filtering/filtering_low_2019.csv")
sum_reads_2019 <- rowSums(otu_table(ps.noncontam_2019))
samples_to_remove_2019<- names(which(sum_reads_2019 <1000))
ps.filtered_2019 <- prune_samples(!sample_names(ps.noncontam_2019) %in% samples_to_remove_2019, ps.noncontam_2019)
write.csv(otu_table(ps.filtered_2019), file = "filtering/supposed_filtered_low_2019.csv")
otu_table(ps.filtered_2019) <- subset(otu_table(ps.filtered_2019),
                                      select = -c(sample083,sample086,sample097,sample098,sample101,sample106,sample109,sample113,
                                                  sample124,sample132,sample137,sample139,sample147,sample152,sample197,sample296,
                                                  sample303,sample336,sample342,sample363,sample374,sample375))
write.csv(otu_table(ps.filtered_2019), file = "filtering/doublechecked_filtered_low_2019.csv")
  #2020
write.csv(otu_table(ps.noncontam_2020), file = "filtering/filtering_low_2020.csv")
sum_reads_2020 <- rowSums(otu_table(ps.noncontam_2020))
samples_to_remove_2020<- names(which(sum_reads_2020 <1000))
ps.filtered_2020 <- prune_samples(!sample_names(ps.noncontam_2020) %in% samples_to_remove_2020, ps.noncontam_2020)
write.csv(otu_table(ps.filtered_2020), file = "filtering/supposed_filtered_low_2020.csv")
otu_table(ps.filtered_2020) <- subset(otu_table(ps.filtered_2020),
                                      select = -c(sample1099,sample629,sample644,sample672,sample695,sample706,
                                                  sample839,sample841,sample850,sample852,sample885,sample886,
                                                  sample963,sample975))
write.csv(otu_table(ps.filtered_2020), file = "filtering/doublechecked_filtered_low_2020.csv")
  #2021
write.csv(otu_table(ps.noncontam_2021), file = "filtering/filtering_low_2021.csv")
sum_reads_2021 <- rowSums(otu_table(ps.noncontam_2021))
samples_to_remove_2021<- names(which(sum_reads_2021 <1000))
ps.filtered_2021 <- prune_samples(!sample_names(ps.noncontam_2021) %in% samples_to_remove_2021, ps.noncontam_2021)
write.csv(otu_table(ps.filtered_2021), file = "filtering/supposed_filtered_low_2021.csv")
otu_table(ps.filtered_2021) <- subset(otu_table(ps.filtered_2021),
                                      select = -c(sample1248,sample1362,sample1436,sample1490))
write.csv(otu_table(ps.filtered_2021), file = "filtering/doublechecked_filtered_low_2021.csv")


# checking read distributions to see how sequencing depth is distributed
  #2018
sample_data(ps.filtered_2018)
sums_2018 <- data.frame(colSums(otu_table(ps.filtered_2018)))
colnames(sums_2018) <- "Sample_totalSeqs_2018"
sums_2018$Sample <- row.names(sums_2018)
sums_2018
ggplot(sums_2018, aes(x=Sample_totalSeqs_2018)) +
  geom_histogram(binwidth=500, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(Sample_totalSeqs_2018, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", linewidth=1)
  #2019
sample_data(ps.filtered_2019)
sums_2019 <- data.frame(colSums(otu_table(ps.filtered_2019)))
colnames(sums_2019) <- "Sample_totalSeqs_2019"
sums_2019$Sample <- row.names(sums_2019)
sums_2019
ggplot(sums_2019, aes(x=Sample_totalSeqs_2019)) +
  geom_histogram(binwidth=500, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(Sample_totalSeqs_2019, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", linewidth=1)
  #2020
sample_data(ps.filtered_2020)
sums_2020 <- data.frame(colSums(otu_table(ps.filtered_2020)))
colnames(sums_2020) <- "Sample_totalSeqs_2020"
sums_2020$Sample <- row.names(sums_2020)
sums_2020
ggplot(sums_2020, aes(x=Sample_totalSeqs_2020)) +
  geom_histogram(binwidth=500, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(Sample_totalSeqs_2020, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", linewidth=1)
  #2021
sample_data(ps.filtered_2021)
sums_2021 <- data.frame(colSums(otu_table(ps.filtered_2021)))
colnames(sums_2021) <- "Sample_totalSeqs_2021"
sums_2021$Sample <- row.names(sums_2021)
sums_2021
ggplot(sums_2021, aes(x=Sample_totalSeqs_2021)) +
  geom_histogram(binwidth=500, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(Sample_totalSeqs_2021, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", linewidth=1)

#Rarefaction Curves
  #2018
otu_rare_2018 <- as.data.frame(otu_table(ps.filtered_2018))
metadata_2018_rare <- as.data.frame(sample_data(ps.filtered_2018))
rarecurve_2018 <- rarecurve(t(otu_rare_2018), col = as.factor(metadata_2018_rare$Origin), label = FALSE, 
                            step = 50,
                            main="Fungi_2018", ylab = "Number of ASVs", xlab = "Number of DNA reads")-> rare_2018_fungi
legend("bottomright", legend=c("leaf","root","stem", "soil"),
       col=c("black", "green", "red","blue"), lty=1, cex=0.8, box.lty=1)
rarecurve_2018
  #2019
otu_rare_2019 <- as.data.frame(otu_table(ps.filtered_2019))
metadata_2019_rare <- as.data.frame(sample_data(ps.filtered_2019))
rarecurve_2019 <- rarecurve(t(otu_rare_2019), col = as.factor(metadata_2019_rare$Origin), label = FALSE, 
                            step = 50,
                            main="Fungi_2019", ylab = "Number of ASVs", xlab = "Number of DNA reads")-> rare_2019_fungi
legend("bottomright", legend=c("stem_leaf","root", "soil"),
       col=c("black", "green", "red"), lty=1, cex=0.8, box.lty=1)
rarecurve_2019
  #2020
otu_rare_2020 <- as.data.frame(otu_table(ps.filtered_2020))
metadata_2020_rare <- as.data.frame(sample_data(ps.filtered_2020))
rarecurve_2020 <- rarecurve(t(otu_rare_2020), col = as.factor(metadata_2020_rare$Origin), label = FALSE, 
                            step = 50,
                            main="Fungi_2020", ylab = "Number of ASVs", xlab = "Number of DNA reads")-> rare_2020_fungi
legend("bottomright", legend=c("leaf","root","stem", "soil","Node"),
       col=c("black", "green", "red","blue","purple"), lty=1, cex=0.8, box.lty=1)
rarecurve_2020
  #2021
otu_rare_2021 <- as.data.frame(otu_table(ps.filtered_2021))
metadata_2021_rare <- as.data.frame(sample_data(ps.filtered_2021))
rarecurve_2021 <- rarecurve(t(otu_rare_2021), col = as.factor(metadata_2021_rare$Origin), label = FALSE, 
                            step = 50,
                            main="Fungi_2021", ylab = "Number of ASVs", xlab = "Number of DNA reads")-> rare_2021_fungi
legend("bottomright", legend=c("leaf","root","stem", "soil"),
       col=c("black", "green", "red","blue"), lty=1, cex=0.8, box.lty=1)
rarecurve_2021


#---------------------------------------------------------------------------------------------------------------------------------------------



#PCOA by Year for origin as a check
#FACTORING AND COLOR MAPPING
management_colors<- c("Conventional" = "#85614F", "No-Till" = "#F3B342", "Organic" = "#427C85")
Management_Order<- c("Conventional", "No-Till", "Organic")
Growth_Stage_Order<- c("C2", "V2", "V5", "C3", "VT", "R2", "R3", "C4", "R6", "R4")
shape_mapping <- scale_shape_manual(values = c("C2" = 8,
                                               "V2" = 8,
                                               "V5" = 8,
                                               "C3" = 5,
                                               "VT" = 5,
                                               "R2" = 5,
                                               "R3" = 5,
                                               "C4" = 16,
                                               "R6" = 16,
                                               "R4" = 16))
#2018
ordinate_2018<- ordinate(ps.filtered_2018, method = "PCoA", distance = "bray")
plot_ord_2018<- plot_ordination(ps.filtered_2018, ordinate_2018, color = "Origin", title = "2018 by Origin")
plot_ord_2018
plot_ord_2018_management<- plot_ordination(ps.filtered_2018, ordinate_2018, color = "Management", title = "2018 by Management")
plot_ord_2018_management
#2019
ordinate_2019<- ordinate(ps.filtered_2019, method = "PCoA", distance = "bray")
plot_ord_2019<- plot_ordination(ps.filtered_2019, ordinate_2019, color = "Origin", title = "2019 by Origin")
plot_ord_2019
plot_ord_2019_management<- plot_ordination(ps.filtered_2019, ordinate_2019, color = "Management", title = "2019 by Management")
plot_ord_2019_management
#2020
ordinate_2020<- ordinate(ps.filtered_2020, method = "PCoA", distance = "bray")
plot_ord_2020<- plot_ordination(ps.filtered_2020, ordinate_2020, color = "Origin", title = "2020 by Origin")
plot_ord_2020
plot_ord_2020_management<- plot_ordination(ps.filtered_2020, ordinate_2020, color = "Management", title = "2020 by Management")
plot_ord_2020_management
#2021
ordinate_2021<- ordinate(ps.filtered_2021, method = "PCoA", distance = "bray")
plot_ord_2021<- plot_ordination(ps.filtered_2021, ordinate_2021, color = "Origin", title = "2021 by Origin")
plot_ord_2021
plot_ord_2021_management<- plot_ordination(ps.filtered_2021, ordinate_2021, color = "Management", title = "2021 by Management")
plot_ord_2021_management
#COMBINED
PCOA_YearxOrigin<- ggarrange(plot_ord_2018,plot_ord_2019,plot_ord_2020,plot_ord_2021, nrow=2, ncol=2)
PCOA_YearxOrigin
PCOA_YearxManagement<-ggarrange(plot_ord_2018_management,plot_ord_2019_management,plot_ord_2020_management,plot_ord_2021_management, nrow=2, ncol=2)
PCOA_YearxManagement
PCOA_YearxOrigin_YearxManagement<- ggarrange(plot_ord_2018, plot_ord_2018_management,
                                             plot_ord_2019, plot_ord_2019_management,
                                             plot_ord_2020, plot_ord_2020_management,
                                             plot_ord_2021, plot_ord_2021_management,
                                             nrow=4, ncol = 2)
PCOA_YearxOrigin_YearxManagement

#By Compartment (above vs belowground)
shape_mapping <- scale_shape_manual(values = c("C2" = 8,
                                               "V2" = 8,
                                               "V5" = 8,
                                               "C3" = 5,
                                               "VT" = 5,
                                               "R2" = 5,
                                               "R3" = 5,
                                               "C4" = 16,
                                               "R6" = 16,
                                               "R4" = 16))
Growth_Stage_Order<- c("C2", "V2", "V5", "C3", "VT", "R2", "R3", "C4", "R6", "R4")
#2018
ps.filtered_2018_df<- as.data.frame(sample_data(ps.filtered_2018))
ps.filtered_2018_df$Management<- factor(ps.filtered_2018_df$Management,levels = Management_Order)
sample_data(ps.filtered_2018) <-sample_data(ps.filtered_2018_df)
print(levels(sample_data(ps.filtered_2018)$Management))
ps.filtered_2018_df<- as.data.frame(sample_data(ps.filtered_2018))
ps.filtered_2018_df$Growth_stage<- factor(ps.filtered_2018_df$Growth_stage,levels = Growth_Stage_Order)
sample_data(ps.filtered_2018) <-sample_data(ps.filtered_2018_df)
print(levels(sample_data(ps.filtered_2018)$Growth_stage))
ps.filtered_2018_above<- subset_samples(ps.filtered_2018, Compartment == "above-ground")
ps.filtered_2018_below<- subset_samples(ps.filtered_2018, Compartment == "below-ground")
#2019
ps.filtered_2019_df<- as.data.frame(sample_data(ps.filtered_2019))
ps.filtered_2019_df$Management<- factor(ps.filtered_2019_df$Management,levels = Management_Order)
sample_data(ps.filtered_2019) <-sample_data(ps.filtered_2019_df)
print(levels(sample_data(ps.filtered_2019)$Management))
ps.filtered_2019_df<- as.data.frame(sample_data(ps.filtered_2019))
ps.filtered_2019_df$Growth_stage<- factor(ps.filtered_2019_df$Growth_stage,levels = Growth_Stage_Order)
sample_data(ps.filtered_2019) <-sample_data(ps.filtered_2019_df)
print(levels(sample_data(ps.filtered_2019)$Growth_stage))
ps.filtered_2019_above<- subset_samples(ps.filtered_2019, Compartment == "above-ground")
ps_2019_above_filtered<- subset_samples(ps.filtered_2019_above, Growth_stage %in% c("C2", "C3", "C4") & Management %in% c("Conventional", "No-Till", "Organic"))
ps.filtered_2019_below<- subset_samples(ps.filtered_2019, Compartment == "below-ground")
ps_2019_below_filtered<- subset_samples(ps.filtered_2019_below, Growth_stage %in% c("C2", "C3", "C4") & Management %in% c("Conventional", "No-Till", "Organic"))
#2020
ps.filtered_2020_df<- as.data.frame(sample_data(ps.filtered_2020))
ps.filtered_2020_df$Management<- factor(ps.filtered_2020_df$Management,levels = Management_Order)
sample_data(ps.filtered_2020) <-sample_data(ps.filtered_2020_df)
print(levels(sample_data(ps.filtered_2020)$Management))
ps.filtered_2020_df<- as.data.frame(sample_data(ps.filtered_2020))
ps.filtered_2020_df$Growth_stage<- factor(ps.filtered_2020_df$Growth_stage,levels = Growth_Stage_Order)
sample_data(ps.filtered_2020) <-sample_data(ps.filtered_2020_df)
print(levels(sample_data(ps.filtered_2020)$Growth_stage))
ps.filtered_2020_above<- subset_samples(ps.filtered_2020, Compartment == "above-ground")
ps.filtered_2020_below<- subset_samples(ps.filtered_2020, Compartment == "below-ground")
#2021
ps.filtered_2021_df<- as.data.frame(sample_data(ps.filtered_2021))
ps.filtered_2021_df$Management<- factor(ps.filtered_2021_df$Management,levels = Management_Order)
sample_data(ps.filtered_2021) <-sample_data(ps.filtered_2021_df)
print(levels(sample_data(ps.filtered_2021)$Management))
ps.filtered_2021_df<- as.data.frame(sample_data(ps.filtered_2021))
ps.filtered_2021_df$Growth_stage<- factor(ps.filtered_2021_df$Growth_stage,levels = Growth_Stage_Order)
sample_data(ps.filtered_2021) <-sample_data(ps.filtered_2021_df)
print(levels(sample_data(ps.filtered_2021)$Growth_stage))
ps.filtered_2021_above<- subset_samples(ps.filtered_2021, Compartment == "above-ground")
ps.filtered_2021_below<- subset_samples(ps.filtered_2021, Compartment == "below-ground")

#PCOA Compartment x Management
#2018
above_ord_2018<- ordinate(ps.filtered_2018_above, method= "PCoA", distance = "bray")
below_ord_2018<- ordinate(ps.filtered_2018_below, method= "PCoA", distance = "bray")
plot_above_ord_2018<- plot_ordination(ps.filtered_2018_above, above_ord_2018,
                                      color="Management", shape = "Growth_stage",
                                      title = "Above Ground") + theme_classic() + geom_point(size = 3, stroke = 1) +
                                      shape_mapping + theme(plot.title = element_text(hjust = 0.5)) +
                                      scale_color_manual(values = management_colors) + stat_ellipse(aes(group = Management, color = Management), level = 0.95, size = 1)
plot_above_ord_2018
plot_below_ord_2018<- plot_ordination(ps.filtered_2018_below, below_ord_2018,
                                      color="Management", shape = "Growth_stage",
                                      title = "Below-Ground") + theme_classic() + geom_point(size = 3) +
                                      shape_mapping + theme(plot.title = element_text(hjust = 0.5)) +
                                      scale_color_manual(values = management_colors) + stat_ellipse(aes(group = Management, color = Management), level = 0.95, size = 1)
                                      plot_below_ord_2018
compartment_ord_2018 <- ggarrange(plot_above_ord_2018, plot_below_ord_2018,
                                                                        nrow = 1, ncol = 2,
                                                                        common.legend = TRUE, legend = "top")
compartment_ord_2018
#2019 (below still has a control?)
above_ord_2019<- ordinate(ps_2019_above_filtered, method= "PCoA", distance = "bray")
below_ord_2019<- ordinate(ps_2019_below_filtered, method= "PCoA", distance = "bray")
plot_above_ord_2019<- plot_ordination(ps.filtered_2019_above, above_ord_2019,
                                      color="Management", shape = "Growth_stage",
                                      title = "Above Ground") + theme_classic() + geom_point(size = 3, stroke = 1) +
                                      shape_mapping + theme(plot.title = element_text(hjust = 0.5),
                                                            legend.position = "none") +  # Remove the legend
                                      scale_color_manual(values = management_colors) + stat_ellipse(aes(group = Management, color = Management), level = 0.95, size = 1)
plot_above_ord_2019
plot_below_ord_2019<- plot_ordination(ps.filtered_2019_below, below_ord_2019,
                                      color="Management", shape = "Growth_stage",
                                      title = "Below-Ground") + theme_classic() + geom_point(size = 3) +
                                      shape_mapping + theme(plot.title = element_text(hjust = 0.5),
                                                            legend.position = "none") +  # Remove the legend
                                      scale_color_manual(values = management_colors) + stat_ellipse(aes(group = Management, color = Management), level = 0.95, size = 1)
plot_below_ord_2019
compartment_ord_2019<-ggarrange(plot_above_ord_2019,plot_below_ord_2019, nrow=1, ncol=2)
compartment_ord_2019
#2020
above_ord_2020<- ordinate(ps.filtered_2020_above, method= "PCoA", distance = "bray")
below_ord_2020<- ordinate(ps.filtered_2020_below, method= "PCoA", distance = "bray")
plot_above_ord_2020<- plot_ordination(ps.filtered_2020_above, above_ord_2020,
                                      color="Management", shape = "Growth_stage",
                                      title = "Above-Ground") + theme_classic() + geom_point(size = 3, stroke = 1) +
                                      shape_mapping + theme(plot.title = element_text(hjust = 0.5),
                                                            legend.position = "none") +  # Remove the legend
                                      scale_color_manual(values = management_colors) + stat_ellipse(aes(group = Management, color = Management), level = 0.95, size = 1)
plot_above_ord_2020
plot_below_ord_2020<- plot_ordination(ps.filtered_2020_below, below_ord_2020,
                                      color="Management", shape = "Growth_stage",
                                      title = "Below-Ground") + theme_classic() + geom_point(size = 3) +
                                      shape_mapping + theme(plot.title = element_text(hjust = 0.5),
                                                            legend.position = "none") +  # Remove the legend
                                      scale_color_manual(values = management_colors) + stat_ellipse(aes(group = Management, color = Management), level = 0.95, size = 1)
plot_below_ord_2020
compartment_ord_2020<-ggarrange(plot_above_ord_2020,plot_below_ord_2020, nrow=1, ncol=2)
compartment_ord_2020
#2021
above_ord_2021<- ordinate(ps.filtered_2021_above, method= "PCoA", distance = "bray")
below_ord_2021<- ordinate(ps.filtered_2021_below, method= "PCoA", distance = "bray")
plot_above_ord_2021<- plot_ordination(ps.filtered_2021_above, above_ord_2021,
                                      color="Management", shape = "Growth_stage",
                                      title = "Above-Ground") + theme_classic() + geom_point(size = 3, stroke = 1) +
                                      shape_mapping + theme(plot.title = element_text(hjust = 0.5),
                                                            legend.position = "none") +  # Remove the legend
                                      scale_color_manual(values = management_colors) + stat_ellipse(aes(group = Management, color = Management), level = 0.95, size = 1)
plot_above_ord_2021
plot_below_ord_2021<- plot_ordination(ps.filtered_2021_below, below_ord_2021,
                                      color="Management", shape = "Growth_stage",
                                      title = "Below-Ground") + theme_classic() + geom_point(size = 3) +
                                      shape_mapping + theme(plot.title = element_text(hjust = 0.5),
                                                            legend.position = "none") +  # Remove the legend
                                      scale_color_manual(values = management_colors) + stat_ellipse(aes(group = Management, color = Management), level = 0.95, size = 1)
plot_below_ord_2021
compartment_ord_2021<-ggarrange(plot_above_ord_2021,plot_below_ord_2021, nrow=1, ncol=2)
compartment_ord_2021
#Combine
above_below_xManagement_ord<- ggarrange(compartment_ord_2018,compartment_ord_2019,compartment_ord_2020,compartment_ord_2021, nrow = 4, ncol=1, labels = c("A) 2018", "B) 2019", "C) 2020", "D) 2021"))
above_below_xManagement_ord

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#PERMANOVA AND BETADISPERSION
  #2018
    #above ground
perma_18_above_otu<- as.data.frame(otu_table(ps.filtered_2018_above))
perma_18_above_taxa<- as.data.frame(as.matrix(tax_table(ps.filtered_2018_above)))
perma_18_above_meta<- as.data.frame(as.matrix(sample_data(ps.filtered_2018_above)))
model.matrix(~ Growth_stage * Management, data=perma_18_above_meta)
model.matrix(~ Growth_stage + Management + Growth_stage : Management, data=perma_18_above_meta)
adonis(t(perma_18_above_otu) ~ Growth_stage * Management, data=perma_18_above_meta, permutations=9999) # by = "margin"
above_18_perm<-adonis(t(perma_18_above_otu) ~ Growth_stage + Management + Growth_stage : Management, data=perma_18_above_meta, permutations=9999)
adonis(t(perma_18_above_otu) ~ Growth_stage + Management + Growth_stage : Management, data=perma_18_above_meta, permutations=9999)
vegan::vegdist(t(perma_18_above_otu), method="bray") -> dist_18_above
permdisp_otu_18_above_Man <- betadisper(dist_18_above, perma_18_above_meta$Management)
permdisp_otu_18_above_GS<- betadisper(dist_18_above, perma_18_above_meta$Growth_stage)
anova(permdisp_otu_18_above_Man, permutations = 9999)
permutest(permdisp_otu_18_above_Man, permutations = 9999, pairwise = T)
plot(permdisp_otu_18_above_Man)
plot(TukeyHSD(permdisp_otu_18_above_Man), las=1)
boxplot(permdisp_otu_18_above_Man)
anova(permdisp_otu_18_above_GS, permutations = 9999)
permutest(permdisp_otu_18_above_GS, permutations = 9999, pairwise = T)
plot(permdisp_otu_18_above_GS)
plot(TukeyHSD(permdisp_otu_18_above_GS), las=1)
boxplot(permdisp_otu_18_above_GS)
    #below ground
perma_18_below_otu<- as.data.frame(otu_table(ps.filtered_2018_below))
perma_18_below_taxa<- as.data.frame(as.matrix(tax_table(ps.filtered_2018_below)))
perma_18_below_meta<- as.data.frame(as.matrix(sample_data(ps.filtered_2018_below)))
model.matrix(~ Growth_stage * Management, data=perma_18_below_meta)
model.matrix(~ Growth_stage + Management + Growth_stage : Management, data=perma_18_below_meta)
adonis(t(perma_18_below_otu) ~ Growth_stage * Management, data=perma_18_below_meta, permutations=9999) # by = "margin"
below_18_perm<-adonis(t(perma_18_below_otu) ~ Growth_stage + Management + Growth_stage : Management, data=perma_18_below_meta, permutations=9999)
adonis(t(perma_18_below_otu) ~ Growth_stage + Management + Growth_stage : Management, data=perma_18_below_meta, permutations=9999)
vegan::vegdist(t(perma_18_below_otu), method="bray") -> dist_18below
permdisp_otu_18_below_Man <- betadisper(dist_18below, perma_18_below_meta$Management)
permdisp_otu_18_below_GS<- betadisper(dist_18below, perma_18_below_meta$Growth_stage)
anova(permdisp_otu_18_below_Man, permutations = 9999)
permutest(permdisp_otu_18_below_Man, permutations = 9999, pairwise = T)
plot(permdisp_otu_18_below_Man)
plot(TukeyHSD(permdisp_otu_18_below_Man), las=1)
boxplot(permdisp_otu_18_below_Man)
anova(permdisp_otu_18_below_GS, permutations = 9999)
permutest(permdisp_otu_18_below_GS, permutations = 9999, pairwise = T)
plot(permdisp_otu_18_below_GS)
plot(TukeyHSD(permdisp_otu_18_below_GS), las=1)
boxplot(permdisp_otu_18_below_GS)

  #2019
    #above ground
perma_19_above_otu<- as.data.frame(otu_table(ps.filtered_2019_above))
perma_19_above_taxa<- as.data.frame(as.matrix(tax_table(ps.filtered_2019_above)))
perma_19_above_meta<- as.data.frame(as.matrix(sample_data(ps.filtered_2019_above)))
model.matrix(~ Growth_stage * Management, data=perma_19_above_meta)
model.matrix(~ Growth_stage + Management + Growth_stage : Management, data=perma_19_above_meta)
adonis(t(perma_19_above_otu) ~ Growth_stage * Management, data=perma_19_above_meta, permutations=9999) # by = "margin"
above_19_perm<-adonis(t(perma_19_above_otu) ~ Growth_stage + Management + Growth_stage : Management, data=perma_19_above_meta, permutations=9999)
adonis(t(perma_19_above_otu) ~ Growth_stage + Management + Growth_stage : Management, data=perma_19_above_meta, permutations=9999)
vegan::vegdist(t(perma_19_above_otu), method="bray") -> dist_19_above
permdisp_otu_19_above_Man <- betadisper(dist_19_above, perma_19_above_meta$Management)
permdisp_otu_19_above_GS<- betadisper(dist_19_above, perma_19_above_meta$Growth_stage)
anova(permdisp_otu_19_above_Man, permutations = 9999)
permutest(permdisp_otu_19_above_Man, permutations = 9999, pairwise = T)
plot(permdisp_otu_19_above_Man)
plot(TukeyHSD(permdisp_otu_19_above_Man), las=1)
boxplot(permdisp_otu_19_above_Man)
anova(permdisp_otu_19_above_GS, permutations = 9999)
permutest(permdisp_otu_19_above_GS, permutations = 9999, pairwise = T)
plot(permdisp_otu_19_above_GS)
plot(TukeyHSD(permdisp_otu_19_above_GS), las=1)
boxplot(permdisp_otu_19_above_GS)
    #below ground
perma_19_below_otu <- as.data.frame(otu_table(ps.filtered_2019_below))
perma_19_below_taxa <- as.data.frame(as.matrix(tax_table(ps.filtered_2019_below)))
perma_19_below_meta <- as.data.frame(as.matrix(sample_data(ps.filtered_2019_below)))
# Filter metadata to exclude NA values in Growth_stage and Management
perma_19_below_meta <- perma_19_below_meta %>%
  filter(!is.na(Growth_stage) & !is.na(Management))
# Ensure OTU table only contains samples in metadata
common_samples <- intersect(rownames(perma_19_below_meta), colnames(perma_19_below_otu))
perma_19_below_otu <- perma_19_below_otu[, common_samples]
perma_19_below_meta <- perma_19_below_meta[common_samples, ]
model.matrix(~ Growth_stage * Management, data=perma_19_below_meta)
model.matrix(~ Growth_stage + Management + Growth_stage : Management, data=perma_19_below_meta)
adonis(t(perma_19_below_otu) ~ Growth_stage * Management, data=perma_19_below_meta, permutations=9999) # by = "margin"
below_19_perm<-adonis(t(perma_19_below_otu) ~ Growth_stage + Management + Growth_stage : Management, data=perma_19_below_meta, permutations=9999)
adonis(t(perma_19_below_otu) ~ Growth_stage + Management + Growth_stage : Management, data=perma_19_below_meta, permutations=9999)
vegan::vegdist(t(perma_19_below_otu), method="bray") -> dist_19below
permdisp_otu_19_below_Man <- betadisper(dist_19below, perma_19_below_meta$Management)
permdisp_otu_19_below_GS<- betadisper(dist_19below, perma_19_below_meta$Growth_stage)
anova(permdisp_otu_19_below_Man, permutations = 9999)
permutest(permdisp_otu_19_below_Man, permutations = 9999, pairwise = T)
plot(permdisp_otu_19_below_Man)
plot(TukeyHSD(permdisp_otu_19_below_Man), las=1)
boxplot(permdisp_otu_19_below_Man)
anova(permdisp_otu_19_below_GS, permutations = 9999)
permutest(permdisp_otu_19_below_GS, permutations = 9999, pairwise = T)
plot(permdisp_otu_19_below_GS)
plot(TukeyHSD(permdisp_otu_19_below_GS), las=1)
boxplot(permdisp_otu_19_below_GS)

  #2020
    #above ground
perma_20_above_otu<- as.data.frame(otu_table(ps.filtered_2020_above))
perma_20_above_taxa<- as.data.frame(as.matrix(tax_table(ps.filtered_2020_above)))
perma_20_above_meta<- as.data.frame(as.matrix(sample_data(ps.filtered_2020_above)))
model.matrix(~ Growth_stage * Management, data=perma_20_above_meta)
model.matrix(~ Growth_stage + Management + Growth_stage : Management, data=perma_20_above_meta)
adonis(t(perma_20_above_otu) ~ Growth_stage * Management, data=perma_20_above_meta, permutations=9999) # by = "margin"
above_20_perm<-adonis(t(perma_20_above_otu) ~ Growth_stage + Management + Growth_stage : Management, data=perma_20_above_meta, permutations=9999)
adonis(t(perma_20_above_otu) ~ Growth_stage + Management + Growth_stage : Management, data=perma_20_above_meta, permutations=9999)
vegan::vegdist(t(perma_20_above_otu), method="bray") -> dist_20_above
permdisp_otu_20_above_Man <- betadisper(dist_20_above, perma_20_above_meta$Management)
permdisp_otu_20_above_GS<- betadisper(dist_20_above, perma_20_above_meta$Growth_stage)
anova(permdisp_otu_20_above_Man, permutations = 9999)
permutest(permdisp_otu_20_above_Man, permutations = 9999, pairwise = T)
plot(permdisp_otu_20_above_Man)
plot(TukeyHSD(permdisp_otu_20_above_Man), las=1)
boxplot(permdisp_otu_20_above_Man)
anova(permdisp_otu_20_above_GS, permutations = 9999)
permutest(permdisp_otu_20_above_GS, permutations = 9999, pairwise = T)
plot(permdisp_otu_20_above_GS)
plot(TukeyHSD(permdisp_otu_20_above_GS), las=1)
boxplot(permdisp_otu_20_above_GS)
    #below ground
perma_20_below_otu<- as.data.frame(otu_table(ps.filtered_2020_below))
perma_20_below_taxa<- as.data.frame(as.matrix(tax_table(ps.filtered_2020_below)))
perma_20_below_meta<- as.data.frame(as.matrix(sample_data(ps.filtered_2020_below)))
model.matrix(~ Growth_stage * Management, data=perma_20_below_meta)
model.matrix(~ Growth_stage + Management + Growth_stage : Management, data=perma_20_below_meta)
adonis(t(perma_20_below_otu) ~ Growth_stage * Management, data=perma_20_below_meta, permutations=9999) # by = "margin"
below_20_perm<-adonis(t(perma_20_below_otu) ~ Growth_stage + Management + Growth_stage : Management, data=perma_20_below_meta, permutations=9999)
adonis(t(perma_20_below_otu) ~ Growth_stage + Management + Growth_stage : Management, data=perma_20_below_meta, permutations=9999)
vegan::vegdist(t(perma_20_below_otu), method="bray") -> dist_20below
permdisp_otu_20_below_Man <- betadisper(dist_20below, perma_20_below_meta$Management)
permdisp_otu_20_below_GS<- betadisper(dist_20below, perma_20_below_meta$Growth_stage)
anova(permdisp_otu_20_below_Man, permutations = 9999)
permutest(permdisp_otu_20_below_Man, permutations = 9999, pairwise = T)
plot(permdisp_otu_20_below_Man)
plot(TukeyHSD(permdisp_otu_20_below_Man), las=1)
boxplot(permdisp_otu_20_below_Man)
anova(permdisp_otu_20_below_GS, permutations = 9999)
permutest(permdisp_otu_20_below_GS, permutations = 9999, pairwise = T)
plot(permdisp_otu_20_below_GS)
plot(TukeyHSD(permdisp_otu_20_below_GS), las=1)
boxplot(permdisp_otu_20_below_GS)

  #2021
#above ground
perma_21_above_otu<- as.data.frame(otu_table(ps.filtered_2021_above))
perma_21_above_taxa<- as.data.frame(as.matrix(tax_table(ps.filtered_2021_above)))
perma_21_above_meta<- as.data.frame(as.matrix(sample_data(ps.filtered_2021_above)))
model.matrix(~ Growth_stage * Management, data=perma_21_above_meta)
model.matrix(~ Growth_stage + Management + Growth_stage : Management, data=perma_21_above_meta)
adonis(t(perma_21_above_otu) ~ Growth_stage * Management, data=perma_21_above_meta, permutations=9999) # by = "margin"
above_21_perm<-adonis(t(perma_21_above_otu) ~ Growth_stage + Management + Growth_stage : Management, data=perma_21_above_meta, permutations=9999)
adonis(t(perma_21_above_otu) ~ Growth_stage + Management + Growth_stage : Management, data=perma_21_above_meta, permutations=9999)
vegan::vegdist(t(perma_21_above_otu), method="bray") -> dist_21_above
permdisp_otu_21_above_Man <- betadisper(dist_21_above, perma_21_above_meta$Management)
permdisp_otu_21_above_GS<- betadisper(dist_21_above, perma_21_above_meta$Growth_stage)
anova(permdisp_otu_21_above_Man, permutations = 9999)
permutest(permdisp_otu_21_above_Man, permutations = 9999, pairwise = T)
plot(permdisp_otu_21_above_Man)
plot(TukeyHSD(permdisp_otu_21_above_Man), las=1)
boxplot(permdisp_otu_21_above_Man)
anova(permdisp_otu_21_above_GS, permutations = 9999)
permutest(permdisp_otu_21_above_GS, permutations = 9999, pairwise = T)
plot(permdisp_otu_21_above_GS)
plot(TukeyHSD(permdisp_otu_21_above_GS), las=1)
boxplot(permdisp_otu_21_above_GS)
#below ground
perma_21_below_otu<- as.data.frame(otu_table(ps.filtered_2021_below))
perma_21_below_taxa<- as.data.frame(as.matrix(tax_table(ps.filtered_2021_below)))
perma_21_below_meta<- as.data.frame(as.matrix(sample_data(ps.filtered_2021_below)))
model.matrix(~ Growth_stage * Management, data=perma_21_below_meta)
model.matrix(~ Growth_stage + Management + Growth_stage : Management, data=perma_21_below_meta)
adonis(t(perma_21_below_otu) ~ Growth_stage * Management, data=perma_21_below_meta, permutations=9999) # by = "margin"
below_21_perm<-adonis(t(perma_21_below_otu) ~ Growth_stage + Management + Growth_stage : Management, data=perma_21_below_meta, permutations=9999)
adonis(t(perma_21_below_otu) ~ Growth_stage + Management + Growth_stage : Management, data=perma_21_below_meta, permutations=9999)
vegan::vegdist(t(perma_21_below_otu), method="bray") -> dist_21below
permdisp_otu_21_below_Man <- betadisper(dist_21below, perma_21_below_meta$Management)
permdisp_otu_21_below_GS<- betadisper(dist_21below, perma_21_below_meta$Growth_stage)
anova(permdisp_otu_21_below_Man, permutations = 9999)
permutest(permdisp_otu_21_below_Man, permutations = 9999, pairwise = T)
plot(permdisp_otu_21_below_Man)
plot(TukeyHSD(permdisp_otu_21_below_Man), las=1)
boxplot(permdisp_otu_21_below_Man)
anova(permdisp_otu_21_below_GS, permutations = 9999)
permutest(permdisp_otu_21_below_GS, permutations = 9999, pairwise = T)
plot(permdisp_otu_21_below_GS)
plot(TukeyHSD(permdisp_otu_21_below_GS), las=1)
boxplot(permdisp_otu_21_below_GS)

##
data <- read.table("fullcroprotationperma.txt", header = TRUE, sep = "\t")
unique_year_compartment <- unique(data$Year.Compartment)
split_data <- split(data, data$Year.Compartment)
for (name in unique_year_compartment) {
  assign(gsub(" ", "_", name), split_data[[name]])
}

Variable_colors<- c("Growth Stage"="#6DA45C", "Management"="#BCE069","Growth Stage : Management"="#A8AE0F", "Residuals"="#427C85")
#2018
  #above
`2018_Above-Ground`$Variable <- factor(`2018_Below-Ground`$Variable, 
                                       levels = c("Growth Stage", "Management", "Growth Stage : Management", "Residuals"))
above_18_perm_plot <- ggplot(`2018_Above-Ground`, aes(x = Variable, y = R.2, fill = Variable)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = R.2), vjust = -0.5, color = "black") +
  theme_minimal() +
  labs(title = "Above-Ground Permanova",
       x = NULL,
       y = "R^2") +
  scale_fill_manual(values = Variable_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.title.y = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(color = "grey75"),
        panel.grid.minor = element_line(color = "white"),
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))
above_18_perm_plot
    #below
`2018_Below-Ground`$Variable <- factor(`2018_Below-Ground`$Variable, 
                                       levels = c("Growth Stage", "Management", "Growth Stage : Management", "Residuals"))
below_18_perm_plot <- ggplot(`2018_Below-Ground`, aes(x = Variable, y = R.2, fill = Variable)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = R.2), vjust = -0.5, color = "black") +
  theme_minimal() +
  labs(title = "Below-Ground Permanova",
       x = NULL,
       y = "R^2") +
  scale_fill_manual(values = Variable_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.title.y = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(color = "grey75"),
        panel.grid.minor = element_line(color = "white"),
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))
below_18_perm_plot
    #combined
 permbetadisp_18<- ggarrange(above_18_perm_plot, below_18_perm_plot, ncol=2, nrow=1)
permbetadisp_18 
#2019
#above
`2019_Above-Ground`$Variable <- factor(`2019_Above-Ground`$Variable, 
                                       levels = c("Growth Stage", "Management", "Growth Stage : Management", "Residuals"))
above_19_perm_plot <- ggplot(`2019_Above-Ground`, aes(x = Variable, y = R.2, fill = Variable)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = R.2), vjust = -0.5, color = "black") +
  theme_minimal() +
  labs(title = "Above-Ground Permanova",
       x = NULL,
       y = "R^2") +
  scale_fill_manual(values = Variable_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.title.y = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(color = "grey75"),
        panel.grid.minor = element_line(color = "white"),
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))
above_19_perm_plot
#below
`2019_Below-Ground`$Variable <- factor(`2019_Below-Ground`$Variable, 
                                       levels = c("Growth Stage", "Management", "Growth Stage : Management", "Residuals"))
below_19_perm_plot <- ggplot(`2019_Below-Ground`, aes(x = Variable, y = R.2, fill = Variable)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = R.2), vjust = -0.5, color = "black") +
  theme_minimal() +
  labs(title = "Below-Ground Permanova",
       x = NULL,
       y = "R^2") +
  scale_fill_manual(values = Variable_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.title.y = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(color = "grey75"),
        panel.grid.minor = element_line(color = "white"),
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))
below_19_perm_plot
#combined
permbetadisp_19<- ggarrange(above_19_perm_plot, below_19_perm_plot, ncol=2, nrow=1)
permbetadisp_19 
#2020
#above
`2020_Above-Ground`$Variable <- factor(`2020_Above-Ground`$Variable, 
                                       levels = c("Growth Stage", "Management", "Growth Stage : Management", "Residuals"))
above_20_perm_plot <- ggplot(`2020_Above-Ground`, aes(x = Variable, y = R.2, fill = Variable)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = R.2), vjust = -0.5, color = "black") +
  theme_minimal() +
  labs(title = "Above-Ground Permanova",
       x = NULL,
       y = "R^2") +
  scale_fill_manual(values = Variable_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.title.y = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(color = "grey75"),
        panel.grid.minor = element_line(color = "white"),
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))
above_20_perm_plot
#below
`2020_Below-Ground`$Variable <- factor(`2020_Below-Ground`$Variable, 
                                       levels = c("Growth Stage", "Management", "Growth Stage : Management", "Residuals"))
below_20_perm_plot <- ggplot(`2020_Below-Ground`, aes(x = Variable, y = R.2, fill = Variable)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = R.2), vjust = -0.5, color = "black") +
  theme_minimal() +
  labs(title = "Below-Ground Permanova",
       x = NULL,
       y = "R^2") +
  scale_fill_manual(values = Variable_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.title.y = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(color = "grey75"),
        panel.grid.minor = element_line(color = "white"),
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))
below_20_perm_plot
#combined
permbetadisp_20<- ggarrange(above_20_perm_plot, below_20_perm_plot, ncol=2, nrow=1)
permbetadisp_20 
#2021
#above
`2021_Above-Ground`$Variable <- factor(`2021_Above-Ground`$Variable, 
                                       levels = c("Growth Stage", "Management", "Growth Stage : Management", "Residuals"))
above_21_perm_plot <- ggplot(`2021_Above-Ground`, aes(x = Variable, y = R.2, fill = Variable)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = R.2), vjust = -0.5, color = "black") +
  theme_minimal() +
  labs(title = "Above-Ground Permanova",
       x = NULL,
       y = "R^2") +
  scale_fill_manual(values = Variable_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.title.y = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(color = "grey75"),
        panel.grid.minor = element_line(color = "white"),
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))
above_21_perm_plot
#below
`2021_Below-Ground`$Variable <- factor(`2021_Below-Ground`$Variable, 
                                       levels = c("Growth Stage", "Management", "Growth Stage : Management", "Residuals"))
below_21_perm_plot <- ggplot(`2021_Below-Ground`, aes(x = Variable, y = R.2, fill = Variable)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = R.2), vjust = -0.5, color = "black") +
  theme_minimal() +
  labs(title = "Below-Ground Permanova",
       x = NULL,
       y = "R^2") +
  scale_fill_manual(values = Variable_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.title.y = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(color = "grey75"),
        panel.grid.minor = element_line(color = "white"),
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))
below_21_perm_plot
#combined
permbetadisp_21<- ggarrange(above_21_perm_plot, below_21_perm_plot, ncol=2, nrow=1)
permbetadisp_21 

#------------------------------Figure1-----------------------------------
pcoa_perma_18<- ggarrange(compartment_ord_2018, permbetadisp_18, ncol=2,nrow=1,
                          labels = c("a) 2018", "b"))
pcoa_perma_18
pcoa_perma_19<- ggarrange(compartment_ord_2019, permbetadisp_19, ncol=2,nrow=1,
                          labels = c("a) 2019", "b"))
pcoa_perma_19
pcoa_perma_20<- ggarrange(compartment_ord_2020, permbetadisp_20, ncol=2,nrow=1,
                          labels = c("a) 2020", "b"))
pcoa_perma_20
pcoa_perma_21<- ggarrange(compartment_ord_2021, permbetadisp_21, ncol=2,nrow=1,
                          labels = c("a) 2021", "b"))
pcoa_perma_21

pcoa_perma_all<- ggarrange(pcoa_perma_18,
                           pcoa_perma_19,
                           pcoa_perma_20,
                           pcoa_perma_21,
                           ncol=1, nrow = 4,
                           labels = c("1", "2", "3","4"))
pcoa_perma_all

#Indicator species
  #2018 above
otu_table_2018_above <- as.data.frame(otu_table(ps.filtered_2018_above))
tax_table_2018_above <- as.data.frame(as.matrix(tax_table(ps.filtered_2018_above)))
sample_data_2018_above <- as.data.frame(as.matrix(sample_data(ps.filtered_2018_above)))
sample_data_2018_above <- sample_data_2018_above[!is.na(sample_data_2018_above$Management), ]
otu_table_2018_above <- otu_table_2018_above[ , rownames(sample_data_2018_above)]
isa_results_2018_above <- multipatt(as.data.frame(t(otu_table_2018_above)), sample_data_2018_above$Management, control=how(nperm=9999))
isa_results_2018_fdr_above <- isa_results_2018_above
isa_results_2018_fdr_above$sign$p.value <- p.adjust(isa_results_2018_fdr_above$sign$p.value, "fdr")
significant_results_2018_above <- isa_results_2018_fdr_above$sign[isa_results_2018_fdr_above$sign$p.value <= 0.04, ]
ps_significant_2018_above <- prune_taxa(rownames(significant_results_2018_above), ps.filtered_2018_above)
ps_significant_2018_above <- transform_sample_counts(ps_significant_2018_above, function(x) 100 * x / sum(x))
otu_significant_2018_above <- as.data.frame(otu_table(ps_significant_2018_above))
metadata_significant_2018_above <- as.data.frame(sample_data(ps_significant_2018_above))
sample_data_2018_above <- sample_data_2018_above[!is.na(sample_data_2018_above$Management), ]
otu_table_2018_above <- otu_table_2018_above[, rownames(sample_data_2018_above)]
print(head(sample_data_2018_above))
print(head(otu_table_2018_above))
agg_otu_2018_above <- otu_significant_2018_above %>%
  rownames_to_column("OTU") %>%
  gather(Sample, Abundance, -OTU) %>%
  left_join(metadata_significant_2018_above %>% rownames_to_column("Sample"), by = "Sample") %>%
  filter(!is.na(Management)) %>%  # Filter out any NA management
  group_by(Management, OTU) %>%
  summarise(Abundance = mean(Abundance), .groups = 'drop') %>%
  spread(Management, Abundance, fill = 0) %>%
  column_to_rownames("OTU")
significant_results_2018_above <- significant_results_2018_above[rownames(significant_results_2018_above) %in% rownames(agg_otu_2018_above), ]
isa_obj_2018_above <- cbind(agg_otu_2018_above[rownames(significant_results_2018_above), ], significant_results_2018_above)
isa_obj_2018_above$readNo <- rowSums(otu_table_2018_above[rownames(significant_results_2018_above), ])
isa_obj_2018_above$relAb <- (isa_obj_2018_above$readNo / sum(colSums(otu_table_2018_above))) * 100
isa_obj_2018_above$logAb <- log(isa_obj_2018_above$readNo)
isa_obj_2018_above$sqrtAb <- sqrt(isa_obj_2018_above$readNo)
isa_obj_2018_above <- isa_obj_2018_above[order(isa_obj_2018_above$relAb, decreasing = TRUE), ]
isa_obj_2018_above <- isa_obj_2018_above[1:50, ]
clean_matrix <- function(matrix_data) {
  matrix_data[is.na(matrix_data)] <- 0 # Replace NAs with zeros
  matrix_data[matrix_data < 0] <- 0    # Replace negative values if needed
  return(matrix_data)
}
cleaned_isa_matrix_2018_above <- clean_matrix(as.matrix(sqrt(isa_obj_2018_above[, 1:ncol(agg_otu_2018_above)] * 10)))
get_highest_rank <- function(tax_row) {
  ranks <- c("Genus", "Family", "Order", "Class", "Phylum", "Kingdom")
  highest_rank <- NA
  for (rank in ranks) {
    if (!is.na(tax_row[rank]) && tax_row[rank] != "") {
      highest_rank <- tax_row[rank]
      break
    }
  }
  if (is.na(highest_rank) || highest_rank == "") {
    highest_rank <- "Unclassified"
  }
  return(highest_rank)
}
tax_table_2018_above <- as.data.frame(as.matrix(tax_table(ps.filtered_2018_above)))
tax_table_2018_above$HighestRank <- apply(tax_table_2018_above, 1, get_highest_rank)
print(table(tax_table_2018_above$HighestRank))
otu_to_highest_rank <- tax_table_2018_above %>%
  rownames_to_column("OTU") %>%
  select(OTU, HighestRank) %>%
  column_to_rownames("OTU")
isa_obj_2018_above$HighestRank <- otu_to_highest_rank[rownames(isa_obj_2018_above), "HighestRank"]
rownames(cleaned_isa_matrix_2018_above) <- isa_obj_2018_above$HighestRank
print(head(rownames(cleaned_isa_matrix_2018_above)))
ht_2018_above <- Heatmap(cleaned_isa_matrix_2018_above, 
                         col = colorRamp2(c(0, 5), c("white", "#6DA45C")), 
                         cluster_rows = FALSE, cluster_columns = FALSE, 
                         name = "Abundance",
                         row_names_gp = gpar(fontsize = 8), 
                         column_names_gp = gpar(fontsize = 8),
                         show_heatmap_legend = FALSE)

ha_bar_2018_above <- HeatmapAnnotation("Rel. Abundance" = row_anno_barplot(isa_obj_2018_above$relAb, axis = FALSE, width = unit(2, "mm")), 
                                       which = "row", annotation_width = unit(1.75, "cm"), 
                                       show_annotation_name = TRUE, 
                                       annotation_name_gp = gpar(fontsize = 8), 
                                       annotation_name_offset = unit(.75, "cm"),
                                       annotation_name_rot = c(0))

ISA_HM_2018_above <- ha_bar_2018_above + ht_2018_above
ISA_HM_2018_above
  #2018 below
otu_table_2018_below <- as.data.frame(otu_table(ps.filtered_2018_below))
tax_table_2018_below <- as.data.frame(as.matrix(tax_table(ps.filtered_2018_below)))
sample_data_2018_below <- as.data.frame(as.matrix(sample_data(ps.filtered_2018_below)))
sample_data_2018_below <- sample_data_2018_below[!is.na(sample_data_2018_below$Management), ]
otu_table_2018_below <- otu_table_2018_below[ , rownames(sample_data_2018_below)]
isa_results_2018_below <- multipatt(as.data.frame(t(otu_table_2018_below)), sample_data_2018_below$Management, control=how(nperm=9999))
isa_results_2018_fdr_below <- isa_results_2018_below
isa_results_2018_fdr_below$sign$p.value <- p.adjust(isa_results_2018_fdr_below$sign$p.value, "fdr")
significant_results_2018_below <- isa_results_2018_fdr_below$sign[isa_results_2018_fdr_below$sign$p.value <= 0.04, ]
ps_significant_2018_below <- prune_taxa(rownames(significant_results_2018_below), ps.filtered_2018_below)
ps_significant_2018_below <- transform_sample_counts(ps_significant_2018_below, function(x) 100 * x / sum(x))
otu_significant_2018_below <- as.data.frame(otu_table(ps_significant_2018_below))
metadata_significant_2018_below <- as.data.frame(sample_data(ps_significant_2018_below))
sample_data_2018_below <- sample_data_2018_below[!is.na(sample_data_2018_below$Management), ]
otu_table_2018_below <- otu_table_2018_below[, rownames(sample_data_2018_below)]
print(head(sample_data_2018_below))
print(head(otu_table_2018_below))
agg_otu_2018_below <- otu_significant_2018_below %>%
  rownames_to_column("OTU") %>%
  gather(Sample, Abundance, -OTU) %>%
  left_join(metadata_significant_2018_below %>% rownames_to_column("Sample"), by = "Sample") %>%
  filter(!is.na(Management)) %>%  # Filter out any NA management
  group_by(Management, OTU) %>%
  summarise(Abundance = mean(Abundance), .groups = 'drop') %>%
  spread(Management, Abundance, fill = 0) %>%
  column_to_rownames("OTU")
significant_results_2018_below <- significant_results_2018_below[rownames(significant_results_2018_below) %in% rownames(agg_otu_2018_below), ]
isa_obj_2018_below <- cbind(agg_otu_2018_below[rownames(significant_results_2018_below), ], significant_results_2018_below)
isa_obj_2018_below$readNo <- rowSums(otu_table_2018_below[rownames(significant_results_2018_below), ])
isa_obj_2018_below$relAb <- (isa_obj_2018_below$readNo / sum(colSums(otu_table_2018_below))) * 100
isa_obj_2018_below$logAb <- log(isa_obj_2018_below$readNo)
isa_obj_2018_below$sqrtAb <- sqrt(isa_obj_2018_below$readNo)
isa_obj_2018_below <- isa_obj_2018_below[order(isa_obj_2018_below$relAb, decreasing = TRUE), ]
isa_obj_2018_below <- isa_obj_2018_below[1:50, ]
clean_matrix <- function(matrix_data) {
  matrix_data[is.na(matrix_data)] <- 0 # Replace NAs with zeros
  matrix_data[matrix_data < 0] <- 0    # Replace negative values if needed
  return(matrix_data)
}
cleaned_isa_matrix_2018_below <- clean_matrix(as.matrix(sqrt(isa_obj_2018_below[, 1:ncol(agg_otu_2018_below)] * 10)))
get_highest_rank <- function(tax_row) {
  ranks <- c("Genus", "Family", "Order", "Class", "Phylum", "Kingdom")
  highest_rank <- NA
  for (rank in ranks) {
    if (!is.na(tax_row[rank]) && tax_row[rank] != "") {
      highest_rank <- tax_row[rank]
      break
    }
  }
  if (is.na(highest_rank) || highest_rank == "") {
    highest_rank <- "Unclassified"
  }
  return(highest_rank)
}
tax_table_2018_below <- as.data.frame(as.matrix(tax_table(ps.filtered_2018_below)))
tax_table_2018_below$HighestRank <- apply(tax_table_2018_below, 1, get_highest_rank)
print(table(tax_table_2018_below$HighestRank))
otu_to_highest_rank <- tax_table_2018_below %>%
  rownames_to_column("OTU") %>%
  select(OTU, HighestRank) %>%
  column_to_rownames("OTU")
isa_obj_2018_below$HighestRank <- otu_to_highest_rank[rownames(isa_obj_2018_below), "HighestRank"]
rownames(cleaned_isa_matrix_2018_below) <- isa_obj_2018_below$HighestRank
print(head(rownames(cleaned_isa_matrix_2018_below)))
ht_2018_below <- Heatmap(cleaned_isa_matrix_2018_below, 
                         col = colorRamp2(c(0, 5), c("white", "#6DA45C")), 
                         cluster_rows = FALSE, cluster_columns = FALSE, 
                         name = "Abundance",
                         row_names_gp = gpar(fontsize = 8), 
                         column_names_gp = gpar(fontsize = 8),
                         show_heatmap_legend = FALSE)

ha_bar_2018_below <- HeatmapAnnotation("Rel. Abundance" = row_anno_barplot(isa_obj_2018_below$relAb, axis = FALSE, width = unit(2, "mm")), 
                                       which = "row", annotation_width = unit(1.75, "cm"), 
                                       show_annotation_name = TRUE, 
                                       annotation_name_gp = gpar(fontsize = 8), 
                                       annotation_name_offset = unit(.75, "cm"),
                                       annotation_name_rot = c(0))

ISA_HM_2018_below <- ha_bar_2018_below + ht_2018_below
ISA_HM_2018_below



##RF----------------------------------------------------------
#2018 above
ordination_results_18_above <- ordinate(ps.filtered_2018_above, method = "PCoA", distance = "bray")
pcoa_scores_18_above <- ordination_results_18_above$vectors
response_variable_18_above <- as.numeric(pcoa_scores_18_above[, 1])
sample_data_df_18_above <- as.data.frame(sample_data(ps.filtered_2018_above))
correct_columns <- c(
  "Growth_stage",
  "Management",
  "Precipitation",
  "Ambient.Temperature",
  "Nitrate",
  "Ammonium",
  "pH",
  "Phosphorus",
  "Potassium",
  "Calcium",
  "Magnesium",
  "Plot"
)
sample_data_df_18_above <- sample_data_df_18_above[, correct_columns]
sample_data_df_18_above <- data.frame(lapply(sample_data_df_18_above, function(x) {
  if (is.factor(x) || is.character(x)) {
    return(as.factor(x))
  } else {
    return(as.numeric(as.character(x)))
  }
}))
missing_indices <- which(!complete.cases(sample_data_df_18_above))
if(length(missing_indices) > 0) {
  response_variable_18_above <- response_variable_18_above[-missing_indices]
  sample_data_df_18_above <- sample_data_df_18_above[complete.cases(sample_data_df_18_above), ]
}
if (length(response_variable_18_above) == nrow(sample_data_df_18_above)) {
  rf_model_18_above <- randomForest(response_variable_18_above ~ ., data = sample_data_df_18_above, importance = TRUE, ntree = 99999)
  importance_matrix_18_above <- importance(rf_model_18_above)
  print(importance_matrix_18_above)
  varImpPlot(rf_model_18_above, main = "2018 Above Ground")
} else {
  stop("Mismatch in number of observations between response variable and predictor data.")
}

#2018 below
ordination_results_18_below <- ordinate(ps.filtered_2018_below, method = "PCoA", distance = "bray")
pcoa_scores_18_below <- ordination_results_18_below$vectors
response_variable_18_below <- as.numeric(pcoa_scores_18_below[, 1])
sample_data_df_18_below <- as.data.frame(sample_data(ps.filtered_2018_below))
correct_columns <- c(
  "Growth_stage",
  "Management",
  "Precipitation",
  "Ambient.Temperature",
  "Nitrate",
  "Ammonium",
  "pH",
  "Phosphorus",
  "Potassium",
  "Calcium",
  "Magnesium"
)
sample_data_df_18_below <- sample_data_df_18_below[, correct_columns]
sample_data_df_18_below <- data.frame(lapply(sample_data_df_18_below, function(x) {
  if (is.factor(x) || is.character(x)) {
    return(as.factor(x))
  } else {
    return(as.numeric(as.character(x)))
  }
}))
missing_indices <- which(!complete.cases(sample_data_df_18_below))
if(length(missing_indices) > 0) {
  response_variable_18_below <- response_variable_18_below[-missing_indices]
  sample_data_df_18_below <- sample_data_df_18_below[complete.cases(sample_data_df_18_below), ]
}
if (length(response_variable_18_below) == nrow(sample_data_df_18_below)) {
  rf_model_18_below <- randomForest(response_variable_18_below ~ ., data = sample_data_df_18_below, importance = TRUE, ntree = 99999)
  importance_matrix_18_below <- importance(rf_model_18_below)
  print(importance_matrix_18_below)
  varImpPlot(rf_model_18_below, main = "2018 Below Ground")
} else {
  stop("Mismatch in number of observations between response variable and predictor data.")
}

#2019 above
ordination_results_19_above <- ordinate(ps.filtered_2019_above, method = "PCoA", distance = "bray")
pcoa_scores_19_above <- ordination_results_19_above$vectors
response_variable_19_above <- as.numeric(pcoa_scores_19_above[, 1])
sample_data_df_19_above <- as.data.frame(sample_data(ps.filtered_2019_above))
correct_columns <- c(
  "Growth_stage",
  "Management",
  "Precipitation",
  "Ambient.Temperature",
  "Nitrate",
  "Ammonium",
  "pH",
  "Phosphorus",
  "Potassium",
  "Calcium",
  "Magnesium"
)
sample_data_df_19_above <- sample_data_df_19_above[, correct_columns]
sample_data_df_19_above <- data.frame(lapply(sample_data_df_19_above, function(x) {
  if (is.factor(x) || is.character(x)) {
    return(as.factor(x))
  } else {
    return(as.numeric(as.character(x)))
  }
}))
missing_indices <- which(!complete.cases(sample_data_df_19_above))
if(length(missing_indices) > 0) {
  response_variable_19_above <- response_variable_19_above[-missing_indices]
  sample_data_df_19_above <- sample_data_df_19_above[complete.cases(sample_data_df_19_above), ]
}
if (length(response_variable_19_above) == nrow(sample_data_df_19_above)) {
  rf_model_19_above <- randomForest(response_variable_19_above ~ ., data = sample_data_df_19_above, importance = TRUE, ntree = 99999)
  importance_matrix_19_above <- importance(rf_model_19_above)
  print(importance_matrix_19_above)
  varImpPlot(rf_model_19_above, main = "2019 Above Ground")
} else {
  stop("Mismatch in number of observations between response variable and predictor data.")
}

#2019 below
ordination_results_19_below <- ordinate(ps.filtered_2019_below, method = "PCoA", distance = "bray")
pcoa_scores_19_below <- ordination_results_19_below$vectors
response_variable_19_below <- as.numeric(pcoa_scores_19_below[, 1])
sample_data_df_19_below <- as.data.frame(sample_data(ps.filtered_2019_below))
correct_columns <- c(
  "Growth_stage",
  "Management",
  "Precipitation",
  "Ambient.Temperature",
  "Nitrate",
  "Ammonium",
  "pH",
  "Phosphorus",
  "Potassium",
  "Calcium",
  "Magnesium"
)
sample_data_df_19_below <- sample_data_df_19_below[, correct_columns]
sample_data_df_19_below <- data.frame(lapply(sample_data_df_19_below, function(x) {
  if (is.factor(x) || is.character(x)) {
    return(as.factor(x))
  } else {
    return(as.numeric(as.character(x)))
  }
}))
missing_indices <- which(!complete.cases(sample_data_df_19_below))
if(length(missing_indices) > 0) {
  response_variable_19_below <- response_variable_19_below[-missing_indices]
  sample_data_df_19_below <- sample_data_df_19_below[complete.cases(sample_data_df_19_below), ]
}
if (length(response_variable_19_below) == nrow(sample_data_df_19_below)) {
  rf_model_19_below <- randomForest(response_variable_19_below ~ ., data = sample_data_df_19_below, importance = TRUE, ntree = 99999)
  importance_matrix_19_below <- importance(rf_model_19_below)
  print(importance_matrix_19_below)
  varImpPlot(rf_model_19_below, main = "2019 Below Ground")
} else {
  stop("Mismatch in number of observations between response variable and predictor data.")
}

#2020 above
ordination_results_20_above <- ordinate(ps.filtered_2020_above, method = "PCoA", distance = "bray")
pcoa_scores_20_above <- ordination_results_20_above$vectors
response_variable_20_above <- as.numeric(pcoa_scores_20_above[, 1])
sample_data_df_20_above <- as.data.frame(sample_data(ps.filtered_2020_above))
correct_columns <- c(
  "Growth_stage",
  "Management",
  "Precipitation",
  "Ambient.Temperature",
  "Nitrate",
  "Ammonium",
  "pH",
  "Phosphorus",
  "Potassium",
  "Calcium",
  "Magnesium"
)
sample_data_df_20_above <- sample_data_df_20_above[, correct_columns]
sample_data_df_20_above <- data.frame(lapply(sample_data_df_20_above, function(x) {
  if (is.factor(x) || is.character(x)) {
    return(as.factor(x))
  } else {
    return(as.numeric(as.character(x)))
  }
}))
missing_indices <- which(!complete.cases(sample_data_df_20_above))
if(length(missing_indices) > 0) {
  response_variable_20_above <- response_variable_20_above[-missing_indices]
  sample_data_df_20_above <- sample_data_df_20_above[complete.cases(sample_data_df_20_above), ]
}
if (length(response_variable_20_above) == nrow(sample_data_df_20_above)) {
  rf_model_20_above <- randomForest(response_variable_20_above ~ ., data = sample_data_df_20_above, importance = TRUE, ntree = 99999)
  importance_matrix_20_above <- importance(rf_model_20_above)
  print(importance_matrix_20_above)
  varImpPlot(rf_model_20_above, main = "2020 Above Ground")
} else {
  stop("Mismatch in number of observations between response variable and predictor data.")
}

#2020 below
ordination_results_20_below <- ordinate(ps.filtered_2020_below, method = "PCoA", distance = "bray")
pcoa_scores_20_below <- ordination_results_20_below$vectors
response_variable_20_below <- as.numeric(pcoa_scores_20_below[, 1])
sample_data_df_20_below <- as.data.frame(sample_data(ps.filtered_2020_below))
correct_columns <- c(
  "Growth_stage",
  "Management",
  "Precipitation",
  "Ambient.Temperature",
  "Nitrate",
  "Ammonium",
  "pH",
  "Phosphorus",
  "Potassium",
  "Calcium",
  "Magnesium"
)
sample_data_df_20_below <- sample_data_df_20_below[, correct_columns]
sample_data_df_20_below <- data.frame(lapply(sample_data_df_20_below, function(x) {
  if (is.factor(x) || is.character(x)) {
    return(as.factor(x))
  } else {
    return(as.numeric(as.character(x)))
  }
}))
missing_indices <- which(!complete.cases(sample_data_df_20_below))
if(length(missing_indices) > 0) {
  response_variable_20_below <- response_variable_20_below[-missing_indices]
  sample_data_df_20_below <- sample_data_df_20_below[complete.cases(sample_data_df_20_below), ]
}
if (length(response_variable_20_below) == nrow(sample_data_df_20_below)) {
  rf_model_20_below <- randomForest(response_variable_20_below ~ ., data = sample_data_df_20_below, importance = TRUE, ntree = 99999)
  importance_matrix_20_below <- importance(rf_model_20_below)
  print(importance_matrix_20_below)
  varImpPlot(rf_model_20_below, main = "2020 Below Ground")
} else {
  stop("Mismatch in number of observations between response variable and predictor data.")
}

#2021 above
ordination_results_21_above <- ordinate(ps.filtered_2021_above, method = "PCoA", distance = "bray")
pcoa_scores_21_above <- ordination_results_21_above$vectors
response_variable_21_above <- as.numeric(pcoa_scores_21_above[, 1])
sample_data_df_21_above <- as.data.frame(sample_data(ps.filtered_2021_above))
correct_columns <- c(
  "Growth_stage",
  "Management",
  "Precipitation",
  "Ambient.Temperature",
  "Nitrate",
  "Ammonium",
  "pH",
  "Phosphorus",
  "Potassium",
  "Calcium",
  "Magnesium"
)
sample_data_df_21_above <- sample_data_df_21_above[, correct_columns]
sample_data_df_21_above <- data.frame(lapply(sample_data_df_21_above, function(x) {
  if (is.factor(x) || is.character(x)) {
    return(as.factor(x))
  } else {
    return(as.numeric(as.character(x)))
  }
}))
missing_indices <- which(!complete.cases(sample_data_df_21_above))
if(length(missing_indices) > 0) {
  response_variable_21_above <- response_variable_21_above[-missing_indices]
  sample_data_df_21_above <- sample_data_df_21_above[complete.cases(sample_data_df_21_above), ]
}
if (length(response_variable_21_above) == nrow(sample_data_df_21_above)) {
  rf_model_21_above <- randomForest(response_variable_21_above ~ ., data = sample_data_df_21_above, importance = TRUE, ntree = 99999)
  importance_matrix_21_above <- importance(rf_model_21_above)
  print(importance_matrix_21_above)
  varImpPlot(rf_model_21_above, main = "2021 Above Ground")
} else {
  stop("Mismatch in number of observations between response variable and predictor data.")
}

#2021 below
ordination_results_21_below <- ordinate(ps.filtered_2021_below, method = "PCoA", distance = "bray")
pcoa_scores_21_below <- ordination_results_21_below$vectors
response_variable_21_below <- as.numeric(pcoa_scores_21_below[, 1])
sample_data_df_21_below <- as.data.frame(sample_data(ps.filtered_2021_below))
correct_columns <- c(
  "Growth_stage",
  "Management",
  "Precipitation",
  "Ambient.Temperature",
  "Nitrate",
  "Ammonium",
  "pH",
  "Phosphorus",
  "Potassium",
  "Calcium",
  "Magnesium"
)
sample_data_df_21_below <- sample_data_df_21_below[, correct_columns]
sample_data_df_21_below <- data.frame(lapply(sample_data_df_21_below, function(x) {
  if (is.factor(x) || is.character(x)) {
    return(as.factor(x))
  } else {
    return(as.numeric(as.character(x)))
  }
}))
missing_indices <- which(!complete.cases(sample_data_df_21_below))
if(length(missing_indices) > 0) {
  response_variable_21_below <- response_variable_21_below[-missing_indices]
  sample_data_df_21_below <- sample_data_df_21_below[complete.cases(sample_data_df_21_below), ]
}
if (length(response_variable_21_below) == nrow(sample_data_df_21_below)) {
  rf_model_21_below <- randomForest(response_variable_21_below ~ ., data = sample_data_df_21_below, importance = TRUE, ntree = 99999)
  importance_matrix_21_below <- importance(rf_model_21_below)
  print(importance_matrix_21_below)
  varImpPlot(rf_model_21_below, main = "2021 Below Ground")
} else {
  stop("Mismatch in number of observations between response variable and predictor data.")
}


#2018
ordination_results_18 <- ordinate(ps.filtered_2018, method = "PCoA", distance = "bray")
pcoa_scores_18 <- ordination_results_18$vectors
response_variable_18 <- as.numeric(pcoa_scores_18[, 1])
sample_data_df_18 <- as.data.frame(sample_data(ps.filtered_2018))
correct_columns <- c(
  "Growth_stage",
  "Management",
  "Precipitation",
  "Ambient.Temperature",
  "Nitrate",
  "Ammonium",
  "pH",
  "Phosphorus",
  "Potassium",
  "Calcium",
  "Magnesium"
)
sample_data_df_18 <- sample_data_df_18[, correct_columns]
sample_data_df_18 <- data.frame(lapply(sample_data_df_18, function(x) {
  if (is.factor(x) || is.character(x)) {
    return(as.factor(x))
  } else {
    return(as.numeric(as.character(x)))
  }
}))
missing_indices <- which(!complete.cases(sample_data_df_18))
if(length(missing_indices) > 0) {
  response_variable_18 <- response_variable_18[-missing_indices]
  sample_data_df_18 <- sample_data_df_18[complete.cases(sample_data_df_18), ]
}
if (length(response_variable_18) == nrow(sample_data_df_18)) {
  rf_model_18 <- randomForest(response_variable_18 ~ ., data = sample_data_df_18, importance = TRUE, ntree = 99999)
  importance_matrix_18 <- importance(rf_model_18)
  print(importance_matrix_18)
  varImpPlot(rf_model_18, main = "2021 Below Ground")
} else {
  stop("Mismatch in number of observations between response variable and predictor data.")
}
print(importance_matrix_18)
importance_df_18 <- data.frame(
  Feature = rownames(importance_matrix_18),
  IncMSE = importance_matrix_18[, 1],
  IncNodePurity = importance_matrix_18[, 2]
)
importance_df_18$TotalImportance <- importance_df_18$IncMSE + importance_df_18$IncNodePurity
importance_plot_18<- ggplot(importance_df_18, aes(x = reorder(Feature, TotalImportance), y = TotalImportance)) +
  geom_bar(stat = "identity", fill = "#6DA45C", color = "black") +
  coord_flip() +
  labs(title = "Feature Importance - 2018 Soybean", 
       x = "Features", 
       y = "Total Importance") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))
importance_plot_18

#2019
ordination_results_19 <- ordinate(ps.filtered_2019, method = "PCoA", distance = "bray")
pcoa_scores_19 <- ordination_results_19$vectors
response_variable_19 <- as.numeric(pcoa_scores_19[, 1])
sample_data_df_19 <- as.data.frame(sample_data(ps.filtered_2019))
correct_columns <- c(
  "Growth_stage",
  "Management",
  "Precipitation",
  "Ambient.Temperature",
  "Nitrate",
  "Ammonium",
  "pH",
  "Phosphorus",
  "Potassium",
  "Calcium",
  "Magnesium"
)
sample_data_df_19 <- sample_data_df_19[, correct_columns]
sample_data_df_19 <- data.frame(lapply(sample_data_df_19, function(x) {
  if (is.factor(x) || is.character(x)) {
    return(as.factor(x))
  } else {
    return(as.numeric(as.character(x)))
  }
}))
missing_indices <- which(!complete.cases(sample_data_df_19))
if(length(missing_indices) > 0) {
  response_variable_19 <- response_variable_19[-missing_indices]
  sample_data_df_19 <- sample_data_df_19[complete.cases(sample_data_df_19), ]
}
if (length(response_variable_19) == nrow(sample_data_df_19)) {
  rf_model_19 <- randomForest(response_variable_19 ~ ., data = sample_data_df_19, importance = TRUE, ntree = 99999)
  importance_matrix_19 <- importance(rf_model_19)
  print(importance_matrix_19)
  varImpPlot(rf_model_19, main = "2021 Below Ground")
} else {
  stop("Mismatch in number of observations between response variable and predictor data.")
}
print(importance_matrix_19)
importance_df_19 <- data.frame(
  Feature = rownames(importance_matrix_19),
  IncMSE = importance_matrix_19[, 1],
  IncNodePurity = importance_matrix_19[, 2]
)
importance_df_19$TotalImportance <- importance_df_19$IncMSE + importance_df_19$IncNodePurity
importance_plot_19<- ggplot(importance_df_19, aes(x = reorder(Feature, TotalImportance), y = TotalImportance)) +
  geom_bar(stat = "identity", fill = "#85614F", color = "black") +
  coord_flip() +
  labs(title = "Feature Importance - 2019 Wheat", 
       x = "Features", 
       y = "Total Importance") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))
importance_plot_19

#2020
ordination_results_20 <- ordinate(ps.filtered_2020, method = "PCoA", distance = "bray")
pcoa_scores_20 <- ordination_results_20$vectors
response_variable_20 <- as.numeric(pcoa_scores_20[, 1])
sample_data_df_20 <- as.data.frame(sample_data(ps.filtered_2020))
correct_columns <- c(
  "Growth_stage",
  "Management",
  "Precipitation",
  "Ambient.Temperature",
  "Nitrate",
  "Ammonium",
  "pH",
  "Phosphorus",
  "Potassium",
  "Calcium",
  "Magnesium"
)
sample_data_df_20 <- sample_data_df_20[, correct_columns]
sample_data_df_20 <- data.frame(lapply(sample_data_df_20, function(x) {
  if (is.factor(x) || is.character(x)) {
    return(as.factor(x))
  } else {
    return(as.numeric(as.character(x)))
  }
}))
missing_indices <- which(!complete.cases(sample_data_df_20))
if(length(missing_indices) > 0) {
  response_variable_20 <- response_variable_20[-missing_indices]
  sample_data_df_20 <- sample_data_df_20[complete.cases(sample_data_df_20), ]
}
if (length(response_variable_20) == nrow(sample_data_df_20)) {
  rf_model_20 <- randomForest(response_variable_20 ~ ., data = sample_data_df_20, importance = TRUE, ntree = 99999)
  importance_matrix_20 <- importance(rf_model_20)
  print(importance_matrix_20)
  varImpPlot(rf_model_20, main = "2020 Below Ground")
} else {
  stop("Mismatch in number of observations between response variable and predictor data.")
}
print(importance_matrix_20)
importance_df_20 <- data.frame(
  Feature = rownames(importance_matrix_20),
  IncMSE = importance_matrix_20[, 1],
  IncNodePurity = importance_matrix_20[, 2]
)
importance_df_20$TotalImportance <- importance_df_20$IncMSE + importance_df_20$IncNodePurity
importance_plot_20<- ggplot(importance_df_20, aes(x = reorder(Feature, TotalImportance), y = TotalImportance)) +
  geom_bar(stat = "identity", fill = "#F3B342", color = "black") +
  coord_flip() +
  labs(title = "Feature Importance - 2020 Maize", 
       x = "Features", 
       y = "Total Importance") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))
importance_plot_20


#2021
ordination_results_21 <- ordinate(ps.filtered_2021, method = "PCoA", distance = "bray")
pcoa_scores_21 <- ordination_results_21$vectors
response_variable_21 <- as.numeric(pcoa_scores_21[, 1])
sample_data_df_21 <- as.data.frame(sample_data(ps.filtered_2021))
correct_columns <- c(
  "Growth_stage",
  "Management",
  "Precipitation",
  "Ambient.Temperature",
  "Nitrate",
  "Ammonium",
  "pH",
  "Phosphorus",
  "Potassium",
  "Calcium",
  "Magnesium"
)
sample_data_df_21 <- sample_data_df_21[, correct_columns]
sample_data_df_21 <- data.frame(lapply(sample_data_df_21, function(x) {
  if (is.factor(x) || is.character(x)) {
    return(as.factor(x))
  } else {
    return(as.numeric(as.character(x)))
  }
}))
missing_indices <- which(!complete.cases(sample_data_df_21))
if(length(missing_indices) > 0) {
  response_variable_21 <- response_variable_21[-missing_indices]
  sample_data_df_21 <- sample_data_df_21[complete.cases(sample_data_df_21), ]
}
if (length(response_variable_21) == nrow(sample_data_df_21)) {
  rf_model_21 <- randomForest(response_variable_21 ~ ., data = sample_data_df_21, importance = TRUE, ntree = 99999)
  importance_matrix_21 <- importance(rf_model_21)
  print(importance_matrix_21)
  varImpPlot(rf_model_21, main = "2021 Below Ground")
} else {
  stop("Mismatch in number of observations between response variable and predictor data.")
}
print(importance_matrix_21)
importance_df_21 <- data.frame(
  Feature = rownames(importance_matrix_21),
  IncMSE = importance_matrix_21[, 1],
  IncNodePurity = importance_matrix_21[, 2]
)
importance_df_21$TotalImportance <- importance_df_21$IncMSE + importance_df_21$IncNodePurity
importance_plot_21<- ggplot(importance_df_21, aes(x = reorder(Feature, TotalImportance), y = TotalImportance)) +
  geom_bar(stat = "identity", fill = "#6DA45C", color = "black") +
  coord_flip() +
  labs(title = "Feature Importance - 2021 Soybean", 
       x = "Features", 
       y = "Total Importance") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))
importance_plot_21

all_rd_plots<- ggarrange(importance_plot_18, importance_plot_19, importance_plot_20, importance_plot_21,
                         nrow= 1, ncol=4)
all_rd_plots

#barplots-------------------------------------------------------------------------------------------------------------------------------------
genus_colors <- c(
  "Alternaria" = "#1B9E77",
  "Apodus" = "#6E8243",
  "Athelia" = "#C2660F",
  "Berkeleyomyces" = "wheat",
  "Blumeria" = "#8C6C89",
  "Bullera" = "blueviolet",
  "Cercospora" = "#BF4198",
  "Chaetomium" = "cadetblue",
  "Colletotrichum" = "#A26B50",
  "Conlarium" = "#69A320",
  "Corynespora" = "#9BA812",
  "Dioszegia" = "aquamarine",
  "Exophiala" = "#D29A0A",
  "Funneliformis" = "hotpink1",
  "Gigaspora" = "pink",
  "Glomus" = "#7D6B4A",
  "Hannaella" = "yellowgreen",
  "Humicola" = "darkgreen",
  "Knufia" = "#A3C9DD",
  "Leucosporidium" = "#70ABD0",
  "Linnemannia" = "steelblue4",
  "Macrophomina" = "seagreen",
  "Metacordyceps" = "#88C295",
  "Minimedusa polyspora" = "tomato4",
  "Mortierella" = "firebrick3",
  "Myrmecridium" = "darkslateblue",
  "Mycoleptodiscus" = "orchid1",
  "Neosetophoma" = "#ED9A91",
  "Oedocephalum" = "#F26A6A",
  "Paecilomyces" = "#E73133",
  "Paraglomus" = "#E94431",
  "Peziza" = "steelblue2",
  "Phaeosphaeria" = "#FDB65F",
  "Phallus" = "#FE992E",
  "Rhizophlyctis" = "#FD8004",
  "Sarocladium" = "midnightblue",
  "Schizothecium" = "#CEADC2",
  "Septoglomus" = "#A889C1",
  "Sporobolomyces" = "mediumseagreen",
  "Taphrina" = "#8C6A99",
  "Tilletiopsis" = "chocolate4",
  "Ustilago" = "#F5EB8B",
  "Vishniacozyma" = "#D3A259",
  "Other" = "#6b6b6b"
)


ps.merged_above<- merge_phyloseq(ps.filtered_2018_above, ps.filtered_2019_above, ps.filtered_2020_above,ps.filtered_2021_above)
ps.merged_above<- subset_samples(ps.merged_above, Management %in% c("Conventional", "No-Till", "Organic"))
BP_above<- ps.merged_above %>%
  merge_samples("Label_2")

BP_meta_above <- as.data.frame(as.matrix(BP_above@sam_data)) %>%
  dplyr::select(-Label_2) %>%
  tibble::rownames_to_column("Label_2") %>%
  dplyr::select("Label_2") %>%
  tidyr::separate(Label_2, c("Year", "Management","Compartment"), remove=FALSE) %>%
  dplyr::mutate(trick = Label_2) %>% 
  tibble::column_to_rownames("trick")

BP_meta_above

BP_above@sam_data <- sample_data(BP_meta_above)
BP_above@sam_data

BP_above_mapping_G <- BP_above %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by Family
BP_above_dt_G <- data.table(BP_above_mapping_G)
BP_above_dt_G
BP_above_dt_G[(Abundance <= 0.04), Genus:= "Other"]
BP_above_dt_G

combined_fungal_list <- c("Alternaria", "Apodus", "Athelia", "Berkeleyomyces", "Blumeria", "Bullera",
                          "Cercospora", "Chaetomium", "Colletotrichum", "Conlarium", "Corynespora",
                          "Dioszegia", "Exophiala", "Funneliformis", "Gigaspora", "Glomus",
                          "Hannaella", "Humicola", "Knufia", "Leucosporidium", "Linnemannia",
                          "Macrophomina", "Metacordyceps", "Minimedusa polyspora", "Mortierella",
                          "Myrmecridium", "Mycoleptodiscus", "Neosetophoma", "Oedocephalum", "Paecilomyces",
                          "Paraglomus", "Peziza", "Phaeosphaeria", "Phallus", "Rhizophlyctis",
                          "Sarocladium", "Schizothecium", "Septoglomus", "Sporobolomyces", "Taphrina",
                          "Tilletiopsis", "Ustilago", "Vishniacozyma", "Other")

BP_above_dt_G$Genus <- as.character(BP_above_dt_G$Genus)
BP_above_dt_G$Genus <- factor(BP_above_dt_G$Genus, levels = combined_fungal_list)

Bar_above_G= ggplot(BP_above_dt_G, aes(x = Label_2, y = Abundance, fill = Genus)) + 
  facet_wrap(~Year, strip.position = "bottom", scales = "free_x", nrow=1) +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = genus_colors)+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 0, face = "bold")) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  ggtitle("Above Ground fungal communities")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = .3)) +
  ylab("Relative Abundance (Genera > 4%) \n") +
  xlab("")
plot(Bar_above_G)

ps.merged_below<- merge_phyloseq(ps.filtered_2018_below, ps.filtered_2019_below, ps.filtered_2020_below,ps.filtered_2021_below)
ps.merged_below<- subset_samples(ps.merged_below, Management %in% c("Conventional", "No-Till", "Organic"))
BP_below<- ps.merged_below %>%
  merge_samples("Label_2")

BP_meta_below <- as.data.frame(as.matrix(BP_below@sam_data)) %>%
  dplyr::select(-Label_2) %>%
  tibble::rownames_to_column("Label_2") %>%
  dplyr::select("Label_2") %>%
  tidyr::separate(Label_2, c("Year", "Management","Compartment"), remove=FALSE) %>%
  dplyr::mutate(trick = Label_2) %>% 
  tibble::column_to_rownames("trick")

BP_meta_below

BP_below@sam_data <- sample_data(BP_meta_below)
BP_below@sam_data

BP_below_mapping_G <- BP_below %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by Family
BP_below_dt_G <- data.table(BP_below_mapping_G)
BP_below_dt_G
BP_below_dt_G[(Abundance <= 0.04), Genus:= "Other"]
BP_below_dt_G

combined_fungal_list <- c("Alternaria", "Apodus", "Athelia", "Berkeleyomyces", "Blumeria", "Bullera",
                          "Cercospora", "Chaetomium", "Colletotrichum", "Conlarium", "Corynespora",
                          "Dioszegia", "Exophiala", "Funneliformis", "Gigaspora", "Glomus",
                          "Hannaella", "Humicola", "Knufia", "Leucosporidium", "Linnemannia",
                          "Macrophomina", "Metacordyceps", "Minimedusa polyspora", "Mortierella",
                          "Myrmecridium", "Mycoleptodiscus", "Neosetophoma", "Oedocephalum", "Paecilomyces",
                          "Paraglomus", "Peziza", "Phaeosphaeria", "Phallus", "Rhizophlyctis",
                          "Sarocladium", "Schizothecium", "Septoglomus", "Sporobolomyces", "Taphrina",
                          "Tilletiopsis", "Ustilago", "Vishniacozyma", "Other")

BP_below_dt_G$Genus <- as.character(BP_below_dt_G$Genus)
BP_below_dt_G$Genus <- factor(BP_below_dt_G$Genus, levels = combined_fungal_list)

Bar_below_G= ggplot(BP_below_dt_G, aes(x = Label_2, y = Abundance, fill = Genus)) + 
  facet_wrap(~Year, strip.position = "bottom", scales = "free_x", nrow=1) +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = genus_colors)+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 0, face = "bold")) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  ggtitle("Below Ground fungal communities")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = .3)) +
  ylab("Relative Abundance (Genera > 4%) \n") +
  xlab("")
plot(Bar_below_G)

all_bar_plot<- ggarrange(Bar_above_G, Bar_below_G, nrow = 1)
all_bar_plot

