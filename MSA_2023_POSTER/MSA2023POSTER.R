###MSA 2023 poster, compiled code
###Ashlynn Morin

#this code is a mixture of Dr. Reid Longley's from his
#github (https://github.com/longleyr/Management-of-Soybean-Code-and-Files), and my code
#updating some of his code, and adding some of my own.


#STARTING
##load in
library(stringi)
library(rhdf5)
library(zlibbioc)
library(S4Vectors)
library(phyloseq)
library(Biostrings)
library(yaml)
library(colorspace)
library(ggplot2)
library(indicspecies)
library(vegan)
library(phyloseq)
library(decontam)
library(devtools)
library(processx)
install.packages("psych")
library("psych")
install.packages("agricolae")
library("agricolae")
library(dplyr)
library(data.table)
library(readxl)
library(ggpubbr)
library(tidyverse)
#setworking directory
setwd("~/Desktop/Microbiome/MSA23")

###Create a phyloseq object
#otu_table
ITS_otus<- read.delim("combined_otu_table_R1.txt",
                      row.names=1) 
head(ITS_otus)
ITS_otus_phy <-otu_table(ITS_otus,
                         taxa_are_rows = TRUE)
ITS_otus_phy
#metadata
ITS_metadata <-read.delim("Combined_Map.txt",
                          row.names=1)
ITS_metadata
ITS_metadata_phy <-sample_data(ITS_metadata)
#Taxonomy
ITS_taxonomy<-read.delim("constax_taxonomy.txt",
                         header=TRUE, 
                         row.names=1)
ITS_taxonomy
ITS_taxonomy_phy <- tax_table(as.matrix(ITS_taxonomy))
#Sequences
ITS_sequences <- readDNAStringSet("otus_R1.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
ITS_sequences
#the object
physeq_object <- phyloseq(ITS_otus_phy,
                          ITS_metadata_phy,
                          ITS_taxonomy_phy,
                          ITS_sequences)
physeq_object
tax_table(physeq_object)
sample_data(physeq_object)

###Format taxonomy
colnames(tax_table(physeq_object)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
tax_table(physeq_object)
tax_table(physeq_object)[, "Kingdom"] <- gsub("_1", "", tax_table(physeq_object)[, "Kingdom"])
tax_table(physeq_object)[, "Phylum"] <- gsub("_1", "", tax_table(physeq_object)[, "Phylum"])
tax_table(physeq_object)[, "Class"] <- gsub("_1", "", tax_table(physeq_object)[, "Class"])
tax_table(physeq_object)[, "Order"] <- gsub("_1", "", tax_table(physeq_object)[, "Order"])
tax_table(physeq_object)[, "Family"] <- gsub("_1", "", tax_table(physeq_object)[, "Family"])
tax_table(physeq_object)[, "Genus"] <- gsub("_1", "", tax_table(physeq_object)[, "Genus"])
tax_table(physeq_object)[, "Species"] <- gsub("_1", "", tax_table(physeq_object)[, "Species"])
tax_table(physeq_object)
physeq_object <- subset_taxa(physeq_object, Phylum!="Chloroplast")
physeq_object <- subset_taxa(physeq_object, Class!="Chloroplast")
physeq_object <- subset_taxa(physeq_object, Order!="Chloroplast")
physeq_object <- subset_taxa(physeq_object, Family!="Chloroplast")
physeq_object <- subset_taxa(physeq_object, Genus!="Chloroplast")
tax_table(physeq_object)
physeq_object <- subset_taxa(physeq_object, Phylum!="Mitochondria")
physeq_object <- subset_taxa(physeq_object, Class!="Mitochondria")
physeq_object <- subset_taxa(physeq_object, Order!="Mitochondria")
physeq_object <- subset_taxa(physeq_object, Family!="Mitochondria")
physeq_object <- subset_taxa(physeq_object, Genus!="Mitochondria")
tax_table(physeq_object)
physeq_object <- subset_taxa(physeq_object, Kingdom!="Anthophyta")
physeq_object <- subset_taxa(physeq_object, Kingdom!="Alveolata")
physeq_object <- subset_taxa(physeq_object, Kingdom!="Ichthyosporia")
physeq_object <- subset_taxa(physeq_object, Kingdom!="Protista")
physeq_object <- subset_taxa(physeq_object, Kingdom!="Metazoa")
physeq_object <- subset_taxa(physeq_object, Kingdom!="Rhizaria")
physeq_object <- subset_taxa(physeq_object, Kingdom!="Viridiplantae")
sort(unique(as.data.frame(tax_table(physeq_object))$Kingdom))
physeq_object

###Break apart by origin (all years combined)
physeq_L<- subset_samples(physeq_object, ORIGIN%in%c("leaf"))
sample_data(physeq_L)
physeq_R<- subset_samples(physeq_object, ORIGIN%in%c("root"))
sample_data(physeq_R)
physeq_So<- subset_samples(physeq_object, ORIGIN%in%c("soil"))
sample_data(physeq_So)
physeq_St<- subset_samples(physeq_object, ORIGIN%in%c("stem"))
sample_data(physeq_St)


###Check library distribution
#all leaves
df_physeq_L <- as.data.frame(sample_data(physeq_L))
df_physeq_L$LibrarySize_leaf <- sample_sums(physeq_L)
df_physeq_L <- df_physeq_L[order(df_physeq_L$LibrarySize_leaf),]
df_physeq_L$Index <- seq(nrow(df_physeq_L))
ggplot(data=df_physeq_L, aes(x=Index, y=LibrarySize_leaf, color=SAMPLE_OR_CONTROL)) + geom_point()

#all roots
df_physeq_R <- as.data.frame(sample_data(physeq_R))
df_physeq_R$LibrarySize_root <- sample_sums(physeq_R)
df_physeq_R <- df_physeq_R[order(df_physeq_R$LibrarySize_root),]
df_physeq_R$Index <- seq(nrow(df_physeq_R))
ggplot(data=df_physeq_R, aes(x=Index, y=LibrarySize_root, color=SAMPLE_OR_CONTROL)) + geom_point()

#all soils
df_physeq_So <- as.data.frame(sample_data(physeq_So))
df_physeq_So$LibrarySize_soil <- sample_sums(physeq_So)
df_physeq_So <- df_physeq_So[order(df_physeq_So$LibrarySize_soil),]
df_physeq_So$Index <- seq(nrow(df_physeq_So))
ggplot(data=df_physeq_So, aes(x=Index, y=LibrarySize_soil, color=SAMPLE_OR_CONTROL)) + geom_point()

#all stems
df_physeq_St <- as.data.frame(sample_data(physeq_St))
df_physeq_St$LibrarySize_stem <- sample_sums(physeq_St)
df_physeq_St <- df_physeq_St[order(df_physeq_St$LibrarySize_stem),]
df_physeq_St$Index <- seq(nrow(df_physeq_St))
ggplot(data=df_physeq_St, aes(x=Index, y=LibrarySize_stem, color=SAMPLE_OR_CONTROL)) + geom_point()

###Filter by prevalence
#leaf
sample_data(physeq_L)$is.neg <- sample_data(physeq_L)$SAMPLE_OR_CONTROL == "Control"
contamdf.prev_L <- isContaminant(physeq_L, method="prevalence", neg="is.neg")
table(contamdf.prev_L$contaminant)

#root
sample_data(physeq_R)$is.neg <- sample_data(physeq_R)$SAMPLE_OR_CONTROL == "Control"
contamdf.prev_R <- isContaminant(physeq_R, method="prevalence", neg="is.neg")
table(contamdf.prev_R$contaminant)

#soil
sample_data(physeq_So)$is.neg <- sample_data(physeq_So)$SAMPLE_OR_CONTROL == "Control"
contamdf.prev_So <- isContaminant(physeq_So, method="prevalence", neg="is.neg")
table(contamdf.prev_So$contaminant)

#stem
sample_data(physeq_St)$is.neg <- sample_data(physeq_St)$SAMPLE_OR_CONTROL == "Control"
contamdf.prev_St <- isContaminant(physeq_St, method="prevalence", neg="is.neg")
table(contamdf.prev_St$contaminant)

###Remove contaminants and negative controls
#leaf
ps.noncontam_L <- prune_taxa(!contamdf.prev_L$contaminant, physeq_L) #remove contaminants
ps.noncontam_L<-subset_samples(ps.noncontam_L,SAMPLE_OR_CONTROL%in%c("Sample")) #remove negative control
ps.noncontam_L
otu_table(ps.noncontam_L) <- otu_table(ps.noncontam_L)[which(rowSums(otu_table(ps.noncontam_L)) >= 1),]
ps.noncontam_L
otu_table(ps.noncontam_L)
sample_data(ps.noncontam_L)
ps.noncontam_L

#root
ps.noncontam_R<- prune_taxa(!contamdf.prev_R$contaminant, physeq_R)#remove contaminants
ps.noncontam_R<- subset_samples(ps.noncontam_R,SAMPLE_OR_CONTROL%in%c("Sample"))#remove negative control
ps.noncontam_R
otu_table(ps.noncontam_R) <- otu_table(ps.noncontam_R)[which(rowSums(otu_table(ps.noncontam_R)) >= 1),]
ps.noncontam_R
otu_table(ps.noncontam_R)
sample_data(ps.noncontam_R)
ps.noncontam_R

#soil
ps.noncontam_So<-prune_taxa(!contamdf.prev_So$contaminant, physeq_So)#remove contaminants
ps.noncontam_So<-subset_samples(ps.noncontam_So, SAMPLE_OR_CONTROL%in%c("Sample"))#remove negative control
ps.noncontam_So
otu_table(ps.noncontam_So) <- otu_table(ps.noncontam_So)[which(rowSums(otu_table(ps.noncontam_So)) >= 1),]
ps.noncontam_So
otu_table(ps.noncontam_So)
sample_data(ps.noncontam_So)
ps.noncontam_So
#stems
ps.noncontam_St<- prune_taxa(!contamdf.prev_St$contaminant, physeq_St)#remove contaminants
ps.noncontam_St<-subset_samples(ps.noncontam_St, SAMPLE_OR_CONTROL%in%c("Sample"))#remove negative control
ps.noncontam_St
otu_table(ps.noncontam_St) <- otu_table(ps.noncontam_St)[which(rowSums(otu_table(ps.noncontam_St)) >= 1),]
ps.noncontam_St
otu_table(ps.noncontam_St)
sample_data(ps.noncontam_St)
ps.noncontam_St

merged_phylo<-merge_phyloseq(ps.noncontam_L,ps.noncontam_R,ps.noncontam_So,ps.noncontam_St)

###RAREFACTION CURVES
#leaf
otu_rare_L <- as.data.frame(otu_table(ps.noncontam_L))
metadata_L_rare <- as.data.frame(sample_data(ps.noncontam_L))
L_rarecurve <- rarecurve(t(otu_rare_L), col = as.factor(metadata_L_rare$MANAGEMENT), label = FALSE, 
                         step = 50,
                         main="Fungi_leaf", ylab = "Number of OTUs", xlab = "Number of DNA reads")-> rare_L_fungi
legend("bottomright", legend=c("No-Till", "Organic", "Conventional"),
       col=c("black", "green", "red"), lty=1, cex=0.8, box.lty=1)

#roots
otu_rare_R <- as.data.frame(otu_table(ps.noncontam_R))
metadata_R_rare <- as.data.frame(sample_data(ps.noncontam_R))
R_rarecurve <- rarecurve(t(otu_rare_R), col = as.factor(metadata_R_rare$MANAGEMENT), label = FALSE, 
                         step = 50,
                         main="Fungi_Root", ylab = "Number of OTUs", xlab = "Number of DNA reads")-> rare_R_fungi
legend("bottomright", legend=c("No-Till", "Organic", "Conventional"),
       col=c("black", "green", "red"), lty=1, cex=0.8, box.lty=1)

#soil
otu_rare_So <- as.data.frame(otu_table(ps.noncontam_So))
metadata_So_rare <- as.data.frame(sample_data(ps.noncontam_So))
So_rarecurve <- rarecurve(t(otu_rare_So), col = as.factor(metadata_So_rare$MANAGEMENT), label = FALSE, 
                          step = 50,
                          main="Fungi_soil", ylab = "Number of OTUs", xlab = "Number of DNA reads")-> rare_So_fungi
legend("bottomright", legend=c("No-Till", "Organic", "Conventional"),
       col=c("black", "green", "red"), lty=1, cex=0.8, box.lty=1)

#Stem
otu_rare_St <- as.data.frame(otu_table(ps.noncontam_St))
metadata_St_rare <- as.data.frame(sample_data(ps.noncontam_St))
St_rarecurve <- rarecurve(t(otu_rare_St), col = as.factor(metadata_St_rare$MANAGEMENT), label = FALSE, 
                          step = 50,
                          main="Fungi_stem", ylab = "Number of OTUs", xlab = "Number of DNA reads")-> rare_St_fungi
legend("bottomright", legend=c("No-Till", "Organic", "Conventional"),
       col=c("black", "green", "red"), lty=1, cex=0.8, box.lty=1)


###ALPHA DIVERSITY
label_names <- c(Observed="Richness", Shannon="Shannon")
label_names
###leaf
ps.noncontam_alpha_L <- ps.noncontam_L
sample_data(ps.noncontam_alpha_L)$BAR_LABEL <- factor(sample_data(ps.noncontam_alpha_L)$BAR_LABEL,
                                                      level=c("2018 V2 Conventional",
                                                              "2018 R2 Conventional",
                                                              "2018 R6 Conventional",
                                                              "2018 V2 No-Till",
                                                              "2018 R2 No-Till",
                                                              "2018 R6 No-Till",
                                                              "2018 V2 Organic",
                                                              "2018 R2 Organic",
                                                              "2018 R6 Organic",
                                                              "2021 V2 Conventional",
                                                              "2021 R3 Conventional",
                                                              "2021 R6 Conventional",
                                                              "2021 V2 No-Till",
                                                              "2021 R3 No-Till",
                                                              "2021 R6 No-Till",
                                                              "2021 V2 Organic",
                                                              "2021 R3 Organic",
                                                              "2021 R6 Organic"))

alpha_L <- estimate_richness(ps.noncontam_alpha_L, split = TRUE, measures = "Observed")
alpha_L_fungi = plot_richness(ps.noncontam_alpha_L, x= "BAR_LABEL", 
                              color="GROWTH_STAGE", measures = c("Observed")) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 2, shape = 16) +
  labs(title="Leaf", y = "Observed OTUs") +
  scale_colour_manual("GROWTH_STAGE",breaks = c("V2","R2", "R3","R6"),
                      values = c("V2" = "#56B4E9", "R2" = "#FFB5B8", "R3" = "#009E73", "R6" = "#E69F00")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 
  theme_set(theme_classic())+
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 11, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  theme(legend.position="none") +
  ylim(0,900)+
  theme(legend.title=element_blank())
plot(alpha_L_fungi)


alpha_L_fungi_shan = plot_richness(ps.noncontam_alpha_L, x= "BAR_LABEL", 
                                   color="GROWTH_STAGE", measures = c( "Shannon")) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 2, shape = 16) +
  labs(title="", x="", y = "Shannon Diversity") +
  scale_colour_manual("GROWTH_STAGE",breaks = c("V2","R2", "R3","R6"),
                      values = c("V2" = "#56B4E9", "R2" = "#FFB5B8", "R3" = "#009E73", "R6" = "#E69F00")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 
  theme_set(theme_classic())+
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 11, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  ylim(0,6)+
  theme(legend.position="none") 
plot(alpha_L_fungi_shan)

alpha_L<- ggarrange(alpha_L_fungi, alpha_L_fungi_shan)
alpha_L
#ROOT
ps.noncontam_alpha_R <- ps.noncontam_R
sample_data(ps.noncontam_alpha_R)$BAR_LABEL <- factor(sample_data(ps.noncontam_alpha_R)$BAR_LABEL,
                                                      level=c("2018 V2 Conventional",
                                                              "2018 R2 Conventional",
                                                              "2018 R6 Conventional",
                                                              "2018 V2 No-Till",
                                                              "2018 R2 No-Till",
                                                              "2018 R6 No-Till",
                                                              "2018 V2 Organic",
                                                              "2018 R2 Organic",
                                                              "2018 R6 Organic",
                                                              "2021 V2 Conventional",
                                                              "2021 R3 Conventional",
                                                              "2021 R6 Conventional",
                                                              "2021 V2 No-Till",
                                                              "2021 R3 No-Till",
                                                              "2021 R6 No-Till",
                                                              "2021 V2 Organic",
                                                              "2021 R3 Organic",
                                                              "2021 R6 Organic"))

alpha_R <- estimate_richness(ps.noncontam_alpha_R, split = TRUE, measures = "Observed")
alpha_R_fungi = plot_richness(ps.noncontam_alpha_R, x= "BAR_LABEL", 
                              color="GROWTH_STAGE", measures = c( "Observed")) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 2, shape = 16) +
  labs(title="Root", y = "Observed OTUs") +
  scale_colour_manual("GROWTH_STAGE",breaks = c("V2","R2", "R3","R6"),
                      values = c("V2" = "#56B4E9", "R2" = "#FFB5B8", "R3" = "#009E73", "R6" = "#E69F00")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 
  theme_set(theme_classic())+
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 11, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  theme(legend.position="none") +
  ylim(0,900)+
  theme(legend.title=element_blank())
plot(alpha_R_fungi)


alpha_R_fungi_shan = plot_richness(ps.noncontam_alpha_R, x= "BAR_LABEL", 
                                   color="GROWTH_STAGE", measures = c( "Shannon")) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 2, shape = 16) +
  labs(title="", x="", y = "Shannon Diversity") +
  scale_colour_manual("GROWTH_STAGE",breaks = c("V2","R2", "R3","R6"),
                      values = c("V2" = "#56B4E9", "R2" = "#FFB5B8", "R3" = "#009E73", "R6" = "#E69F00")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 
  theme_set(theme_classic())+
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 11, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  ylim(0,6)+
  theme(legend.position="none") 
plot(alpha_R_fungi_shan)

alpha_R<- ggarrange(alpha_R_fungi, alpha_R_fungi_shan)
alpha_R
#soil
ps.noncontam_alpha_So <- ps.noncontam_So
sample_data(ps.noncontam_alpha_So)$BAR_LABEL <- factor(sample_data(ps.noncontam_alpha_So)$BAR_LABEL,
                                                       level=c("2018 V2 Conventional",
                                                               "2018 R2 Conventional",
                                                               "2018 R6 Conventional",
                                                               "2018 V2 No-Till",
                                                               "2018 R2 No-Till",
                                                               "2018 R6 No-Till",
                                                               "2018 V2 Organic",
                                                               "2018 R2 Organic",
                                                               "2018 R6 Organic",
                                                               "2021 V2 Conventional",
                                                               "2021 R3 Conventional",
                                                               "2021 R6 Conventional",
                                                               "2021 V2 No-Till",
                                                               "2021 R3 No-Till",
                                                               "2021 R6 No-Till",
                                                               "2021 V2 Organic",
                                                               "2021 R3 Organic",
                                                               "2021 R6 Organic"))

alpha_So <- estimate_richness(ps.noncontam_alpha_So, split = TRUE, measures = "Observed")
alpha_So_fungi = plot_richness(ps.noncontam_alpha_So, x= "BAR_LABEL", 
                               color="GROWTH_STAGE", measures = c( "Observed")) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 2, shape = 16) +
  labs(title="Soil", y = "Observed OTUs") +
  scale_colour_manual("GROWTH_STAGE",breaks = c("V2","R2", "R3","R6"),
                      values = c("V2" = "#56B4E9", "R2" = "#FFB5B8", "R3" = "#009E73", "R6" = "#E69F00")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 
  theme_set(theme_classic())+
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 11, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  theme(legend.position="none") +
  ylim(0,900)+
  theme(legend.title=element_blank())
plot(alpha_So_fungi)


alpha_So_fungi_shan = plot_richness(ps.noncontam_alpha_So, x= "BAR_LABEL", 
                                    color="GROWTH_STAGE", measures = c( "Shannon")) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 2, shape = 16) +
  labs(title="", x="", y = "Shannon Diversity") +
  scale_colour_manual("GROWTH_STAGE",breaks = c("V2","R2", "R3","R6"),
                      values = c("V2" = "#56B4E9", "R2" = "#FFB5B8", "R3" = "#009E73", "R6" = "#E69F00")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 
  theme_set(theme_classic())+
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 11, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  ylim(0,6)+
  theme(legend.position="none") 
plot(alpha_So_fungi_shan)

alpha_So<- ggarrange(alpha_So_fungi, alpha_So_fungi_shan)
#stem
ps.noncontam_alpha_St <- ps.noncontam_St
sample_data(ps.noncontam_alpha_St)$BAR_LABEL <- factor(sample_data(ps.noncontam_alpha_St)$BAR_LABEL,
                                                       level=c("2018 V2 Conventional",
                                                               "2018 R2 Conventional",
                                                               "2018 R6 Conventional",
                                                               "2018 V2 No-Till",
                                                               "2018 R2 No-Till",
                                                               "2018 R6 No-Till",
                                                               "2018 V2 Organic",
                                                               "2018 R2 Organic",
                                                               "2018 R6 Organic",
                                                               "2021 V2 Conventional",
                                                               "2021 R3 Conventional",
                                                               "2021 R6 Conventional",
                                                               "2021 V2 No-Till",
                                                               "2021 R3 No-Till",
                                                               "2021 R6 No-Till",
                                                               "2021 V2 Organic",
                                                               "2021 R3 Organic",
                                                               "2021 R6 Organic"))

alpha_St <- estimate_richness(ps.noncontam_alpha_St, split = TRUE, measures = "Observed")
alpha_St_fungi = plot_richness(ps.noncontam_alpha_St, x= "BAR_LABEL", 
                               color="GROWTH_STAGE", measures = c( "Observed")) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 2, shape = 16) +
  labs(title="Stem", y = "Observed OTUs") +
  scale_colour_manual("GROWTH_STAGE",breaks = c("V2","R2", "R3","R6"),
                      values = c("V2" = "#56B4E9", "R2" = "#FFB5B8", "R3" = "#009E73", "R6" = "#E69F00")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 
  theme_set(theme_classic())+
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 11, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  theme(legend.position="none") +
  ylim(0,900)+
  theme(legend.title=element_blank())
plot(alpha_St_fungi)


alpha_St_fungi_shan = plot_richness(ps.noncontam_alpha_St, x= "BAR_LABEL", 
                                    color="GROWTH_STAGE", measures = c( "Shannon")) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 2, shape = 16) +
  labs(title="", x="", y = "Shannon Diversity") +
  scale_colour_manual("GROWTH_STAGE",breaks = c("V2","R2", "R3","R6"),
                      values = c("V2" = "#56B4E9", "R2" = "#FFB5B8", "R3" = "#009E73", "R6" = "#E69F00")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 
  theme_set(theme_classic())+
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 11, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  ylim(0,6)+
  theme(legend.position="none") 
plot(alpha_St_fungi_shan)

alpha_St<- ggarrange(alpha_St_fungi, alpha_St_fungi_shan)

alpha_overall<- ggarrange(alpha_L, alpha_R, alpha_So, alpha_St,
                          labels = c("A", "B", "C", "D"))
alpha_overall

###rarify
set.seed(219)
ps.noncontam_Rare_L<- rarefy_even_depth(ps.noncontam_L, sample.size = min(sample_sums(ps.noncontam_L)), rngseed = 219)
ps.noncontam_Rare_R<- rarefy_even_depth(ps.noncontam_R, sample.size = min(sample_sums(ps.noncontam_R)), rngseed = 219)
ps.noncontam_Rare_So<- rarefy_even_depth(ps.noncontam_So, sample.size = min(sample_sums(ps.noncontam_So)), rngseed = 219)
ps.noncontam_Rare_St<- rarefy_even_depth(ps.noncontam_St, sample.size = min(sample_sums(ps.noncontam_St)), rngseed = 219)
merged_rare<- merge_phyloseq(ps.noncontam_Rare_L,ps.noncontam_Rare_R,ps.noncontam_Rare_So,ps.noncontam_Rare_St)

#PLOTTING STACKED_BAR LEAF
ps.noncontam_Rare_meta_L<- sample_data(ps.noncontam_Rare_L)
HM_L<- ps.noncontam_Rare_L %>%
  merge_samples("BAR_LABEL")

HM_meta_L <- as.data.frame(as.matrix(HM_L@sam_data)) %>%
  select(-BAR_LABEL) %>%
  rownames_to_column("BAR_LABEL") %>% 
  select("BAR_LABEL") %>% 
  separate(BAR_LABEL, c("Year", "Growth_Stage","Management"), remove=FALSE) %>%
  mutate(trick = BAR_LABEL) %>% 
  column_to_rownames("trick")

HM_meta_L

HM_L@sam_data <- sample_data(HM_meta_L)
HM_L@sam_data

HM_L_mapping <- HM_L %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by Genus
HM_L_dt <- data.table(HM_L_mapping)
HM_L_dt
HM_L_dt[(Abundance <= 0.04), Genus:= "Other"]
HM_L_dt
HM_L_dt$Genus <- factor(HM_L_dt$Genus, levels= c('Alternaria',
                                                 'Bullera',
                                                 'Cladosporium',
                                                 'Coniothyrium',
                                                 'Didymella',
                                                 "Dioszegia",
                                                 'Fusarium',
                                                 'Hannaella',
                                                 'Neoascochyta',
                                                 'Neosetophoma',
                                                 "Ophiosphaerella",
                                                 "Phaeosphaeria",
                                                 "Phomopsis",
                                                 "Sarocladium",
                                                 "Sporobolomyces",
                                                 'Symmetrospora',
                                                 'Tilletiopsis',
                                                 'Other'))


Bar_L= ggplot(HM_L_dt, aes(x = BAR_LABEL, y = Abundance, fill = Genus)) + 
  #facet_wrap(~Stage, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_x_discrete(limits= c("2018 V2 Conventional",
                             "2018 R2 Conventional",
                             "2018 R6 Conventional",
                             "2018 V2 No-Till",
                             "2018 R2 No-Till",
                             "2018 R6 No-Till",
                             "2018 V2 Organic",
                             "2018 R2 Organic",
                             "2018 R6 Organic",
                             "2021 V2 Conventional",
                             "2021 R3 Conventional",
                             "2021 R6 Conventional",
                             "2021 V2 No-Till",
                             "2021 R3 No-Till",
                             "2021 R6 No-Till",
                             "2021 V2 Organic",
                             "2021 R3 Organic",
                             "2021 R6 Organic"))+
  scale_fill_manual(values = c("Albifimbria" = "#8B0606",
                               "Alternaria" = "#BF0D0D",
                               "Apodus" = "#F50B0B",
                               "Arnium" = "#F49595",
                               "Berkeleyomyces" = "#C06F4C",
                               "Bullera" = "#E24400",
                               "Chloridium" = "#F0711B",
                               "Cladosporium" = "#8E5803",
                               "Coniothyrium" = "#F09300",
                               "Conlarium" = "#FCAF36",
                               "Corynespora" = "#FCCF36",
                               "Dictyochaeta" = "#FFF224",
                               "Didymella" = "#FFEBA6",
                               "Dioszegia" = "#A07C00",
                               "Edenia" = "#685A2A",
                               "Emericellopsis" = "#707B40",
                               "Fusarium" = "#89A705",
                               "Fusidium" = "#B8DF09",
                               "Gaeumannomyces" = "#DDFF45",
                               "Gigaspora" = "#7DFF45",
                               "Glomus" = "#40CB05",
                               "Hannaella" = "#309A04",
                               "Herpotrichia" = "#3A7621",
                               "Knufia" = "#739775",
                               "Macrophomina" = "#6CB99D",
                               "Minimedusa" = "#46DEA6",
                               "Mortierella" = "#0DF5A0",
                               "Myrothecium" = "#07A76C",
                               "Neoascochyta" = "#036642",
                               "Neonectria" = "#0E8C86",
                               "Neosetophoma" = "#0DBCB4",
                               "Oliveonia" = "#13ECE2",
                               "Ophiosphaerella" = "#92F8F3",
                               "Papiliotrema" = "#8DC0BD",
                               "Paraphoma" = "#6F8C9B",
                               "Penicillium" = "#2A5A73",
                               "Phaeosphaeria" = "#1374A6",
                               "Phallus" = "#079AE5",
                               "Phomopsis" = "#53C5FF",
                               "Plectosphaerella" = "#528CDA",
                               "Protocreopsis" = "#1857AD",
                               "Rhizophyctis" = "#0A3D84",
                               "Sarocladium" = "#253A58",
                               "Schizothecium" = "#463278",
                               "Septoglomus" = "#24007C",
                               "Sistotrema" = "#3A00C7",
                               "Sporobolomyces" = "#4A00FF",
                               "straticonidium" = "#9064FC",
                               "Symmetrospora" = "#CD93FE",
                               "Tilletiopsis" = "#FAB5FF",
                               "Trechispora" = "#F56FFE",
                               "Trichoderma" = "#EF01FF",
                               "Tuber" = "#9700A2",
                               "Vishniacozyma" = "#5E0265",
                               "Xenasmatella" = "#593350",
                               "Other" = "grey26"))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 0, face = "bold")) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  ggtitle("leaf")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = .3)) +
  ylab("Relative Abundance (Genera > 4%) \n") +
  xlab("")
plot(Bar_L)

#PLOTTING STACKED_BAR ROOT
ps.noncontam_Rare_meta_R<- sample_data(ps.noncontam_Rare_R)
HM_R<- ps.noncontam_Rare_R %>%
  merge_samples("BAR_LABEL")

HM_meta_R <- as.data.frame(as.matrix(HM_R@sam_data)) %>%
  select(-BAR_LABEL) %>%
  rownames_to_column("BAR_LABEL") %>% 
  select("BAR_LABEL") %>% 
  separate(BAR_LABEL, c("Year", "Growth_Stage","Management"), remove=FALSE) %>%
  mutate(trick = BAR_LABEL) %>% 
  column_to_rownames("trick")

HM_meta_R

HM_R@sam_data <- sample_data(HM_meta_R)
HM_R@sam_data

HM_R_mapping <- HM_R %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by Genus
HM_R_dt <- data.table(HM_R_mapping)
HM_R_dt
HM_R_dt[(Abundance <= 0.04), Genus:= "Other"]
HM_R_dt
HM_R_dt$Genus <- factor(HM_R_dt$Genus, levels= c("Alternaria",
                                                 "Arnium",
                                                 "Berkeleyomyces",
                                                 "Conlarium",
                                                 "Corynespora",
                                                 "Dictyochaeta",
                                                 "Fusarium",
                                                 "Fusidium",
                                                 "Gigaspora",
                                                 "Glomus",
                                                 "Knufia",
                                                 "Macrophomina",
                                                 "Mortierella",
                                                 "Neonectria",
                                                 "Oliveonia",
                                                 "Paraphoma",
                                                 "Penicillium",
                                                 "Plectosphaerella",
                                                 "Septoglomus",
                                                 "Striaticonidium",
                                                 "Trechispora",
                                                 "Trichoderma",
                                                 "Tuber",
                                                 "Other"))

Bar_R= ggplot(HM_R_dt, aes(x = BAR_LABEL, y = Abundance, fill = Genus)) + 
  #facet_wrap(~Stage, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_x_discrete(limits= c("2018 V2 Conventional",
                             "2018 R2 Conventional",
                             "2018 R6 Conventional",
                             "2018 V2 No-Till",
                             "2018 R2 No-Till",
                             "2018 R6 No-Till",
                             "2018 V2 Organic",
                             "2018 R2 Organic",
                             "2018 R6 Organic",
                             "2021 V2 Conventional",
                             "2021 R3 Conventional",
                             "2021 R6 Conventional",
                             "2021 V2 No-Till",
                             "2021 R3 No-Till",
                             "2021 R6 No-Till",
                             "2021 V2 Organic",
                             "2021 R3 Organic",
                             "2021 R6 Organic"))+
  scale_fill_manual(values = c("Albifimbria" = "#8B0606",
                               "Alternaria" = "#BF0D0D",
                               "Apodus" = "#F50B0B",
                               "Arnium" = "#F49595",
                               "Berkeleyomyces" = "#C06F4C",
                               "Bullera" = "#E24400",
                               "Chloridium" = "#F0711B",
                               "Cladosporium" = "#8E5803",
                               "Coniothyrium" = "#F09300",
                               "Conlarium" = "#FCAF36",
                               "Corynespora" = "#FCCF36",
                               "Dictyochaeta" = "#FFF224",
                               "Didymella" = "#FFEBA6",
                               "Dioszegia" = "#A07C00",
                               "Edenia" = "#685A2A",
                               "Emericellopsis" = "#707B40",
                               "Fusarium" = "#89A705",
                               "Fusidium" = "#B8DF09",
                               "Gaeumannomyces" = "#DDFF45",
                               "Gigaspora" = "#7DFF45",
                               "Glomus" = "#40CB05",
                               "Hannaella" = "#309A04",
                               "Herpotrichia" = "#3A7621",
                               "Knufia" = "#739775",
                               "Macrophomina" = "#6CB99D",
                               "Minimedusa" = "#46DEA6",
                               "Mortierella" = "#0DF5A0",
                               "Myrothecium" = "#07A76C",
                               "Neoascochyta" = "#036642",
                               "Neonectria" = "#0E8C86",
                               "Neosetophoma" = "#0DBCB4",
                               "Oliveonia" = "#13ECE2",
                               "Ophiosphaerella" = "#92F8F3",
                               "Papiliotrema" = "#8DC0BD",
                               "Paraphoma" = "#6F8C9D",
                               "Penicillium" = "#2A5A73",
                               "Phaeosphaeria" = "#1374A6",
                               "Phallus" = "#079AE5",
                               "Phomopsis" = "#53C5FF",
                               "Plectosphaerella" = "#528CDA",
                               "Protocreopsis" = "#1857AD",
                               "Rhizophyctis" = "#0A3D84",
                               "Sarocladium" = "#253A58",
                               "Schizothecium" = "#463278",
                               "Septoglomus" = "#24007C",
                               "Sistotrema" = "#3A00C7",
                               "Sporobolomyces" = "#4A00FF",
                               "Striaticonidium" = "#9064FC",
                               "Symmetrospora" = "#CD93FE",
                               "Tilletiopsis" = "#FAB5FF",
                               "Trechispora" = "#F56FFE",
                               "Trichoderma" = "#EF01FF",
                               "Tuber" = "#9700A2",
                               "Vishniacozyma" = "#5E0265",
                               "Xenasmatella" = "#593350",
                               "Other" = "grey26"))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 0, face = "bold")) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  ggtitle("root")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = .5)) +
  ylab("Relative Abundance (Genera > 4%) \n") +
  xlab("")
plot(Bar_R)

#PLOTTING STACKED_BAR SOIL
ps.noncontam_Rare_meta_So<- sample_data(ps.noncontam_Rare_So)
HM_So<- ps.noncontam_Rare_So %>%
  merge_samples("BAR_LABEL")

HM_meta_So <- as.data.frame(as.matrix(HM_So@sam_data)) %>%
  select(-BAR_LABEL) %>%
  rownames_to_column("BAR_LABEL") %>% 
  select("BAR_LABEL") %>% 
  separate(BAR_LABEL, c("Year", "Growth_Stage","Management"), remove=FALSE) %>%
  mutate(trick = BAR_LABEL) %>% 
  column_to_rownames("trick")

HM_meta_So

HM_So@sam_data <- sample_data(HM_meta_So)
HM_So@sam_data

HM_So_mapping <- HM_So %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by Genus
HM_So_dt <- data.table(HM_So_mapping)
HM_So_dt
HM_So_dt[(Abundance <= 0.04), Genus:= "Other"]
HM_So_dt
HM_So_dt$Genus<- factor(HM_So_dt$Genus, levels= c("Albifimbria",
                                                  "Alternaria",
                                                  "Apodus",
                                                  "Arnium",
                                                  "Chloridium",
                                                  "Cladosporium",
                                                  "Coniothyrium",
                                                  "Didymella",
                                                  "Emericellopsis",
                                                  "Fusarium",
                                                  "Gaeumannomyces",
                                                  "Herpotrichia",
                                                  "Minimedusa",
                                                  "Mortierella",
                                                  "Myrothecium",
                                                  "Neosetophoma",
                                                  "Phallus",
                                                  "Phomopsis",
                                                  "Protocreopsis",
                                                  "Rhizophlyctis",
                                                  "Schizothecium",
                                                  "Sistotrema",
                                                  "Xenasmatella",
                                                  "Other"))

Bar_So= ggplot(HM_So_dt, aes(x = BAR_LABEL, y = Abundance, fill = Genus)) + 
  #facet_wrap(~Stage, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_x_discrete(limits= c("2018 V2 Conventional",
                             "2018 R2 Conventional",
                             "2018 R6 Conventional",
                             "2018 V2 No-Till",
                             "2018 R2 No-Till",
                             "2018 R6 No-Till",
                             "2018 V2 Organic",
                             "2018 R2 Organic",
                             "2018 R6 Organic",
                             "2021 V2 Conventional",
                             "2021 R3 Conventional",
                             "2021 R6 Conventional",
                             "2021 V2 No-Till",
                             "2021 R3 No-Till",
                             "2021 R6 No-Till",
                             "2021 V2 Organic",
                             "2021 R3 Organic",
                             "2021 R6 Organic"))+
  scale_fill_manual(values = c("Albifimbria" = "#8B0606",
                               "Alternaria" = "#BF0D0D",
                               "Apodus" = "#F50B0B",
                               "Arnium" = "#F49595",
                               "Berkeleyomyces" = "#C06F4C",
                               "Bullera" = "#E24400",
                               "Chloridium" = "#F0711B",
                               "Cladosporium" = "#8E5803",
                               "Coniothyrium" = "#F09300",
                               "Conlarium" = "#FCAF36",
                               "Corynespora" = "#FCCF36",
                               "Dictyochaeta" = "#FFF224",
                               "Didymella" = "#FFEBA6",
                               "Dioszegia" = "#A07C00",
                               "Edenia" = "#685A2A",
                               "Emericellopsis" = "#707B40",
                               "Fusarium" = "#89A705",
                               "Fusidium" = "#B8DF09",
                               "Gaeumannomyces" = "#DDFF45",
                               "Gigaspora" = "#7DFF45",
                               "Glomus" = "#40CB05",
                               "Hannaella" = "#309A04",
                               "Herpotrichia" = "#3A7621",
                               "Knufia" = "#739775",
                               "Macrophomina" = "#6CB99D",
                               "Minimedusa" = "#46DEA6",
                               "Mortierella" = "#0DF5A0",
                               "Myrothecium" = "#07A76C",
                               "Neoascochyta" = "#036642",
                               "Neonectria" = "#0E8C86",
                               "Neosetophoma" = "#0DBCB4",
                               "Oliveonia" = "#13ECE2",
                               "Ophiosphaerella" = "#92F8F3",
                               "Papiliotrema" = "#8DC0BD",
                               "Paraphoma" = "#6F8C9D",
                               "Penicillium" = "#2A5A73",
                               "Phaeosphaeria" = "#1374A6",
                               "Phallus" = "#079AE5",
                               "Phomopsis" = "#53C5FF",
                               "Plectosphaerella" = "#528CDA",
                               "Protocreopsis" = "#1857AD",
                               "Rhizophlyctis" = "#0A3D84",
                               "Sarocladium" = "#253A58",
                               "Schizothecium" = "#463278",
                               "Septoglomus" = "#24007C",
                               "Sistotrema" = "#3A00C7",
                               "Sporobolomyces" = "#4A00FF",
                               "Striaticonidium" = "#9064FC",
                               "Symmetrospora" = "#CD93FE",
                               "Tilletiopsis" = "#FAB5FF",
                               "Trechispora" = "#F56FFE",
                               "Trichoderma" = "#EF01FF",
                               "Tuber" = "#9700A2",
                               "Vishniacozyma" = "#5E0265",
                               "Xenasmatella" = "#593350",
                               "Other" = "grey26"))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 0, face = "bold")) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  ggtitle("soil")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = .5)) +
  ylab("Relative Abundance (Genera > 4%) \n") +
  xlab("")
plot(Bar_So)

#stem
#PLOTTING STACKED_BAR STEM
ps.noncontam_Rare_meta_St<- sample_data(ps.noncontam_Rare_St)
HM_St<- ps.noncontam_Rare_St %>%
  merge_samples("BAR_LABEL")

HM_meta_St <- as.data.frame(as.matrix(HM_So@sam_data)) %>%
  select(-BAR_LABEL) %>%
  rownames_to_column("BAR_LABEL") %>% 
  select("BAR_LABEL") %>% 
  separate(BAR_LABEL, c("Year", "Growth_Stage","Management"), remove=FALSE) %>%
  mutate(trick = BAR_LABEL) %>% 
  column_to_rownames("trick")

HM_meta_St

HM_St@sam_data <- sample_data(HM_meta_St)
HM_St@sam_data

HM_St_mapping <- HM_St %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by Genus
HM_St_dt <- data.table(HM_St_mapping)
HM_St_dt
HM_St_dt[(Abundance <= 0.04), Genus:= "Other"]
HM_St_dt
HM_St_dt$Genus<- factor(HM_St_dt$Genus, levels= c("Alternaria",
                                                  "Bullera",
                                                  "Cladosporium",
                                                  "Didymella",
                                                  "Edenia",
                                                  "Fusarium",
                                                  "Hannaella",
                                                  "Neoascochyta",
                                                  "Neosetophoma",
                                                  "Papiliotrema",
                                                  "Phomopsis",
                                                  "Plectosphaerella",
                                                  "Sarocladium",
                                                  "Symmetrospora",
                                                  "Tilletiopsis",
                                                  "Vishniacozyma",
                                                  "Other"))

Bar_St= ggplot(HM_St_dt, aes(x = BAR_LABEL, y = Abundance, fill = Genus)) + 
  #facet_wrap(~Stage, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_x_discrete(limits= c("2018 V2 Conventional",
                             "2018 R2 Conventional",
                             "2018 R6 Conventional",
                             "2018 V2 No-Till",
                             "2018 R2 No-Till",
                             "2018 R6 No-Till",
                             "2018 V2 Organic",
                             "2018 R2 Organic",
                             "2018 R6 Organic",
                             "2021 V2 Conventional",
                             "2021 R3 Conventional",
                             "2021 R6 Conventional",
                             "2021 V2 No-Till",
                             "2021 R3 No-Till",
                             "2021 R6 No-Till",
                             "2021 V2 Organic",
                             "2021 R3 Organic",
                             "2021 R6 Organic"))+
  scale_fill_manual(values = c("Albifimbria" = "#8B0606",
                               "Alternaria" = "#BF0D0D",
                               "Apodus" = "#F50B0B",
                               "Arnium" = "#F49595",
                               "Berkeleyomyces" = "#C06F4C",
                               "Bullera" = "#E24400",
                               "Chloridium" = "#F0711B",
                               "Cladosporium" = "#8E5803",
                               "Coniothyrium" = "#F09300",
                               "Conlarium" = "#FCAF36",
                               "Corynespora" = "#FCCF36",
                               "Dictyochaeta" = "#FFF224",
                               "Didymella" = "#FFEBA6",
                               "Dioszegia" = "#A07C00",
                               "Edenia" = "#685A2A",
                               "Emericellopsis" = "#707B40",
                               "Fusarium" = "#89A705",
                               "Fusidium" = "#B8DF09",
                               "Gaeumannomyces" = "#DDFF45",
                               "Gigaspora" = "#7DFF45",
                               "Glomus" = "#40CB05",
                               "Hannaella" = "#309A04",
                               "Herpotrichia" = "#3A7621",
                               "Knufia" = "#739775",
                               "Macrophomina" = "#6CB99D",
                               "Minimedusa" = "#46DEA6",
                               "Mortierella" = "#0DF5A0",
                               "Myrothecium" = "#07A76C",
                               "Neoascochyta" = "#036642",
                               "Neonectria" = "#0E8C86",
                               "Neosetophoma" = "#0DBCB4",
                               "Oliveonia" = "#13ECE2",
                               "Ophiosphaerella" = "#92F8F3",
                               "Papiliotrema" = "#8DC0BD",
                               "Paraphoma" = "#6F8C9D",
                               "Penicillium" = "#2A5A73",
                               "Phaeosphaeria" = "#1374A6",
                               "Phallus" = "#079AE5",
                               "Phomopsis" = "#53C5FF",
                               "Plectosphaerella" = "#528CDA",
                               "Protocreopsis" = "#1857AD",
                               "Rhizophlyctis" = "#0A3D84",
                               "Sarocladium" = "#253A58",
                               "Schizothecium" = "#463278",
                               "Septoglomus" = "#24007C",
                               "Sistotrema" = "#3A00C7",
                               "Sporobolomyces" = "#4A00FF",
                               "Striaticonidium" = "#9064FC",
                               "Symmetrospora" = "#CD93FE",
                               "Tilletiopsis" = "#FAB5FF",
                               "Trechispora" = "#F56FFE",
                               "Trichoderma" = "#EF01FF",
                               "Tuber" = "#9700A2",
                               "Vishniacozyma" = "#5E0265",
                               "Xenasmatella" = "#593350",
                               "Other" = "grey26"))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 0, face = "bold")) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  ggtitle("stem")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = .5)) +
  ylab("Relative Abundance (Genera > 4%) \n") +
  xlab("")
plot(Bar_St)

ggarrange(Bar_L, Bar_R, Bar_So, Bar_St,
          labels= c("A", "B", "C", "D"),
          align = "h", ncol = 2, nrow = 2)

###Beta Diversity
merged_phyloseq<- merge_phyloseq(ps.noncontam_L,ps.noncontam_R,ps.noncontam_So, ps.noncontam_St)
merged_phyloseq<- merge_phyloseq(ps.noncontam_L,ps.noncontam_R,ps.noncontam_So, ps.noncontam_St)
merged_phyloseq_norm <- transform_sample_counts(merged_phyloseq, function(x) x / sum(x))

management_colors<- c("Conventional" = "#feffbe", "No-Till"="#feff40", "Organic"="#06a651")
growth_stage_colros<- c("V2" = "#56B4E9", "R2" = "#FFB5B8", "R3" = "#009E73", "R6" = "#E69F00")
collection_year_colors<-c("2018"="darkslateblue", "2021"="grey65")

pcoa_man <- ordinate(merged_phyloseq_norm, method = "PCoA", distance = "bray")
p_man <- plot_ordination(merged_phyloseq, pcoa, aes(color = MANAGEMENT)) +
  geom_point(size = 3, shape = 21, aes(fill = MANAGEMENT)) +
  scale_color_manual(values = management_colors) +
  scale_fill_manual(values = management_colors) +
  theme_minimal() +
  theme(text = element_text(size = 18),
        legend.position = "right") +
  labs(title = "PCoA plot with Bray-Curtis distance, Management",x = "PCoA 1",y = "PCoA 2")
print(p_man)

p_origin <- plot_ordination(merged_phyloseq, pcoa, aes(color = ORIGIN)) +
  geom_point(size = 3, shape = 21, aes(fill = ORIGIN)) +
  theme_minimal() +
  theme(text = element_text(size = 18),
        legend.position = "right") +
  labs(title = "PCoA plot with Bray-Curtis distance, Origin",x = "PCoA 1",y = "PCoA 2")
print(p_origin)

p_gs <- plot_ordination(merged_phyloseq, pcoa, aes(color = GROWTH_STAGE)) +
  geom_point(size = 3, shape = 21, aes(fill = GROWTH_STAGE)) +
  scale_color_manual(values = growth_stage_colros) +
  scale_fill_manual(values = growth_stage_colros) +
  theme_minimal() +
  theme(text = element_text(size = 18),
        legend.position = "right") +
  labs(title = "PCoA plot with Bray-Curtis distance, Management",x = "PCoA 1",y = "PCoA 2")
print(p_gs)

p_cy <- plot_ordination(merged_phyloseq, pcoa, aes(color = COLLECTION_YEAR)) +
  geom_point(size = 3, shape = 21, aes(fill = COLLECTION_YEAR)) +
  scale_color_manual(values = collection_year_colors) +
  scale_fill_manual(values = collection_year_colors) +
  theme_minimal() +
  theme(text = element_text(size = 18),
        legend.position = "right") +
  labs(title = "PCoA plot with Bray-Curtis distance, Collection_Year",x = "PCoA 1",y = "PCoA 2")
print(p_cy)

PCOA_Combined<- ggarrange(p_man,p_origin,p_gs,p_cy)

###FUNGAL TRAITS
setwd("~/Desktop/Microbiome/MSA23/fungaltraits")

data<- read.delim("Overall_fungaltraits.txt", header= TRUE)

data_leaf<- filter(data, Origin == "leaf")
data_leaf$Primary_guild_unique <- reorder(data_leaf$Primary_guild_unique, -data_leaf$Primary_guild_percent)
data_leaf <- filter(data_leaf, Primary_guild_percent >= 1)
data_root<- filter(data, Origin == "root")
data_root$Primary_guild_unique <- reorder(data_root$Primary_guild_unique, -data_root$Primary_guild_percent)
data_root <- filter(data_root, Primary_guild_percent >= 1)
data_soil<- filter(data, Origin == "soil")
data_soil$Primary_guild_unique <- reorder(data_soil$Primary_guild_unique, -data_soil$Primary_guild_percent)
data_soil <- filter(data_soil, Primary_guild_percent >= 1)
data_stem<- filter(data, Origin == "stem")
data_stem$Primary_guild_unique <- reorder(data_stem$Primary_guild_unique, -data_stem$Primary_guild_percent)
data_stem <- filter(data_stem, Primary_guild_percent >= 1)
library(ggplot2)

# Create the bar plot
leaf_bar<- ggplot(data_leaf, aes(x = Primary_guild_unique, y = Primary_guild_percent, fill = Management)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("Organic" = "#06a651", "Conventional" = "#feffbe", "No-Till" = "#feff40")) +
  theme_minimal() +
  xlab("Primary Guild") +
  ylab("Percentage") +
  ggtitle("Leaves")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(breaks = seq(0, 40, by = 2))+
  guides(fill = FALSE)
leaf_bar

root_bar<- ggplot(data_root, aes(x = Primary_guild_unique, y = Primary_guild_percent, fill = Management)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("Organic" = "#06a651", "Conventional" = "#feffbe", "No-Till" = "#feff40")) +
  theme_minimal() +
  xlab("Primary Guild") +
  ylab("Percentage") +
  ggtitle("Roots")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(breaks = seq(0, 50, by = 2))+
  guides(fill = FALSE)
root_bar

soil_bar<- ggplot(data_soil, aes(x = Primary_guild_unique, y = Primary_guild_percent, fill = Management)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("Organic" = "#06a651", "Conventional" = "#feffbe", "No-Till" = "#feff40")) +
  theme_minimal() +
  xlab("Primary Guild") +
  ylab("Percentage") +
  ggtitle("Soil")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(breaks = seq(0, 55, by = 2))+
  guides(fill = FALSE)
soil_bar

stem_bar<- ggplot(data_stem, aes(x = Primary_guild_unique, y = Primary_guild_percent, fill = Management)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("Organic" = "#06a651", "Conventional" = "#feffbe", "No-Till" = "#feff40")) +
  theme_minimal() +
  xlab("Primary Guild") +
  ylab("Percentage") +
  ggtitle("Stems")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(breaks = seq(0, 40, by = 2))
stem_bar

library(ggpubr)
ggarrange(leaf_bar,root_bar,soil_bar,stem_bar,
          labels = c("A", "B","C","D"),
          align = c("h"),
          ncol = 4, nrow = 1)
