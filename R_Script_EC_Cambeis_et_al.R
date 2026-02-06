#This is the Script for the analysis of 16S rRNA data in the manuscript: 
#Bacterial quorum sensing molecule functions as postbiotics and modifies barley physiology in a microbial community-dependent manner
#By Cambeis et al. 2026
#Submitted to Enivormental Microbiome

#Path to data have to be adjusted!!


# DADA2 -------------------------------------------------------------------
library(dada2)
packageVersion("dada2")

# Set Path
data_path <- "~/CleanData_fq"
list.files(data_path)

# List all forward and reverse read files
fnFs <- sort(list.files(data_path, pattern="_1.fq", full.names = TRUE))
fnRs <- sort(list.files(data_path, pattern="_2.fq", full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_1.fq"), `[`, 1)

# Inspect read quality profiles
plotQualityProfile(fnFs[1:4])  # Adjust the number of plots as needed
plotQualityProfile(fnRs[1:4])

# Define filtered output file paths
filtFs <- file.path(data_path, "filtered", paste0(sample.names, "_1_filt.fastq.gz"))
filtRs <- file.path(data_path, "filtered", paste0(sample.names, "_2_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

# Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.table(track, "16S_ReadProcessing_Tracking.txt", sep="\t", quote=FALSE)

# Assign taxonomy using SILVA database
taxa_silva <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)
write.table(taxa_silva, "16S_Taxonomy_Silva.txt", sep="\t", quote=FALSE)

# Save ASV table
asvtable_count <- t(seqtab.nochim)
head(asvtable_count)
write.table(asvtable_count, "16S_ASVcount.txt", sep="\t", quote=FALSE)


asvtable_count <- read.table("16S_ASVcount.txt", header=T, sep="\t", dec=".", row.names=1 ) # samples in columns, no tax info; check.names false to accept my sample names
dim(asvtable_count)
ncol(asvtable_count)
nrow(asvtable_count)
colnames(asvtable_count)
rownames(asvtable_count)[1:5]  # Show first 5 ASV IDs
sum(is.na(asvtable_count))




#Phyloseq----
###create phyloseq objects from files----
pacman::p_load(
  broom,
  conflicted,
  here,
  janitor,
  writexl, 
  naniar,
  readxl,
  tidyverse,
  corrr,
  MetBrewer,
  pheatmap,
  ggpattern,
  psych,
  conflicted,
  desplot,
  emmeans,
  ggtext,
  MetBrewer,
  multcomp,
  multcompView,
  conflicted,
  here,
  naniar,
  readxl,
  rstatix,
  ggpubr,
  ggtext,
  rstatix, 
  corrr,
  vegan,
  ggpmisc,
  factoextra,
  ggfortify, 
  openxlsx2, 
  knitr,
  reshape, 
  viridis,
  microbiome,
  phyloseq
  
)


# handle function conflicts
conflict_prefer("filter", "dplyr") 
conflict_prefer("select", "dplyr")

out_dir<-"/Out" # Adjust!


##### phyloseq object ---------------------------------------------------------

#ASV-count-table
asv_count <- read.table("16S_ASVcount.txt", header=T, sep="\t", dec=".", row.names=1 ) ##incl.ASV-names AND sample-names
#Check dimensions
dim(asv_count)
ncol(asv_count)
nrow(asv_count)
#Check names of columns and row --> have to match to metadata and tax table
colnames(asv_count)
rownames(asv_count)[1:5]  # Show first 5 ASV IDs
sum(is.na(asv_count))

#Taxonomy-table with asv-names matching asv_count
tax <- read.table("16S_Taxonomy_Silva.txt", header=T, sep="\t", dec=".", row.names=1 ) ##incl.ASV-names
#Check if import was successfull and dataframe is correct
ncol(tax)
nrow(tax)
colnames(tax)
rownames(tax)[1:5] 

#Metadata
metadata <- read.table("Metadata.txt", header=T, sep="\t", dec=".", row.names=1, fileEncoding="ISO-8859-1")
#Check if import was successfull and dataframe is formated as the others
ncol(metadata)
nrow(metadata)
colnames(metadata)
rownames(metadata)[1:5] 

#combine to phyloseq object
ASV = otu_table(asv_count, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax))
META = sample_data(metadata)

ASV
TAX
META

#    ##   MAIN phyloseq object
physeq = phyloseq(ASV,TAX,META)
physeq


library(phyloseqCompanion)
physeq_numbered <- numbered.ASVs(physeq)
head(taxa_names(physeq_numbered))
head(taxa_names(physeq))


ntaxa(physeq)
ntaxa(physeq_numbered)

nsamples(physeq)
nsamples(physeq_numbered)


sum(otu_table(physeq))
sum(otu_table(physeq_numbered))

dim(tax_table(physeq))
dim(tax_table(physeq_numbered))

str(physeq)
str(physeq_numbered)
#überschreiben der originalen phyloseq datei
physeq<- numbered.ASVs(physeq)

# Clean phyloseq ----------------------------------------------------------

##Basis of the script is by Hauschild et. al https://github.com/krishauschi/ORDIAmur_ms2/tree/main/scripts/microbiome/16S_rRNA_amplicon
#Modification by MattCambeis 

#Visualize data
sample_names(physeq)
rank_names(physeq)
sample_variables(physeq)


#clean dataset and remove ASVs identified as Chloroplast or Mitochondria 

##Create table with number of features for each phylum
##Remove respective artefacts on a specific taxonomic level (Keeping Archaea)
#Kingdom
table(tax_table(physeq)[,"Kingdom"], exclude = NULL)
physeq <- subset_taxa(physeq, !is.na(Kingdom) & !Kingdom %in% c("Unassigned", "Eukaryota", "NA", ""))
table(tax_table(physeq)[,"Kingdom"], exclude = NULL)

#Phylum
table(tax_table(physeq)[,"Phylum"], exclude = NULL)
physeq <- subset_taxa(physeq, !Phylum %in% c("Bacteria", "Archaea"))
table(tax_table(physeq)[,"Phylum"], exclude = NULL)

#Order
table(tax_table(physeq)[,"Order"], exclude = NULL)
physeq <- subset_taxa(physeq, !is.na(Order) & !Order %in% c("Chloroplast"))
table(tax_table(physeq)[,"Order"], exclude = NULL)

#Family
table(tax_table(physeq)[,"Family"], exclude = NULL)
physeq <- subset_taxa(physeq, !is.na(Family) & !Family %in% c("Mitochondria"))
table(tax_table(physeq)[,"Family"], exclude = NULL)

#Genus
#table(tax_table(psO_WP3_16S)[,"Genus"], exclude = NULL)

#Annotation
#table(tax_table(psO_WP3_16S)[,"Annotation"], exclude = NULL)

##Observe psO after clean spurious taxa
physeq

#Compute prevalence of each feature and store as data.frame
prevdf = apply(X = otu_table(physeq),
               MARGIN = ifelse(taxa_are_rows(physeq), yes = 1, no = 2),
               FUN = function (x) {sum(x>0)})

#Add taxonomy and total reads counts to the data frame
prevdf = data.frame(Prevalence=prevdf, TotalAbundance = taxa_sums(physeq), tax_table(physeq), otu_table(physeq))
#Write data frame in csv format
write.csv(prevdf, "prevdf_physeq.csv")


cat("Number of ASVs:", ntaxa(physeq), "\n")
cat("Number of samples:", nsamples(physeq), "\n")
cat("Average reads per sample:", mean(sample_sums(physeq)), "\n")


# Rename_NAs --------------------------------------------------------------


#now we have a phyloseq object that was cleaned from missassigned taxa. As next step, go to "16S_rename_NAs" to fill missing taxonimic information.

#input is the phyloseq object that was created using the script "16S_create_phyloseq_object.R"
#Load packages
library("phyloseq")
library("stringr")


# Export tax table
tax <- data.frame(phyloseq::tax_table(physeq))

# Change NA to empty string
tax[is.na(tax)] <- ""

# Ensure all columns are characters
tax.clean <- as.data.frame(lapply(tax, as.character), stringsAsFactors = FALSE)
rownames(tax.clean) <- rownames(tax)

# Fill missing taxonomy hierarchically
for (i in 1:nrow(tax.clean)) {
  if (tax.clean[i, "Phylum"] == "") {
    tax.clean[i, 2:7] <- tax.clean[i, "Kingdom"]
  } else if (tax.clean[i, "Class"] == "") {
    tax.clean[i, 3:7] <- tax.clean[i, "Phylum"]
  } else if (tax.clean[i, "Order"] == "") {
    tax.clean[i, 4:7] <- tax.clean[i, "Class"]
  } else if (tax.clean[i, "Family"] == "") {
    tax.clean[i, 5:7] <- tax.clean[i, "Order"]
  } else if (tax.clean[i, "Genus"] == "") {
    tax.clean[i, 6:7] <- tax.clean[i, "Family"]
  } else if (tax.clean[i, "Species"] == "") {
    tax.clean[i, "Species"] <- tax.clean[i, "Genus"]
  }
}


head(tax.clean, 10) 
tax.clean[tax.clean$Species == "", ]

# View the cleaned taxonomy table as a well-formatted table
library(knitr)
kable(head(tax.clean, 20))  # Shows the first 20 rows of the cleaned taxonomy


#Return data.frame to a phyloseq object
phyloseq::tax_table(physeq) <- as.matrix(tax.clean)
head(phyloseq::tax_table(physeq))
tail(phyloseq::tax_table(physeq))
physeq

#this is the phyloseq object with taxonomy information for all ranks. It includes raw count data and will be the input object for differential abundance analysis. For this, go to script "16S_ANCOM-BC2.R"


# Rarefaction -------------------------------------------------------------


#Create rarefied data and rarefaction curve using the phyloseq object from script "ITS_rename_NAs.R"  

library("phyloseq")
library("dplyr")
library("tidyverse")
library("vegan")

#rarefaction
otu.rare = otu_table(physeq)
otu.rare = as.data.frame(t(otu.rare))
sample_names = rownames(otu.rare)
mindata<- min(rowSums(otu.rare))
mindata # smallest sample size = 3534 <- kommt von Bulk soil t0 --> wird rausgeschmissen!
#rarefaction curve
otu.rarecurve = rarecurve(otu.rare, step = 1000, label = T, xlim=c(0, 60000))


##Check if kick out 0 days samples; mit CHatGPT gemacht
total_reads <- sum(otu.rare)  # Replace `otu_table` with your OTU/ASV matrix
total_reads
n_asvs <- sum(colSums(otu.rare) > 0)
n_asvs
median_reads <- median(rowSums(otu.rare))
median_reads

library(vegan)

rarecurve(otu.rare, step = 1000, cex = 0.6, label = TRUE)

reads_per_sample <- rowSums(otu.rare)
reads_per_sample_sorted <- sort(reads_per_sample, decreasing = TRUE)
# 
# Plot for a quick look
hist(reads_per_sample,  col = "lightgreen", main = "reads per Sample", xlab = "Number of Reads")



low_samples <- reads_per_sample[reads_per_sample < 10000]  # or 10,000, depending on context

# View sample names and read counts
low_samples #Due to low sequencing quality the samples 81 - 88 were removed from the data set. 

# > low_samples
# RemovePrimer_Final.EC_81 RemovePrimer_Final.EC_82 RemovePrimer_Final.EC_83 RemovePrimer_Final.EC_84 RemovePrimer_Final.EC_85 
# 6910                     7164                     6872                     7069                     3562 
# RemovePrimer_Final.EC_86 RemovePrimer_Final.EC_87 RemovePrimer_Final.EC_88 
# 3744                     3534                     3772 
# 

asvs_per_sample <- rowSums(otu.rare > 0)

# Plot for a quick look
hist(asvs_per_sample, breaks = 30, col = "lightgreen", main = "ASVs per Sample", xlab = "Number of ASVs")
###CHECK 0 days Bulk soil kicked out

# rarefy without replacement
set.seed(2608)
physeq_raref = rarefy_even_depth(physeq, rngseed=1, sample.size=21707 
                                 , replace=F) #rarefaction on min. sample size =21707 reads, ist die niedrigste Probe wenn bulk soil t0 rausgeschmissen!

###export data frame
df_physeq_raref <- data.frame(tax_table(physeq_raref),otu_table(physeq_raref)) ### Create tables
write.csv(df_physeq_raref, "df_physeq_raref.csv")

df_physeq_raref
#this is the phyloseq object containing count data that will be used for all subsequent analyses, except for differential abundance testing.
#For the next steps, go to the diversity scripts "16S_shannon.Rmd" for alpha-diversity and "16S_ordination.R" for beta-diversity.



# MDS -----------------------------------------------------

library(patchwork)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(forcats)


#JKI
# Subset rhizosphere samples
physeq_raref_JKI <- subset_samples(physeq_raref, Soil == "Field soil")
physeq_raref_JKI_wo_bulk <- subset_samples(physeq_raref_JKI, Microhabitat == "RZ")

# Remove taxa that are zero *within the subset*
physeq_raref_JKI_wo_bulk <- prune_taxa(
  taxa_sums(physeq_raref_JKI_wo_bulk) > 0,
  physeq_raref_JKI_wo_bulk
)



# Extract sample_data as a data.frame
samp_df <- data.frame(sample_data(physeq_raref_JKI_wo_bulk))

# Make sure the relevant columns are factors
samp_df$Treatment <- as.factor(samp_df$Treatment)
samp_df$Genotype  <- as.factor(samp_df$Genotype)
samp_df$Duration  <- as.factor(samp_df$Duration)

# Reorder factor levels
samp_df <- samp_df %>%
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Duration = fct_relevel(Duration, "Three-week-old", "18-week-old")
  )

# Put it back into phyloseq object
sample_data(physeq_raref_JKI_wo_bulk) <- sample_data(samp_df)



# Define colors for each Genotype × Treatment combination
combo_colors <- c(
  "Golden_Promise.mock" = "#A67C00",
  "Golden_Promise.C6-HSL" = "#ffcf40",
  "BCC436.mock" = "#4F7942",
  "BCC436.C6-HSL" = "#4aba91"
)

# Define shapes for genotypes
genotype_shapes <- c(
  "BCC436" = 16,
  "Golden_Promise" = 17
)

# Get unique durations
durations <- unique(sample_data(physeq_raref_JKI_wo_bulk)$Duration)

# Initialize list to store plots
mds_plots <- list()

for (d in durations) {
  # Subset by duration
  physeq_sub <- subset_samples(physeq_raref_JKI_wo_bulk, Duration == d)
  
  # Create combination factor for color mapping
  sample_data(physeq_sub)$Genotype_Treatment <- 
    paste(sample_data(physeq_sub)$Genotype,
          sample_data(physeq_sub)$Treatment,
          sep = ".")
  
  
  sample_data(physeq_sub)$Genotype_Treatment <- 
    factor(sample_data(physeq_sub)$Genotype_Treatment,
           levels = c("Golden_Promise.mock", 
                      "Golden_Promise.C6-HSL",
                      "BCC436.mock", 
                      "BCC436.C6-HSL"))
  
  # Run MDS
  ord <- ordinate(physeq_sub, method = "NMDS", distance = "bray")
  
  # # Optional: percent variation explained
  # eig <- ord$values$Eigenvalues
  # percent_explained <- round(eig / sum(eig) * 100, 1)
  # xlab <- paste0("Axis 1 (", percent_explained[1], "%)")
  # ylab <- paste0("Axis 2 (", percent_explained[2], "%)")
  
  # Plot
  p <- plot_ordination(physeq_sub, ord, type = "samples", 
                       shape = "Genotype", color = "Genotype_Treatment") +
    geom_point(size = 4, alpha = 0.8) +
    stat_ellipse(aes(group = Genotype_Treatment, fill = Genotype_Treatment), 
                 geom = "polygon", alpha = 0.2) +
    scale_color_manual(values = combo_colors,
                       labels = c(
                         "Golden_Promise.mock" = "Golden Promise + mock",
                         "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
                         "BCC436.mock" = "BCC436 + mock",
                         "BCC436.C6-HSL" = "BCC436 + C6-HSL"
                       )) +
    scale_fill_manual(values = combo_colors,
                      labels = c(
                        "Golden_Promise.mock" = "Golden Promise + mock",
                        "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
                        "BCC436.mock" = "BCC436 + mock",
                        "BCC436.C6-HSL" = "BCC436 + C6-HSL"
                      )) +
    scale_shape_manual(values = genotype_shapes,labels = c(
      "Golden_Promise" = "Golden Promise",
      "BCC436" = "BCC436"
    )) +
    labs(
      title = paste("MDS Plot - Rhizosphere - Duration:", d),
      x = xlab,
      y = ylab,
      color = "Genotype & Treatment",
      fill = "Genotype & Treatment",
      shape = "Genotype"
    ) +
    my_theme
  # Store plot in list
  mds_plots[[as.character(d)]] <- p
}



combined_plot <- wrap_plots(mds_plots, ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom"))

combined_plot





# Initialize list to store plots
mds_plots <- list()

for (d in durations) {
  # Subset by duration
  physeq_sub <- subset_samples(physeq_raref_JKI_wo_bulk, Duration == d)
  
  # Create combination factor for color mapping
  sample_data(physeq_sub)$Genotype_Treatment <- 
    paste(sample_data(physeq_sub)$Genotype,
          sample_data(physeq_sub)$Treatment,
          sep = ".")
  
  sample_data(physeq_sub)$Genotype_Treatment <- 
    factor(sample_data(physeq_sub)$Genotype_Treatment,
           levels = c("Golden_Promise.mock", 
                      "Golden_Promise.C6-HSL",
                      "BCC436.mock", 
                      "BCC436.C6-HSL"))
  
  # Run NMDS
  ord <- ordinate(physeq_sub, method = "NMDS", distance = "bray")
  
  # Extract stress
  stress_val <- round(ord$stress, 3)
  stress_label <- paste0("Stress = ", stress_val)
  
  # Plot
  p <- plot_ordination(physeq_sub, ord, type = "samples", 
                       shape = "Genotype", color = "Genotype_Treatment") +
    geom_point(size = 4, alpha = 0.8) +
    stat_ellipse(aes(group = Genotype_Treatment, fill = Genotype_Treatment), 
                 geom = "polygon", alpha = 0.2) +
    scale_color_manual(values = combo_colors,
                       labels = c(
                         "Golden_Promise.mock" = "Golden Promise + mock",
                         "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
                         "BCC436.mock" = "BCC436 + mock",
                         "BCC436.C6-HSL" = "BCC436 + C6-HSL"
                       )) +
    scale_fill_manual(values = combo_colors,
                      labels = c(
                        "Golden_Promise.mock" = "Golden Promise + mock",
                        "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
                        "BCC436.mock" = "BCC436 + mock",
                        "BCC436.C6-HSL" = "BCC436 + C6-HSL"
                      )) +
    scale_shape_manual(values = genotype_shapes,
                       labels = c(
                         "Golden_Promise" = "Golden Promise",
                         "BCC436" = "BCC436"
                       )) +
    labs(
      title = paste("MDS Plot - Rhizosphere - Duration:", d),
      subtitle = stress_label,
      x = xlab,
      y = ylab,
      color = "Genotype & Treatment",
      fill = "Genotype & Treatment",
      shape = "Genotype"
    ) +
    my_theme 
  
  # Store plot in list
  mds_plots[[as.character(d)]] <- p
}



combined_plot <- wrap_plots(mds_plots, ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom"))

combined_plot



file_name <- paste0("MDS/MDS_JKI.png")

ggsave(
  filename = file_name,
  plot = last_plot(),
  width = 22,       # width in cm
  height = 15,      # height in cm
  units = "cm",     # specify centimeters
  dpi = 300,
  path = out_dir# resolution (300 dpi = publication quality)
)




# FE

# Subset rhizosphere samples
physeq_raref_FE <- subset_samples(physeq_raref, Soil == "Substrate")
physeq_raref_FE_wo_bulk <- subset_samples(physeq_raref_FE, Microhabitat == "RZ")

# Remove taxa that are zero *within the subset*
physeq_raref_FE_wo_bulk <- prune_taxa(
  taxa_sums(physeq_raref_FE_wo_bulk) > 0,
  physeq_raref_FE_wo_bulk
)


# Extract sample_data as a data.frame
samp_df <- data.frame(sample_data(physeq_raref_FE_wo_bulk))

# Make sure the relevant columns are factors
samp_df$Treatment <- as.factor(samp_df$Treatment)
samp_df$Genotype  <- as.factor(samp_df$Genotype)
samp_df$Duration  <- as.factor(samp_df$Duration)

# Reorder factor levels
samp_df <- samp_df %>%
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Duration = fct_relevel(Duration, "Two-week-old", "16-week-old")
  )

# Put it back into phyloseq object
sample_data(physeq_raref_FE_wo_bulk) <- sample_data(samp_df)



# Define colors for each Genotype × Treatment combination
combo_colors <- c(
  "Golden_Promise.mock" = "#A67C00",
  "Golden_Promise.C6-HSL" = "#ffcf40",
  "BCC436.mock" = "#4F7942",
  "BCC436.C6-HSL" = "#4aba91"
)

# Define shapes for genotypes
genotype_shapes <- c(
  "BCC436" = 16,
  "Golden_Promise" = 17
)

# Get unique durations
durations <- unique(sample_data(physeq_raref_FE_wo_bulk)$Duration)

# Initialize list to store plots
mds_plots <- list()

for (d in durations) {
  # Subset by duration
  physeq_sub <- subset_samples(physeq_raref_FE_wo_bulk, Duration == d)
  
  # Create combination factor for color mapping
  sample_data(physeq_sub)$Genotype_Treatment <- 
    paste(sample_data(physeq_sub)$Genotype,
          sample_data(physeq_sub)$Treatment,
          sep = ".")
  
  
  sample_data(physeq_sub)$Genotype_Treatment <- 
    factor(sample_data(physeq_sub)$Genotype_Treatment,
           levels = c("Golden_Promise.mock", 
                      "Golden_Promise.C6-HSL",
                      "BCC436.mock", 
                      "BCC436.C6-HSL"))
  
  # Run MDS
  ord <- ordinate(physeq_sub, method = "NMDS", distance = "bray")
  
  # # Optional: percent variation explained
  # eig <- ord$values$Eigenvalues
  # percent_explained <- round(eig / sum(eig) * 100, 1)
  # xlab <- paste0("Axis 1 (", percent_explained[1], "%)")
  # ylab <- paste0("Axis 2 (", percent_explained[2], "%)")
  # 
  # Plot
  p <- plot_ordination(physeq_sub, ord, type = "samples", 
                       shape = "Genotype", color = "Genotype_Treatment") +
    geom_point(size = 4, alpha = 0.8) +
    stat_ellipse(aes(group = Genotype_Treatment, fill = Genotype_Treatment), 
                 geom = "polygon", alpha = 0.2) +
    scale_color_manual(values = combo_colors,
                       labels = c(
                         "Golden_Promise.mock" = "Golden Promise + mock",
                         "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
                         "BCC436.mock" = "BCC436 + mock",
                         "BCC436.C6-HSL" = "BCC436 + C6-HSL"
                       )) +
    scale_fill_manual(values = combo_colors,
                      labels = c(
                        "Golden_Promise.mock" = "Golden Promise + mock",
                        "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
                        "BCC436.mock" = "BCC436 + mock",
                        "BCC436.C6-HSL" = "BCC436 + C6-HSL"
                      )) +
    scale_shape_manual(values = genotype_shapes,labels = c(
      "Golden_Promise" = "Golden Promise",
      "BCC436" = "BCC436"
    )) +
    labs(
      title = paste("MDS Plot - Rhizosphere - Duration:", d),
      x = xlab,
      y = ylab,
      color = "Genotype & Treatment",
      fill = "Genotype & Treatment",
      shape = "Genotype"
    ) +
    my_theme
  # Store plot in list
  mds_plots[[as.character(d)]] <- p
}



combined_plot <- wrap_plots(mds_plots, ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom"))

combined_plot




# Initialize list to store plots
mds_plots <- list()

for (d in durations) {
  # Subset by duration
  physeq_sub <- subset_samples(physeq_raref_FE_wo_bulk, Duration == d)
  
  # Create combination factor for color mapping
  sample_data(physeq_sub)$Genotype_Treatment <- 
    paste(sample_data(physeq_sub)$Genotype,
          sample_data(physeq_sub)$Treatment,
          sep = ".")
  
  
  sample_data(physeq_sub)$Genotype_Treatment <- 
    factor(sample_data(physeq_sub)$Genotype_Treatment,
           levels = c("Golden_Promise.mock", 
                      "Golden_Promise.C6-HSL",
                      "BCC436.mock", 
                      "BCC436.C6-HSL"))
  
  # Run MDS
  ord <- ordinate(physeq_sub, method = "NMDS", distance = "bray")
  
  # Extract stress
  stress_val <- round(ord$stress, 3)
  stress_label <- paste0("Stress = ", stress_val)
  
  # Plot
  p <- plot_ordination(physeq_sub, ord, type = "samples", 
                       shape = "Genotype", color = "Genotype_Treatment") +
    geom_point(size = 4, alpha = 0.8) +
    stat_ellipse(aes(group = Genotype_Treatment, fill = Genotype_Treatment), 
                 geom = "polygon", alpha = 0.2) +
    scale_color_manual(values = combo_colors,
                       labels = c(
                         "Golden_Promise.mock" = "Golden Promise + mock",
                         "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
                         "BCC436.mock" = "BCC436 + mock",
                         "BCC436.C6-HSL" = "BCC436 + C6-HSL"
                       )) +
    scale_fill_manual(values = combo_colors,
                      labels = c(
                        "Golden_Promise.mock" = "Golden Promise + mock",
                        "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
                        "BCC436.mock" = "BCC436 + mock",
                        "BCC436.C6-HSL" = "BCC436 + C6-HSL"
                      )) +
    scale_shape_manual(values = genotype_shapes,labels = c(
      "Golden_Promise" = "Golden Promise",
      "BCC436" = "BCC436"
    )) +
    labs(
      title = paste("MDS Plot - Rhizosphere - Duration:", d),
      subtitle = stress_label, 
      x = xlab,
      y = ylab,
      color = "Genotype & Treatment",
      fill = "Genotype & Treatment",
      shape = "Genotype"
    ) +
    my_theme
  # Store plot in list
  mds_plots[[as.character(d)]] <- p
}



combined_plot <- wrap_plots(mds_plots, ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom"))

combined_plot


file_name <- paste0("MDS/MDS_FE.png")

ggsave(
  filename = file_name,
  plot = last_plot(),
  width = 22,       # width in cm
  height = 15,      # height in cm
  units = "cm",     # specify centimeters
  dpi = 300,
  path = out_dir# resolution (300 dpi = publication quality)
)

library(cowplot)   # or library(ggpubr)

legendy <- cowplot::get_legend(combined_plot)   # extract legend

legend_plot <- cowplot::ggdraw(legendy)

ggsave(
  filename = "legend_MDS_FE.png",
  plot = legend_plot,
  width = 15,      # cm
  height = 10,     # cm
  units = "cm",
  dpi = 300,
  path = out_dir# resolution (300 dpi = publication quality)
)



# PERMANOVA ---------------------------------------------------------------


library(vegan)
library(pairwiseAdonis)
library(microbiome)


#Substrate
# permanova on sqrt-count data
physeq_raref_FE_wo_bulk <- subset_samples(physeq_raref_FE, Microhabitat == "RZ")

# Metadata extrahieren
metadata <- as(sample_data(physeq_raref_JKI_wo_bulk), "data.frame")

# Ergebnisse vorbereiten
global_results <- list()
pairwise_results <- list()

for (d in unique(metadata$Duration)) {
  
  message(">>> Duration: ", d)
  
  # Subset Metadata
  metadata_d <- metadata[metadata$Duration == d, ]
  
  # Subset Phyloseq
  physeq_d <- prune_samples(rownames(metadata_d), physeq_raref_JKI_wo_bulk)
  
  # Distanzmatrix (z. B. Bray-Curtis)
  bray_dist <- phyloseq::distance(physeq_d, method = "bray")
  
  # --- Globale PERMANOVA ---
  adonis_res <- adonis2(bray_dist ~ Treatment * Genotype,
                        data = metadata_d,
                        permutations = 10000,
                        by = "terms",
                        strata = metadata_d$strata_combined)
  
  global_results[[d]] <- adonis_res
  
  # --- Pairwise Post-hoc ---
  # Treatment
  pairwise_results[[paste0("Treatment_", d)]] <- 
    pairwise.adonis(as.matrix(bray_dist), metadata_d$Treatment, 
                    perm = 9999, p.adjust.m = "fdr")
  
  # Genotype
  pairwise_results[[paste0("Genotype_", d)]] <- 
    pairwise.adonis(as.matrix(bray_dist), metadata_d$Genotype, 
                    perm = 9999, p.adjust.m = "fdr")
}

# Ergebnisse als Dataframes zusammenfassen
global_df <- do.call(rbind, lapply(names(global_results), function(d) {
  res <- as.data.frame(global_results[[d]])
  res$Duration <- d
  res$Factor <- rownames(res)
  res
}))

pairwise_df <- do.call(rbind, lapply(names(pairwise_results), function(name) {
  res <- pairwise_results[[name]]
  res$ComparisonType <- name
  res
}))






##JKI

library(vegan)
library(pairwiseAdonis)
library(microbiome)


physeq_raref_JKI_wo_bulk <- subset_samples(physeq_raref_JKI, Microhabitat == "RZ")

# Metadata extrahieren
metadata <- as(sample_data(physeq_raref_JKI_wo_bulk), "data.frame")

# Ergebnisse vorbereiten
global_results <- list()
pairwise_results <- list()
treatment_per_genotype_results <- list()

for (d in unique(metadata$Duration)) {
  
  message(">>> Duration: ", d)
  
  # Subset Metadata
  metadata_d <- metadata[metadata$Duration == d, ]
  
  # Subset Phyloseq
  physeq_d <- prune_samples(rownames(metadata_d), physeq_raref_JKI_wo_bulk)
  
  # Distanzmatrix (z. B. Bray-Curtis)
  bray_dist <- phyloseq::distance(physeq_d, method = "bray")
  
  # --- Globale PERMANOVA ---
  adonis_res <- adonis2(bray_dist ~ Treatment * Genotype,
                        data = metadata_d,
                        permutations = 10000,
                        by = "terms",
                        strata = metadata_d$strata_combined)
  
  global_results[[d]] <- adonis_res
  
  # --- Pairwise Post-hoc ---
  # Treatment global
  pairwise_results[[paste0("Treatment_", d)]] <- pairwise.adonis(
    as.matrix(bray_dist),
    metadata_d$Treatment,
    perm = 9999,
    p.adjust.m = "fdr"
  )
  
  # Genotype global
  pairwise_results[[paste0("Genotype_", d)]] <- pairwise.adonis(
    as.matrix(bray_dist),
    metadata_d$Genotype,
    perm = 9999,
    p.adjust.m = "fdr"
  )
  
  # --- Treatment innerhalb jeder Genotype ---
  for(gen in unique(metadata_d$Genotype)){
    
    # Metadata Subset für diese Genotype
    meta_gen <- subset(metadata_d, Genotype == gen)
    
    if(length(unique(meta_gen$Treatment)) > 1){  # nur, wenn mehr als 1 Treatment
      # Distanzmatrix für Subgruppe
      dist_gen <- as.matrix(phyloseq::distance(prune_samples(rownames(meta_gen), physeq_d), method = "bray"))
      
      treatment_per_genotype_results[[paste0(gen, "_", d)]] <- pairwise.adonis(
        dist_gen,
        meta_gen$Treatment,
        perm = 9999,
        p.adjust.m = "fdr"
      )
    }
  }
}

# --- Ergebnisse zusammenfassen ---
# Global PERMANOVA
global_df <- do.call(rbind, lapply(names(global_results), function(d) {
  res <- as.data.frame(global_results[[d]])
  res$Duration <- d
  res$Factor <- rownames(res)
  res
}))

# Pairwise Post-hoc (global)
pairwise_df <- do.call(rbind, lapply(names(pairwise_results), function(name) {
  res <- pairwise_results[[name]]
  res$ComparisonType <- name
  res
}))

# Pairwise Treatment pro Genotype
treatment_per_genotype_df <- do.call(rbind, lapply(names(treatment_per_genotype_results), function(name) {
  res <- treatment_per_genotype_results[[name]]
  res$Genotype_Duration <- name
  res
}))

# Ergebnis-Listen anschauen
global_df
pairwise_df
treatment_per_genotype_df


# Installiere ggf. einmal das Paket
# install.packages("openxlsx")
library(openxlsx)

wb <- createWorkbook()

addWorksheet(wb, "Global_PERMANOVA")
writeData(wb, "Global_PERMANOVA", global_df)

addWorksheet(wb, "Pairwise_Global")
writeData(wb, "Pairwise_Global", pairwise_df)

addWorksheet(wb, "Treatment_per_Genotype")
writeData(wb, "Treatment_per_Genotype", treatment_per_genotype_df)

# Speichern
saveWorkbook(wb, "PERMANOVA_results.xlsx", overwrite = TRUE)



##FE

library(vegan)
library(pairwiseAdonis)
library(microbiome)

# Metadata extrahieren
metadata <- as(sample_data(physeq_raref_FE_wo_bulk), "data.frame")

# Ergebnisse vorbereiten
global_results <- list()
pairwise_results <- list()
treatment_per_genotype_results <- list()

for (d in unique(metadata$Duration)) {
  
  message(">>> Duration: ", d)
  
  # Subset Metadata
  metadata_d <- metadata[metadata$Duration == d, ]
  
  # Subset Phyloseq
  physeq_d <- prune_samples(rownames(metadata_d), physeq_raref_FE_wo_bulk)
  
  # Distanzmatrix (z. B. Bray-Curtis)
  bray_dist <- phyloseq::distance(physeq_d, method = "bray")
  
  # --- Globale PERMANOVA ---
  adonis_res <- adonis2(bray_dist ~ Treatment * Genotype,
                        data = metadata_d,
                        permutations = 10000,
                        by = "terms",
                        strata = metadata_d$strata_combined)
  
  global_results[[d]] <- adonis_res
  
  # --- Pairwise Post-hoc ---
  # Treatment global
  pairwise_results[[paste0("Treatment_", d)]] <- pairwise.adonis(
    as.matrix(bray_dist),
    metadata_d$Treatment,
    perm = 9999,
    p.adjust.m = "fdr"
  )
  
  # Genotype global
  pairwise_results[[paste0("Genotype_", d)]] <- pairwise.adonis(
    as.matrix(bray_dist),
    metadata_d$Genotype,
    perm = 9999,
    p.adjust.m = "fdr"
  )
  
  # --- Treatment innerhalb jeder Genotype ---
  for(gen in unique(metadata_d$Genotype)){
    
    # Metadata Subset für diese Genotype
    meta_gen <- subset(metadata_d, Genotype == gen)
    
    if(length(unique(meta_gen$Treatment)) > 1){  # nur, wenn mehr als 1 Treatment
      # Distanzmatrix für Subgruppe
      dist_gen <- as.matrix(phyloseq::distance(prune_samples(rownames(meta_gen), physeq_d), method = "bray"))
      
      treatment_per_genotype_results[[paste0(gen, "_", d)]] <- pairwise.adonis(
        dist_gen,
        meta_gen$Treatment,
        perm = 9999,
        p.adjust.m = "fdr"
      )
    }
  }
}

# --- Ergebnisse zusammenfassen ---
# Global PERMANOVA
global_df <- do.call(rbind, lapply(names(global_results), function(d) {
  res <- as.data.frame(global_results[[d]])
  res$Duration <- d
  res$Factor <- rownames(res)
  res
}))

# Pairwise Post-hoc (global)
pairwise_df <- do.call(rbind, lapply(names(pairwise_results), function(name) {
  res <- pairwise_results[[name]]
  res$ComparisonType <- name
  res
}))

# Pairwise Treatment pro Genotype
treatment_per_genotype_df <- do.call(rbind, lapply(names(treatment_per_genotype_results), function(name) {
  res <- treatment_per_genotype_results[[name]]
  res$Genotype_Duration <- name
  res
}))

# Ergebnis-Listen anschauen
global_df
pairwise_df
treatment_per_genotype_df


# Installiere ggf. einmal das Paket
# install.packages("openxlsx")
library(openxlsx)

wb <- createWorkbook()

addWorksheet(wb, "Global_PERMANOVA")
writeData(wb, "Global_PERMANOVA", global_df)

addWorksheet(wb, "Pairwise_Global")
writeData(wb, "Pairwise_Global", pairwise_df)

addWorksheet(wb, "Treatment_per_Genotype")
writeData(wb, "Treatment_per_Genotype", treatment_per_genotype_df)

# Speichern
saveWorkbook(wb, "PERMANOVA_results_FE.xlsx", overwrite = TRUE)


# global

# Subset rhizosphere samples



library(vegan)
library(pairwiseAdonis)
library(microbiome)


physeq_raref_wo_bulk <- subset_samples(physeq_raref, Microhabitat == "RZ")

# Metadata extrahieren
metadata <- as(sample_data(physeq_raref_wo_bulk), "data.frame")

# Subset Phyloseq
physeq_d <- prune_samples(rownames(metadata), physeq_raref_wo_bulk)

# Distanzmatrix (z. B. Bray-Curtis)
bray_dist <- phyloseq::distance(physeq_d, method = "bray")

# --- Globale PERMANOVA ---
adonis_res <- adonis2(bray_dist ~ Treatment * Genotype * Soil * Duration,
                      data = metadata,
                      permutations = 10000,
                      by = "terms",
                      strata = metadata$strata_combined)

adonis_res

adonis_df <- as.data.frame(adonis_res)
adonis_df$Term <- rownames(adonis_df)

adonis_df <- adonis_df |>
  dplyr::filter(!Term %in% c("Residual", "Total")) |>
  dplyr::mutate(
    Significance = case_when(
      `Pr(>F)` < 0.001 ~ "***",
      `Pr(>F)` < 0.01  ~ "**",
      `Pr(>F)` < 0.05  ~ "*",
      `Pr(>F)` < 0.1   ~ ".",
      TRUE             ~ ""
    )
  )


# Installiere ggf. einmal das Paket
# install.packages("openxlsx")
library(openxlsx)

wb <- createWorkbook()

addWorksheet(wb, "Global_PERMANOVA")
writeData(wb, "Global_PERMANOVA", adonis_df)


# Speichern
saveWorkbook(wb, "global_PERMANOVA.xlsx", overwrite = TRUE)


# Relative abundance ------------------------------------------------------
library(ggh4x)
#relative abundance on phylum level
PS_phyl <- tax_glom(physeq_raref, "Phylum")  ##phylum level
PS_phyl_rel <- transform_sample_counts(PS_phyl, function(x) x/sum(x))  ##transform to relabu
PS_phyl_rel <-prune_taxa(taxa_sums(PS_phyl_rel) >0 , PS_phyl_rel) 
PS_phyl_rel <- filter_taxa(PS_phyl_rel, 
                           function(x) mean(x) >= 0.01, 
                           prune = TRUE)

#relative abundance on class level
PS_cl <- tax_glom(physeq_raref, "Class")  ##Class level
PS_cl_rel <- transform_sample_counts(PS_cl, function(x) x/sum(x))  ##transform to relabu
PS_cl_rel <-prune_taxa(taxa_sums(PS_cl_rel) >0 , PS_cl_rel) 

PS_cl_rel <- filter_taxa(PS_cl_rel, 
                         function(x) mean(x) >= 0.01, 
                         prune = TRUE)  


# Relative Abundance Phylum ----------------------------------------------------------------


#create data.frame for plotting
relabu_phyla <- psmelt(PS_phyl_rel) #create data.frame for plotting
relabu_phyla$Abundance <- relabu_phyla$Abundance * 100


# Order phyla by total abundance
phyla_order <- relabu_phyla %>%
  group_by(Phylum) %>%
  summarise(total = sum(Abundance, na.rm=TRUE)) %>%
  arrange(desc(total)) %>%
  pull(Phylum)

# Set factor levels
relabu_phyla$Phylum <- factor(relabu_phyla$Phylum, levels = phyla_order)

ggplot(relabu_phyla, aes(x = Number, y = Abundance, fill = Phylum, group = Number)) +
  theme_bw() +
  labs(title = "Relative Abundance", x = "", y = "Relative Abundance (%)", fill="") +
  geom_bar(position="fill", stat = "identity") +
  facet_nested(~ Soil + Duration +Genotype +  Microhabitat + Treatment,
               scales = "free_x", space = "free") + 
  scale_fill_viridis_d()+
  scale_y_continuous(labels = function(x) x * 100, expand = c(0,0))+
  theme(
    plot.title=element_text(color="black",size = 12, hjust=0, vjust = 2, face= "bold"),
    plot.title.position = "panel",
    strip.text = element_text(size = 10, face = "plain"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.placement = "outside",
    axis.text.x = element_blank(), 
    axis.title.y = element_text(size = 10, face = "plain"),  
    legend.text = element_text(size = 10, face = "plain", colour = "black"),
    axis.text.y = element_text(colour = "black", size = 7, face = "plain"),
    axis.title.x =element_text(size = 10, face = "plain"))


#nur JKI soil + Ungeordnet
#mit 0 days und bulk
PS_JKIS_phyl_rel <- subset_samples(PS_phyl_rel, Soil == "Field soil")

relabu_phyla_JKIS <- psmelt(PS_JKIS_phyl_rel) #create data.frame for plotting



relabu_phyla_JKIS$Abundance <- relabu_phyla_JKIS$Abundance * 100

ggplot(relabu_phyla_JKIS, aes(x = Number, y = Abundance, fill = Phylum, group = Number)) +
  theme_bw() +
  labs(title = "Relative Abundance", x = "", y = "Relative Abundance (%)", fill="") +
  geom_bar(position="fill", stat = "identity") +
  facet_nested(~  Duration +Genotype +  Microhabitat + Treatment,
               scales = "free_x", space = "free") + 
  scale_fill_viridis_d()+
  scale_y_continuous(labels = function(x) x * 100, expand = c(0,0))+
  theme(
    plot.title=element_text(color="black",size = 12, hjust=0, vjust = 2, face= "bold"),
    plot.title.position = "panel",
    strip.text = element_text(size = 10, face = "plain"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.placement = "outside",
    axis.text.x = element_blank(), 
    axis.title.y = element_text(size = 10, face = "plain"),  
    legend.text = element_text(size = 10, face = "plain", colour = "black"),
    axis.text.y = element_text(colour = "black", size = 7, face = "plain"),
    axis.title.x =element_text(size = 10, face = "plain"))


#nur JKI soil + GEORDNET
#ohne 0 days  bulk


relabu_phyla_JKIS_wo_0_bulk <- relabu_phyla_JKIS %>%
  filter(Duration != "Zero-day")

relabu_phyla_JKIS_wo_0_bulk<- relabu_phyla_JKIS_wo_0_bulk %>%   
  mutate(Treatment = as_factor(Treatment)) %>% 
  mutate(Genotype = as_factor(Genotype)) %>% 
  mutate(Duration   = as_factor(Duration  )) %>% 
  mutate(Treatment = fct_relevel(Treatment, "none", "mock", "C6-HSL")) %>% 
  mutate(Genotype = fct_relevel(Genotype, "none","Golden_Promise", "BCC436")) %>% 
  mutate(Duration   = fct_relevel(Duration  , "Three-week-old","18-week-old"))


ggplot(relabu_phyla_JKIS_wo_0_bulk, aes(x = Number, y = Abundance, fill = Phylum, group = Number)) +
  theme_bw() +
  labs(title = "Relative Abundance", x = "", y = "Relative Abundance (%)", fill="") +
  geom_bar(position="fill", stat = "identity") +
  facet_nested(~  Duration +Genotype +  Microhabitat + Treatment,
               scales = "free_x", space = "free") + 
  scale_fill_viridis_d()+
  scale_y_continuous(labels = function(x) x * 100, expand = c(0,0))+
  theme(
    plot.title=element_text(color="black",size = 12, hjust=0, vjust = 2, face= "bold"),
    plot.title.position = "panel",
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(colour = "white"),
    #strip.background = element_rect(fill = "white", color = "black"),
    strip.placement = "outside",
    axis.text.x = element_blank(), 
    axis.title.y = element_text(size = 10, face = "plain"),  
    legend.text = element_text(size = 10, face = "plain", colour = "black"),
    axis.text.y = element_text(colour = "black", size = 7, face = "plain"),
    axis.title.x =element_text(size = 10, face = "plain"))



#nur FE soil + GEORDNET
#ohne 0 days  bulk

PS_FE_phyl_rel <- subset_samples(PS_phyl_rel, Soil == "Substrate")

relabu_phyla_FE <- psmelt(PS_FE_phyl_rel) #create data.frame for plotting
relabu_phyla_FE$Abundance <- relabu_phyla_FE$Abundance * 100

relabu_phyla_FE_wo_0_bulk <- relabu_phyla_FE %>%
  filter(Duration != "Zero-day")

relabu_phyla_FE_wo_0_bulk<- relabu_phyla_FE_wo_0_bulk %>%   
  mutate(Treatment = as_factor(Treatment)) %>% 
  mutate(Genotype = as_factor(Genotype)) %>% 
  mutate(Duration   = as_factor(Duration  )) %>% 
  mutate(Treatment = fct_relevel(Treatment, "none", "mock", "C6-HSL")) %>% 
  mutate(Genotype = fct_relevel(Genotype, "none","Golden_Promise", "BCC436")) %>% 
  mutate(Duration   = fct_relevel(Duration  , "Two-week-old","16-week-old"))

# Order phyla by total abundance
phyla_order <- relabu_phyla_FE_wo_0_bulk %>%
  group_by(Phylum) %>%
  summarise(total = sum(Abundance, na.rm=TRUE)) %>%
  arrange(desc(total)) %>%
  pull(Phylum)

# Set factor levels
relabu_phyla_FE_wo_0_bulk$Phylum <- factor(relabu_phyla_FE_wo_0_bulk$Phylum, levels = phyla_order)

relabu_phyla_FE_wo_0_bulk <- relabu_phyla_FE_wo_0_bulk %>%
  arrange(Number, factor(Phylum, levels = phyla_order))




ggplot(relabu_phyla_FE_wo_0_bulk, aes(x = Number, y = Abundance, fill = Phylum, group = Number)) +
  theme_bw() +
  labs(title = "Relative Abundance", x = "", y = "Relative Abundance (%)", fill="") +
  geom_bar(position="fill", stat = "identity") +
  facet_nested(~  Duration +Microhabitat  +Genotype +   Treatment,
               scales = "free_x", space = "free") + 
  scale_fill_viridis_d(option = "cividis")+
  scale_y_continuous(labels = function(x) x * 100, expand = c(0,0))+
  theme(
    plot.title=element_text(color="black",size = 10, hjust=0, vjust = 2, face= "bold"),
    plot.title.position = "panel",
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
    strip.placement = "outside",
    axis.text.x = element_blank(), 
    axis.title.y = element_text(size = 10, face = "plain"),  
    legend.text = element_text(size = 10, face = "plain", colour = "black"),
    axis.text.y = element_text(colour = "black", size = 7, face = "plain"),
    axis.title.x =element_text(size = 10, face = "plain"))


file_name <- paste0("RA/RA_Phylum_FE.png")

ggsave(
  filename = file_name,
  plot = last_plot(),
  width = 22,       # width in cm
  height = 15,      # height in cm
  units = "cm",     # specify centimeters
  dpi = 300,
  path = out_dir# resolution (300 dpi = publication quality)
)




# Relative Abundance Class -------------------------------------------------------------------


#create data.frame for plotting
relabu_class <- psmelt(PS_cl_rel) #create data.frame for plotting
relabu_class$Abundance <- relabu_class$Abundance * 100

#nur JKI soil + GEORDNET
#ohne 0 days  bulk
PS_JKIS_phyl_class <- subset_samples(PS_cl_rel, Soil == "Field soil")

relabu_class_JKIS <- psmelt(PS_JKIS_phyl_class) #create data.frame for plotting

relabu_class_JKIS$Abundance <- relabu_class_JKIS$Abundance * 100


relabu_class_JKIS_wo_0_bulk <- relabu_class_JKIS %>%
  filter(Duration != "Zero-day")

relabu_class_JKIS_wo_0_bulk<- relabu_class_JKIS_wo_0_bulk %>%   
  mutate(Treatment = as_factor(Treatment)) %>% 
  mutate(Genotype = as_factor(Genotype)) %>% 
  mutate(Duration   = as_factor(Duration  )) %>% 
  mutate(Treatment = fct_relevel(Treatment, "none", "mock", "C6-HSL")) %>% 
  mutate(Genotype = fct_relevel(Genotype, "none","Golden_Promise", "BCC436")) %>% 
  mutate(Duration   = fct_relevel(Duration  , "Three-week-old","18-week-old"))


ggplot(relabu_class_JKIS_wo_0_bulk, aes(x = Number, y = Abundance, fill = Class, group = Number)) +
  theme_bw() +
  labs(title = "Relative Abundance", x = "", y = "Relative Abundance (%)", fill="") +
  geom_bar(position="fill", stat = "identity") +
  facet_nested(~  Duration +Genotype +  Microhabitat + Treatment,
               scales = "free_x", space = "free") + 
  scale_fill_viridis_d(option = "plasma")+
  scale_y_continuous(labels = function(x) x * 100, expand = c(0,0))+
  theme(
    plot.title=element_text(color="black",size = 12, hjust=0, vjust = 2, face= "bold"),
    plot.title.position = "panel",
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(colour = "white"),
    #strip.background = element_rect(fill = "white", color = "black"),
    strip.placement = "outside",
    axis.text.x = element_blank(), 
    axis.title.y = element_text(size = 10, face = "plain"),  
    legend.text = element_text(size = 10, face = "plain", colour = "black"),
    axis.text.y = element_text(colour = "black", size = 7, face = "plain"),
    axis.title.x =element_text(size = 10, face = "plain"))



#nur FE soil + GEORDNET
#ohne 0 days  bulk
PS_FE_phyl_class <- subset_samples(PS_cl_rel, Soil == "Substrate")

relabu_class_FE <- psmelt(PS_FE_phyl_class) #create data.frame for plotting

relabu_class_FE$Abundance <- relabu_class_FE$Abundance * 100


relabu_class_FE_wo_0_bulk <- relabu_class_FE %>%
  filter(Duration != "Zero-day")

relabu_class_FE_wo_0_bulk<- relabu_class_FE_wo_0_bulk %>%   
  mutate(Treatment = as_factor(Treatment)) %>% 
  mutate(Genotype = as_factor(Genotype)) %>% 
  mutate(Duration   = as_factor(Duration  )) %>% 
  mutate(Treatment = fct_relevel(Treatment, "none", "mock", "C6-HSL")) %>% 
  mutate(Genotype = fct_relevel(Genotype, "none","Golden_Promise", "BCC436")) %>% 
  mutate(Duration   = fct_relevel(Duration  , "Two-week-old","16-week-old"))

# Order phyla by total abundance
class_order <- relabu_class_FE_wo_0_bulk %>%
  group_by(Class) %>%
  summarise(total = sum(Abundance, na.rm=TRUE)) %>%
  arrange(desc(total)) %>%
  pull(Class)

# Set factor levels
relabu_class_FE_wo_0_bulk$Class <- factor(relabu_class_FE_wo_0_bulk$Class, levels = class_order)

relabu_class_FE_wo_0_bulk <- relabu_class_FE_wo_0_bulk %>%
  arrange(Number, factor(Class, levels = class_order))



ggplot(relabu_class_FE_wo_0_bulk, aes(x = Number, y = Abundance, fill = Class, group = Number)) +
  theme_bw() +
  labs(title = "Relative Abundance", x = "", y = "Relative Abundance (%)", fill="") +
  geom_bar(position="fill", stat = "identity") +
  facet_nested(~  Duration+Microhabitat +Genotype + Treatment,
               scales = "free_x", space = "free") + 
  scale_fill_viridis_d(option = "plasma")+
  scale_y_continuous(labels = function(x) x * 100, expand = c(0,0))+
  theme(
    plot.title=element_text(color="black",size = 10, hjust=0, vjust = 2, face= "bold"),
    plot.title.position = "panel",
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
    strip.placement = "outside",
    axis.text.x = element_blank(), 
    axis.title.y = element_text(size = 10, face = "plain"),  
    legend.text = element_text(size = 10, face = "plain", colour = "black"),
    axis.text.y = element_text(colour = "black", size = 7, face = "plain"),
    axis.title.x =element_text(size = 10, face = "plain"))



file_name <- paste0("RA/RA_Class_FE.png")

ggsave(
  filename = file_name,
  plot = last_plot(),
  width = 22,       # width in cm
  height = 15,      # height in cm
  units = "cm",     # specify centimeters
  dpi = 300,
  path = out_dir# resolution (300 dpi = publication quality)
)


# TOP 10 most abundance ---------------------------------------------------


library(phyloseq)
library(dplyr)
library(tibble)

# Agglomerate to Genus level
ps_spec <- tax_glom(physeq_raref_FE_wo_bulk, "Species")

# Transform to relative abundance
ps_spec_rel <- transform_sample_counts(ps_spec, function(x) x / sum(x))
ps_spec_rel <-prune_taxa(taxa_sums(ps_spec_rel) >0 , ps_spec_rel) # -> Entfernt **alle Taxa**, die in der gesamten Tabelle nur Nullen haben. 
# (Also „tote“/nie beobachtete Taxa, Summen = 0).

ps_spec_rel <- filter_taxa(ps_spec_rel, 
                           function(x) mean(x) >= 0.01, 
                           prune = TRUE)  # Diese Zeile behält nur die Taxa, deren mittlere relative Abundanz 
# über alle Samples hinweg ≥ 1 % ist.
# -> Alles darunter wird entfernt.


# Melt phyloseq to long format
melted <- psmelt(ps_spec_rel)  # Already at genus level, relative abundance
melted$Abundance <- melted$Abundance * 100


top10_per_group <- melted %>%
  group_by(Duration,Genotype,  Treatment, Species) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  arrange(Genotype, Duration, Treatment, desc(Abundance)) %>%
  group_by(Genotype, Duration, Treatment) %>%
  slice_max(order_by = Abundance, n = 10, with_ties = FALSE)

# Liste aller relevanten Species behalten
species_keep <- unique(top10_per_group$Species)

# Filtere dein melted-Objekt auf diese Species
melted_top <- melted %>%
  filter(Species %in% species_keep)

melted_top<- melted_top %>%   
  mutate(Treatment = as_factor(Treatment)) %>% 
  mutate(Genotype = as_factor(Genotype)) %>% 
  mutate(Duration   = as_factor(Duration  )) %>% 
  mutate(Treatment = fct_relevel(Treatment,  "mock", "C6-HSL")) %>% 
  mutate(Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436")) %>% 
  mutate(Duration   = fct_relevel(Duration  , "Two-week-old","16-week-old"))



ggplot(melted_top, aes(
  x = Sample, 
  y = OTU, 
  fill = interaction(Genotype, Treatment, sep = "."),
  alpha = Abundance
)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    name = "Genotype × Treatment",
    values = c(
      "Golden_Promise.mock"   = "#A67C00",
      "Golden_Promise.C6-HSL" = "#ffcf40",
      "BCC436.mock"           = "#4F7942",
      "BCC436.C6-HSL"         = "#4aba91"
    ),
    labels = c(
      "Golden_Promise.mock" = "Golden Promise + mock",
      "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
      "BCC436.mock" = "BCC436 + mock",
      "BCC436.C6-HSL" = "BCC436 + C6-HSL"
    ),
    breaks =   c(
      "Golden_Promise.mock" ,
      "Golden_Promise.C6-HSL" ,
      "BCC436.mock" ,
      "BCC436.C6-HSL" 
    )
  ) +
  scale_alpha_continuous(
    name = "Rel. Abundance [%]",
    range = c(0.1, 1)   # 0.1 = sehr blass, 1 = volle Sättigung
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))+
  facet_nested(~ Duration + Genotype + Treatment, scales = "free_x") +
  theme(
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Top 10 Most Abundant ASV per Group",
    x = "Sample",
    y = "ASV"
  )




#Tabelle für unten drunter
tab1<-melted_top %>% 
  select(OTU, Phylum,Order,Family,Genus,Species) %>%
  distinct()  # entfernt Duplikate





top10_per_group <- melted %>%
  group_by(Duration, Genotype, Treatment, OTU) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  arrange(Genotype, Duration, Treatment, desc(Abundance)) %>%
  group_by(Genotype, Duration, Treatment) %>%
  slice_max(order_by = Abundance, n = 10, with_ties = FALSE)

OTU_keep <- unique(top10_per_group$OTU)

melted_top <- melted %>%
  filter(OTU %in% OTU_keep)


melted_top<- melted_top %>%   
  mutate(Treatment = as_factor(Treatment)) %>% 
  mutate(Genotype = as_factor(Genotype)) %>% 
  mutate(Duration   = as_factor(Duration  )) %>% 
  mutate(Treatment = fct_relevel(Treatment,  "mock", "C6-HSL")) %>% 
  mutate(Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436")) %>% 
  mutate(Duration   = fct_relevel(Duration  , "Two-week-old","16-week-old"))


ggplot(melted_top, aes(
  x = Sample, 
  y = OTU, 
  fill = interaction(Genotype, Treatment, sep = "."),
  alpha = Abundance
)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    name = "Genotype & Treatment",
    values = c(
      "Golden_Promise.mock"   = "#A67C00",
      "Golden_Promise.C6-HSL" = "#ffcf40",
      "BCC436.mock"           = "#4F7942",
      "BCC436.C6-HSL"         = "#4aba91"
    ),
    labels = c(
      "Golden_Promise.mock" = "Golden Promise + mock",
      "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
      "BCC436.mock" = "BCC436 + mock",
      "BCC436.C6-HSL" = "BCC436 + C6-HSL"
    ),
    breaks =   c(
      "Golden_Promise.mock" ,
      "Golden_Promise.C6-HSL" ,
      "BCC436.mock" ,
      "BCC436.C6-HSL"
    )
  ) +
  scale_alpha_continuous(
    name = "Rel. Abundance [%]",
    range = c(0.1, 1)   # 0.1 = sehr blass, 1 = volle Sättigung
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))+
  facet_nested(~ Duration + Genotype + Treatment, scales = "free_x") +
  theme(
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Top 10 Most Abundant ASV per Group",
    x = "Sample",
    y = "ASV"
  )+ theme(
    plot.title=element_text(color="black",size = 10, hjust=0, vjust = 2, face= "bold"),
    plot.title.position = "panel",
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
    strip.placement = "outside",
    axis.text.x = element_blank(), 
    axis.title.y = element_text(size = 10, face = "plain"),  
    legend.text = element_text(size = 10, face = "plain", colour = "black"),
    #legend.position = "bottom",
    axis.text.y = element_text(colour = "black", size = 7, face = "plain"),
    axis.title.x =element_text(size = 10, face = "plain"))


file_name <- paste0("TOP10/TOP10_Phylum_FE.png")



ggsave(
  filename = file_name,
  plot = last_plot(),
  width = 25,      # cm
  height = 10,     # cm
  units = "cm",     # specify centimeters
  dpi = 300,
  path = out_dir# resolution (300 dpi = publication quality)
)



# Alpha diversity for FE without bulk -------------------------------------
# ###
# Observed & Chao1 → tell you how many species are present.
# 
# Shannon & Simpson → tell you how diverse the community is in terms of both species number and relative abundances.
####
my_theme <-  theme(
  legend.text = element_text(size = 14, colour = "black"),
  axis.title.x = element_text(colour = "black", size = 14),
  axis.title.y = element_text( size = 14, colour = "black"),
  axis.text.x = element_text(vjust = 0.5, size = 14, colour = "black"),
  axis.text.y = element_text(vjust = 0.5, size = 14, colour = "black"),
  axis.line = element_line(colour = "black"),
  panel.background = element_blank(),
  plot.title = element_text(size = 14, color = "#000", hjust = 0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.title = element_text(size = 14, face = "bold", colour = "black")  # Legend title
)


physeq_raref_FE <- subset_samples(physeq_raref, Soil == "Substrate")
physeq_raref_FE_wo_bulk <- subset_samples(physeq_raref_FE, Microhabitat == "RZ")


alpha_div <- estimate_richness(physeq_raref_FE_wo_bulk, measures=c("Shannon", "Simpson", "Chao1", "Observed"))
head(alpha_div)
# Add metadata to the alpha diversity table
alpha_div$SampleID <- rownames(alpha_div)
alpha_div <- merge(alpha_div, metadata, by="row.names")
#alpha_div$Soil<-as.factor(alpha_div$Soil)

ggplot(alpha_div, aes(x=Duration, y=Simpson, fill=Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title="Alpha Diversity (Simpson Index)", x="Treatment", y="Simpson Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~Genotype, scales = "free_y")+
  stat_compare_means(method = "kruskal.test", label = "p.signif")


# Optional: reorder factor levels if needed
alpha_div <- alpha_div %>%
  mutate(Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
         Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
         Duration = fct_relevel(Duration, "Two-week-old", "16-week-old"))



mod_shannon <- aov(Shannon ~ Genotype *Treatment * Duration, data = alpha_div)
anova(mod_shannon)

# für jede Duration seperat
mean_comp <- mod_shannon %>% 
  emmeans(specs = ~  Genotype : Treatment| Duration ) %>% # adj. mean per cultivar
  cld(adjust = "Tukey", Letters = letters) # compact letter display (CLD)

# über die Durations hinaus
mean_comp <- mod_shannon %>% 
  emmeans(specs = ~  Genotype *Treatment * Duration ) %>% # adj. mean per cultivar
  cld(adjust = "Tukey", Letters = letters) # compact letter display (CLD)

# Rename .group to something clearer
mean_comp <- mean_comp %>% 
  rename(letter = .group)
mean_comp$letter <- trimws(mean_comp$letter)

alpha_div <- alpha_div %>%
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Duration = fct_relevel(Duration, "Two-week-old", "16-week-old"),
    
    # FIX: Create interaction factor
    Genotype_Treatment = interaction(Genotype, Treatment, sep = "."),
    
    # FIX: Force the correct order
    Genotype_Treatment = factor(
      Genotype_Treatment,
      levels = c(
        "Golden_Promise.mock",
        "Golden_Promise.C6-HSL",
        "BCC436.mock",
        "BCC436.C6-HSL"
      )
    )
  )

mean_comp<- mean_comp %>% 
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Duration = fct_relevel(Duration, "Two-week-old", "16-week-old"),
    
    # FIX: Create interaction factor
    Genotype_Treatment = interaction(Genotype, Treatment, sep = "."),
    
    # FIX: Force the correct order
    Genotype_Treatment = factor(
      Genotype_Treatment,
      levels = c(
        "Golden_Promise.mock",
        "Golden_Promise.C6-HSL",
        "BCC436.mock",
        "BCC436.C6-HSL"
      )
    )
  )




plot_Shannon_FE<-ggplot(alpha_div, aes(x=Duration, y=Shannon, fill=Treatment)) +
  geom_boxplot(aes(fill = Genotype_Treatment, color = Genotype_Treatment),
               width = 0.6,
               outlier.color = "white",
               alpha = 0.3,
               lwd = 0.6,
               fatten = 1.5) +
  geom_jitter(
    aes(color = Genotype_Treatment),
    size = 1.8,
    shape = 9,
    alpha = 1,
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.6)
  ) +
  geom_text(
    data = mean_comp,
    aes(
      x = Duration,
      y = upper.CL,
      label = letter,
      group = Genotype_Treatment,
      color = Genotype_Treatment
    ),
    color = "black",
    position = position_dodge(width = 0.59),
    inherit.aes = FALSE,
    size = 5, hjust = 0.5
  ) +
  labs(
    x = "",
    y = "Shannon Index",
    fill = "Genotype + Treatment",
    color = "Genotype + Treatment"
  ) +
  
  # Colors for both box and jitter
  scale_fill_manual(
    values = c(
      "Golden_Promise.mock" = "#A67C00",
      "Golden_Promise.C6-HSL" = "#ffcf40",
      "BCC436.mock" = "#4F7942",
      "BCC436.C6-HSL" = "#4aba91"
    ),
    labels = c(
      "Golden_Promise.mock" = "Golden Promise + mock",
      "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
      "BCC436.mock" = "BCC436 + mock",
      "BCC436.C6-HSL" = "BCC436 + C6-HSL"
    ),
    breaks =   c(
      "Golden_Promise.mock" ,
      "Golden_Promise.C6-HSL" ,
      "BCC436.mock" ,
      "BCC436.C6-HSL" 
    )
  )+
  scale_color_manual(  values = c(
    "Golden_Promise.mock" = "#A67C00",
    "Golden_Promise.C6-HSL" = "#ffcf40",
    "BCC436.mock" = "#4F7942",
    "BCC436.C6-HSL" = "#4aba91"
  ),
  labels = c(
    "Golden_Promise.mock" = "Golden Promise + mock",
    "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
    "BCC436.mock" = "BCC436 + mock",
    "BCC436.C6-HSL" = "BCC436 + C6-HSL"
  ),
  breaks  =   c(
    "Golden_Promise.mock" ,
    "Golden_Promise.C6-HSL" ,
    "BCC436.mock" ,
    "BCC436.C6-HSL" 
  )
  ) +
  my_theme+
  theme(axis.text.x = element_blank())




#Für alle normalverteilten ANOVA
mod_Chao1 <- aov(Chao1 ~ Genotype *Treatment * Duration, data = alpha_div)
anova(mod_Chao1)

# für jede Duration seperat
mean_comp <- mod_Chao1 %>% 
  emmeans(specs = ~  Genotype : Treatment| Duration ) %>% # adj. mean per cultivar
  cld(adjust = "Tukey", Letters = letters) # compact letter display (CLD)

# über die Durations hinaus
mean_comp <- mod_Chao1 %>% 
  emmeans(specs = ~  Genotype *Treatment * Duration ) %>% # adj. mean per cultivar
  cld(adjust = "Tukey", 
      Letters = letters
  ) # compact letter display (CLD)

# Rename .group to something clearer
mean_comp <- mean_comp %>% 
  rename(letter = .group)
mean_comp$letter <- trimws(mean_comp$letter)

alpha_div <- alpha_div %>%
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Duration = fct_relevel(Duration, "Two-week-old", "16-week-old"),
    
    # FIX: Create interaction factor
    Genotype_Treatment = interaction(Genotype, Treatment, sep = "."),
    
    # FIX: Force the correct order
    Genotype_Treatment = factor(
      Genotype_Treatment,
      levels = c(
        "Golden_Promise.mock",
        "Golden_Promise.C6-HSL",
        "BCC436.mock",
        "BCC436.C6-HSL"
      )
    )
  )

mean_comp<- mean_comp %>% 
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Duration = fct_relevel(Duration, "Two-week-old", "16-week-old"),
    
    # FIX: Create interaction factor
    Genotype_Treatment = interaction(Genotype, Treatment, sep = "."),
    
    # FIX: Force the correct order
    Genotype_Treatment = factor(
      Genotype_Treatment,
      levels = c(
        "Golden_Promise.mock",
        "Golden_Promise.C6-HSL",
        "BCC436.mock",
        "BCC436.C6-HSL"
      )
    )
  )

plot_Chao1_FE<-ggplot(alpha_div, aes(x=Duration, y=Chao1, fill=Treatment)) +
  geom_boxplot(aes(fill = Genotype_Treatment, color = Genotype_Treatment),
               width = 0.6,
               outlier.color = "white",
               alpha = 0.3,
               lwd = 0.6,
               fatten = 1.5) +
  geom_jitter(
    aes(color = Genotype_Treatment),
    size = 1.8,
    shape = 9,
    alpha = 1,
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.6)
  ) +
  geom_text(
    data = mean_comp,
    aes(
      x = Duration,
      y = upper.CL,
      label = letter,
      group = Genotype_Treatment,
      color = Genotype_Treatment
    ),
    color = "black",
    position = position_dodge(width = 0.59),
    inherit.aes = FALSE,
    size = 5, hjust = 0.5
  ) +
  labs(
    x = "",
    y = "Chao1 Index",
    fill = "Genotype + Treatment",
    color = "Genotype + Treatment"
  ) +
  
  # Colors for both box and jitter
  scale_fill_manual(
    values = c(
      "Golden_Promise.mock" = "#A67C00",
      "Golden_Promise.C6-HSL" = "#ffcf40",
      "BCC436.mock" = "#4F7942",
      "BCC436.C6-HSL" = "#4aba91"
    ),
    labels = c(
      "Golden_Promise.mock" = "Golden Promise + mock",
      "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
      "BCC436.mock" = "BCC436 + mock",
      "BCC436.C6-HSL" = "BCC436 + C6-HSL"
    ),
    breaks =   c(
      "Golden_Promise.mock" ,
      "Golden_Promise.C6-HSL" ,
      "BCC436.mock" ,
      "BCC436.C6-HSL" 
    )
  )+
  scale_color_manual(  values = c(
    "Golden_Promise.mock" = "#A67C00",
    "Golden_Promise.C6-HSL" = "#ffcf40",
    "BCC436.mock" = "#4F7942",
    "BCC436.C6-HSL" = "#4aba91"
  ),
  labels = c(
    "Golden_Promise.mock" = "Golden Promise + mock",
    "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
    "BCC436.mock" = "BCC436 + mock",
    "BCC436.C6-HSL" = "BCC436 + C6-HSL"
  ),
  breaks  =   c(
    "Golden_Promise.mock" ,
    "Golden_Promise.C6-HSL" ,
    "BCC436.mock" ,
    "BCC436.C6-HSL" 
  )
  ) +
  my_theme



#Für alle normalverteilten ANOVA
mod_Observed <- aov(Observed ~ Genotype *Treatment * Duration, data = alpha_div)
anova(mod_Observed)

# für jede Duration seperat
mean_comp <- mod_Observed %>% 
  emmeans(specs = ~  Genotype : Treatment| Duration ) %>% # adj. mean per cultivar
  cld(adjust = "Tukey", Letters = letters) # compact letter display (CLD)

# über die Durations hinaus
mean_comp <- mod_Observed %>% 
  emmeans(specs = ~  Genotype *Treatment * Duration ) %>% # adj. mean per cultivar
  cld(adjust = "Tukey", Letters = letters) # compact letter display (CLD)

# Rename .group to something clearer
mean_comp <- mean_comp %>% 
  rename(letter = .group)
mean_comp$letter <- trimws(mean_comp$letter)

alpha_div <- alpha_div %>%
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Duration = fct_relevel(Duration, "Two-week-old", "16-week-old"),
    
    # FIX: Create interaction factor
    Genotype_Treatment = interaction(Genotype, Treatment, sep = "."),
    
    # FIX: Force the correct order
    Genotype_Treatment = factor(
      Genotype_Treatment,
      levels = c(
        "Golden_Promise.mock",
        "Golden_Promise.C6-HSL",
        "BCC436.mock",
        "BCC436.C6-HSL"
      )
    )
  )

mean_comp<- mean_comp %>% 
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Duration = fct_relevel(Duration, "Two-week-old", "16-week-old"),
    
    # FIX: Create interaction factor
    Genotype_Treatment = interaction(Genotype, Treatment, sep = "."),
    
    # FIX: Force the correct order
    Genotype_Treatment = factor(
      Genotype_Treatment,
      levels = c(
        "Golden_Promise.mock",
        "Golden_Promise.C6-HSL",
        "BCC436.mock",
        "BCC436.C6-HSL"
      )
    )
  )

plot_Observed_FE<-ggplot(alpha_div, aes(x=Duration, y=Observed, fill=Treatment)) +
  geom_boxplot(aes(fill = Genotype_Treatment, color = Genotype_Treatment),
               width = 0.6,
               outlier.color = "white",
               alpha = 0.3,
               lwd = 0.6,
               fatten = 1.5) +
  geom_jitter(
    aes(color = Genotype_Treatment),
    size = 1.8,
    shape = 9,
    alpha = 1,
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.6)
  ) +
  geom_text(
    data = mean_comp,
    aes(
      x = Duration,
      y = upper.CL,
      label = letter,
      group = Genotype_Treatment,
      color = Genotype_Treatment
    ),
    color = "black",
    position = position_dodge(width = 0.59),
    inherit.aes = FALSE,
    size = 5, hjust = 0.5
  ) +
  labs(
    x = "",
    y = "Observed Index",
    fill = "Genotype + Treatment",
    color = "Genotype + Treatment"
  ) +
  
  # Colors for both box and jitter
  scale_fill_manual(
    values = c(
      "Golden_Promise.mock" = "#A67C00",
      "Golden_Promise.C6-HSL" = "#ffcf40",
      "BCC436.mock" = "#4F7942",
      "BCC436.C6-HSL" = "#4aba91"
    ),
    labels = c(
      "Golden_Promise.mock" = "Golden Promise + mock",
      "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
      "BCC436.mock" = "BCC436 + mock",
      "BCC436.C6-HSL" = "BCC436 + C6-HSL"
    ),
    breaks =   c(
      "Golden_Promise.mock" ,
      "Golden_Promise.C6-HSL" ,
      "BCC436.mock" ,
      "BCC436.C6-HSL" 
    )
  )+
  scale_color_manual(  values = c(
    "Golden_Promise.mock" = "#A67C00",
    "Golden_Promise.C6-HSL" = "#ffcf40",
    "BCC436.mock" = "#4F7942",
    "BCC436.C6-HSL" = "#4aba91"
  ),
  labels = c(
    "Golden_Promise.mock" = "Golden Promise + mock",
    "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
    "BCC436.mock" = "BCC436 + mock",
    "BCC436.C6-HSL" = "BCC436 + C6-HSL"
  ),
  breaks  =   c(
    "Golden_Promise.mock" ,
    "Golden_Promise.C6-HSL" ,
    "BCC436.mock" ,
    "BCC436.C6-HSL" 
  )
  ) +
  my_theme



#Für alle normalverteilten ANOVA
mod_Simpson <- aov(Simpson ~ Genotype *Treatment * Duration, data = alpha_div)
anova(mod_Simpson)

# für jede Duration seperat
mean_comp <- mod_Simpson %>% 
  emmeans(specs = ~  Genotype : Treatment| Duration ) %>% # adj. mean per cultivar
  cld(adjust = "Tukey", Letters = letters) # compact letter display (CLD)

# über die Durations hinaus
mean_comp <- mod_Simpson %>% 
  emmeans(specs = ~  Genotype *Treatment * Duration ) %>% # adj. mean per cultivar
  cld(adjust = "Tukey", Letters = letters) # compact letter display (CLD)

# Rename .group to something clearer
mean_comp <- mean_comp %>% 
  rename(letter = .group)
mean_comp$letter <- trimws(mean_comp$letter)


alpha_div <- alpha_div %>%
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Duration = fct_relevel(Duration, "Two-week-old", "16-week-old"),
    
    # FIX: Create interaction factor
    Genotype_Treatment = interaction(Genotype, Treatment, sep = "."),
    
    # FIX: Force the correct order
    Genotype_Treatment = factor(
      Genotype_Treatment,
      levels = c(
        "Golden_Promise.mock",
        "Golden_Promise.C6-HSL",
        "BCC436.mock",
        "BCC436.C6-HSL"
      )
    )
  )

mean_comp<- mean_comp %>% 
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Duration = fct_relevel(Duration, "Two-week-old", "16-week-old"),
    
    # FIX: Create interaction factor
    Genotype_Treatment = interaction(Genotype, Treatment, sep = "."),
    
    # FIX: Force the correct order
    Genotype_Treatment = factor(
      Genotype_Treatment,
      levels = c(
        "Golden_Promise.mock",
        "Golden_Promise.C6-HSL",
        "BCC436.mock",
        "BCC436.C6-HSL"
      )
    )
  )

plot_Simpson_FE<-ggplot(alpha_div, aes(x=Duration, y=Simpson, fill=Treatment)) +
  geom_boxplot(aes(fill = Genotype_Treatment, color = Genotype_Treatment),
               width = 0.6,
               outlier.color = "white",
               alpha = 0.3,
               lwd = 0.6,
               fatten = 1.5) +
  geom_jitter(
    aes(color = Genotype_Treatment),
    size = 1.8,
    shape = 9,
    alpha = 1,
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.6)
  ) +
  geom_text(
    data = mean_comp,
    aes(
      x = Duration,
      y = upper.CL,
      label = letter,
      group = Genotype_Treatment,
      color = Genotype_Treatment
    ),
    color = "black",
    position = position_dodge(width = 0.59),
    inherit.aes = FALSE,
    size = 5, hjust = 0.5
  ) +
  labs(
    x = "",
    y = "Simpson Index",
    fill = "Genotype + Treatment",
    color = "Genotype + Treatment"
  ) +
  
  # Colors for both box and jitter
  scale_fill_manual(
    values = c(
      "Golden_Promise.mock" = "#A67C00",
      "Golden_Promise.C6-HSL" = "#ffcf40",
      "BCC436.mock" = "#4F7942",
      "BCC436.C6-HSL" = "#4aba91"
    ),
    labels = c(
      "Golden_Promise.mock" = "Golden Promise + mock",
      "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
      "BCC436.mock" = "BCC436 + mock",
      "BCC436.C6-HSL" = "BCC436 + C6-HSL"
    ),
    breaks =   c(
      "Golden_Promise.mock" ,
      "Golden_Promise.C6-HSL" ,
      "BCC436.mock" ,
      "BCC436.C6-HSL" 
    )
  )+
  scale_color_manual(  values = c(
    "Golden_Promise.mock" = "#A67C00",
    "Golden_Promise.C6-HSL" = "#ffcf40",
    "BCC436.mock" = "#4F7942",
    "BCC436.C6-HSL" = "#4aba91"
  ),
  labels = c(
    "Golden_Promise.mock" = "Golden Promise + mock",
    "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
    "BCC436.mock" = "BCC436 + mock",
    "BCC436.C6-HSL" = "BCC436 + C6-HSL"
  ),
  breaks  =   c(
    "Golden_Promise.mock" ,
    "Golden_Promise.C6-HSL" ,
    "BCC436.mock" ,
    "BCC436.C6-HSL" 
  )
  ) +
  my_theme+
  theme(axis.text.x = element_blank())



plot_Shannon_FE
plot_Chao1_FE
plot_Observed_FE
plot_Simpson_FE


ggarrange(plot_Shannon_FE,plot_Simpson_FE   ,
          plot_Chao1_FE,plot_Observed_FE ,  
          common.legend = T,
          legend = "bottom",
          align = "hv")

file_name <- paste0("Diversity/Diversity_FE.png")

ggsave(
  filename = file_name,
  plot = last_plot(),
  width = 22,       # width in cm
  height = 15,      # height in cm
  units = "cm",     # specify centimeters
  dpi = 300,
  path = out_dir# resolution (300 dpi = publication quality)
)



# Alpha diversity for JKI without bulk -------------------------------------

physeq_raref_JKI <- subset_samples(physeq_raref, Soil == "Field soil")

physeq_raref_JKI_wo_bulk <- subset_samples(physeq_raref_JKI, Microhabitat == "RZ")


alpha_div <- estimate_richness(physeq_raref_JKI_wo_bulk, measures=c("Shannon", "Simpson", "Chao1", "Observed"))
head(alpha_div)
# Add metadata to the alpha diversity table
alpha_div$SampleID <- rownames(alpha_div)
alpha_div <- merge(alpha_div, metadata, by="row.names")
#alpha_div$Soil<-as.factor(alpha_div$Soil)

ggplot(alpha_div, aes(x=Duration, y=Simpson, fill=Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title="Alpha Diversity (Simpson Index)", x="Treatment", y="Simpson Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~Genotype, scales = "free_y")+
  stat_compare_means(method = "kruskal.test", label = "p.signif")


# Optional: reorder factor levels if needed
alpha_div <- alpha_div %>%
  mutate(Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
         Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
         Duration = fct_relevel(Duration, "Three-week-old", "18-week-old"))



mod_shannon <- aov(Shannon ~ Genotype *Treatment * Duration, data = alpha_div)
anova(mod_shannon)

# für jede Duration seperat
mean_comp <- mod_shannon %>% 
  emmeans(specs = ~  Genotype : Treatment| Duration ) %>% # adj. mean per cultivar
  cld(adjust = "Tukey", Letters = letters) # compact letter display (CLD)

# über die Durations hinaus
mean_comp <- mod_shannon %>% 
  emmeans(specs = ~  Genotype *Treatment * Duration ) %>% # adj. mean per cultivar
  cld(adjust = "Tukey", Letters = letters) # compact letter display (CLD)

# Rename .group to something clearer
mean_comp <- mean_comp %>% 
  rename(letter = .group)
mean_comp$letter <- trimws(mean_comp$letter)

alpha_div <- alpha_div %>%
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Duration = fct_relevel(Duration, "Two-week-old", "16-week-old"),
    
    # FIX: Create interaction factor
    Genotype_Treatment = interaction(Genotype, Treatment, sep = "."),
    
    # FIX: Force the correct order
    Genotype_Treatment = factor(
      Genotype_Treatment,
      levels = c(
        "Golden_Promise.mock",
        "Golden_Promise.C6-HSL",
        "BCC436.mock",
        "BCC436.C6-HSL"
      )
    )
  )

mean_comp<- mean_comp %>% 
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Duration = fct_relevel(Duration, "Two-week-old", "16-week-old"),
    
    # FIX: Create interaction factor
    Genotype_Treatment = interaction(Genotype, Treatment, sep = "."),
    
    # FIX: Force the correct order
    Genotype_Treatment = factor(
      Genotype_Treatment,
      levels = c(
        "Golden_Promise.mock",
        "Golden_Promise.C6-HSL",
        "BCC436.mock",
        "BCC436.C6-HSL"
      )
    )
  )


plot_Shannon_JKI<-ggplot(alpha_div, aes(x=Duration, y=Shannon, fill=Treatment)) +
  geom_boxplot(aes(fill = Genotype_Treatment, color = Genotype_Treatment),
               width = 0.6,
               outlier.color = "white",
               alpha = 0.3,
               lwd = 0.6,
               fatten = 1.5) +
  geom_jitter(
    aes(color = Genotype_Treatment),
    size = 1.8,
    shape = 9,
    alpha = 1,
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.6)
  ) +
  geom_text(
    data = mean_comp,
    aes(
      x = Duration,
      y = upper.CL,
      label = letter,
      group = Genotype_Treatment,
      color = Genotype_Treatment
    ),
    color = "black",
    position = position_dodge(width = 0.59),
    inherit.aes = FALSE,
    size = 5, hjust = 0.5
  ) +
  labs(
    x = "",
    y = "Shannon Index",
    fill = "Genotype + Treatment",
    color = "Genotype + Treatment"
  ) +
  
  # Colors for both box and jitter
  scale_fill_manual(
    values = c(
      "Golden_Promise.mock" = "#A67C00",
      "Golden_Promise.C6-HSL" = "#ffcf40",
      "BCC436.mock" = "#4F7942",
      "BCC436.C6-HSL" = "#4aba91"
    ),
    labels = c(
      "Golden_Promise.mock" = "Golden Promise + mock",
      "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
      "BCC436.mock" = "BCC436 + mock",
      "BCC436.C6-HSL" = "BCC436 + C6-HSL"
    ),
    breaks =   c(
      "Golden_Promise.mock" ,
      "Golden_Promise.C6-HSL" ,
      "BCC436.mock" ,
      "BCC436.C6-HSL" 
    )
  )+
  scale_color_manual(  values = c(
    "Golden_Promise.mock" = "#A67C00",
    "Golden_Promise.C6-HSL" = "#ffcf40",
    "BCC436.mock" = "#4F7942",
    "BCC436.C6-HSL" = "#4aba91"
  ),
  labels = c(
    "Golden_Promise.mock" = "Golden Promise + mock",
    "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
    "BCC436.mock" = "BCC436 + mock",
    "BCC436.C6-HSL" = "BCC436 + C6-HSL"
  ),
  breaks  =   c(
    "Golden_Promise.mock" ,
    "Golden_Promise.C6-HSL" ,
    "BCC436.mock" ,
    "BCC436.C6-HSL" 
  )
  ) +
  my_theme+
  theme(axis.text.x = element_blank())




mod_Chao1 <- aov(Chao1 ~ Genotype *Treatment * Duration, data = alpha_div)
anova(mod_Chao1)

# für jede Duration seperat
mean_comp <- mod_Chao1 %>% 
  emmeans(specs = ~  Genotype : Treatment| Duration ) %>% # adj. mean per cultivar
  cld(adjust = "Tukey", Letters = letters) # compact letter display (CLD)

# über die Durations hinaus
mean_comp <- mod_Chao1 %>% 
  emmeans(specs = ~  Genotype *Treatment * Duration ) %>% # adj. mean per cultivar
  cld(adjust = "Tukey", 
      Letters = letters
  ) # compact letter display (CLD)

# Rename .group to something clearer
mean_comp <- mean_comp %>% 
  rename(letter = .group)
mean_comp$letter <- trimws(mean_comp$letter)

alpha_div <- alpha_div %>%
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Duration = fct_relevel(Duration, "Two-week-old", "16-week-old"),
    
    # FIX: Create interaction factor
    Genotype_Treatment = interaction(Genotype, Treatment, sep = "."),
    
    # FIX: Force the correct order
    Genotype_Treatment = factor(
      Genotype_Treatment,
      levels = c(
        "Golden_Promise.mock",
        "Golden_Promise.C6-HSL",
        "BCC436.mock",
        "BCC436.C6-HSL"
      )
    )
  )

mean_comp<- mean_comp %>% 
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Duration = fct_relevel(Duration, "Two-week-old", "16-week-old"),
    
    # FIX: Create interaction factor
    Genotype_Treatment = interaction(Genotype, Treatment, sep = "."),
    
    # FIX: Force the correct order
    Genotype_Treatment = factor(
      Genotype_Treatment,
      levels = c(
        "Golden_Promise.mock",
        "Golden_Promise.C6-HSL",
        "BCC436.mock",
        "BCC436.C6-HSL"
      )
    )
  )


plot_Chao1_JKI<-ggplot(alpha_div, aes(x=Duration, y=Chao1, fill=Treatment)) +
  geom_boxplot(aes(fill = Genotype_Treatment, color = Genotype_Treatment),
               width = 0.6,
               outlier.color = "white",
               alpha = 0.3,
               lwd = 0.6,
               fatten = 1.5) +
  geom_jitter(
    aes(color = Genotype_Treatment),
    size = 1.8,
    shape = 9,
    alpha = 1,
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.6)
  ) +
  geom_text(
    data = mean_comp,
    aes(
      x = Duration,
      y = upper.CL,
      label = letter,
      group = Genotype_Treatment,
      color = Genotype_Treatment
    ),
    color = "black",
    position = position_dodge(width = 0.59),
    inherit.aes = FALSE,
    size = 5, hjust = 0.5
  ) +
  labs(
    x = "",
    y = "Chao1 Index",
    fill = "Genotype + Treatment",
    color = "Genotype + Treatment"
  ) +
  
  # Colors for both box and jitter
  scale_fill_manual(
    values = c(
      "Golden_Promise.mock" = "#A67C00",
      "Golden_Promise.C6-HSL" = "#ffcf40",
      "BCC436.mock" = "#4F7942",
      "BCC436.C6-HSL" = "#4aba91"
    ),
    labels = c(
      "Golden_Promise.mock" = "Golden Promise + mock",
      "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
      "BCC436.mock" = "BCC436 + mock",
      "BCC436.C6-HSL" = "BCC436 + C6-HSL"
    ),
    breaks =   c(
      "Golden_Promise.mock" ,
      "Golden_Promise.C6-HSL" ,
      "BCC436.mock" ,
      "BCC436.C6-HSL" 
    )
  )+
  scale_color_manual(  values = c(
    "Golden_Promise.mock" = "#A67C00",
    "Golden_Promise.C6-HSL" = "#ffcf40",
    "BCC436.mock" = "#4F7942",
    "BCC436.C6-HSL" = "#4aba91"
  ),
  labels = c(
    "Golden_Promise.mock" = "Golden Promise + mock",
    "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
    "BCC436.mock" = "BCC436 + mock",
    "BCC436.C6-HSL" = "BCC436 + C6-HSL"
  ),
  breaks  =   c(
    "Golden_Promise.mock" ,
    "Golden_Promise.C6-HSL" ,
    "BCC436.mock" ,
    "BCC436.C6-HSL" 
  )
  ) +
  my_theme



mod_Observed <- aov(Observed ~ Genotype *Treatment * Duration, data = alpha_div)
anova(mod_Observed)

# für jede Duration seperat
mean_comp <- mod_Observed %>% 
  emmeans(specs = ~  Genotype : Treatment| Duration ) %>% # adj. mean per cultivar
  cld(adjust = "Tukey", Letters = letters) # compact letter display (CLD)

# über die Durations hinaus
mean_comp <- mod_Observed %>% 
  emmeans(specs = ~  Genotype *Treatment * Duration ) %>% # adj. mean per cultivar
  cld(adjust = "Tukey", Letters = letters) # compact letter display (CLD)

# Rename .group to something clearer
mean_comp <- mean_comp %>% 
  rename(letter = .group)
mean_comp$letter <- trimws(mean_comp$letter)

alpha_div <- alpha_div %>%
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Duration = fct_relevel(Duration, "Two-week-old", "16-week-old"),
    
    # FIX: Create interaction factor
    Genotype_Treatment = interaction(Genotype, Treatment, sep = "."),
    
    # FIX: Force the correct order
    Genotype_Treatment = factor(
      Genotype_Treatment,
      levels = c(
        "Golden_Promise.mock",
        "Golden_Promise.C6-HSL",
        "BCC436.mock",
        "BCC436.C6-HSL"
      )
    )
  )

mean_comp<- mean_comp %>% 
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Duration = fct_relevel(Duration, "Two-week-old", "16-week-old"),
    
    # FIX: Create interaction factor
    Genotype_Treatment = interaction(Genotype, Treatment, sep = "."),
    
    # FIX: Force the correct order
    Genotype_Treatment = factor(
      Genotype_Treatment,
      levels = c(
        "Golden_Promise.mock",
        "Golden_Promise.C6-HSL",
        "BCC436.mock",
        "BCC436.C6-HSL"
      )
    )
  )


plot_Observed_JKI<-ggplot(alpha_div, aes(x=Duration, y=Observed, fill=Treatment)) +
  geom_boxplot(aes(fill = Genotype_Treatment, color = Genotype_Treatment),
               width = 0.6,
               outlier.color = "white",
               alpha = 0.3,
               lwd = 0.6,
               fatten = 1.5) +
  geom_jitter(
    aes(color = Genotype_Treatment),
    size = 1.8,
    shape = 9,
    alpha = 1,
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.6)
  ) +
  geom_text(
    data = mean_comp,
    aes(
      x = Duration,
      y = upper.CL,
      label = letter,
      group = Genotype_Treatment,
      color = Genotype_Treatment
    ),
    color = "black",
    position = position_dodge(width = 0.59),
    inherit.aes = FALSE,
    size = 5, hjust = 0.5
  ) +
  labs(
    x = "",
    y = "Observed Index",
    fill = "Genotype + Treatment",
    color = "Genotype + Treatment"
  ) +
  
  # Colors for both box and jitter
  scale_fill_manual(
    values = c(
      "Golden_Promise.mock" = "#A67C00",
      "Golden_Promise.C6-HSL" = "#ffcf40",
      "BCC436.mock" = "#4F7942",
      "BCC436.C6-HSL" = "#4aba91"
    ),
    labels = c(
      "Golden_Promise.mock" = "Golden Promise + mock",
      "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
      "BCC436.mock" = "BCC436 + mock",
      "BCC436.C6-HSL" = "BCC436 + C6-HSL"
    ),
    breaks =   c(
      "Golden_Promise.mock" ,
      "Golden_Promise.C6-HSL" ,
      "BCC436.mock" ,
      "BCC436.C6-HSL" 
    )
  )+
  scale_color_manual(  values = c(
    "Golden_Promise.mock" = "#A67C00",
    "Golden_Promise.C6-HSL" = "#ffcf40",
    "BCC436.mock" = "#4F7942",
    "BCC436.C6-HSL" = "#4aba91"
  ),
  labels = c(
    "Golden_Promise.mock" = "Golden Promise + mock",
    "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
    "BCC436.mock" = "BCC436 + mock",
    "BCC436.C6-HSL" = "BCC436 + C6-HSL"
  ),
  breaks  =   c(
    "Golden_Promise.mock" ,
    "Golden_Promise.C6-HSL" ,
    "BCC436.mock" ,
    "BCC436.C6-HSL" 
  )
  ) +
  my_theme


mod_Simpson <- aov(Simpson ~ Genotype *Treatment * Duration, data = alpha_div)
anova(mod_Simpson)

# für jede Duration seperat
mean_comp <- mod_Simpson %>% 
  emmeans(specs = ~  Genotype : Treatment| Duration ) %>% # adj. mean per cultivar
  cld(adjust = "Tukey", Letters = letters) # compact letter display (CLD)

# über die Durations hinaus
mean_comp <- mod_Simpson %>% 
  emmeans(specs = ~  Genotype *Treatment * Duration ) %>% # adj. mean per cultivar
  cld(adjust = "Tukey", Letters = letters) # compact letter display (CLD)

# Rename .group to something clearer
mean_comp <- mean_comp %>% 
  rename(letter = .group)
mean_comp$letter <- trimws(mean_comp$letter)


alpha_div <- alpha_div %>%
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Duration = fct_relevel(Duration, "Two-week-old", "16-week-old"),
    
    # FIX: Create interaction factor
    Genotype_Treatment = interaction(Genotype, Treatment, sep = "."),
    
    # FIX: Force the correct order
    Genotype_Treatment = factor(
      Genotype_Treatment,
      levels = c(
        "Golden_Promise.mock",
        "Golden_Promise.C6-HSL",
        "BCC436.mock",
        "BCC436.C6-HSL"
      )
    )
  )

mean_comp<- mean_comp %>% 
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Duration = fct_relevel(Duration, "Two-week-old", "16-week-old"),
    
    # FIX: Create interaction factor
    Genotype_Treatment = interaction(Genotype, Treatment, sep = "."),
    
    # FIX: Force the correct order
    Genotype_Treatment = factor(
      Genotype_Treatment,
      levels = c(
        "Golden_Promise.mock",
        "Golden_Promise.C6-HSL",
        "BCC436.mock",
        "BCC436.C6-HSL"
      )
    )
  )


plot_Simpson_JKI<-ggplot(alpha_div, aes(x=Duration, y=Simpson, fill=Treatment)) +
  geom_boxplot(aes(fill = Genotype_Treatment, color = Genotype_Treatment),
               width = 0.6,
               outlier.color = "white",
               alpha = 0.3,
               lwd = 0.6,
               fatten = 1.5) +
  geom_jitter(
    aes(color = Genotype_Treatment),
    size = 1.8,
    shape = 9,
    alpha = 1,
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.6)
  ) +
  geom_text(
    data = mean_comp,
    aes(
      x = Duration,
      y = upper.CL,
      label = letter,
      group = Genotype_Treatment,
      color = Genotype_Treatment
    ),
    color = "black",
    position = position_dodge(width = 0.59),
    inherit.aes = FALSE,
    size = 5, hjust = 0.5
  ) +
  labs(
    x = "",
    y = "Simpson Index",
    fill = "Genotype + Treatment",
    color = "Genotype + Treatment"
  ) +
  
  # Colors for both box and jitter
  scale_fill_manual(
    values = c(
      "Golden_Promise.mock" = "#A67C00",
      "Golden_Promise.C6-HSL" = "#ffcf40",
      "BCC436.mock" = "#4F7942",
      "BCC436.C6-HSL" = "#4aba91"
    ),
    labels = c(
      "Golden_Promise.mock" = "Golden Promise + mock",
      "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
      "BCC436.mock" = "BCC436 + mock",
      "BCC436.C6-HSL" = "BCC436 + C6-HSL"
    ),
    breaks =   c(
      "Golden_Promise.mock" ,
      "Golden_Promise.C6-HSL" ,
      "BCC436.mock" ,
      "BCC436.C6-HSL" 
    )
  )+
  scale_color_manual(  values = c(
    "Golden_Promise.mock" = "#A67C00",
    "Golden_Promise.C6-HSL" = "#ffcf40",
    "BCC436.mock" = "#4F7942",
    "BCC436.C6-HSL" = "#4aba91"
  ),
  labels = c(
    "Golden_Promise.mock" = "Golden Promise + mock",
    "Golden_Promise.C6-HSL" = "Golden Promise + C6-HSL",
    "BCC436.mock" = "BCC436 + mock",
    "BCC436.C6-HSL" = "BCC436 + C6-HSL"
  ),
  breaks  =   c(
    "Golden_Promise.mock" ,
    "Golden_Promise.C6-HSL" ,
    "BCC436.mock" ,
    "BCC436.C6-HSL" 
  )
  ) +
  my_theme+
  theme(axis.text.x = element_blank())



plot_Shannon_JKI
plot_Chao1_JKI
plot_Observed_JKI
plot_Simpson_JKI


ggarrange(plot_Shannon_JKI,plot_Simpson_JKI   ,
          plot_Chao1_JKI,plot_Observed_JKI ,  
          common.legend = T,
          legend = "bottom",
          align = "hv")


file_name <- paste0("Diversity/Diversity_JKI.png")

ggsave(
  filename = file_name,
  plot = last_plot(),
  width = 22,       # width in cm
  height = 15,      # height in cm
  units = "cm",     # specify centimeters
  dpi = 300,
  path = out_dir# resolution (300 dpi = publication quality)
)





# Deseq clean -------------------------------------------------------------

# PIPELINE 1 — Soil effect on the microbiome
# 
# Question:
#   Does Substrate differ from Field soil in microbial composition?


library(phyloseq)
library(DESeq2)
library(dplyr)

physeq_soil <- subset_samples(
  physeq_raref,
  Soil %in% c("Substrate", "Field soil")
)

# Remove low-abundance taxa
physeq_soil <- prune_taxa(taxa_sums(physeq_soil) > 10, physeq_soil)

# Ensure Soil is a factor with correct reference level
sample_data(physeq_soil)$Soil <- factor(
  sample_data(physeq_soil)$Soil,
  levels = c("Field soil", "Substrate")  # reference = Field soil
)

# This ensures:
#   
# log2FC > 0 → enriched in Substrate
# log2FC < 0 → enriched in Field soil

dds_soil <- phyloseq_to_deseq2(physeq_soil, ~ Soil)

dds_soil <- DESeq(
  dds_soil,
  test = "Wald",
  fitType = "parametric",
  sfType = "poscounts"
)

resultsNames(dds_soil)

res_soil <- results(
  dds_soil,
  name = "Soil_Substrate_vs_Field.soil",
  cooksCutoff = FALSE
)


res_soil_sig <- res_soil %>%
  as.data.frame() %>%
  filter(!is.na(padj) & padj < 0.05)


res_soil_sig <- cbind(
  res_soil_sig,
  as(tax_table(physeq_soil)[rownames(res_soil_sig), ], "matrix")
)

str(res_soil_sig)
# Result interpretation
# 
# log2FoldChange > 0
# → taxon enriched in Substrate
# 
# log2FoldChange < 0
# → taxon enriched in Field soil

table(sample_data(physeq_soil)$Soil)

summary(res_soil)




##Zeigt jetzt alle verschiedene ASV von allen Genotypen und alle Behandlungen
ggplot(res_soil_sig, aes(x = reorder(Genus, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE), 
                width = 0.2, color = "black") +
  coord_flip() +
  labs(title = paste("Differentially Abundant Taxa:"),
       x = "Genus", y = "Log2 Fold Change") +
  scale_fill_manual(values = c("red", "blue"), labels = c("Decreased", "Increased")) +
  theme_minimal() +
  theme(legend.position = "none")


# PIPELINE 2 — Genotype effect on the microbiome
# 
# Question:
#   Do Golden Promise and BCC436 differ in their microbiome under identical conditions?


SOIL_SYSTEM <- "Substrate"          # or "JKI"
FIXED_TREATMENT <- "mock"
FIXED_DURATION <- "Two-week-old"

physeq_geno <- physeq_raref_FE_wo_bulk

physeq_geno <- subset_samples(
  physeq_geno,
  Treatment == FIXED_TREATMENT &
    Duration  == FIXED_DURATION &
    Genotype %in% c("Golden_Promise", "BCC436")
)

physeq_geno <- prune_taxa(taxa_sums(physeq_geno) > 10, physeq_geno)

sample_data(physeq_geno)$Genotype <- factor(
  sample_data(physeq_geno)$Genotype,
  levels = c("BCC436", "Golden_Promise")  # reference = BCC436
)

dds_geno <- phyloseq_to_deseq2(physeq_geno, ~ Genotype)

dds_geno <- DESeq(
  dds_geno,
  test = "Wald",
  fitType = "parametric",
  sfType = "poscounts"
)

resultsNames(dds_geno)

res_geno <- results(
  dds_geno,
  name = "Genotype_Golden_Promise_vs_BCC436",
  cooksCutoff = FALSE
)

res_geno_sig <- res_geno %>%
  as.data.frame() %>%
  filter(!is.na(padj) & padj < 0.05)

res_geno_sig <- cbind(
  res_geno_sig,
  as(tax_table(physeq_geno)[rownames(res_geno_sig), ], "matrix")
)

# Interpretation
# 
# log2FoldChange > 0
# → taxon enriched in Golden Promise
# 
# log2FoldChange < 0
# → taxon enriched in BCC436

table(sample_data(physeq_geno)$Genotype)
table(sample_data(physeq_geno)$Treatment)
table(sample_data(physeq_geno)$Duration)



library(tibble)

res_geno_sig_wo_rownames <- res_geno_sig %>%
  rownames_to_column(var = "ASV")


ggplot(res_geno_sig_wo_rownames, aes(x = reorder(ASV, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE), 
                width = 0.2, color = "black") +
  coord_flip() +
  labs(title = paste("Differentially Abundant Taxa:"),
       x = "Genus", y = "Log2 Fold Change") +
  scale_fill_manual(values = c("red", "blue"), labels = c("Decreased", "Increased")) +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(res_geno_sig_wo_rownames, aes(x = reorder(ASV, log2FoldChange), y = log2FoldChange, fill = Phylum)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE), 
                width = 0.2, color = "black") +
  coord_flip() +
  labs(title = paste("Differentially Abundant Taxa:"),
       x = "Genus", y = "Log2 Fold Change") +
  #scale_fill_manual(values = c("red", "blue"), labels = c("Decreased", "Increased")) +
  theme_minimal() 

res_geno_sig_wo_rownames <- res_geno_sig_wo_rownames %>% 
  mutate(ASV_Phylum = paste(ASV, Phylum, sep = "_"))



ggplot(res_geno_sig_wo_rownames, aes(x = reorder(ASV_Phylum, log2FoldChange), y = log2FoldChange, fill = Phylum)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE), 
                width = 0.2, color = "black") +
  coord_flip() +
  labs(title = paste("Differentially Abundant Taxa:"),
       x = "Genus", y = "Log2 Fold Change") +
  #scale_fill_manual(values = c("red", "blue"), labels = c("Decreased", "Increased")) +
  theme_minimal() 






dif<-ggplot(res_geno_sig_wo_rownames, aes(log2FoldChange, reorder(ASV_Phylum, log2FoldChange), color = Phylum ))+
  geom_point() +theme_bw() +
  geom_segment( aes(x=0, xend=log2FoldChange, y=ASV_Phylum, yend=ASV_Phylum, color = Phylum)) +
  geom_vline(xintercept = 0, size=0.3) +xlab("log2FoldChange") +ylab(NULL)


PS_asv_rel <- transform_sample_counts(
  physeq_geno,   # SAME object used for DESeq2
  function(x) x / sum(x)
)

relabu_asv <- psmelt(PS_asv_rel)
relabu_asv$Abundance <- relabu_asv$Abundance * 100
length(unique(relabu_asv$OTU))
# hundreds/thousands (correct)



str(res_geno_sig)


sig_asvs <- rownames(res_geno_sig)
length(sig_asvs)
head(sig_asvs)



SOIL_SYSTEM <- "Substrate"          # or "JKI"
FIXED_TREATMENT <- "mock"
FIXED_DURATION <- "Two-week-old"


relabu_sig <- relabu_sig %>%
  filter(
    Soil == SOIL_SYSTEM,
    Duration == FIXED_DURATION,
    Treatment == FIXED_TREATMENT
    # add Genotype / Treatment filters if needed
  )


relabu_sig <- relabu_asv %>%
  filter(OTU %in% sig_asvs)

unique(relabu_sig$OTU)

library(dplyr)

res_geno_sig2 <- res_geno_sig %>%
  rownames_to_column(var = "ASV")  # adds rownames as a column

asv_order <- res_geno_sig2 %>%
  arrange(log2FoldChange) %>%
  pull(ASV)

relabu_sig$OTU <- factor(relabu_sig$OTU, levels = asv_order)



ab<-ggplot(relabu_sig, aes(Abundance, OTU, fill = Genotype)) +
  geom_boxplot(outlier.size = 0.5) + theme_bw() + xlab("Relative abundance") +
  theme( axis.text.y=element_blank()) + ylab(NULL)



ggarrange(dif,ab, ncol = 2, nrow = 1,align = "h", common.legend = TRUE, widths = c(1.85, 0.9))


#### GP 2 weeks ---------------------------
# PIPELINE 2 — Genotype effect on the microbiome 
# 
# Question:
#   Do Golden Promise and BCC436 differ in their microbiome under identical conditions?


SOIL_SYSTEM <- "FE"          # or "JKI"
FIXED_GENOTYPE <- "Golden_Promise"
FIXED_DURATION <- "Two-week-old"

physeq_treat_GP <- physeq_raref_FE_wo_bulk

physeq_treat_GP <- subset_samples(
  physeq_treat_GP,
  Genotype == FIXED_GENOTYPE &
    Duration  == FIXED_DURATION &
    Treatment %in% c("mock", "C6-HSL")
)


otu <- as.data.frame(otu_table(physeq_treat_GP))

sample_data(physeq_treat_GP)$Treatment <- factor(
  sample_data(physeq_treat_GP)$Treatment,
  levels = c("mock", "C6-HSL")
)

meta <- as.data.frame(sample_data(physeq_treat_GP))

ntaxa(physeq_treat_GP)


if (!taxa_are_rows(physeq_treat_GP)) {
  otu <- t(otu)
}

prev_by_treat <- sapply(levels(meta$Treatment), function(tr) {
  samples_tr <- rownames(meta)[meta$Treatment == tr]
  rowSums(otu[, samples_tr] > 0)
})

keep_asvs <- apply(prev_by_treat, 1, function(x) any(x >= 3))

keep_asvs_final <- keep_asvs & (taxa_sums(physeq_treat_GP) > 10)
physeq_treat_GP <- prune_taxa(keep_asvs_final, physeq_treat_GP)




colnames(sample_data(physeq_treat_GP))


dds_treat_GP <- phyloseq_to_deseq2(physeq_treat_GP, ~ Treatment)

dds_treat_GP <- DESeq(
  dds_treat_GP,
  test = "Wald",
  fitType = "parametric",
  sfType = "poscounts"
)

resultsNames(dds_treat_GP)

res_treat_GP <- results(
  dds_treat_GP,
  name = "Treatment_C6.HSL_vs_mock",
  cooksCutoff = FALSE
)


##speichern für unten
res_treat_GP_2weeks<-res_treat_GP


res_treat_GP_sig <- res_treat_GP %>%
  as.data.frame() %>%
  filter(!is.na(padj) & padj < 0.05)

res_treat_GP_sig <- cbind(
  res_treat_GP_sig,
  as(tax_table(physeq_treat_GP)[rownames(res_treat_GP_sig), ], "matrix")
)

# Interpretation
# 
# log2FoldChange > 0 → enriched in C6-HSL
# log2FoldChange < 0 → enriched in mock


table(sample_data(physeq_treat_GP)$Genotype)
table(sample_data(physeq_treat_GP)$Treatment)
table(sample_data(physeq_treat_GP)$Duration)



library(tibble)

res_treat_GP_sig_wo_rownames <- res_treat_GP_sig %>%
  rownames_to_column(var = "ASV")


ggplot(res_treat_GP_sig_wo_rownames, aes(x = reorder(ASV, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE), 
                width = 0.2, color = "black") +
  coord_flip() +
  labs(title = paste("Differentially Abundant Taxa:"),
       x = "ASV", y = "Log2 Fold Change") +
  scale_fill_manual(values = c("red", "blue"), labels = c("Decreased", "Increased")) +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(res_treat_GP_sig_wo_rownames, aes(x = reorder(ASV, log2FoldChange), y = log2FoldChange, fill = Phylum)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE), 
                width = 0.2, color = "black") +
  coord_flip() +
  labs(title = paste("Differentially Abundant Taxa:"),
       x = "Genus", y = "Log2 Fold Change") +
  #scale_fill_manual(values = c("red", "blue"), labels = c("Decreased", "Increased")) +
  theme_minimal() 

res_treat_GP_sig_wo_rownames <- res_treat_GP_sig_wo_rownames %>% 
  mutate(ASV_Phylum = paste(ASV, Phylum, sep = "_"))



ggplot(res_treat_GP_sig_wo_rownames, aes(x = reorder(ASV_Phylum, log2FoldChange), y = log2FoldChange, fill = Phylum)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE), 
                width = 0.2, color = "black") +
  coord_flip() +
  labs(title = paste("Differentially Abundant Taxa:"),
       x = "Genus", y = "Log2 Fold Change") +
  #scale_fill_manual(values = c("red", "blue"), labels = c("Decreased", "Increased")) +
  theme_minimal() 



##speichern für unten
res_treat_GP_sig_wo_rownames_2weeks<-res_treat_GP_sig_wo_rownames



palette_24  <- c(
  "#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD",
  "#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF",
  "#393B79","#637939","#8C6D31","#843C39","#7B4173",
  "#3182BD","#31A354","#756BB1","#636363","#E6550D",
  "#969696","#CEDB9C","#9EDAE5","#FDD0A2","#A6CEE3",
  "#B2DF8A","#FB9A99","#FDBF6F","#CAB2D6","#FFFF99",
  "#6A3D9A","#33A02C","#E31A1C","#FF7F00","#1B9E77",
  "#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
  "#A6761D","#666666","#8DD3C7","#FCCDE5","#B3DE69",
  "#BC80BD","#CCEBC5","#FFED6F","#4E79A7","#59A14F",
  "#9C755F","#EDC948","#AF7AA1","#76B7B2"
)


dif<-ggplot(res_treat_GP_sig_wo_rownames, aes(log2FoldChange, reorder(ASV, log2FoldChange), color = Genus ))+
  geom_point() +theme_bw() +
  geom_segment( aes(x=0, xend=log2FoldChange, y=ASV, yend=ASV, color = Genus)) +
  geom_vline(xintercept = 0, size=0.3) +xlab("log2FoldChange") +
  ylab(NULL)+
  scale_color_manual(values = palette_24)


PS_asv_rel <- transform_sample_counts(
  physeq_treat_GP,   # SAME object used for DESeq2
  function(x) x / sum(x)
)

relabu_asv <- psmelt(PS_asv_rel)
relabu_asv$Abundance <- relabu_asv$Abundance * 100
length(unique(relabu_asv$OTU))
# hundreds/thousands (correct)



str(res_treat_GP_sig)


sig_asvs <- rownames(res_treat_GP_sig)
length(sig_asvs)
head(sig_asvs)



SOIL_SYSTEM <- "FE"          # or "JKI"
FIXED_GENOTYPE <- "Golden_Promise"
FIXED_DURATION <- "Two-week-old"


relabu_sig <- relabu_sig %>%
  filter(
    Soil == SOIL_SYSTEM,
    Duration == FIXED_DURATION,
    Genotype == FIXED_GENOTYPE
    # add Genotype / Treatment filters if needed
  )


relabu_sig <- relabu_asv %>%
  filter(OTU %in% sig_asvs)

unique(relabu_sig$OTU)

library(dplyr)

res_treat_GP_sig2 <- res_treat_GP_sig %>%
  rownames_to_column(var = "ASV")  # adds rownames as a column

asv_order <- res_treat_GP_sig2 %>%
  arrange(log2FoldChange) %>%
  pull(ASV)

relabu_sig$OTU <- factor(relabu_sig$OTU, levels = asv_order)



ab<-ggplot(relabu_sig, aes(Abundance, OTU, fill = Treatment)) +
  geom_boxplot(outlier.size = 0.5) + theme_bw() + xlab("Relative abundance") +
  ggtitle ("GP two weeks")+
  theme( axis.text.y=element_blank()) + ylab(NULL)



ggarrange(dif,ab, ncol = 2, nrow = 1,align = "h", common.legend = F, widths = c(1.85, 0.9),  legend = "bottom"  )




####GP 16 weeks ---------------------------------------------------

SOIL_SYSTEM <- "FE"          # or "JKI"
FIXED_GENOTYPE <- "Golden_Promise"
FIXED_DURATION <- "16-week-old"


physeq_treat_GP <- physeq_raref_FE_wo_bulk

physeq_treat_GP <- subset_samples(
  physeq_treat_GP,
  Genotype == FIXED_GENOTYPE &
    Duration  == FIXED_DURATION &
    Treatment %in% c("mock", "C6-HSL")
)


otu <- as.data.frame(otu_table(physeq_treat_GP))

sample_data(physeq_treat_GP)$Treatment <- factor(
  sample_data(physeq_treat_GP)$Treatment,
  levels = c("mock", "C6-HSL")
)

meta <- as.data.frame(sample_data(physeq_treat_GP))

ntaxa(physeq_treat_GP)


if (!taxa_are_rows(physeq_treat_GP)) {
  otu <- t(otu)
}

prev_by_treat <- sapply(levels(meta$Treatment), function(tr) {
  samples_tr <- rownames(meta)[meta$Treatment == tr]
  rowSums(otu[, samples_tr] > 0)
})

keep_asvs <- apply(prev_by_treat, 1, function(x) any(x >= 3))

keep_asvs_final <- keep_asvs & (taxa_sums(physeq_treat_GP) > 10)
physeq_treat_GP <- prune_taxa(keep_asvs_final, physeq_treat_GP)




colnames(sample_data(physeq_treat_GP))


dds_treat_GP <- phyloseq_to_deseq2(physeq_treat_GP, ~ Treatment)


dds_treat_GP <- DESeq(
  dds_treat_GP,
  test = "Wald",
  fitType = "parametric",
  sfType = "poscounts"
)

resultsNames(dds_treat_GP)

res_treat_GP <- results(
  dds_treat_GP,
  name = "Treatment_C6.HSL_vs_mock",
  cooksCutoff = FALSE
)



##speichern für unten
res_treat_GP_16weeks<-res_treat_GP

res_treat_GP_sig <- res_treat_GP %>%
  as.data.frame() %>%
  filter(!is.na(padj) & padj < 0.05)

res_treat_GP_sig <- cbind(
  res_treat_GP_sig,
  as(tax_table(physeq_treat_GP)[rownames(res_treat_GP_sig), ], "matrix")
)

# Interpretation
# Positiver log2FoldChange
# log2FC > 0
# 
# 
# → ASV ist angereichert unter C6-HSL
# → höhere relative (normalisierte) Abundanz im C6-HSL-Treatment
# → „C6-HSL fördert diese ASV“
# 
# Negativer log2FoldChange
# log2FC < 0
# 
# 
# → ASV ist reduziert unter C6-HSL
# → höhere Abundanz im mock
# → „C6-HSL unterdrückt / verdrängt diese ASV“



table(sample_data(physeq_treat_GP)$Genotype)
table(sample_data(physeq_treat_GP)$Treatment)
table(sample_data(physeq_treat_GP)$Duration)



library(tibble)

res_treat_GP_sig_wo_rownames <- res_treat_GP_sig %>%
  rownames_to_column(var = "ASV")


ggplot(res_treat_GP_sig_wo_rownames, aes(x = reorder(ASV, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE), 
                width = 0.2, color = "black") +
  coord_flip() +
  labs(title = paste("Differentially Abundant Taxa:"),
       x = "Genus", y = "Log2 Fold Change") +
  scale_fill_manual(values = c("red", "blue"), labels = c("Decreased", "Increased")) +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(res_treat_GP_sig_wo_rownames, aes(x = reorder(ASV, log2FoldChange), y = log2FoldChange, fill = Phylum)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE), 
                width = 0.2, color = "black") +
  coord_flip() +
  labs(title = paste("Differentially Abundant Taxa:"),
       x = "Genus", y = "Log2 Fold Change") +
  #scale_fill_manual(values = c("red", "blue"), labels = c("Decreased", "Increased")) +
  theme_minimal() 

res_treat_GP_sig_wo_rownames <- res_treat_GP_sig_wo_rownames %>% 
  mutate(ASV_Phylum = paste(ASV, Phylum, sep = "_"))



ggplot(res_treat_GP_sig_wo_rownames, aes(x = reorder(ASV_Phylum, log2FoldChange), y = log2FoldChange, fill = Phylum)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE), 
                width = 0.2, color = "black") +
  coord_flip() +
  labs(title = paste("Differentially Abundant Taxa:"),
       x = "Genus", y = "Log2 Fold Change") +
  #scale_fill_manual(values = c("red", "blue"), labels = c("Decreased", "Increased")) +
  theme_minimal() 



##speichern für unten
res_treat_GP_sig_wo_rownames_16weeks<-res_treat_GP_sig_wo_rownames

dif<-ggplot(res_treat_GP_sig_wo_rownames, aes(log2FoldChange, reorder(ASV, log2FoldChange), color = Genus ))+
  geom_point() +theme_bw() +
  geom_segment( aes(x=0, xend=log2FoldChange, y=ASV, yend=ASV, color = Genus)) +
  geom_vline(xintercept = 0, size=0.3) +xlab("log2FoldChange") +
  ylab(NULL)+
  scale_color_manual(values = palette_24)


PS_asv_rel <- transform_sample_counts(
  physeq_treat_GP,   # SAME object used for DESeq2
  function(x) x / sum(x)
)

relabu_asv <- psmelt(PS_asv_rel)
relabu_asv$Abundance <- relabu_asv$Abundance * 100
length(unique(relabu_asv$OTU))
# hundreds/thousands (correct)



str(res_treat_GP_sig)


sig_asvs <- rownames(res_treat_GP_sig)
length(sig_asvs)
head(sig_asvs)



SOIL_SYSTEM <- "FE"          # or "JKI"
FIXED_GENOTYPE <- "Golden_Promise"
FIXED_DURATION <- "16-week-old"


relabu_sig <- relabu_sig %>%
  filter(
    Soil == SOIL_SYSTEM,
    Duration == FIXED_DURATION,
    Genotype == FIXED_GENOTYPE
    # add Genotype / Treatment filters if needed
  )


relabu_sig <- relabu_asv %>%
  filter(OTU %in% sig_asvs)

unique(relabu_sig$OTU)

library(dplyr)

res_treat_GP_sig2 <- res_treat_GP_sig %>%
  rownames_to_column(var = "ASV")  # adds rownames as a column

asv_order <- res_treat_GP_sig2 %>%
  arrange(log2FoldChange) %>%
  pull(ASV)

relabu_sig$OTU <- factor(relabu_sig$OTU, levels = asv_order)



ab<-ggplot(relabu_sig, aes(Abundance, OTU, fill = Treatment)) +
  geom_boxplot(outlier.size = 0.5) + theme_bw() + xlab("Relative abundance") +
  ggtitle ("GP 16 weeks")+
  theme( axis.text.y=element_blank()) + ylab(NULL)



ggarrange(dif,ab, ncol = 2, nrow = 1,align = "h", common.legend = F, widths = c(1.85, 0.9),  legend = "bottom"  )


#### BCC 2 weeks ---------------------------
# PIPELINE 2 — Genotype effect on the microbiome 
# 
# Question:
#   Do Golden Promise and BCC436 differ in their microbiome under identical conditions?


SOIL_SYSTEM <- "FE"          # or "JKI"
FIXED_GENOTYPE <- "BCC436"
FIXED_DURATION <- "Two-week-old"

physeq_treat_BCC <- physeq_raref_FE_wo_bulk

physeq_treat_BCC <- subset_samples(
  physeq_treat_BCC,
  Genotype == FIXED_GENOTYPE &
    Duration  == FIXED_DURATION &
    Treatment %in% c("mock", "C6-HSL")
)

otu <- as.data.frame(otu_table(physeq_treat_BCC))

sample_data(physeq_treat_BCC)$Treatment <- factor(
  sample_data(physeq_treat_BCC)$Treatment,
  levels = c("mock", "C6-HSL")
)

meta <- as.data.frame(sample_data(physeq_treat_BCC))

ntaxa(physeq_treat_BCC)


if (!taxa_are_rows(physeq_treat_BCC)) {
  otu <- t(otu)
}

prev_by_treat <- sapply(levels(meta$Treatment), function(tr) {
  samples_tr <- rownames(meta)[meta$Treatment == tr]
  rowSums(otu[, samples_tr] > 0)
})

keep_asvs <- apply(prev_by_treat, 1, function(x) any(x >= 3))

keep_asvs_final <- keep_asvs & (taxa_sums(physeq_treat_BCC) > 10)
physeq_treat_BCC <- prune_taxa(keep_asvs_final, physeq_treat_BCC)




colnames(sample_data(physeq_treat_BCC))


dds_treat_BCC <- phyloseq_to_deseq2(physeq_treat_BCC, ~ Treatment)
dds_treat_BCC <- DESeq(
  dds_treat_BCC,
  test = "Wald",
  fitType = "parametric",
  sfType = "poscounts"
)

resultsNames(dds_treat_BCC)

res_treat_BCC <- results(
  dds_treat_BCC,
  name = "Treatment_C6.HSL_vs_mock",
  cooksCutoff = FALSE
)

##speichern für unten
res_treat_BCC_2weeks<-res_treat_BCC


res_treat_BCC_sig <- res_treat_BCC %>%
  as.data.frame() %>%
  filter(!is.na(padj) & padj < 0.05)

res_treat_BCC_sig <- cbind(
  res_treat_BCC_sig,
  as(tax_table(physeq_treat_BCC)[rownames(res_treat_BCC_sig), ], "matrix")
)

# Interpretation
# 
# log2FoldChange > 0
# → taxon enriched in Golden Promise
# 
# log2FoldChange < 0
# → taxon enriched in BCC436

table(sample_data(physeq_treat_BCC)$Genotype)
table(sample_data(physeq_treat_BCC)$Treatment)
table(sample_data(physeq_treat_BCC)$Duration)



library(tibble)

res_treat_BCC_sig_wo_rownames <- res_treat_BCC_sig %>%
  rownames_to_column(var = "ASV")


ggplot(res_treat_BCC_sig_wo_rownames, aes(x = reorder(ASV, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE), 
                width = 0.2, color = "black") +
  coord_flip() +
  labs(title = paste("Differentially Abundant Taxa:"),
       x = "Genus", y = "Log2 Fold Change") +
  scale_fill_manual(values = c("red", "blue"), labels = c("Decreased", "Increased")) +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(res_treat_BCC_sig_wo_rownames, aes(x = reorder(ASV, log2FoldChange), y = log2FoldChange, fill = Phylum)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE), 
                width = 0.2, color = "black") +
  coord_flip() +
  labs(title = paste("Differentially Abundant Taxa:"),
       x = "Genus", y = "Log2 Fold Change") +
  #scale_fill_manual(values = c("red", "blue"), labels = c("Decreased", "Increased")) +
  theme_minimal() 

res_treat_BCC_sig_wo_rownames <- res_treat_BCC_sig_wo_rownames %>% 
  mutate(ASV_Phylum = paste(ASV, Phylum, sep = "_"))



ggplot(res_treat_BCC_sig_wo_rownames, aes(x = reorder(ASV_Phylum, log2FoldChange), y = log2FoldChange, fill = Phylum)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE), 
                width = 0.2, color = "black") +
  coord_flip() +
  labs(title = paste("Differentially Abundant Taxa:"),
       x = "Genus", y = "Log2 Fold Change") +
  #scale_fill_manual(values = c("red", "blue"), labels = c("Decreased", "Increased")) +
  theme_minimal() 



##speichern für unten
res_treat_BCC_sig_wo_rownames_2weeks<-res_treat_BCC_sig_wo_rownames


dif<-ggplot(res_treat_BCC_sig_wo_rownames, aes(log2FoldChange, reorder(ASV, log2FoldChange), color = Genus ))+
  geom_point() +theme_bw() +
  geom_segment( aes(x=0, xend=log2FoldChange, y=ASV, yend=ASV, color = Genus)) +
  geom_vline(xintercept = 0, size=0.3) +xlab("log2FoldChange") +
  ylab(NULL)+
  scale_color_manual(values = palette_24)

PS_asv_rel <- transform_sample_counts(
  physeq_treat_BCC,   # SAME object used for DESeq2
  function(x) x / sum(x)
)

relabu_asv <- psmelt(PS_asv_rel)
relabu_asv$Abundance <- relabu_asv$Abundance * 100
length(unique(relabu_asv$OTU))
# hundreds/thousands (correct)



str(res_treat_BCC_sig)


sig_asvs <- rownames(res_treat_BCC_sig)
length(sig_asvs)
head(sig_asvs)



SOIL_SYSTEM <- "FE"          # or "JKI"
FIXED_GENOTYPE <- "Golden_Promise"
FIXED_DURATION <- "Two-week-old"


relabu_sig <- relabu_sig %>%
  filter(
    Soil == SOIL_SYSTEM,
    Duration == FIXED_DURATION,
    Genotype == FIXED_GENOTYPE
    # add Genotype / Treatment filters if needed
  )


relabu_sig <- relabu_asv %>%
  filter(OTU %in% sig_asvs)

unique(relabu_sig$OTU)

library(dplyr)

res_treat_BCC_sig2 <- res_treat_BCC_sig %>%
  rownames_to_column(var = "ASV")  # adds rownames as a column

asv_order <- res_treat_BCC_sig2 %>%
  arrange(log2FoldChange) %>%
  pull(ASV)

relabu_sig$OTU <- factor(relabu_sig$OTU, levels = asv_order)



ab<-ggplot(relabu_sig, aes(Abundance, OTU, fill = Treatment)) +
  geom_boxplot(outlier.size = 0.5) + theme_bw() + xlab("Relative abundance") +
  ggtitle ("BCC two weeks")+
  theme( axis.text.y=element_blank()) + ylab(NULL)



ggarrange(dif,ab, ncol = 2, nrow = 1,align = "h", common.legend = F, widths = c(1.85, 0.9),  legend = "bottom"  )




####BCC 16 weeks ---------------------------------------------------

SOIL_SYSTEM <- "FE"          # or "JKI"
FIXED_GENOTYPE <- "BCC436"
FIXED_DURATION <- "16-week-old"

physeq_treat_BCC <- physeq_raref_FE_wo_bulk

physeq_treat_BCC <- subset_samples(
  physeq_treat_BCC,
  Genotype == FIXED_GENOTYPE &
    Duration  == FIXED_DURATION &
    Treatment %in% c("mock", "C6-HSL")
)

otu <- as.data.frame(otu_table(physeq_treat_BCC))

sample_data(physeq_treat_BCC)$Treatment <- factor(
  sample_data(physeq_treat_BCC)$Treatment,
  levels = c("mock", "C6-HSL")
)

meta <- as.data.frame(sample_data(physeq_treat_BCC))

ntaxa(physeq_treat_BCC)


if (!taxa_are_rows(physeq_treat_BCC)) {
  otu <- t(otu)
}

prev_by_treat <- sapply(levels(meta$Treatment), function(tr) {
  samples_tr <- rownames(meta)[meta$Treatment == tr]
  rowSums(otu[, samples_tr] > 0)
})

keep_asvs <- apply(prev_by_treat, 1, function(x) any(x >= 3))

keep_asvs_final <- keep_asvs & (taxa_sums(physeq_treat_BCC) > 10)
physeq_treat_BCC <- prune_taxa(keep_asvs_final, physeq_treat_BCC)




colnames(sample_data(physeq_treat_BCC))


dds_treat_BCC <- phyloseq_to_deseq2(physeq_treat_BCC, ~ Treatment)

dds_treat_BCC <- DESeq(
  dds_treat_BCC,
  test = "Wald",
  fitType = "parametric",
  sfType = "poscounts"
)

resultsNames(dds_treat_BCC)

res_treat_BCC <- results(
  dds_treat_BCC,
  name = "Treatment_C6.HSL_vs_mock",
  cooksCutoff = FALSE
)

##speichern für unten
res_treat_BCC_16weeks<-res_treat_BCC


res_treat_BCC_sig <- res_treat_BCC %>%
  as.data.frame() %>%
  filter(!is.na(padj) & padj < 0.05)

res_treat_BCC_sig <- cbind(
  res_treat_BCC_sig,
  as(tax_table(physeq_treat_BCC)[rownames(res_treat_BCC_sig), ], "matrix")
)

# Interpretation
# Positiver log2FoldChange
# log2FC > 0
# 
# 
# → ASV ist angereichert unter C6-HSL
# → höhere relative (normalisierte) Abundanz im C6-HSL-Treatment
# → „C6-HSL fördert diese ASV“
# 
# Negativer log2FoldChange
# log2FC < 0
# 
# 
# → ASV ist reduziert unter C6-HSL
# → höhere Abundanz im mock
# → „C6-HSL unterdrückt / verdrängt diese ASV“



table(sample_data(physeq_treat_BCC)$Genotype)
table(sample_data(physeq_treat_BCC)$Treatment)
table(sample_data(physeq_treat_BCC)$Duration)



library(tibble)

res_treat_BCC_sig_wo_rownames <- res_treat_BCC_sig %>%
  rownames_to_column(var = "ASV")


ggplot(res_treat_BCC_sig_wo_rownames, aes(x = reorder(ASV, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE), 
                width = 0.2, color = "black") +
  coord_flip() +
  labs(title = paste("Differentially Abundant Taxa:"),
       x = "Genus", y = "Log2 Fold Change") +
  scale_fill_manual(values = c("red", "blue"), labels = c("Decreased", "Increased")) +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(res_treat_BCC_sig_wo_rownames, aes(x = reorder(ASV, log2FoldChange), y = log2FoldChange, fill = Phylum)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE), 
                width = 0.2, color = "black") +
  coord_flip() +
  labs(title = paste("Differentially Abundant Taxa:"),
       x = "Genus", y = "Log2 Fold Change") +
  #scale_fill_manual(values = c("red", "blue"), labels = c("Decreased", "Increased")) +
  theme_minimal() 

res_treat_BCC_sig_wo_rownames <- res_treat_BCC_sig_wo_rownames %>% 
  mutate(ASV_Phylum = paste(ASV, Phylum, sep = "_"))



ggplot(res_treat_BCC_sig_wo_rownames, aes(x = reorder(ASV_Phylum, log2FoldChange), y = log2FoldChange, fill = Phylum)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE), 
                width = 0.2, color = "black") +
  coord_flip() +
  labs(title = paste("Differentially Abundant Taxa:"),
       x = "Genus", y = "Log2 Fold Change") +
  #scale_fill_manual(values = c("red", "blue"), labels = c("Decreased", "Increased")) +
  theme_minimal() 



##speichern für unten
res_treat_BCC_sig_wo_rownames_16weeks<-res_treat_BCC_sig_wo_rownames

dif<-ggplot(res_treat_BCC_sig_wo_rownames, aes(log2FoldChange, reorder(ASV, log2FoldChange), color = Genus ))+
  geom_point() +theme_bw() +
  geom_segment( aes(x=0, xend=log2FoldChange, y=ASV, yend=ASV, color = Genus)) +
  geom_vline(xintercept = 0, size=0.3) +xlab("log2FoldChange") +
  ylab(NULL)+
  scale_color_manual(values = palette_24)


PS_asv_rel <- transform_sample_counts(
  physeq_treat_BCC,   # SAME object used for DESeq2
  function(x) x / sum(x)
)

relabu_asv <- psmelt(PS_asv_rel)
relabu_asv$Abundance <- relabu_asv$Abundance * 100
length(unique(relabu_asv$OTU))
# hundreds/thousands (correct)



str(res_treat_BCC_sig)


sig_asvs <- rownames(res_treat_BCC_sig)
length(sig_asvs)
head(sig_asvs)



SOIL_SYSTEM <- "FE"          # or "JKI"
FIXED_GENOTYPE <- "Golden_Promise"
FIXED_DURATION <- "16-week-old"


relabu_sig <- relabu_sig %>%
  filter(
    Soil == SOIL_SYSTEM,
    Duration == FIXED_DURATION,
    Genotype == FIXED_GENOTYPE
    # add Genotype / Treatment filters if needed
  )


relabu_sig <- relabu_asv %>%
  filter(OTU %in% sig_asvs)

unique(relabu_sig$OTU)

library(dplyr)

res_treat_BCC_sig2 <- res_treat_BCC_sig %>%
  rownames_to_column(var = "ASV")  # adds rownames as a column

asv_order <- res_treat_BCC_sig2 %>%
  arrange(log2FoldChange) %>%
  pull(ASV)

relabu_sig$OTU <- factor(relabu_sig$OTU, levels = asv_order)



ab<-ggplot(relabu_sig, aes(Abundance, OTU, fill = Treatment)) +
  geom_boxplot(outlier.size = 0.5) + theme_bw() + xlab("Relative abundance") +
  ggtitle ("BCC 16 weeks")+
  theme( axis.text.y=element_blank()) + ylab(NULL)



ggarrange(dif,ab, ncol = 2, nrow = 1,align = "h", common.legend = F, widths = c(1.85, 0.9),  legend = "bottom"  )




# Volcano plots -----------------------------------------------------------


library(ggplot2)
library(dplyr)
library(ggrepel)


res_treat_GP_2weeks
res_treat_GP_16weeks

res_treat_BCC_2weeks
res_treat_BCC_16weeks

##GP 2 weeks

# 1. DESeq2 Ergebnis in data.frame umwandeln
res_df <- as.data.frame(res_treat_GP_2weeks)

# 2. Nur vollständige Werte verwenden
res_df <- res_df %>%
  filter(!is.na(log2FoldChange) & !is.na(padj) & 
           is.finite(log2FoldChange) & is.finite(padj))

# 3. Status nach Signifikanz und Richtung festlegen
res_df <- res_df %>%
  mutate(
    status = case_when(
      padj >= 0.05 | abs(log2FoldChange) < 1 ~ "ns",   # nicht signifikant oder kleine Änderung
      log2FoldChange > 0 & padj < 0.05 & log2FoldChange >= 1 ~ "up",
      log2FoldChange < 0 & padj < 0.05 & log2FoldChange <= -1 ~ "down"
    )
  )

# 4. Top 3 positive und negative Veränderungen auswählen
top_up <- res_df %>%
  filter(status == "up") %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 3)

top_down <- res_df %>%
  filter(status == "down") %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 3)

top_genes_GP_2wk <- bind_rows(top_up, top_down)

# 5. Volcano-Plot erstellen
plot_GP_2wk<-ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = status), alpha = 0.6) +
  scale_color_manual(values = c("ns" = "grey", "down" = "red", "up" = "blue")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = top_genes_GP_2wk,
    aes(label = rownames(top_genes_GP_2wk)),
    nudge_x = ifelse(top_genes_GP_2wk$log2FoldChange > 0, 0.8, -0.8),
    direction = "y",
    force = 2,
    box.padding = 0.6,
    point.padding = 0.4,
    max.overlaps = Inf,
    segment.color = "black",
    segment.size = 0.3,
    min.segment.length = 0,
    size = 3
  )+
  theme_bw() +
  theme(
    legend.text = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(colour = "black", size = 10),
    axis.title.y = element_text( size = 10, colour = "black"),
    axis.text.x = element_text(vjust = 0.5, size = 10, colour = "black"),
    axis.text.y = element_text(vjust = 0.5, size = 10, colour = "black"),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank(),
    plot.title = element_text(size = 10, color = "#000", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 12, face = "bold", colour = "black")  # Legend title
  )+
  labs(
    x = "log2 fold change (C6-HSL vs mock)",
    y = "-log10(adjusted p-value)",
    title = "Golden Promise 2 weeks-old"
  )



##GP 16 weeks

# 1. DESeq2 Ergebnis in data.frame umwandeln
res_df <- as.data.frame(res_treat_GP_16weeks)

# 2. Nur vollständige Werte verwenden
res_df <- res_df %>%
  filter(!is.na(log2FoldChange) & !is.na(padj) & 
           is.finite(log2FoldChange) & is.finite(padj))

# 3. Status nach Signifikanz und Richtung festlegen
res_df <- res_df %>%
  mutate(
    status = case_when(
      padj >= 0.05 | abs(log2FoldChange) < 1 ~ "ns",   # nicht signifikant oder kleine Änderung
      log2FoldChange > 0 & padj < 0.05 & log2FoldChange >= 1 ~ "up",
      log2FoldChange < 0 & padj < 0.05 & log2FoldChange <= -1 ~ "down"
    )
  )

# 4. Top 3 positive und negative Veränderungen auswählen
top_up <- res_df %>%
  filter(status == "up") %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 3)

top_down <- res_df %>%
  filter(status == "down") %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 3)

top_genes_GP_16wk <- bind_rows(top_up, top_down)

# 5. Volcano-Plot erstellen
plot_GP_16wk<-ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = status), alpha = 0.6) +
  scale_color_manual(values = c("ns" = "grey", "down" = "red", "up" = "blue")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = top_genes_GP_16wk,
    aes(label = rownames(top_genes_GP_16wk)),
    nudge_x = ifelse(top_genes_GP_16wk$log2FoldChange > 0, 0.8, -0.8),
    direction = "y",
    force = 2,
    box.padding = 0.6,
    point.padding = 0.4,
    max.overlaps = Inf,
    segment.color = "black",
    segment.size = 0.3,
    min.segment.length = 0,
    size = 3
  )+
  theme_bw() +
  theme(
    legend.text = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(colour = "black", size = 10),
    axis.title.y = element_text( size = 10, colour = "black"),
    axis.text.x = element_text(vjust = 0.5, size = 10, colour = "black"),
    axis.text.y = element_text(vjust = 0.5, size = 10, colour = "black"),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank(),
    plot.title = element_text(size = 10, color = "#000", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 12, face = "bold", colour = "black")  # Legend title
  )+
  labs(
    x = "log2 fold change (C6-HSL vs mock)",
    y = "-log10(adjusted p-value)",
    title = "Golden Promise 16 weeks-old"
  )


##BCC436 2 weeks

# 1. DESeq2 Ergebnis in data.frame umwandeln
res_df <- as.data.frame(res_treat_BCC_2weeks)

# 2. Nur vollständige Werte verwenden
res_df <- res_df %>%
  filter(!is.na(log2FoldChange) & !is.na(padj) & 
           is.finite(log2FoldChange) & is.finite(padj))

# 3. Status nach Signifikanz und Richtung festlegen
res_df <- res_df %>%
  mutate(
    status = case_when(
      padj >= 0.05 | abs(log2FoldChange) < 1 ~ "ns",   # nicht signifikant oder kleine Änderung
      log2FoldChange > 0 & padj < 0.05 & log2FoldChange >= 1 ~ "up",
      log2FoldChange < 0 & padj < 0.05 & log2FoldChange <= -1 ~ "down"
    )
  )

# 4. Top 3 positive und negative Veränderungen auswählen
top_up <- res_df %>%
  filter(status == "up") %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 3)

top_down <- res_df %>%
  filter(status == "down") %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 3)

top_genes_BCC_2wk <- bind_rows(top_up, top_down)

# 5. Volcano-Plot erstellen
plot_BCC_2wk<-ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = status), alpha = 0.6) +
  scale_color_manual(values = c("ns" = "grey", "down" = "red", "up" = "blue")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = top_genes_BCC_2wk,
    aes(label = rownames(top_genes_BCC_2wk)),
    nudge_x = ifelse(top_genes_BCC_2wk$log2FoldChange > 0, 0.8, -0.8),
    direction = "y",
    force = 2,
    box.padding = 0.6,
    point.padding = 0.4,
    max.overlaps = Inf,
    segment.color = "black",
    segment.size = 0.3,
    min.segment.length = 0,
    size = 3
  )+
  theme_bw() +
  theme(
    legend.text = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(colour = "black", size = 10),
    axis.title.y = element_text( size = 10, colour = "black"),
    axis.text.x = element_text(vjust = 0.5, size = 10, colour = "black"),
    axis.text.y = element_text(vjust = 0.5, size = 10, colour = "black"),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank(),
    plot.title = element_text(size = 10, color = "#000", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 12, face = "bold", colour = "black")  # Legend title
  )+
  labs(
    x = "log2 fold change (C6-HSL vs mock)",
    y = "-log10(adjusted p-value)",
    title = "BCC436 2 weeks-old"
  )



##BCC 16 weeks

# 1. DESeq2 Ergebnis in data.frame umwandeln
res_df <- as.data.frame(res_treat_BCC_16weeks)

# 2. Nur vollständige Werte verwenden
res_df <- res_df %>%
  filter(!is.na(log2FoldChange) & !is.na(padj) & 
           is.finite(log2FoldChange) & is.finite(padj))

# 3. Status nach Signifikanz und Richtung festlegen
res_df <- res_df %>%
  mutate(
    status = case_when(
      padj >= 0.05 | abs(log2FoldChange) < 1 ~ "ns",   # nicht signifikant oder kleine Änderung
      log2FoldChange > 0 & padj < 0.05 & log2FoldChange >= 1 ~ "up",
      log2FoldChange < 0 & padj < 0.05 & log2FoldChange <= -1 ~ "down"
    )
  )

# 4. Top 3 positive und negative Veränderungen auswählen
top_up <- res_df %>%
  filter(status == "up") %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 3)

top_down <- res_df %>%
  filter(status == "down") %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 3)

top_genes_BCC_16wk <- bind_rows(top_up, top_down)

# 5. Volcano-Plot erstellen
plot_BCC_16wk<-ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = status), alpha = 0.6) +
  scale_color_manual(values = c("ns" = "grey", "down" = "red", "up" = "blue")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = top_genes_BCC_16wk,
    aes(label = rownames(top_genes_BCC_16wk)),
    nudge_x = ifelse(top_genes_BCC_16wk$log2FoldChange > 0, 0.8, -0.8),
    direction = "y",
    force = 2,
    box.padding = 0.6,
    point.padding = 0.4,
    max.overlaps = Inf,
    segment.color = "black",
    segment.size = 0.3,
    min.segment.length = 0,
    size = 3
  ) +
  theme_bw() +
  theme(
    legend.text = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(colour = "black", size = 10),
    axis.title.y = element_text( size = 10, colour = "black"),
    axis.text.x = element_text(vjust = 0.5, size = 10, colour = "black"),
    axis.text.y = element_text(vjust = 0.5, size = 10, colour = "black"),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank(),
    plot.title = element_text(size = 10, color = "#000", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 12, face = "bold", colour = "black")  # Legend title
  )+
  labs(
    x = "log2 fold change (C6-HSL vs mock)",
    y = "-log10(adjusted p-value)",
    title = "BCC436 16 weeks-old"
  )



plot_GP_2wk
plot_GP_16wk
plot_BCC_2wk
plot_BCC_16wk



ggarrange(plot_GP_2wk,plot_GP_16wk,plot_BCC_2wk,plot_BCC_16wk ,
          ncol = 2, nrow = 2,align = "hv", common.legend = T,
          legend = "bottom"  )



file_name <- paste0("volcano/volcano_deseq_FE.png")

ggsave(
  filename = file_name,
  plot = last_plot(),
  width = 22,       # width in cm
  height = 15,      # height in cm
  units = "cm",     # specify centimeters
  dpi = 300,
  path = out_dir# resolution (300 dpi = publication quality)
)





# MDS Control different soils ---------------------------------------

library(patchwork)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(forcats)





# Subset rhizosphere samples
physeq_raref_GP<- subset_samples(physeq_raref, Genotype == "Golden_Promise")
physeq_raref_GP <- subset_samples(physeq_raref_GP, Microhabitat == "RZ")
physeq_raref_GP <- subset_samples(physeq_raref_GP, Treatment == "mock")

# Remove taxa that are zero *within the subset*
physeq_raref_GP <- prune_taxa(
  taxa_sums(physeq_raref_GP) > 0,
  physeq_raref_GP
)


durations <- unique(sample_data(physeq_raref_GP)$Duration)


# Run MDS
ord <- ordinate(physeq_raref_GP, method = "MDS", distance = "bray")

# Plot
plot_ordination(physeq_raref_GP, ord, type = "Sample", shape = "Soil", color = "Duration") + 
  geom_point(size = 3) +
  stat_ellipse(aes(color = Duration, group = Duration)) + 
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 16)) +
  scale_color_manual(values = c("#D6C0AE", "#B48E6A", "#A76D60", "#601700")) +
  ggtitle(paste("Golden Promise - mock treatments"))


# Subset rhizosphere samples
physeq_raref_BCC<- subset_samples(physeq_raref, Genotype == "BCC436")
physeq_raref_BCC <- subset_samples(physeq_raref_BCC, Microhabitat == "RZ")
physeq_raref_BCC <- subset_samples(physeq_raref_BCC, Treatment == "mock")
# Remove taxa that are zero *within the subset*
physeq_raref_BCC <- prune_taxa(
  taxa_sums(physeq_raref_BCC) > 0,
  physeq_raref_BCC
)


durations <- unique(sample_data(physeq_raref_BCC)$Duration)


# Run MDS
ord <- ordinate(physeq_raref_BCC, method = "MDS", distance = "bray")

# Plot
plot_ordination(physeq_raref_BCC, ord, type = "Sample", shape = "Soil", color = "Duration") + 
  geom_point(size = 3) +
  stat_ellipse(aes(color = Duration, group = Duration)) + 
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 16)) +
  scale_color_manual(values = c("#84DCC6", "#D6EDFF", "#ACD7EC", "#8B95C9")) +
  ggtitle(paste("BCC436 - mock treatments"))

####

# Subset rhizosphere samples
physeq_raref_GP<- subset_samples(physeq_raref, Genotype == "Golden_Promise")
physeq_raref_GP <- subset_samples(physeq_raref_GP, Microhabitat == "RZ")
physeq_raref_GP <- subset_samples(physeq_raref_GP, Treatment == "C6-HSL")
# Remove taxa that are zero *within the subset*
physeq_raref_GP <- prune_taxa(
  taxa_sums(physeq_raref_GP) > 0,
  physeq_raref_GP
)


durations <- unique(sample_data(physeq_raref_GP)$Duration)


# Run MDS
ord <- ordinate(physeq_raref_GP, method = "MDS", distance = "bray")

# Plot
plot_ordination(physeq_raref_GP, ord, type = "Sample", shape = "Soil", color = "Duration") + 
  geom_point(size = 3) +
  stat_ellipse(aes(color = Duration, group = Duration)) + 
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 16)) +
  scale_color_manual(values = c("#D6C0AE", "#B48E6A", "#A76D60", "#601700")) +
  ggtitle(paste("Golden Promise - C6-HSL treatments"))


# Subset rhizosphere samples
physeq_raref_BCC<- subset_samples(physeq_raref, Genotype == "BCC436")
physeq_raref_BCC <- subset_samples(physeq_raref_BCC, Microhabitat == "RZ")
physeq_raref_BCC <- subset_samples(physeq_raref_BCC, Treatment == "C6-HSL")
# Remove taxa that are zero *within the subset*
physeq_raref_BCC <- prune_taxa(
  taxa_sums(physeq_raref_BCC) > 0,
  physeq_raref_BCC
)


durations <- unique(sample_data(physeq_raref_BCC)$Duration)


# Run MDS
ord <- ordinate(physeq_raref_BCC, method = "MDS", distance = "bray")

# Plot
plot_ordination(physeq_raref_BCC, ord, type = "Sample", shape = "Soil", color = "Duration") + 
  geom_point(size = 3) +
  stat_ellipse(aes(color = Duration, group = Duration)) + 
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 16)) +
  scale_color_manual(values = c("#84DCC6", "#D6EDFF", "#ACD7EC", "#8B95C9")) +
  ggtitle(paste("BCC436 - C6-HSL treatments"))


###
### Abschlussgrafik
###


# Subset rhizosphere samples
physeq_raref_GP <- subset_samples(physeq_raref, Microhabitat == "RZ")

# Remove taxa that are zero *within the subset*
physeq_raref_GP <- prune_taxa(
  taxa_sums(physeq_raref_GP) > 0,
  physeq_raref_GP
)


# Extract sample_data as a data.frame
samp_df <- data.frame(sample_data(physeq_raref_GP))

# Make sure the relevant columns are factors
samp_df$Treatment <- as.factor(samp_df$Treatment)
samp_df$Genotype  <- as.factor(samp_df$Genotype)
samp_df$Duration  <- as.factor(samp_df$Duration)
samp_df$Soil  <- as.factor(samp_df$Soil)

# Reorder factor levels
samp_df <- samp_df %>%
  mutate(
    Treatment = fct_relevel(Treatment, "mock", "C6-HSL"),
    Genotype = fct_relevel(Genotype, "Golden_Promise", "BCC436"),
    Soil = fct_relevel(Soil, "Substrate", "Field soil"),
    Duration = fct_relevel(Duration, "Two-week-old", "16-week-old", "Three-week-old","18-week-old")
  )



samp_df <- samp_df %>%
  mutate(
    Timepoint = case_when(
      Duration %in% c("Two-week-old", "Three-week-old")   ~ "Sampling Time Point 1",
      Duration %in% c("16-week-old", "18-week-old") ~ "Sampling Time Point 2",
      TRUE                                             ~ NA_character_
    )
  )

# Put it back into phyloseq object
sample_data(physeq_raref_GP) <- sample_data(samp_df)
durations <- unique(sample_data(physeq_raref_GP)$Duration)


# Run MDS
ord <- ordinate(physeq_raref_GP, method = "MDS", distance = "bray")

# Plot
plot_ordination(physeq_raref_GP, ord, type = "Sample", shape = "Genotype", color = "Treatment",) + 
  geom_point(size = 5.5, alpha = 0.5, stroke = 0) +
  #stat_ellipse(aes(color = Duration, group = Duration)) + 
  theme_bw() +
  theme(axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        strip.text = element_text(
          size = 13,
          colour = "black",
          face = "bold"
        )) +
  scale_color_manual(values = c( "#B48E6A",  "#601700")) +
  facet_grid(Soil ~ Timepoint)



file_name <- paste0("MDS/MDS_Abschluss.png")

ggsave(
  filename = file_name,
  plot = last_plot(),
  width = 22,       # width in cm
  height = 15,      # height in cm
  units = "cm",     # specify centimeters
  dpi = 300,
  path = out_dir# resolution (300 dpi = publication quality)
)



