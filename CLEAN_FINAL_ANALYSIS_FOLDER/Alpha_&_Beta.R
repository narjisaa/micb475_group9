#Final Script 

#### Load All Libraries #### 
library(phyloseq)
library(ape) 
library(indicspecies)
library(dplyr)
library(tidyverse)
library(vegan)
library(picante)
library(permute)
library(ggplot2)
library(patchwork)

#### Load All Data from Exported Qiime2 Files #### 
metafp <- "qiime2_export/parkinsons_metadata.txt"
meta <- read_delim(metafp, delim="\t")

otufp <- "qiime2_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "qiime2_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "qiime2_export/tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format Qiime2 Files ####
##### Format OTU table #####
#save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])

#Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`

#Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

##### Format Sample Metadata #####
#Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])

#Make sampleids the rownames
rownames(samp_df) <- meta[["#SampleID"]]

#Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)

##### Format Taxonomy #####
#Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% 
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix

#Save everything except feature IDs 
tax_mat <- tax_mat[,-1]

# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`

# Make taxa table
TAX <- tax_table(tax_mat)

#### Create Phyloseq Object ####
# Merge all into a phyloseq object
parkinsons_phyloseq <- phyloseq(OTU, SAMP, TAX, phylotree)

#save phyloseq object
save(parkinsons_phyloseq, file = "phyloseq_object.RData")

##### Subset for PD Samples #####
#select samples with PD 
parkinsons_pd <- subset_samples(
  parkinsons_phyloseq,
  Disease == "PD") 

##### Rarefaction #####
#set.seed for reproducability 
set.seed(123)

#rarefy the phyloseq object
parkinsons_pd_rare <- rarefy_even_depth(
  parkinsons_pd,
  sample.size = 8000,
  rngseed = 123,
  replace = FALSE, #done without replacement 
  verbose = FALSE)

##### MetaData Editing from Phyloseq Object #####
#Extract metadata
meta <- data.frame(sample_data(parkinsons_pd_rare))

#Keep only relevant columns
meta <- meta %>%
  select(Disease, updrs3total, moca_total, Sex)

#Keep only PD samples
meta <- meta %>%
  filter(Disease == "PD")

#Bin UPDRS scores 
meta <- meta %>%
  mutate(
    UPDRS_bin = case_when(
      updrs3total <= 33 ~ "Mild",
      updrs3total <= 59 ~ "Moderate",
      updrs3total > 59 ~ "Severe"))

#Bin MOCA scores
meta <- meta %>%
  mutate(
    MOCA_bin = case_when(
      moca_total < 26 ~ "Impaired",
      moca_total >= 26 ~ "Normal"))

#Combine Sex + UPDRS
meta <- meta %>%
  mutate(
    Sex_UPDRSScore = paste(Sex, UPDRS_bin, sep = "_"))

#Combine Sex + MOCA
meta <- meta %>%
  mutate(
    Sex_MOCAScore = paste(Sex, MOCA_bin, sep = "_"))

#Put back into phyloseq
sample_data(parkinsons_phyloseq) <- sample_data(meta)

#save the new phyloseq object
save(parkinsons_phyloseq, file = "phyloseq_object_updated.RData")

#### AIM 1: Alpha Diversity Anaylsis ####
##### UPDRS Faith's Phylogenetic Diversity (Figure 1) #####
#subset rarefied phyloseq object by sex 
updrs_faith_ps_female <- subset_samples(parkinsons_pd_rare, Sex == "Female")
updrs_faith_ps_male   <- subset_samples(parkinsons_pd_rare, Sex == "Male")

#remove missing metadata
updrs_faith_ps_female <- subset_samples(updrs_faith_ps_female, !is.na(updrstotal))
updrs_faith_ps_male   <- subset_samples(updrs_faith_ps_male, !is.na(updrstotal))

#prune taxa 
updrs_faith_ps_female <- prune_taxa(taxa_sums(updrs_faith_ps_female) > 0, updrs_faith_ps_female)
updrs_faith_ps_male   <- prune_taxa(taxa_sums(updrs_faith_ps_male) > 0, updrs_faith_ps_male)

#Extract and Convert otu Data into a Matrix 
updrs_otu_female <- as(otu_table(updrs_faith_ps_female), "matrix")
updrs_otu_male   <- as(otu_table(updrs_faith_ps_male), "matrix")

#Matrix Orientation: Taxa are Rows
if (taxa_are_rows(updrs_faith_ps_female)) {
  updrs_otu_female <- t(updrs_otu_female)}

if (taxa_are_rows(updrs_faith_ps_male)) {
  updrs_otu_male <- t(updrs_otu_male)}

#extract tree
updrs_tree_female <- phy_tree(updrs_faith_ps_female)
updrs_tree_male   <- phy_tree(updrs_faith_ps_male)

#calculate faith's PD 
updrs_female_faith_pd <- pd(updrs_otu_female,
                            updrs_tree_female,
                            include.root = FALSE)

updrs_male_faith_pd <- pd(updrs_otu_male,
                          updrs_tree_male,
                          include.root = FALSE)

#Add PD Values to Metadata
updrs_female_metadata <- data.frame(sample_data(updrs_faith_ps_female))
updrs_male_metadata   <- data.frame(sample_data(updrs_faith_ps_male))

updrs_female_metadata$Faith_PD <- updrs_female_faith_pd$PD
updrs_male_metadata$Faith_PD   <- updrs_male_faith_pd$PD


#Spearman Correlations for Faith's PD 
#Correlation for UPDRS (motor) 'female'
female_cor_updrs <- cor.test(~ Faith_PD + updrs3total, data = updrs_female_metadata, method = "spearman", exact = FALSE)
female_cor_updrs

#Correlation for UPDRS (motor) 'male'
male_cor_updrs <- cor.test(~ Faith_PD + updrs3total, data = updrs_male_metadata, method = "spearman", exact = FALSE)
male_cor_updrs

#Extract Spearman Correlation Values 
#Spearman Correlation for UPDRS (motor) 'female'
female_rho_updrs <- round(female_cor_updrs$estimate, 2)
female_p_updrs   <- signif(female_cor_updrs$p.value, 2)

#Spearman Correlation for UPDRS (motor) 'male'
male_rho_updrs <- round(male_cor_updrs$estimate, 2)
male_p_updrs   <- signif(male_cor_updrs$p.value, 2)

#Make Faith's PD Correlation Plots with RIBBON
#Ribbon is confidence interval
female_updrs_faith_plot_ribbon <- ggplot(updrs_female_metadata, aes(x = updrs3total, y = Faith_PD)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, colour = '#8B1A00', fill = "red", alpha = 0.4) +
  theme_classic() +
  labs(x = "Female UPDRS 3 Total",
       y = "Faith's PD") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman rho = ", female_rho_updrs,
                          "\np = ", female_p_updrs),
           hjust = 1.1, vjust = 1.5)
female_updrs_faith_plot_ribbon

#UPDRS Faith's PD Plot 'males' with RIBBON
male_updrs_faith_plot_ribbon <- ggplot(updrs_male_metadata, aes(x = updrs3total, y = Faith_PD)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, colour = 'blue4', fill = 'turquoise4', alpha = 0.4) +
  theme_classic() +
  labs(x = "Male UPDRS 3 Total",
       y = "Faith's PD") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman rho = ", male_rho_updrs,
                          "\np = ", male_p_updrs),
           hjust = 1.1, vjust = 1.5)
male_updrs_faith_plot_ribbon

#Combine All Plots to Save Only One File 
updrs_faith_combined_plot_ribboned <- female_updrs_faith_plot_ribbon + male_updrs_faith_plot_ribbon

#Save the File 
ggsave("UPDRS_Sex_Faith.png",
       updrs_faith_combined_plot_ribboned,
       width = 12,
       height = 4,
       dpi = 300)


##### MoCA Faith's Phylogenetic Diversity (Figure 1) #####
moca_faith_ps_female <- subset_samples(parkinsons_pd_rare, Sex == "Female")
moca_faith_ps_male   <- subset_samples(parkinsons_pd_rare, Sex == "Male")

#remove missing metadata
moca_faith_ps_female <- subset_samples(moca_faith_ps_female, !is.na(moca_total))
moca_faith_ps_male   <- subset_samples(moca_faith_ps_male, !is.na(moca_total))

#prune taxa 
moca_faith_ps_female <- prune_taxa(taxa_sums(moca_faith_ps_female) > 0, moca_faith_ps_female)
moca_faith_ps_male   <- prune_taxa(taxa_sums(moca_faith_ps_male) > 0, moca_faith_ps_male)

#Extract and Convert otu Data into a Matrix 
moca_otu_female <- as(otu_table(moca_faith_ps_female), "matrix")
moca_otu_male   <- as(otu_table(moca_faith_ps_male), "matrix")

#Matrix Orientation: Taxa are Rows
if (taxa_are_rows(moca_faith_ps_female)) {
  moca_otu_female <- t(moca_otu_female)
}

if (taxa_are_rows(moca_faith_ps_male)) {
  moca_otu_male <- t(moca_otu_male)
}

#extract tree
moca_tree_female <- phy_tree(moca_faith_ps_female)
moca_tree_male   <- phy_tree(moca_faith_ps_male)

#calculate faith's PD 
moca_female_faith_pd <- pd(moca_otu_female,
                           moca_tree_female,
                           include.root = FALSE)

moca_male_faith_pd <- pd(moca_otu_male,
                         moca_tree_male,
                         include.root = FALSE)

#Add PD Values to Metadata
moca_female_metadata <- data.frame(sample_data(moca_faith_ps_female))
moca_male_metadata   <- data.frame(sample_data(moca_faith_ps_male))

moca_female_metadata$Faith_PD <- moca_female_faith_pd$PD
moca_male_metadata$Faith_PD   <- moca_male_faith_pd$PD


#Spearman Correlations for Faith's PD 
#Correlation for MoCA 'female'
female_cor_moca <- cor.test(~ Faith_PD + moca_total,
                            data = moca_female_metadata,
                            method = "spearman",
                            exact = FALSE)
female_cor_moca

#Correlation for MoCA 'male'
male_cor_moca <- cor.test(~ Faith_PD + moca_total,
                          data = moca_male_metadata,
                          method = "spearman",
                          exact = FALSE)
male_cor_moca

#Extract Spearman Correlation Values 
#Spearman Correlation for MoCA 'female'
female_rho_moca <- round(female_cor_moca$estimate, 2)
female_p_moca   <- signif(female_cor_moca$p.value, 2)

#Spearman Correlation for MoCA 'male'
male_rho_moca <- round(male_cor_moca$estimate, 2)
male_p_moca   <- signif(male_cor_moca$p.value, 2)

#Make Faith's PD Correlation Plots with RIBBON
#MoCA Faith's PD Plot 'females' 
female_moca_faith_plot_ribbon <- ggplot(moca_female_metadata, aes(x = moca_total, y = Faith_PD)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, colour = '#8B1A00', fill = 'red', alpha=0.4) +
  theme_classic() +
  labs(x = "Female MoCA Total",
       y = "Faith's PD") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman rho = ", female_rho_moca,
                          "\np = ", female_p_moca),
           hjust = 1.1, vjust = 1.5)
female_moca_faith_plot_ribbon

#MoCA Faith's PD Plot 'males' 
male_moca_faith_plot <- ggplot(moca_male_metadata, aes(x = moca_total, y = Faith_PD)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(x = "Male MoCA Total",
       y = "Faith's PD") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman rho = ", male_rho_moca,
                          "\np = ", male_p_moca),
           hjust = 1.1, vjust = 1.5)
male_moca_faith_plot 

#MoCA Faith's PD Plot 'males' with ribbon
male_moca_faith_plot_ribbon <- ggplot(moca_male_metadata, aes(x = moca_total, y = Faith_PD)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, colour = 'blue4', fill = 'turquoise4', alpha=0.4) +
  theme_classic() +
  labs(x = "Male MoCA Total",
       y = "Faith's PD") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman rho = ", male_rho_moca,
                          "\np = ", male_p_moca),
           hjust = 1.1, vjust = 1.5)
male_moca_faith_plot_ribbon 

#Combine All Plots to Save Only One File 
moca_faith_combined_plot_ribbon <- female_moca_faith_plot_ribbon + male_moca_faith_plot_ribbon


#Save the File 
ggsave("MoCA_Sex_Faith_Ribbon.png",
       moca_faith_combined_plot_ribbon,
       width = 12,
       height = 4,
       dpi = 300)


##### UPDRS Shannon Diversity (Supplemental Analysis) #####
updrs_female_shannon <- subset_samples(parkinsons_pd_rare, Sex == "Female")
updrs_male_shannon <- subset_samples(parkinsons_pd_rare, Sex == "Male")

# Remove missing metadata
updrs_female_shannon <- subset_samples(updrs_female_shannon, !is.na(updrs3total))
updrs_male_shannon   <- subset_samples(updrs_male_shannon, !is.na(updrs3total))

# Prune taxa
updrs_female_shannon <- prune_taxa(taxa_sums(updrs_female_shannon) > 0, updrs_female_shannon)
updrs_male_shannon   <- prune_taxa(taxa_sums(updrs_male_shannon) > 0, updrs_male_shannon)

# Calculate Shannon diversity directly from phyloseq object
updrs_female_shannon_calc <- estimate_richness(updrs_female_shannon, measures = "Shannon")
updrs_male_shannon_calc   <- estimate_richness(updrs_male_shannon, measures = "Shannon")

# Add Shannon Values to Metadata
updrs_female_metadata <- data.frame(sample_data(updrs_female_shannon))
updrs_male_metadata   <- data.frame(sample_data(updrs_male_shannon))
updrs_female_metadata$Shannon <- updrs_female_shannon_calc$Shannon
updrs_male_metadata$Shannon   <- updrs_male_shannon_calc$Shannon

# Spearman Correlations for Shannon Diversity
# Correlation for UPDRS 'female'
female_cor_updrs <- cor.test(~ Shannon + updrs3total,
                             data = updrs_female_metadata,
                             method = "spearman",
                             exact = FALSE)
female_cor_updrs

# Correlation for UPDRS 'male'
male_cor_updrs <- cor.test(~ Shannon + updrs3total,
                           data = updrs_male_metadata,
                           method = "spearman",
                           exact = FALSE)
male_cor_updrs

# Extract Spearman Correlation Values
# Spearman Correlation for UPDRS 'female'
female_rho_updrs <- round(female_cor_updrs$estimate, 2)
female_p_updrs   <- signif(female_cor_updrs$p.value, 2)

# Spearman Correlation for UPDRS 'male'
male_rho_updrs <- round(male_cor_updrs$estimate, 2)
male_p_updrs   <- signif(male_cor_updrs$p.value, 2)

#Make Shannon Correlation Plots with RIBBON
#UPDRS Shannon Plot 'females' 
female_updrs_s_plot_ribbon <- ggplot(updrs_female_metadata, aes(x = updrs3total, y = Shannon)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, colour = '#8B1A00', fill = 'red', alpha=0.4) +
  theme_classic() +
  labs(x = "Female UPDRS 3 Total",
       y = "Shannon Diversity Index") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman rho = ", female_rho_updrs,
                          "\np = ", female_p_updrs),
           hjust = 1.1, vjust = 1.5)
female_updrs_s_plot_ribbon

#Make Shannon Correlation Plots with RIBBON 
#UPDRS Shannon Plot 'males' 
male_updrs_s_plot_ribbon <- ggplot(updrs_male_metadata, aes(x = updrs3total, y = Shannon)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, colour = 'blue4', fill = 'turquoise4', alpha=0.4) +
  theme_classic() +
  labs(x = "Male UPDRS 3 Total",
       y = "Shannon Diversity Index") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman rho = ", female_rho_updrs,
                          "\np = ", female_p_updrs),
           hjust = 1.1, vjust = 1.0)
male_updrs_s_plot_ribbon

#Combine All Plots to Save Only One File 
updrs_shannon_combined_plot_ribbon <- female_updrs_s_plot_ribbon + male_updrs_s_plot_ribbon

#Save the File 
ggsave("UPDRS_Sex_Shannon_Ribbon.png",
       updrs_shannon_combined_plot_ribbon,
       width = 12,
       height = 4,
       dpi = 300)


##### MOCA Shannon Diversity (Supplemental Analysis) #####
moca_female_shannon <- subset_samples(parkinsons_pd_rare, Sex == "Female")
moca_male_shannon <- subset_samples(parkinsons_pd_rare, Sex == "Male")

# Remove missing metadata
moca_female_shannon <- subset_samples(moca_female_shannon, !is.na(updrstotal))
moca_male_shannon   <- subset_samples(moca_male_shannon, !is.na(updrstotal))

# Prune taxa
moca_female_shannon <- prune_taxa(taxa_sums(moca_female_shannon) > 0, moca_female_shannon)
moca_male_shannon   <- prune_taxa(taxa_sums(moca_male_shannon) > 0, moca_male_shannon)

# Calculate Shannon diversity directly from phyloseq object
moca_female_shannon_calc <- estimate_richness(moca_female_shannon, measures = "Shannon")
moca_male_shannon_calc   <- estimate_richness(moca_male_shannon, measures = "Shannon")

# Add Shannon Values to Metadata
moca_female_metadata_shannon <- data.frame(sample_data(moca_female_shannon))
moca_male_metadata_shannon   <- data.frame(sample_data(moca_male_shannon))
moca_female_metadata_shannon$Shannon <- moca_female_shannon_calc$Shannon
moca_male_metadata_shannon$Shannon   <- moca_male_shannon_calc$Shannon

# Spearman Correlations for Shannon Diversity
# Correlation for MoCA 'female'
female_cor_moca_shannon <- cor.test(~ Shannon + moca_total,
                                    data = moca_female_metadata_shannon,
                                    method = "spearman",
                                    exact = FALSE)
female_cor_moca_shannon

# Correlation for MoCA 'male'
male_cor_moca_shannon <- cor.test(~ Shannon + moca_total,
                                  data = moca_male_metadata_shannon,
                                  method = "spearman",
                                  exact = FALSE)
male_cor_moca_shannon

# Extract Spearman Correlation Values
# Spearman Correlation for MoCa 'female'
female_rho_moca_s <- round(female_cor_moca_shannon$estimate, 2)
female_p_moca_s   <- signif(female_cor_moca_shannon$p.value, 2)

# Spearman Correlation for MoCA 'male'
male_rho_moca_s <- round(male_cor_moca_shannon$estimate, 2)
male_p_moca_s   <- signif(male_cor_moca_shannon$p.value, 2)

#Make Shannon Correlation Plots with RIBBON
#MoCA Shannon Plot 'females' 
female_moca_s_plot_ribbon <- ggplot(moca_female_metadata_shannon, aes(x = moca_total, y = Shannon)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, colour = '#8B1A00', fill = 'red', alpha=0.4) +
  theme_classic() +
  labs(x = "Female MoCA Total",
       y = "Shannon Diversity Index") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman rho = ", female_rho_moca_s,
                          "\np = ", female_p_moca_s),
           hjust = 1.1, vjust = 1.5)
female_moca_s_plot_ribbon

#Make Shannon Correlation Plots with RIBBON 
#MoCA Shannon Plot 'males' 
male_moca_s_plot_ribbon <- ggplot(moca_male_metadata_shannon, aes(x = moca_total, y = Shannon)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, colour = 'blue4', fill = 'turquoise4', alpha=0.4) +
  theme_classic() +
  labs(x = "Male MoCA Total",
       y = "Shannon Diversity Index") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman rho = ", male_rho_moca_s,
                          "\np = ", male_p_moca_s),
           hjust = 1.1, vjust = 1.0)
male_moca_s_plot_ribbon

#Combine All Plots to Save Only One File 
moca_shannon_combined_plot_ribbon <- female_moca_s_plot_ribbon + male_moca_s_plot_ribbon


#Save the File 
ggsave("MoCA_Sex_Shannon_Ribbon.png",
       moca_shannon_combined_plot_ribbon,
       width = 12,
       height = 4,
       dpi = 300)

#### AIM 2: Beta Diversity ####
##### UPDRS Beta Diversity : UnWeighted UniFrac (Figure 2) ######
#subset rarefied phyloseq object by sex 
ps_female <- subset_samples(parkinsons_pd_rare, Sex == "Female")
ps_male   <- subset_samples(parkinsons_pd_rare, Sex == "Male")

#remove na and taxa that dissapeared for 'females'
ps_female <- subset_samples(ps_female, !is.na(updrs3total))
ps_female <- prune_taxa(taxa_sums(ps_female) > 0, ps_female)
#remove na and taxa that dissapeared for 'males'
ps_male <- subset_samples(ps_male, !is.na(updrs3total))
ps_male <- prune_taxa(taxa_sums(ps_male) > 0, ps_male)

#calculate unweighted unifrac for 'females' 
updrs_unifrac_female <- distance(ps_female, method="unifrac", weighted=FALSE)
#calculate unweighted unifrac for 'males' 
updrs_unifrac_male <- distance(ps_male, method="unifrac", weighted=FALSE)

#run a PERMANOVA for 'females' 
updrs_permanova_female <- adonis2(updrs_unifrac_female ~ updrs3total,
                                  data = as(sample_data(ps_female), "data.frame"))
#run a PERMANOVA for 'males'
updrs_permanova_male <- adonis2(updrs_unifrac_male ~ updrs3total,
                                data = as(sample_data(ps_male), "data.frame"))

#extract values for 'females' plot 
updrs_female_r2 <- updrs_permanova_female["Model", "R2"]
updrs_female_p <- updrs_permanova_female["Model", "Pr(>F)"]
#extract values for 'males' plot 
updrs_male_r2 <- updrs_permanova_male["Model", "R2"]
updrs_male_p <- updrs_permanova_male["Model", "Pr(>F)"]

#format the values for the 'females' plot 
updrs_female_label_text <- paste0("UPDRS Female:\nR² = ",
                                  round(updrs_female_r2, 3),
                                  "\np = ",
                                  signif(updrs_female_p, 3))
#format the values for the 'males' plot 
updrs_male_label_text <- paste0("UPDRS Male:\nR² = ",
                                round(updrs_male_r2, 3),
                                "\np = ",
                                signif(updrs_male_p, 3))

#PCoA Plot for 'females' 
updrs_female_ordu <- ordinate(ps_female, method = "PCoA", distance = updrs_unifrac_female)

updrs_female_wuni_plot <- plot_ordination(ps_female,
                                          updrs_female_ordu,
                                          color = "updrs3total") +
  geom_point(size = 3) +
  annotate("text",
           x = Inf, y = Inf,
           label = updrs_female_label_text,
           hjust = 1,
           vjust = 1.1,   # moved up
           size = 5) +
  scale_color_viridis_c(option = "rocket", direction = -1) +
  theme_classic()
updrs_female_wuni_plot

#PCoA Plot for 'males' 
updrs_male_ordu <- ordinate(ps_male, method = "PCoA", distance = updrs_unifrac_male)

updrs_male_wuni_plot <- plot_ordination(ps_male,
                                        updrs_male_ordu,
                                        color = "updrs3total") +
  geom_point(size = 3) +
  annotate("text",
           x = Inf, y = Inf,
           label = updrs_male_label_text,
           hjust = 1,
           vjust = 1.1,
           size = 5) +
  scale_color_viridis_c(option = "mako", direction = -1) +
  theme_classic()
updrs_male_wuni_plot

#Combine All Plots to Save Only One File 
updrs_combined_plot <- updrs_female_wuni_plot +
  updrs_male_wuni_plot

#Save the File 
ggsave("UPDRS_Sex_UnWeighted_UniFrac.png",
       updrs_combined_plot,
       width = 12,
       height = 4,
       dpi = 300)


##### MOCA Beta Diversity : UnWeighted UniFrac (Figure 2) ######
#subset rarefied phyloseq object by sex 
moca_ps_female <- subset_samples(parkinsons_pd_rare, Sex == "Female")
moca_ps_male   <- subset_samples(parkinsons_pd_rare, Sex == "Male")

#remove na and taxa that dissapeared for 'females'
moca_ps_female <- subset_samples(moca_ps_female, !is.na(moca_total))
moca_ps_female <- prune_taxa(taxa_sums(moca_ps_female) > 0, moca_ps_female)
#remove na and taxa that dissapeared for 'males'
moca_ps_male <- subset_samples(moca_ps_male, !is.na(moca_total))
moca_ps_male <- prune_taxa(taxa_sums(moca_ps_male) > 0, moca_ps_male)

#calculate unweighted unifrac for 'females' 
moca_unifrac_female <- distance(moca_ps_female, method="unifrac", weighted=FALSE)
#calculate unweighted unifrac for 'males' 
moca_unifrac_male <- distance(moca_ps_male, method="unifrac", weighted=FALSE)

#run a PERMANOVA for 'females' 
moca_permanova_female <- adonis2(moca_unifrac_female ~ moca_total,
                                 data = as(sample_data(moca_ps_female), "data.frame"))
#run a PERMANOVA for 'males'
moca_permanova_male <- adonis2(moca_unifrac_male ~ moca_total,
                               data = as(sample_data(moca_ps_male), "data.frame"))

#extract values for 'females' plot 
moca_female_r2 <- moca_permanova_female["Model", "R2"]
moca_female_p <- moca_permanova_female["Model", "Pr(>F)"]
#extract values for 'males' plot 
moca_male_r2 <- moca_permanova_male["Model", "R2"]
moca_male_p <- moca_permanova_male["Model", "Pr(>F)"]

#format the values for the 'females' plot 
moca_female_label_text <- paste0("MOCA Female:\nR² = ",
                                 round(moca_female_r2, 3),
                                 "\np = ",
                                 signif(moca_female_p, 3))
#format the values for the 'males' plot 
moca_male_label_text <- paste0("MOCA Male:\nR² = ",
                               round(moca_male_r2, 3),
                               "\np = ",
                               signif(moca_male_p, 3))

#PCoA Plot for 'females' 
moca_female_ordu <- ordinate(moca_ps_female, method = "PCoA", distance = moca_unifrac_female)

moca_female_wuni_plot <- plot_ordination(moca_ps_female,
                                         moca_female_ordu,
                                         color = "moca_total") +
  geom_point(size = 3) +
  annotate("text",
           x = Inf, y = Inf,
           label = moca_female_label_text,
           hjust = 1.1,
           vjust = 1.5,
           size = 5) +
  scale_colour_viridis_c(option = "rocket", direction= -1)+
  theme_classic()
moca_female_wuni_plot

#PCoA Plot for 'males' 
moca_male_ordu <- ordinate(moca_ps_male, method = "PCoA", distance = moca_unifrac_male)

moca_male_wuni_plot <- plot_ordination(moca_ps_male,
                                       moca_male_ordu,
                                       color = "moca_total") +
  geom_point(size = 3) +
  annotate("text",
           x = Inf, y = Inf,
           label = moca_male_label_text,
           hjust = 1.1,
           vjust = 0.99,
           size = 5) +
  scale_colour_viridis_c(option = "mako", direction= -1)+
  theme_classic()
moca_male_wuni_plot

#Combine All Plots to Save Only One File 
moca_combined_plot <- moca_female_wuni_plot +
  moca_male_wuni_plot

#Save the File 
ggsave("MOCA_Sex_UnWeighted_UniFrac.png",
       moca_combined_plot,
       width = 12,
       height = 4,
       dpi = 300)



