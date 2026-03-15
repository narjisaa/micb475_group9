# MICB 475 Team 9 Script 

#### Load All Libraries #### 
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)
library(picante)
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

#### Format OTU table ####
#save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])

#Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`

#Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

#### Format Sample Metadata ####
#Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])

#Make sampleids the rownames
rownames(samp_df) <- meta[["#SampleID"]]

#Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)

#### Format Taxonomy ####
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

#### Subset for PD Samples ####
#select samples with PD 
parkinsons_pd <- subset_samples(
  parkinsons_phyloseq,
  Disease == "PD") 

#### Rarefaction ####
#set.seed for reproducability 
set.seed(123)

#rarefy the phyloseq object
parkinsons_pd_rare <- rarefy_even_depth(
  parkinsons_pd,
  sample.size = 8000,
  rngseed = 123,
  replace = FALSE, #done without replacement 
  verbose = FALSE)

#### Faith's Phylogenetic Diversity (Faith's PD) ####
#Extract and Convert otu Data into a Matrix 
otu <- as(otu_table(parkinsons_pd_rare), "matrix")

#Matrix Orientation: Taxa are Rows 
if (taxa_are_rows(parkinsons_pd_rare)) {
  otu <- t(otu)}

#Extract Tree
tree <- phy_tree(parkinsons_pd_rare)

#Calculate Faith's PD
faith_pd <- pd(otu, tree, include.root = FALSE)

#Add PD Values to Metadata
sample_data(parkinsons_pd_rare)$Faith_PD <- faith_pd$PD

#Spearman Correlations for Faith's PD 
#Make Rarefied Phyloseq Object a Data Frame 
metadata <- data.frame(sample_data(parkinsons_pd_rare))

#Correlation for UPDRS (motor) 
cor_updrs <- cor.test(~ Faith_PD + updrstotal, data = metadata, method = "spearman", exact = FALSE)
cor_updrs

#Correlation for MOCA (cognitive)
cor_moca <- cor.test(~ Faith_PD + moca_total, data = metadata, method = "spearman", exact = FALSE)
cor_moca

#Correlation for Duration (Time)
cor_dur <- cor.test(~ Faith_PD + Disease_dur, data = metadata, method = "spearman", exact = FALSE)
cor_dur

#Extract Spearman Correlation Values 
#Spearman Correlation for UPDRS
rho_updrs <- round(cor_updrs$estimate, 2)
p_updrs   <- signif(cor_updrs$p.value, 2)

#Spearman Correlation for MOCA
rho_moca <- round(cor_moca$estimate, 2)
p_moca   <- signif(cor_moca$p.value, 2)

#Spearman Correlation for Duration
rho_dur <- round(cor_dur$estimate, 2)
p_dur   <- signif(cor_dur$p.value, 2)

#Make Faith's PD Correlation Plots
#UPDRS Faith's PD Plot 
updrs_faith_plot <- ggplot(metadata, aes(x = updrstotal, y = Faith_PD)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(x = "UPDRS Total",
       y = "Faith's PD") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman rho = ", rho_updrs,
                          "\np = ", p_updrs),
           hjust = 1.1, vjust = 1.5)

#MOCA Faith's PD Plot
moca_faith_plot <- ggplot(metadata, aes(x = moca_total, y = Faith_PD)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(x = "MoCA Score",
       y = "Faith's PD") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman rho = ", rho_moca,
                          "\np = ", p_moca),
           hjust = 1.1, vjust = 1.5)

#Disease Duration Faith's PD Plot
dur_faith_plot <- ggplot(metadata, aes(x = Disease_dur, y = Faith_PD)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(x = "Disease Duration",
       y = "Faith's PD") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman rho = ", rho_dur,
                          "\np = ", p_dur),
           hjust = 1.1, vjust = 1.5)

#Combine All Plots to Save Only One File 
combined_plot <- updrs_faith_plot +
  moca_faith_plot +
  dur_faith_plot

#Save the File 
ggsave("Faith_PD_correlations.png",
       combined_plot,
       width = 12,
       height = 4,
       dpi = 300)


#### UPDRS Beta Diversity : UnWeighted UniFrac #####
#calculate distance 
updrs_unifrac_w <- phyloseq::distance(parkinsons_pd_rare, 
                                      method = "wunifrac")

#make 'Sex' column a factor for PERMANOVA analysis 
metadata$Sex <- as.factor(metadata$Sex)

#check if data contains NA's
colSums(is.na(metadata[, c("updrstotal", "Age", "Sex", "Age_onset")]))

#since there are NA's present the data should be subset 
updrs_metadata_complete <- metadata[
  complete.cases(metadata[, 
                          c("updrstotal", "Age", "Sex", "Age_onset")]), ]

#subset the distance matrix too
updrs_samples_keep <- rownames(updrs_metadata_complete)

updrs_unifrac_subset <- as.dist(
  as.matrix(updrs_unifrac_w)[updrs_samples_keep, updrs_samples_keep])

#Now the model is using matched data and not dropping samples 

# PERMANOVA with Continuous UPDRS -- NO BINNING 
#controls for age, sex and age onset
updrs_permanova <- adonis2(updrs_unifrac_subset ~ 
                             updrstotal + Age + Sex + Age_onset,
                           data = updrs_metadata_complete,
                           by = "margin")

#extract values for plot 
updrs_r2 <- updrs_permanova["updrstotal", "R2"]
updrs_p <- updrs_permanova["updrstotal", "Pr(>F)"]

#format the values for the plot 
updrs_label_text <- paste0("UPDRS:\nR² = ",
                           round(updrs_r2, 3),
                           "\np = ",
                           signif(updrs_p, 3))

#PCoA Plot with Continuous UPDRS -- NO BINNING
updrs_ps_subset <- prune_samples(updrs_samples_keep, parkinsons_pd_rare)

updrs_ordu_unw <- ordinate(updrs_ps_subset,
                           method = "PCoA",
                           distance = updrs_unifrac_subset)

updrs_wuni_plot <- plot_ordination(updrs_ps_subset,
                                   updrs_ordu_unw,
                                   color = "updrstotal") +
  geom_point(size = 3) +
  annotate("text",
           x = Inf, y = Inf,
           label = updrs_label_text,
           hjust = 1.1,
           vjust = 1.5,
           size = 5) +
  theme_minimal()
updrs_wuni_plot

#Save the File 
ggsave("UPDRS_Weighted_UniFrac.png",
       updrs_wuni_plot,
       width = 12,
       height = 4,
       dpi = 300)

#### MOCA Beta Diversity : UnWeighted UniFrac #####
#calculate distance 
moca_unifrac_w <- phyloseq::distance(parkinsons_pd_rare, 
                                     method = "wunifrac")

#make 'Sex' column a factor for PERMANOVA analysis 
metadata$Sex <- as.factor(metadata$Sex)

#check if data contains NA's
colSums(is.na(metadata[, c("moca_total", "Age", "Sex", "Age_onset")]))

#since there are NA's present the data should be subset 
moca_metadata_complete <- metadata[
  complete.cases(metadata[, 
                          c("moca_total", "Age", "Sex", "Age_onset")]), ]

#subset the distance matrix too
moca_samples_keep <- rownames(moca_metadata_complete)

moca_unifrac_subset <- as.dist(
  as.matrix(moca_unifrac_w)[moca_samples_keep, moca_samples_keep])

#Now the model is using matched data and not dropping samples 

#PERMANOVA with Continuous MOCA -- NO BINNING 
#controls for age, sex and age onset
moca_permanova <- adonis2(moca_unifrac_subset ~ 
                            moca_total + Age + Sex + Age_onset,
                          data = moca_metadata_complete,
                          by = "margin")

#extract values for plot 
moca_r2 <- moca_permanova["moca_total", "R2"]
moca_p <- moca_permanova["moca_total", "Pr(>F)"]

#format the values for the plot 
moca_label_text <- paste0("UPDRS:\nR² = ",
                          round(moca_r2, 3),
                          "\np = ",
                          signif(moca_p, 3))

#PCoA Plot with Continuous MOCA -- NO BINNING 
moca_ps_subset <- prune_samples(moca_samples_keep, parkinsons_pd_rare)

moca_ordu_unw <- ordinate(moca_ps_subset,
                          method = "PCoA",
                          distance = moca_unifrac_subset)

moca_wuni_plot <- plot_ordination(moca_ps_subset,
                                  moca_ordu_unw,
                                  color = "moca_total") +
  geom_point(size = 3) +
  annotate("text",
           x = Inf, y = Inf,
           label = moca_label_text,
           hjust = 1.1,
           vjust = 1.5,
           size = 5) +
  theme_minimal()
moca_wuni_plot

#Save the File 
ggsave("Moca_Weighted_UniFrac.png",
       moca_wuni_plot,
       width = 12,
       height = 4,
       dpi = 300)

#the results above show that 'Sex' is the main confounding variable

#### (Sex Controlled) UPDRS Beta Diversity : UnWeighted UniFrac #####
#subset rarefied phyloseq object by sex 
ps_female <- subset_samples(parkinsons_pd_rare, Sex == "Female")
ps_male   <- subset_samples(parkinsons_pd_rare, Sex == "Male")

#remove na and taxa that dissapeared for 'females'
ps_female <- subset_samples(ps_female, !is.na(updrstotal))
ps_female <- prune_taxa(taxa_sums(ps_female) > 0, ps_female)
#remove na and taxa that dissapeared for 'males'
ps_male <- subset_samples(ps_male, !is.na(updrstotal))
ps_male <- prune_taxa(taxa_sums(ps_male) > 0, ps_male)

#calculate unweighted unifrac for 'females' 
updrs_unifrac_female <- distance(ps_female, method="unifrac", weighted=FALSE)
#calculate unweighted unifrac for 'males' 
updrs_unifrac_male <- distance(ps_male, method="unifrac", weighted=FALSE)

#run a PERMANOVA for 'females' 
updrs_permanova_female <- adonis2(updrs_unifrac_female ~ updrstotal,
                                  data = as(sample_data(ps_female), "data.frame"))
#run a PERMANOVA for 'males'
updrs_permanova_male <- adonis2(updrs_unifrac_male ~ updrstotal,
                                data = as(sample_data(ps_male), "data.frame"))

#extract values for 'females' plot 
updrs_female_r2 <- updrs_permanova_female["Model", "R2"]
updrs_female_p <- updrs_permanova_female["Model", "Pr(>F)"]
#extract values for 'males' plot 
updrs_male_r2 <- updrs_permanova_male["Model", "R2"]
updrs_male_p <- updrs_permanova_male["Model", "Pr(>F)"]

#format the values for the 'females' plot 
updrs_female_label_text <- paste0("UPDRS:\nR² = ",
                           round(updrs_female_r2, 3),
                           "\np = ",
                           signif(updrs_female_p, 3))
#format the values for the 'males' plot 
updrs_male_label_text <- paste0("UPDRS:\nR² = ",
                           round(updrs_male_r2, 3),
                           "\np = ",
                           signif(updrs_male_p, 3))

#PCoA Plot for 'females' 
updrs_female_ordu <- ordinate(ps_female, method = "PCoA", distance = updrs_unifrac_female)

updrs_female_wuni_plot <- plot_ordination(ps_female,
                                   updrs_female_ordu,
                                   color = "updrstotal") +
  geom_point(size = 3) +
  annotate("text",
           x = Inf, y = Inf,
           label = updrs_female_label_text,
           hjust = 1.1,
           vjust = 1.5,
           size = 5) +
  theme_minimal()
updrs_female_wuni_plot

#PCoA Plot for 'males' 
updrs_male_ordu <- ordinate(ps_male, method = "PCoA", distance = updrs_unifrac_male)

updrs_male_wuni_plot <- plot_ordination(ps_male,
                                           updrs_male_ordu,
                                           color = "updrstotal") +
  geom_point(size = 3) +
  annotate("text",
           x = Inf, y = Inf,
           label = updrs_male_label_text,
           hjust = 1.1,
           vjust = 1.5,
           size = 5) +
  theme_minimal()
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


#### (Sex Controlled) MOCA Beta Diversity : UnWeighted UniFrac #####
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
moca_female_label_text <- paste0("MOCA:\nR² = ",
                                 round(moca_female_r2, 3),
                                 "\np = ",
                                 signif(moca_female_p, 3))
#format the values for the 'males' plot 
moca_male_label_text <- paste0("MOCA:\nR² = ",
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
  theme_minimal()
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
           vjust = 1.5,
           size = 5) +
  theme_minimal()
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


#### (Sex Controlled) DURATION Beta Diversity : UnWeighted UniFrac #####
#subset rarefied phyloseq object by sex 
dur_ps_female <- subset_samples(parkinsons_pd_rare, Sex == "Female")
dur_ps_male   <- subset_samples(parkinsons_pd_rare, Sex == "Male")

#remove na and taxa that dissapeared for 'females'
dur_ps_female <- subset_samples(dur_ps_female, !is.na(Disease_dur))
dur_ps_female <- prune_taxa(taxa_sums(dur_ps_female) > 0, dur_ps_female)
#remove na and taxa that dissapeared for 'males'
dur_ps_male <- subset_samples(dur_ps_male, !is.na(Disease_dur))
dur_ps_male <- prune_taxa(taxa_sums(dur_ps_male) > 0, dur_ps_male)

#calculate unweighted unifrac for 'females' 
dur_unifrac_female <- distance(dur_ps_female, method="unifrac", weighted=FALSE)
#calculate unweighted unifrac for 'males' 
dur_unifrac_male <- distance(dur_ps_male, method="unifrac", weighted=FALSE)

#run a PERMANOVA for 'females' 
dur_permanova_female <- adonis2(dur_unifrac_female ~ Disease_dur,
                                data = as(sample_data(dur_ps_female), "data.frame"))
#run a PERMANOVA for 'males'
dur_permanova_male <- adonis2(dur_unifrac_male ~ Disease_dur,
                              data = as(sample_data(dur_ps_male), "data.frame"))

#extract values for 'females' plot 
dur_female_r2 <- dur_permanova_female["Model", "R2"]
dur_female_p <- dur_permanova_female["Model", "Pr(>F)"]
#extract values for 'males' plot 
dur_male_r2 <- dur_permanova_male["Model", "R2"]
dur_male_p <- dur_permanova_male["Model", "Pr(>F)"]

#format the values for the 'females' plot 
dur_female_label_text <- paste0("DURATION:\nR² = ",
                                round(dur_female_r2, 3),
                                "\np = ",
                                signif(dur_female_p, 3))
#format the values for the 'males' plot 
dur_male_label_text <- paste0("DURATION:\nR² = ",
                              round(dur_male_r2, 3),
                              "\np = ",
                              signif(dur_male_p, 3))

#PCoA Plot for 'females' 
dur_female_ordu <- ordinate(dur_ps_female, method = "PCoA", distance = dur_unifrac_female)

dur_female_wuni_plot <- plot_ordination(dur_ps_female,
                                        dur_female_ordu,
                                        color = "Disease_dur") +
  geom_point(size = 3) +
  annotate("text",
           x = Inf, y = Inf,
           label = dur_female_label_text,
           hjust = 1.1,
           vjust = 1.5,
           size = 5) +
  theme_minimal()
dur_female_wuni_plot

#PCoA Plot for 'males' 
dur_male_ordu <- ordinate(dur_ps_male, method = "PCoA", distance = dur_unifrac_male)

dur_male_wuni_plot <- plot_ordination(dur_ps_male,
                                      dur_male_ordu,
                                      color = "Disease_dur") +
  geom_point(size = 3) +
  annotate("text",
           x = Inf, y = Inf,
           label = dur_male_label_text,
           hjust = 1.1,
           vjust = 1.5,
           size = 5) +
  theme_minimal()
dur_male_wuni_plot

#Combine All Plots to Save Only One File 
dur_combined_plot <- dur_female_wuni_plot +
  dur_male_wuni_plot

#Save the File 
ggsave("DURATION_Sex_UnWeighted_UniFrac.png",
       dur_combined_plot,
       width = 12,
       height = 4,
       dpi = 300)



#### (Sex Controlled) UPDRS Faith's Phylogenetic Diversity (Faith's PD) ####

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
female_cor_updrs <- cor.test(~ Faith_PD + updrstotal, data = updrs_female_metadata, method = "spearman", exact = FALSE)
female_cor_updrs

#Correlation for UPDRS (motor) 'male'
male_cor_updrs <- cor.test(~ Faith_PD + updrstotal, data = updrs_male_metadata, method = "spearman", exact = FALSE)
male_cor_updrs

#Extract Spearman Correlation Values 
#Spearman Correlation for UPDRS (motor) 'female'
female_rho_updrs <- round(female_cor_updrs$estimate, 2)
female_p_updrs   <- signif(female_cor_updrs$p.value, 2)

#Spearman Correlation for UPDRS (motor) 'male'
male_rho_updrs <- round(male_cor_updrs$estimate, 2)
male_p_updrs   <- signif(male_cor_updrs$p.value, 2)

#Make Faith's PD Correlation Plots
#UPDRS Faith's PD Plot 'females' 
female_updrs_faith_plot <- ggplot(updrs_female_metadata, aes(x = updrstotal, y = Faith_PD)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(x = "Female UPDRS Total",
       y = "Faith's PD") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman rho = ", female_rho_updrs,
                          "\np = ", female_p_updrs),
           hjust = 1.1, vjust = 1.5)
female_updrs_faith_plot

#UPDRS Faith's PD Plot 'males' 
male_updrs_faith_plot <- ggplot(updrs_male_metadata, aes(x = updrstotal, y = Faith_PD)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(x = "Male UPDRS Total",
       y = "Faith's PD") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman rho = ", male_rho_updrs,
                          "\np = ", male_p_updrs),
           hjust = 1.1, vjust = 1.5)
male_updrs_faith_plot 

#Combine All Plots to Save Only One File 
updrs_faith_combined_plot <- female_updrs_faith_plot + male_updrs_faith_plot

#Save the File 
ggsave("UPDRS_Sex_Faith.png",
       updrs_faith_combined_plot,
       width = 12,
       height = 4,
       dpi = 300)


#### (Sex Controlled) MoCA Faith's Phylogenetic Diversity (Faith's PD) ####

#subset rarefied phyloseq object by sex 
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

#Make Faith's PD Correlation Plots
#MoCA Faith's PD Plot 'females' 
female_moca_faith_plot <- ggplot(moca_female_metadata, aes(x = moca_total, y = Faith_PD)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(x = "Female MoCA Total",
       y = "Faith's PD") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman rho = ", female_rho_moca,
                          "\np = ", female_p_moca),
           hjust = 1.1, vjust = 1.5)
female_moca_faith_plot

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

#Combine All Plots to Save Only One File 
moca_faith_combined_plot <- female_moca_faith_plot + male_moca_faith_plot

#Save the File 
ggsave("MoCA_Sex_Faith.png",
       moca_faith_combined_plot,
       width = 12,
       height = 4,
       dpi = 300)


#### (Sex Controlled) Disease Duration Faith's Phylogenetic Diversity (Faith's PD) ####

#subset rarefied phyloseq object by sex 
dur_faith_ps_female <- subset_samples(parkinsons_pd_rare, Sex == "Female")
dur_faith_ps_male   <- subset_samples(parkinsons_pd_rare, Sex == "Male")

#remove missing metadata
dur_faith_ps_female <- subset_samples(dur_faith_ps_female, !is.na(Disease_dur))
dur_faith_ps_male   <- subset_samples(dur_faith_ps_male, !is.na(Disease_dur))

#prune taxa 
dur_faith_ps_female <- prune_taxa(taxa_sums(dur_faith_ps_female) > 0, dur_faith_ps_female)
dur_faith_ps_male   <- prune_taxa(taxa_sums(dur_faith_ps_male) > 0, dur_faith_ps_male)

#Extract and Convert otu Data into a Matrix 
dur_otu_female <- as(otu_table(dur_faith_ps_female), "matrix")
dur_otu_male   <- as(otu_table(dur_faith_ps_male), "matrix")

#Matrix Orientation: Taxa are Rows
if (taxa_are_rows(dur_faith_ps_female)) {
  dur_otu_female <- t(dur_otu_female)
}

if (taxa_are_rows(dur_faith_ps_male)) {
  dur_otu_male <- t(dur_otu_male)
}

#extract tree
dur_tree_female <- phy_tree(dur_faith_ps_female)
dur_tree_male   <- phy_tree(dur_faith_ps_male)

#calculate faith's PD 
dur_female_faith_pd <- pd(dur_otu_female,
                          dur_tree_female,
                          include.root = FALSE)

dur_male_faith_pd <- pd(dur_otu_male,
                        dur_tree_male,
                        include.root = FALSE)

#Add PD Values to Metadata
dur_female_metadata <- data.frame(sample_data(dur_faith_ps_female))
dur_male_metadata   <- data.frame(sample_data(dur_faith_ps_male))

dur_female_metadata$Faith_PD <- dur_female_faith_pd$PD
dur_male_metadata$Faith_PD   <- dur_male_faith_pd$PD


#Spearman Correlations for Faith's PD 
#Correlation for Disease Duration 'female'
female_cor_dur <- cor.test(~ Faith_PD + Disease_dur,
                           data = dur_female_metadata,
                           method = "spearman",
                           exact = FALSE)
female_cor_dur

#Correlation for Disease Duration 'male'
male_cor_dur <- cor.test(~ Faith_PD + Disease_dur,
                         data = dur_male_metadata,
                         method = "spearman",
                         exact = FALSE)
male_cor_dur

#Extract Spearman Correlation Values 
#Spearman Correlation for Disease Duration 'female'
female_rho_dur <- round(female_cor_dur$estimate, 2)
female_p_dur   <- signif(female_cor_dur$p.value, 2)

#Spearman Correlation for Disease Duration 'male'
male_rho_dur <- round(male_cor_dur$estimate, 2)
male_p_dur   <- signif(male_cor_dur$p.value, 2)

#Make Faith's PD Correlation Plots
#Disease Duration Faith's PD Plot 'females' 
female_dur_faith_plot <- ggplot(dur_female_metadata, aes(x = Disease_dur, y = Faith_PD)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(x = "Female Disease Duration",
       y = "Faith's PD") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman rho = ", female_rho_dur,
                          "\np = ", female_p_dur),
           hjust = 1.1, vjust = 1.5)
female_dur_faith_plot

#Disease Duration Faith's PD Plot 'males' 
male_dur_faith_plot <- ggplot(dur_male_metadata, aes(x = Disease_dur, y = Faith_PD)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(x = "Male Disease Duration",
       y = "Faith's PD") +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman rho = ", male_rho_dur,
                          "\np = ", male_p_dur),
           hjust = 1.1, vjust = 1.5)
male_dur_faith_plot 

#Combine All Plots to Save Only One File 
dur_faith_combined_plot <- female_dur_faith_plot + male_dur_faith_plot

#Save the File 
ggsave("DiseaseDuration_Sex_Faith.png",
       dur_faith_combined_plot,
       width = 12,
       height = 4,
       dpi = 300)












