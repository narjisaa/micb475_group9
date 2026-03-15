#Alpha Diversity Analysis 

#load library 
library(phyloseq)

#load phyloseq object
load("parkinsons_phyloseq.RData")

#select samples with PD 
parkinsons_pd <- subset_samples(
  parkinsons_phyloseq,
  Disease == "PD") 

#set.seed for reproducability 
set.seed(123)

#rarefy the phyloseq object
parkinsons_pd_rare <- rarefy_even_depth(
  parkinsons_pd,
  sample.size = 8000,
  rngseed = 123,
  replace = FALSE, #done without replacement 
  verbose = FALSE)

#perform faiths anaylsis

library(picante)

otu <- as(otu_table(parkinsons_pd_rare), "matrix")

if (!taxa_are_rows(parkinsons_pd_rare)) {
  otu <- t(otu)}

tree <- phy_tree(parkinsons_pd_rare)

faith_pd <- pd(otu, tree, include.root = FALSE)


#this code got an error so trying to fix 
sample_data(parkinsons_pd_rare)$Faith_PD <- faith_pd$PD

if (taxa_are_rows(parkinsons_pd_rare)) {
  otu <- t(otu)
}

# Extract tree
tree <- phy_tree(parkinsons_pd_rare)

# Calculate Faith's PD
faith_pd <- pd(otu, tree, include.root = FALSE)

# Add PD values to metadata
sample_data(parkinsons_pd_rare)$Faith_PD <- faith_pd$PD

length(faith_pd$PD)
nsamples(parkinsons_pd_rare)

#fix ended here 

#perfrom spearman correlations 

metadata <- data.frame(sample_data(parkinsons_pd_rare))

cor_updrs <- cor.test(metadata$Faith_PD, metadata$updrstotal, method = "spearman")
cor_moca  <- cor.test(metadata$Faith_PD, metadata$moca_total,  method = "spearman")
cor_dur   <- cor.test(metadata$Faith_PD, metadata$Disease_dur, method = "spearman")

cor_updrs
cor_moca
cor_dur

#make plots 

library(ggplot2)

p1 <- ggplot(metadata, aes(x = updrstotal, y = Faith_PD)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(x = "UPDRS total", y = "Faith's PD")

p2 <- ggplot(metadata, aes(x = moca_total, y = Faith_PD)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(x = "MoCA score", y = "Faith's PD")

p3 <- ggplot(metadata, aes(x = Disease_dur, y = Faith_PD)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(x = "Disease duration", y = "Faith's PD")

p1; p2; p3

#save the plots 
library(ggplot2)

# Save p1
ggsave("UPDRS_vs_Faith_PD.png", plot = p1, width = 6, height = 4, dpi = 300)

# Save p2
ggsave("MoCA_vs_Faith_PD.png", plot = p2, width = 6, height = 4, dpi = 300)

# Save p3
ggsave("DiseaseDur_vs_Faith_PD.png", plot = p3, width = 6, height = 4, dpi = 300)