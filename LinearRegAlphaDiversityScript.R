#Alpha Diversity Analysis Linear Regression

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


# perform linear regressions instead of spearman

metadata <- data.frame(sample_data(parkinsons_pd_rare))

# make sure sex is treated as categorical
metadata$sex <- as.factor(metadata$Sex)

# Linear regression models

model_updrs <- lm(Faith_PD ~ updrstotal + BMI + Sex + Age_onset,
                  data = metadata)

model_moca  <- lm(Faith_PD ~ moca_total + BMI + Sex + Age_onset,
                  data = metadata)

model_dur   <- lm(Faith_PD ~ Disease_dur + BMI + Sex + Age_onset,
                  data = metadata)

# View results
summary(model_updrs)
summary(model_moca)
summary(model_dur)

# Make plots (unadjusted visual trends)

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
ggsave("Linear-Reg_UPDRS_vs_Faith_PD.png", plot = p1, width = 6, height = 4, dpi = 300)

# Save p2
ggsave("Linear-Reg_MoCA_vs_Faith_PD.png", plot = p2, width = 6, height = 4, dpi = 300)

# Save p3
ggsave("Linear-Reg_DiseaseDur_vs_Faith_PD.png", plot = p3, width = 6, height = 4, dpi = 300)