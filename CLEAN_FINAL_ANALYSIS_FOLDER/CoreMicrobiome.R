# Final Script for Core Microbiome Analysis 

#### Load Libraries ####
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(indicspecies)
library(DESeq2)
library(ggplot2)
library(ggtext)

#### Load Data ####
load("phyloseq_object_updated.RData")
sample_data(parkinsons_phyloseq)

#### Inspect group labels ####
cat("=== Sex_UPDRSScore group labels ===\n")
print(table(sample_data(parkinsons_phyloseq)$Sex_UPDRSScore, useNA = "ifany"))

cat("\n=== Sex_MOCAScore group labels ===\n")
print(table(sample_data(parkinsons_phyloseq)$Sex_MOCAScore, useNA = "ifany"))

cat("\n=== All metadata columns ===\n")
print(colnames(sample_data(parkinsons_phyloseq)))

# Exclude samples with _NA bins (missing UPDRS or MoCA data)
ps_updrs <- subset_samples(parkinsons_phyloseq, Sex_UPDRSScore %in% 
                             c("Male_Mild", "Male_Moderate", "Female_Mild", "Female_Moderate", "Male_Severe"))
ps_moca  <- subset_samples(parkinsons_phyloseq, Sex_MOCAScore %in% 
                             c("Male_Normal", "Male_Impaired", "Female_Normal", "Female_Impaired"))

cat("\nSamples retained:\n")
cat("  Sex_UPDRSScore:", nsamples(ps_updrs), "\n")
cat("  Sex_MOCAScore: ", nsamples(ps_moca),  "\n")

#CORE MICROBIOME 
# Groups: Male_Mild, Male_Moderate, Male_Severe, Female_Mild, Female_Moderate
# Key comparisons:
# Within males: Mild vs Moderate motor severity
# Within females: Mild vs Moderate motor severity
# Across sexes: Male vs Female within same UPDRS severity bin

#### Core Microbiome - Sex x UPDRS ####

ps_updrs_RA <- transform_sample_counts(ps_updrs, fun = function(x) x / sum(x))

# Subset by group
updrs_mild_m     <- subset_samples(ps_updrs_RA, Sex_UPDRSScore == "Male_Mild")
updrs_moderate_m <- subset_samples(ps_updrs_RA, Sex_UPDRSScore == "Male_Moderate")
updrs_severe_m <- subset_samples(ps_updrs_RA, Sex_UPDRSScore == "Male_Severe")

updrs_mild_f     <- subset_samples(ps_updrs_RA, Sex_UPDRSScore == "Female_Mild")
updrs_moderate_f <- subset_samples(ps_updrs_RA, Sex_UPDRSScore == "Female_Moderate")


# Core ASVs per group
# detection = 0.001: ASV must be >= 0.1% relative abundance
# prevalence = 0.1:  ASV must appear in >= 10% of samples in that group
updrs_mild_m_ASVs     <- core_members(updrs_mild_m,     detection = 0.001, prevalence = 0.1)
updrs_moderate_m_ASVs <- core_members(updrs_moderate_m, detection = 0.001, prevalence = 0.1)
updrs_severe_m_ASVs <- core_members(updrs_severe_m, detection = 0.001, prevalence = 0.1)

updrs_mild_f_ASVs     <- core_members(updrs_mild_f,     detection = 0.001, prevalence = 0.1)
updrs_moderate_f_ASVs <- core_members(updrs_moderate_f, detection = 0.001, prevalence = 0.1)

cat("\nCore ASV counts - Sex x UPDRS:\n")
cat("  Male Mild:       ", length(updrs_mild_m_ASVs),     "\n")
cat("  Male Moderate:   ", length(updrs_moderate_m_ASVs), "\n")
cat("  Male Severe:     ", length(updrs_severe_m_ASVs), "\n")
cat("  Female Mild:     ", length(updrs_mild_f_ASVs),     "\n")
cat("  Female Moderate: ", length(updrs_moderate_f_ASVs), "\n")


# Venn diagram - Males: Mild vs Moderate
venn_updrs_male <- ggVennDiagram(
  x = list('Male Mild'     = updrs_mild_m_ASVs,
           'Male Moderate' = updrs_moderate_m_ASVs)) +
  scale_fill_gradient(low = "#F4F4F4", high = "#5B9BD5") +
  labs(title = "Core Microbiome (UPDRS): Male Mild vs Moderate", fill = "Count") +
  coord_fixed(clip = "off") +
  theme(plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text  = element_text(size = 10),
        plot.margin  = margin(10, 10, 10, 40, unit = "pt"))
print(venn_updrs_male)
ggsave("G9_A3_Core_UPDRS_Male_Venn.png", venn_updrs_male, bg = "transparent",
       width = 8, height = 7, dpi = 300)

# Venn diagram - Males: Mild vs Moderate vs Severe
venn_updrs_male_2 <- ggVennDiagram(
  x = list('Male Mild'     = updrs_mild_m_ASVs,
           'Male Moderate' = updrs_moderate_m_ASVs,
           'Male Severe' = updrs_severe_m_ASVs)) +
  scale_fill_gradient(low = "#F4F4F4", high = "#5B9BD5") +
  labs(title = "Core Microbiome (UPDRS): Male Mild vs Moderate vs Severe", fill = "Count") +
  coord_fixed(clip = "off") +
  theme(plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text  = element_text(size = 10),
        plot.margin  = margin(10, 10, 10, 40, unit = "pt"))
print(venn_updrs_male_2)
ggsave("G9_A3_Core_UPDRS_Male_Venn_2.png", venn_updrs_male_2, bg = "transparent",
       width = 8, height = 7, dpi = 300)

# Venn diagram - Females: Mild vs Moderate
venn_updrs_female <- ggVennDiagram(
  x = list('Female Mild'     = updrs_mild_f_ASVs,
           'Female Moderate' = updrs_moderate_f_ASVs)) +
  scale_fill_gradient(low = "#F4F4F4", high = "#C0392B") +
  labs(title = "Core Microbiome (UPDRS): Female Mild vs Moderate", fill = "Count") +
  coord_fixed(clip = "off") +
  theme(plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text  = element_text(size = 10),
        plot.margin  = margin(10, 10, 10, 40, unit = "pt"))
print(venn_updrs_female)
ggsave("G9_A3_Core_UPDRS_Female_Venn.png", venn_updrs_female, bg = "transparent",
       width = 8, height = 7, dpi = 300)

# Venn diagram - All four groups together
venn_updrs_all <- ggVennDiagram(
  x = list('Male Mild'       = updrs_mild_m_ASVs,
           'Male Moderate'   = updrs_moderate_m_ASVs,
           'Female Mild'     = updrs_mild_f_ASVs,
           'Female Moderate' = updrs_moderate_f_ASVs)) +
  scale_fill_gradient(low = "#F4F4F4", high = "#8B6BB1", limits = c(0, 155)) +
  labs(title = "Core Microbiome (UPDRS): All Sex x Severity Groups", fill = "Count") +
  coord_fixed(clip = "off") +
  theme(plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text  = element_text(size = 10),
        plot.margin  = margin(10, 10, 10, 50, unit = "pt"))
print(venn_updrs_all)
ggsave("G9_A3_Core_UPDRS_AllGroups_Venn.png", venn_updrs_all, bg = "transparent",
       width = 9, height = 8, dpi = 300)
ggsave("G9_A3_Core_UPDRS_AllGroups_Venn.pdf", venn_updrs_all, bg = "transparent",
       width = 9, height = 8)

# Venn diagram - All four groups together
venn_updrs_all_2 <- ggVennDiagram(
  x = list('Male Mild'       = updrs_mild_m_ASVs,
           'Male Moderate'   = updrs_moderate_m_ASVs,
           'Male Severe' = updrs_severe_m_ASVs,
           'Female Mild'     = updrs_mild_f_ASVs,
           'Female Moderate' = updrs_moderate_f_ASVs)) +
  scale_fill_gradient(low = "#F4F4F4", high = "#8B6BB1") +
  labs(title = "Core Microbiome (UPDRS): All Sex x Severity Groups", fill = "Count") +
  coord_fixed(clip = "off") +
  theme(plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text  = element_text(size = 10),
        plot.margin  = margin(10, 10, 10, 50, unit = "pt"))
print(venn_updrs_all_2)
ggsave("G9_A3_Core_UPDRS_AllGroups_Venn_2.png", venn_updrs_all_2, bg = "transparent",
       width = 10, height = 9, dpi = 300)

# CORE MICROBIOME: Sex x MoCA 
# Groups: Male_Normal, Male_Impaired, Female_Normal, Female_Impaired
# Key comparisons:
#   - Within males:   Normal vs Impaired cognition
#   - Within females: Normal vs Impaired cognition
#   - Across sexes:   Male vs Female within same MoCA cognitive bin


#### Core Microbiome - Sex x MoCA ####

ps_moca_RA <- transform_sample_counts(ps_moca, fun = function(x) x / sum(x))

# Subset by group
moca_normal_m   <- subset_samples(ps_moca_RA, Sex_MOCAScore == "Male_Normal")
moca_impaired_m <- subset_samples(ps_moca_RA, Sex_MOCAScore == "Male_Impaired")
moca_normal_f   <- subset_samples(ps_moca_RA, Sex_MOCAScore == "Female_Normal")
moca_impaired_f <- subset_samples(ps_moca_RA, Sex_MOCAScore == "Female_Impaired")

# Core ASVs per group
moca_normal_m_ASVs   <- core_members(moca_normal_m,   detection = 0.001, prevalence = 0.1)
moca_impaired_m_ASVs <- core_members(moca_impaired_m, detection = 0.001, prevalence = 0.1)
moca_normal_f_ASVs   <- core_members(moca_normal_f,   detection = 0.001, prevalence = 0.1)
moca_impaired_f_ASVs <- core_members(moca_impaired_f, detection = 0.001, prevalence = 0.1)

cat("\nCore ASV counts - Sex x MoCA:\n")
cat("  Male Normal:     ", length(moca_normal_m_ASVs),   "\n")
cat("  Male Impaired:   ", length(moca_impaired_m_ASVs), "\n")
cat("  Female Normal:   ", length(moca_normal_f_ASVs),   "\n")
cat("  Female Impaired: ", length(moca_impaired_f_ASVs), "\n")


# Venn diagram - Males: Normal vs Impaired
venn_moca_male <- ggVennDiagram(x = list('Male Normal'   = moca_normal_m_ASVs,
                                         'Male Impaired' = moca_impaired_m_ASVs)) +
  scale_fill_gradient(low = "#F4F4F4", high = "#5B9BD5") +
  labs(title = "Core Microbiome (MoCA): Male Normal vs Impaired", fill = "Count") +
  coord_fixed(clip = "off") +
  theme(plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text  = element_text(size = 10),
        plot.margin  = margin(10, 10, 10, 40, unit = "pt"))
print(venn_moca_male)
ggsave("G9_A3_Core_MoCA_Male_Venn.png", venn_moca_male, bg = "transparent",
       width = 8, height = 7, dpi = 300)

# Venn diagram - Females: Normal vs Impaired
venn_moca_female <- ggVennDiagram(x = list('Female Normal'   = moca_normal_f_ASVs,
                                           'Female Impaired' = moca_impaired_f_ASVs)) +
  scale_fill_gradient(low = "#F4F4F4", high = "#C0392B") +
  labs(title = "Core Microbiome (MoCA): Female Normal vs Impaired", fill = "Count") +
  coord_fixed(clip = "off") +
  theme(plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text  = element_text(size = 10),
        plot.margin  = margin(10, 10, 10, 40, unit = "pt"))
print(venn_moca_female)
ggsave("G9_A3_Core_MoCA_Female_Venn.png", venn_moca_female, bg = "transparent",
       width = 8, height = 7, dpi = 300)

#Venn diagram - All four groups together
venn_moca_all <- ggVennDiagram(x = list('Male Normal'     = moca_normal_m_ASVs,
                                        'Male Impaired'   = moca_impaired_m_ASVs,
                                        'Female Normal'   = moca_normal_f_ASVs,
                                        'Female Impaired' = moca_impaired_f_ASVs)) +
  scale_fill_gradient(low = "#F4F4F4", high = "#8B6BB1", limits = c(0,155)) +
  labs(title = "Core Microbiome (MoCA): All Sex x Cognitive Groups", fill = "Count") +
  coord_fixed(clip = "off") +
  theme(plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text  = element_text(size = 10),
        plot.margin  = margin(10, 10, 10, 50, unit = "pt"))
print(venn_moca_all)
ggsave("G9_A3_Core_MoCA_AllGroups_Venn.png", venn_moca_all, bg = "transparent",
       width = 9, height = 8, dpi = 300)
ggsave("G9_A3_Core_MoCA_AllGroups_Venn.pdf", venn_moca_all, bg = "transparent",
       width = 9, height = 8)



