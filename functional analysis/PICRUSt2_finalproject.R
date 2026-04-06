### MICB 475 Group 9 Final Project WT2026
## PICRUSt2

# Load data and libraries
# install.packages("devtools")
# devtools::install_github("cafferychen777/ggpicrust2")
# install.packages('MicrobiomeStat')
# BiocManager::install("KEGGREST")
# install.packages("GGally")

library(tidyverse)
library(phyloseq)
library(ggpicrust2)

# Extracted metadata from phyloseq object saved in project as a .RData object
meta <- sample_data(parkinsons_phyloseq) %>% 
  data.frame() %>% 
  rownames_to_column("sample_name")

# Load in combined KO data from server
ko = read.delim('pred_metagenome_unstrat.tsv')

#### B: Running each step separately ####

#  Differential analysis
# First column should be row names when running pipeline manually, so I'm just reloading it as such
ko = read.delim('pred_metagenome_unstrat.tsv',row.names = 1)
rownames(ko) <- gsub("^ko:", "", 
                                 rownames(ko))

# 1a) UPDRS: Male Mild vs Male Moderate

# Filter meta data and ko data to compare two groups (UPDRS male mild vs male moderate)
M_mildvsmod = meta %>% filter(Sex_UPDRSScore %in% c("Male_Mild", "Male_Moderate"))
UPDRSM_ko_filt = ko %>% select(all_of(M_mildvsmod$sample_name))

# Perform pathway differential abundance analysis (DAA) using LinDA method
UPDRSM_daa_results_df = pathway_daa(abundance = UPDRSM_ko_filt,
                                    metadata = M_mildvsmod, 
                                    group = "Sex_UPDRSScore", 
                                    daa_method = "LinDA", 
                                    select = NULL, reference = NULL)

# UPDRSM_daa_annotated_results_df = pathway_annotation(pathway = "KO",
#                                              daa_results_df = UPDRSM_daa_results_df,
#                                              ko_to_kegg = TRUE)

# No significant pathways above so I'm using this instead where ko_to_kegg = FALSE instead
UPDRSM_daa_annotated_results_df = pathway_annotation(pathway = "KO",
                                  daa_results_df = UPDRSM_daa_results_df,
                                  ko_to_kegg = TRUE)

# Save as R object so it can be reopened without needing to rerun the annotation
saveRDS(UPDRSM_daa_annotated_results_df,'UPDRSM_daa_annotated_results_df.rds')

# Plot results
source('ggpicrust2_errorbar_function_fixed.R')

UPDRSM_daa_annotated_results_df %>% filter(p_adjust< 0.1, abs(log2_fold_change) > 2) %>% nrow()

UPDRSM_daa_filtered <- UPDRSM_daa_annotated_results_df %>%
  dplyr::filter(abs(log2_fold_change) >= 2, p_adjust <= 0.07)

UPDRSM_ko_filt_filtered <- UPDRSM_ko_filt[UPDRSM_daa_filtered$feature, ]


## This function does NOT work without sigificant daa results
# UPDRSM_plot <- pathway_errorbar_fixed(
#   abundance = UPDRSM_ko_filt_filtered,
#   daa_results_df = UPDRSM_daa_filtered,
#   Group = M_mildvsmod$Sex_UPDRSScore,
#   wrap_label = TRUE,
#   wraplength = 60,
#   fc_cutoff = 2,
#   order_by_log = FALSE,
#   p_values_threshold = 0.07,
#   order = "description",  # use KO descriptions
#   ko_to_kegg = FALSE,
#   p_value_bar = TRUE,
#   x_lab = "description"
# )

# UPDRSM_plot = pathway_errorbar_fixed(abundance = UPDRSM_ko_filt,
#                              daa_results_df = UPDRSM_daa_annotated_results_df,
#                              Group = M_mildvsmod$Sex_UPDRSScore,
#                              wrap_label = T, wraplength=60,
#                              fc_cutoff = 2, order_by_log = F,
#                              p_values_threshold = 0.07,
#                              order = "description",
#                              ko_to_kegg = FALSE,
#                              p_value_bar = TRUE,
#                              x_lab = "description")


# 1b) UPDRS: Female Mild vs Female Moderate
F_mildvsmod = meta %>% filter(Sex_UPDRSScore %in% c("Female_Mild", "Female_Moderate"))
UPDRSF_ko_filt = ko %>% select(all_of(F_mildvsmod$sample_name))

UPDRSF_daa_results_df = pathway_daa(abundance = UPDRSF_ko_filt,
                                    metadata = F_mildvsmod, 
                                    group = "Sex_UPDRSScore", 
                                    daa_method = "LinDA", 
                                    select = NULL, reference = NULL)

UPDRSF_daa_annotated_results_df = pathway_annotation(pathway = "KO",
                                                     daa_results_df = UPDRSF_daa_results_df,
                                                     ko_to_kegg = TRUE)

saveRDS(UPDRSF_daa_annotated_results_df,'UPDRSF_daa_annotated_results_df.rds')

# 2a) MoCA: Male Normal vs Male Impaired
M_normvsimp = meta %>% filter(Sex_MOCAScore %in% c("Male_Normal", "Male_Impaired"))
MoCAM_ko_filt = ko %>% select(all_of(M_normvsimp$sample_name))

MoCAM_daa_results_df = pathway_daa(abundance = MoCAM_ko_filt,
                                    metadata = M_normvsimp, 
                                    group = "Sex_MOCAScore", 
                                    daa_method = "LinDA", 
                                    select = NULL, reference = NULL)

MoCAM_daa_annotated_results_df = pathway_annotation(pathway = "KO",
                                                     daa_results_df = MoCAM_daa_results_df,
                                                     ko_to_kegg = TRUE)

saveRDS(MoCAM_daa_annotated_results_df,'MoCAM_daa_annotated_results_df.rds')

# 2b) MoCA: Female Mild vs Female Impaired
F_normvsimp = meta %>% filter(Sex_MOCAScore %in% c("Female_Normal", "Female_Impaired"))
MoCAF_ko_filt = ko %>% select(all_of(F_normvsimp$sample_name))

MoCAF_daa_results_df = pathway_daa(abundance = MoCAF_ko_filt,
                                   metadata = F_normvsimp, 
                                   group = "Sex_MOCAScore", 
                                   daa_method = "LinDA", 
                                   select = NULL, reference = NULL)

MoCAF_daa_annotated_results_df = pathway_annotation(pathway = "KO",
                                                    daa_results_df = MoCAF_daa_results_df,
                                                    ko_to_kegg = TRUE)

saveRDS(MoCAF_daa_annotated_results_df,'MoCAF_daa_annotated_results_df.rds')

# 3a) UPDRS: Male Mild vs Female Mild

MmildvsFmile = meta %>% filter(Sex_UPDRSScore %in% c("Male_Mild", "Female_Mild"))
UPDRSMmildvsFmile_ko_filt = ko %>% select(all_of(MmildvsFmile$sample_name))

UPDRSMmildvsFmile_daa_results_df = pathway_daa(abundance = UPDRSMmildvsFmile_ko_filt,
                                    metadata = MmildvsFmile, 
                                    group = "Sex_UPDRSScore", 
                                    daa_method = "LinDA", 
                                    select = NULL, reference = NULL)

UPDRSMmildvsFmile_daa_annotated_results_df = pathway_annotation(pathway = "KO",
                                                     daa_results_df = UPDRSMmildvsFmile_daa_results_df,
                                                     ko_to_kegg = TRUE)

saveRDS(UPDRSMmildvsFmile_daa_annotated_results_df,'UPDRSMmildvsFmile_daa_annotated_results_df.rds')

# 3b) UPDRS: Male Moderate vs Female Moderate
MmodvsFmod = meta %>% filter(Sex_UPDRSScore %in% c("Male_Moderate", "Female_Moderate"))
UPDRSMmodvsFmod_ko_filt = ko %>% select(all_of(MmodvsFmod$sample_name))

UPDRSMmodvsFmod_daa_results_df = pathway_daa(abundance = UPDRSMmodvsFmod_ko_filt,
                                               metadata = MmodvsFmod, 
                                               group = "Sex_UPDRSScore", 
                                               daa_method = "LinDA", 
                                               select = NULL, reference = NULL)

UPDRSMmodvsFmod_daa_annotated_results_df = pathway_annotation(pathway = "KO",
                                                                daa_results_df = UPDRSMmodvsFmod_daa_results_df,
                                                                ko_to_kegg = TRUE)

saveRDS(UPDRSMmodvsFmod_daa_annotated_results_df,'UPDRSMmodvsFmod_daa_annotated_results_df.rds')

# 4a) MoCA: Male Normal vs Female Normal
MnormvsFnorm = meta %>% filter(Sex_MOCAScore %in% c("Female_Normal", "Male_Normal"))
MoCAMnormvsFnorm_ko_filt = ko %>% select(all_of(MnormvsFnorm$sample_name))

MoCAMnormvsFnorm_daa_results_df = pathway_daa(abundance = MoCAMnormvsFnorm_ko_filt,
                                   metadata = MnormvsFnorm, 
                                   group = "Sex_MOCAScore", 
                                   daa_method = "LinDA", 
                                   select = NULL, reference = NULL)

MoCAMnormvsFnorm_daa_annotated_results_df = pathway_annotation(pathway = "KO",
                                                    daa_results_df = MoCAMnormvsFnorm_daa_results_df,
                                                    ko_to_kegg = TRUE)

saveRDS(MoCAMnormvsFnorm_daa_annotated_results_df,'MoCAMnormvsFnorm_daa_annotated_results_df.rds')

# 4b) MoCA: Male Impaired vs Female Impaired
MimpvsFimp = meta %>% filter(Sex_MOCAScore %in% c("Female_Impaired", "Male_Impaired"))
MimpvsFimp_ko_filt = ko %>% select(all_of(MimpvsFimp$sample_name))

MimpvsFimp_daa_results_df = pathway_daa(abundance = MimpvsFimp_ko_filt,
                                              metadata = MimpvsFimp, 
                                              group = "Sex_MOCAScore", 
                                              daa_method = "LinDA", 
                                              select = NULL, reference = NULL)

MimpvsFimp_daa_annotated_results_df = pathway_annotation(pathway = "KO",
                                                               daa_results_df = MimpvsFimp_daa_results_df,
                                                               ko_to_kegg = TRUE)

saveRDS(MimpvsFimp_daa_annotated_results_df,'MimpvsFimp_daa_annotated_results_df.rds')

# To read the .rds file
daa_annotated_results_df = readRDS('file_name.rds')

#### C: Generate Plots ####
# source('ggpicrust2_errorbar_function_fixed.R')


#### Other ####

# Code below needed to generate UPDRS male heatmap (not included as a figure) 

UPDRSM_sig_features = UPDRSM_daa_annotated_results_df %>% 
  filter(p_adjust < 0.07, abs(log2_fold_change) > 2) %>% 
  pull('feature') %>% unique()

UPDRSM_sig_features <- paste0("ko:", UPDRSM_sig_features)

ko_relab = read.delim('pred_metagenome_unstrat.tsv',row.names = 1) %>% 
  apply(2,function(x) x/sum(x)) %>% as.data.frame()
colSums(ko_relab) # Should be all 1

stats <- ko_relab[sig_features, ] %>% 
  t() %>% 
  cor(method = "pearson")

# stats = ko_relab %>% t() %>%  # switch rows and columns
#   as.data.frame() %>% # Lets us use select()
#   select(all_of(sig_features)) %>% 
#   cor(method = 'pearson') #Pearson tests between every feature pair

# Define the color palette
color_palette <- colorRampPalette(c("blue", "white", "red"))(40)

# Define breaks for the color scale
breaks <- seq(-1, 1, length.out = 41)  # 40 colors + 1 for the endpoint

# Create the clustered heatmap with centered color at zero
# install.packages("pheatmap")
library('pheatmap')
pheatmap(stats, 
         clustering_distance_rows = "euclidean",  
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",          
         color = color_palette, 
         breaks = breaks,
         main = '',
         fontsize_row = 10, 
         fontsize_col = 10,
         filename = "UPDRSM Kegg Heatmap.png",
         height= 9, width = 9)

# Generate pathway PCA plots (included as figures)

# Male UPDRS PCA plot
pathway_pca(abundance = UPDRSM_ko_filt,
            metadata = M_mildvsmod, 
            group = "Sex_UPDRSScore")

# Female UPDRS PCA Plot
pathway_pca(abundance = UPDRSF_ko_filt,
            metadata = F_mildvsmod, 
            group = "Sex_UPDRSScore")

