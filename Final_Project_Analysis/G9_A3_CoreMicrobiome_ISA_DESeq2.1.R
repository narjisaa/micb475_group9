# Aim 3: Core Microbiome, Indicator Species Analysis, and DESeq2

#### Load Libraries ####
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(indicspecies)
library(DESeq2)
library(ggplot2)

#### Load Data ####
load("new_phyloseq.RData")
sample_data(parkinsons_phyloseq)

################################ STEP 1: Inspect group labels ################################
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

################################ AIM 3.1 - CORE MICROBIOME: Sex x UPDRS ################################
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
  scale_fill_gradient(low = "#F4F4F4", high = "#8B6BB1") +
  labs(title = "Core Microbiome (UPDRS): All Sex x Severity Groups", fill = "Count") +
  coord_fixed(clip = "off") +
  theme(plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text  = element_text(size = 10),
        plot.margin  = margin(10, 10, 10, 50, unit = "pt"))
print(venn_updrs_all)
ggsave("G9_A3_Core_UPDRS_AllGroups_Venn.png", venn_updrs_all, bg = "transparent",
       width = 9, height = 8, dpi = 300)

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


################################ AIM 3.1 - CORE MICROBIOME: Sex x MoCA ################################ 
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

# Venn diagram - All four groups together
venn_moca_all <- ggVennDiagram(x = list('Male Normal'     = moca_normal_m_ASVs,
                               'Male Impaired'   = moca_impaired_m_ASVs,
                               'Female Normal'   = moca_normal_f_ASVs,
                               'Female Impaired' = moca_impaired_f_ASVs)) +
  scale_fill_gradient(low = "#F4F4F4", high = "#8B6BB1") +
  labs(title = "Core Microbiome (MoCA): All Sex x Cognitive Groups", fill = "Count") +
  coord_fixed(clip = "off") +
  theme(plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text  = element_text(size = 10),
        plot.margin  = margin(10, 10, 10, 50, unit = "pt"))
print(venn_moca_all)
ggsave("G9_A3_Core_MoCA_AllGroups_Venn.png", venn_moca_all, bg = "transparent",
       width = 9, height = 8, dpi = 300)

################################ Helper Functions for ISA and DESeq2 ################################ 
run_isa <- function(ps_obj, group_col, filename_csv, filename_plot,
                    plot_title, group_colours) {

  ps_genus    <- tax_glom(ps_obj, "Genus", NArm = FALSE)
  ps_genus_RA <- transform_sample_counts(ps_genus, fun = function(x) x / sum(x))

  otu_matrix <- t(otu_table(ps_genus_RA))
  group_vec  <- sample_data(ps_genus_RA)[[group_col]]

  set.seed(475)
  isa_res <- multipatt(as.data.frame(otu_matrix),
                       cluster = group_vec,
                       func    = "r.g",
                       control = how(nperm = 999))
  summary(isa_res)

  taxtable <- tax_table(ps_obj) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ASV")

  results_df <- isa_res$sign %>%
    rownames_to_column(var = "ASV") %>%
    left_join(taxtable, by = "ASV") %>%
    filter(p.value < 0.05) %>%
    arrange(p.value)

  cat("\nSignificant indicator taxa (p < 0.05) for", group_col, ":",
      nrow(results_df), "\n")

  write.csv(results_df, filename_csv, row.names = FALSE)

  grp_levels <- sort(unique(as.character(group_vec)))
  s_cols     <- paste0("s.", seq_along(grp_levels))

  if (nrow(results_df) > 0) {

    top_df <- results_df %>%
      slice_max(order_by = stat, n = 20) %>%
      mutate(Genus = ifelse(is.na(Genus) | Genus == "",
                            paste0("Unknown_", Family), Genus),
             Genus = make.unique(Genus),
             Genus = factor(Genus, levels = rev(unique(Genus))))

    top_df$indicator_group <- apply(
      top_df[, intersect(s_cols, colnames(top_df)), drop = FALSE], 1,
      function(row) {
        idx <- which(row == 1)
        if (length(idx) == 1) grp_levels[idx] else "Multi-group"
      })

    isa_plot <- ggplot(top_df,
                       aes(x = stat, y = Genus, fill = indicator_group)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = group_colours, na.value = "#AAAAAA") +
      labs(title = plot_title,
           x = "Indicator Value (r)", y = "Genus", fill = "Group") +
      theme_bw() +
      theme(plot.title      = element_text(hjust = 0.5, face = "bold", size = 14),
            axis.title      = element_text(size = 12),
            axis.text.x     = element_text(size = 11),
            axis.text.y     = element_text(size = 11),
            legend.title    = element_text(size = 11),
            legend.text     = element_text(size = 10),
            legend.position = "right")

    print(isa_plot)
    ggsave(filename_plot, isa_plot, bg = "transparent",
           width = 12, height = max(7, nrow(top_df) * 0.4 + 3), dpi = 300)

  } else {
    cat("No significant indicator taxa to plot for", group_col, "\n")
  }

  invisible(results_df)
}

################################ Helper Functions for DESeq2 ################################ 
run_deseq2 <- function(ps_obj, group_col, contrasts, filename_prefix) {

  ps_plus1  <- transform_sample_counts(ps_obj, function(x) x + 1)
  formula   <- as.formula(paste("~", group_col))
  ps_ds     <- phyloseq_to_deseq2(ps_plus1, formula)

  # Set Healthy Male as reference if present, otherwise first alphabetical level
  lvls      <- levels(ps_ds[[group_col]])
  ref_level <- if (any(grepl("Healthy", lvls))) lvls[grepl("Healthy", lvls)][1] else sort(lvls)[1]
  ps_ds[[group_col]] <- relevel(ps_ds[[group_col]], ref = ref_level)

  deseq_obj <- DESeq(ps_ds)

  taxtable <- tax_table(ps_obj) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ASV")

  result_list <- list()

  for (contrast in contrasts) {

    label     <- paste(contrast[2], "vs", contrast[3])
    safe_name <- gsub("[^A-Za-z0-9]", "_", label)

    res <- results(deseq_obj, tidy = TRUE, contrast = contrast)

    sig <- res %>%
      filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
      dplyr::rename(ASV = row)

    cat("\nSignificant ASVs -", label, ":", nrow(sig), "\n")

    if (nrow(sig) == 0) {
      cat("  No significant ASVs. Try: padj < 0.1 and/or |LFC| > 1.\n")
      result_list[[label]] <- NULL
      next
    }

    sig_vec <- sig %>% pull(ASV)
    ps_filt <- prune_taxa(sig_vec, ps_obj)

    merged <- tax_table(ps_filt) %>%
      as.data.frame() %>%
      rownames_to_column(var = "ASV") %>%
      right_join(sig, by = "ASV") %>%
      arrange(log2FoldChange) %>%
      mutate(
        Genus     = ifelse(is.na(Genus) | Genus == "",
                           paste0("Unknown_", Family), Genus),
        Genus     = make.unique(Genus),
        Genus     = factor(Genus, levels = unique(Genus)),
        direction = ifelse(log2FoldChange > 0, "Enriched", "Depleted")
      )

    bar_plot <- ggplot(merged,
                       aes(x = Genus, y = log2FoldChange, fill = direction)) +
      geom_bar(stat = "identity") +
      geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                        ymax = log2FoldChange + lfcSE),
                    width = 0.3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      scale_fill_manual(values = c("Enriched" = "#C0392B",
                                   "Depleted" = "#4CAF82")) +
      coord_flip() +
      labs(title = paste("Differential Abundance:", label),
           x = "Genus", y = "Log2 Fold Change", fill = "Direction") +
      theme_bw() +
      theme(plot.title      = element_text(hjust = 0.5, face = "bold", size = 14),
            axis.title      = element_text(size = 12),
            axis.text.x     = element_text(size = 11),
            axis.text.y     = element_text(size = 11),
            legend.title    = element_text(size = 11),
            legend.text     = element_text(size = 10),
            legend.position = "right")

    print(bar_plot)
    ggsave(paste0(filename_prefix, "_", safe_name, ".png"),
           bar_plot, bg = "transparent",
           width = 12, height = max(7, nrow(merged) * 0.4 + 3), dpi = 300)

    write.csv(merged,
              paste0(filename_prefix, "_", safe_name, ".csv"),
              row.names = FALSE)

    result_list[[label]] <- list(merged = merged, full_res = res)
  }

################################ Combined Heatmap ################################ 
  all_sig_ASVs <- unique(unlist(lapply(result_list, function(x) {
    if (!is.null(x)) x$merged$ASV else character(0)
  })))

  if (length(all_sig_ASVs) > 0) {

    lfc_rows <- lapply(names(result_list), function(label) {
      if (is.null(result_list[[label]])) return(NULL)
      result_list[[label]]$full_res %>%
        dplyr::rename(ASV = row) %>%
        filter(ASV %in% all_sig_ASVs) %>%
        mutate(Comparison = label)
    })

    lfc_df <- bind_rows(lfc_rows)

    tax_df <- taxtable %>%
      mutate(Genus_label = ifelse(is.na(Genus) | Genus == "",
                                  paste0("Unknown_", Family), Genus))

    lfc_df <- left_join(lfc_df, tax_df[, c("ASV", "Genus_label")], by = "ASV") %>%
      mutate(Genus_label = make.unique(Genus_label))

    heatmap <- ggplot(lfc_df,
                      aes(x = Comparison, y = Genus_label, fill = log2FoldChange)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "#4CAF82", mid = "white", high = "#C0392B",
                           midpoint = 0, name = "Log2FC") +
      labs(title = paste("DESeq2 Heatmap -", group_col),
           x = "Comparison", y = "Genus") +
      theme_bw() +
      theme(plot.title  = element_text(hjust = 0.5, face = "bold", size = 14),
            axis.title  = element_text(size = 12),
            axis.text.x = element_text(angle = 35, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            legend.title = element_text(size = 11),
            legend.text  = element_text(size = 10))

    print(heatmap)
    ggsave(paste0(filename_prefix, "_Heatmap.png"), heatmap, bg = "transparent",
           width = max(10, length(unique(lfc_df$Comparison)) * 1.2 + 4),
           height = max(7, length(all_sig_ASVs) * 0.35 + 3),
           dpi = 300)
  }

  invisible(result_list)
}

################################ AIM 3.2 + 3.3 - ISA & DESeq2: Sex x UPDRS################################ 
# Contrasts:
#   Within-male disease severity:   Male_Moderate vs Male_Mild
#   Within-female disease severity: Female_Moderate vs Female_Mild
#   Sex effect at Mild UPDRS:       Male_Mild vs Female_Mild
#   Sex effect at Moderate UPDRS:   Male_Moderate vs Female_Moderate

cat("\n\n=== ANALYSIS: Sex x UPDRS ===\n")

## ISA - Sex x UPDRS
run_isa(
  ps_obj        = ps_updrs,
  group_col     = "Sex_UPDRSScore",
  filename_csv  = "G9_A3_ISA_UPDRS_Results.csv",
  filename_plot = "G9_A3_ISA_UPDRS_Plot.png",
  plot_title    = "Indicator Taxa - Sex x UPDRS Groups",
  group_colours = c(
    "Male_Mild"       = "#5B9BD5",
    "Male_Moderate"   = "#1A5276",
    "Male_Severe"     = "#003354",
    "Female_Mild"     = "#E8857A",
    "Female_Moderate" = "#C0392B",
    "Multi-group"     = "#8B6BB1"
  )
)

## DESeq2 - Sex x UPDRS
run_deseq2(
  ps_obj          = ps_updrs,
  group_col       = "Sex_UPDRSScore",
  contrasts       = list(
    # Disease severity within each sex
    c("Sex_UPDRSScore", "Male_Severe",   "Male_Mild"),
    c("Sex_UPDRSScore", "Male_Moderate",   "Male_Severe"),
    
    c("Sex_UPDRSScore", "Male_Moderate",   "Male_Mild"),
    c("Sex_UPDRSScore", "Female_Moderate", "Female_Mild"),
    
    # Sex effect within each severity bin
    c("Sex_UPDRSScore", "Male_Severe",       "Female_Mild"),
    c("Sex_UPDRSScore", "Male_Severe",   "Female_Moderate"),
    
    c("Sex_UPDRSScore", "Male_Mild",       "Female_Mild"),
    c("Sex_UPDRSScore", "Male_Moderate",   "Female_Moderate")
    
  ),
  filename_prefix = "G9_A3_DESeq2_UPDRS"
)

################################ AIM 3.2 + 3.3 - ISA & DESeq2: Sex x MoCA ################################ 
# Contrasts:
#   Within-male cognitive status:   Male_Impaired vs Male_Normal
#   Within-female cognitive status: Female_Impaired vs Female_Normal
#   Sex effect in Normal cognition:   Male_Normal vs Female_Normal
#   Sex effect in Impaired cognition: Male_Impaired vs Female_Impaired

cat("\n\n=== ANALYSIS: Sex x MoCA ===\n")

## ISA - Sex x MoCA
run_isa(
  ps_obj        = ps_moca,
  group_col     = "Sex_MOCAScore",
  filename_csv  = "G9_A3_ISA_MoCA_Results.csv",
  filename_plot = "G9_A3_ISA_MoCA_Plot.png",
  plot_title    = "Indicator Taxa - Sex x MoCA Groups",
  group_colours = c(
    "Male_Normal"     = "#5B9BD5",
    "Male_Impaired"   = "#1A5276",
    "Female_Normal"   = "#E8857A",
    "Female_Impaired" = "#C0392B",
    "Multi-group"     = "#8B6BB1"
  )
)

## DESeq2 - Sex x MoCA
run_deseq2(
  ps_obj          = ps_moca,
  group_col       = "Sex_MOCAScore",
  contrasts       = list(
    # Cognitive impairment effect within each sex
    c("Sex_MOCAScore", "Male_Impaired",   "Male_Normal"),
    c("Sex_MOCAScore", "Female_Impaired", "Female_Normal"),
    # Sex effect within each cognitive bin
    c("Sex_MOCAScore", "Male_Normal",     "Female_Normal"),
    c("Sex_MOCAScore", "Male_Impaired",   "Female_Impaired")
  ),
  filename_prefix = "G9_A3_DESeq2_MoCA"
)

