library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(dplyr)
library(picante)
library(ggpubr)
library(rstatix)
library(indicspecies)
library(forcats)
library(tibble)
library(ggplot2)
library(DESeq2)

#### --- ISA / (DESeq2) Differential Abundance / Integrated Analysis ---####
#### --- ISA --- ####

# --- 1. Load phyloseq objects ---
load("covid_mexico_phyloseq_final.RData")  # unrarefied
load("covid_mexico_phyloseq_rare.RData")   # rarefied for ISA

# --- 2. Define sex-severity groups for both phyloseq objects ---
prepare_groups <- function(ps) {
  meta <- data.frame(sample_data(ps))
  if(!"Group_combined" %in% colnames(meta)) {
    meta$Severity <- NA_character_
    meta$Severity[grepl("asymptomatic", tolower(meta$severity))] <- "Control"
    meta$Severity[grepl("ambulatory positive", tolower(meta$severity))] <- "Mild"
    meta$Severity[grepl("hospitalized positive", tolower(meta$severity))] <- "Severe"
    meta$Severity[grepl("deceased", tolower(meta$severity))] <- "Deceased"
    meta$Group_combined <- paste(meta$sex, meta$Severity, sep="_")
    sample_data(ps) <- sample_data(meta)
  }
  return(ps)
}
# Apply the same 'prepare_groups' function to the rarefied phyloseq object
# This ensures that both unrarefied and rarefied datasets have consistent metadata/grouping
covid_mexico_phyloseq_final <- prepare_groups(covid_mexico_phyloseq_final)
covid_mexico_phyloseq_rare <- prepare_groups(covid_mexico_phyloseq_rare)

# --- 3. Collapse to Genus level ---
ps_genus_rare <- tax_glom(covid_mexico_phyloseq_rare, taxrank="Genus", NArm=TRUE)
ps_genus_unrare <- tax_glom(covid_mexico_phyloseq_final, taxrank="Genus", NArm=TRUE)

# --- 4. Indicator Species Analysis (ISA) - rarefied ---
otu_mat <- as.data.frame(t(otu_table(ps_genus_rare)))  # transpose OTU table
meta_df <- data.frame(sample_data(ps_genus_rare)) # sample metadata
meta_clean <- meta_df[!is.na(meta_df$Group_combined), ]  # remove samples with NA
otu_clean <- otu_mat[rownames(meta_clean), ]   # subset OTUs to clean samples
otu_clean[is.na(otu_clean)] <- 0  # replace any NA counts with 0

# Run ISA using IndVal.g
isa_mpt <- multipatt(
  otu_clean,
  meta_clean$Group_combined,
  func = "IndVal.g",
  duleg = TRUE,
  control = how(nperm = 999)
)
#Extract significant OTUs (p <= 0.05)
sig_table <- data.frame(isa_mpt$sign) %>%
  rownames_to_column("OTU") %>%
  filter(p.value <= 0.05)
# Assign clusters
group_names <- colnames(isa_mpt$comb)
sig_table$cluster <- group_names[sig_table$index]
# Add Genus names to significant OTUs
significant_indicators <- sig_table %>%
  select(OTU, stat, p.value, cluster)

# --- 5. Add Genus names to significant OTUs---
tax_df <- as.data.frame(tax_table(ps_genus_rare)) %>%
  rownames_to_column("OTU") %>%
  mutate(Label = Genus) %>%
  filter(!is.na(Label) & Label != "" & Label != "g__")

isa_plot_df <- significant_indicators %>%
  left_join(tax_df %>% select(OTU, Label), by="OTU") %>%
  arrange(cluster, desc(stat)) %>%
  group_by(cluster) %>% slice_max(stat, n=10) %>% ungroup()

# --- 6. Plot ISA ---

isa_plot <- ggplot(isa_plot_df, aes(x=reorder(Label, stat), y=stat, fill=cluster)) +
  geom_col() +
  coord_flip() +
  labs(title="Top Indicator Genera by Sex-Severity Cluster (ISA)",
       x="Genus", y="Indicator Value", fill="Sex-Severity Group") +
  theme_minimal() +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=12, face="bold"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=10))

print(isa_plot)
ggsave("ISA_and_DESeq2_Plots/ISA_top_indicator_genera.png", isa_plot, width=10, height=6, dpi=300)


#### ---- Differential Abundance Analysis ---- ####

# --- 2. Load Data ---
covid_mex_metadata_fp <- "COVID_Mexico_Final_Updated/covid_mex_metadata_edited.tsv"
covid_mex_metadata <- read_delim(covid_mex_metadata_fp, delim="\t")

covid_mex_otu_fp <- "COVID_Mexico_Final_Updated/covid_mexico_final_updated_export/table_export/feature-table.txt"
covid_mex_otu <- read_delim(covid_mex_otu_fp, delim="\t", skip=1)

covid_mex_tax_fp <- "COVID_Mexico_Final_Updated/covid_mexico_final_updated_export/taxonomy_export/taxonomy.tsv"
covid_mex_taxonomy <- read_delim(covid_mex_tax_fp, delim="\t")

covid_mex_phylotree_fp <- "COVID_Mexico_Final_Updated/covid_mexico_final_updated_export/rooted_tree_export/tree.nwk"
covid_mex_phylotree <- read.tree(covid_mex_phylotree_fp)

#load input files
meta <- read.delim(covid_mex_metadata_fp, sep = "\t", stringsAsFactors = FALSE)
otu  <- read.delim(covid_mex_otu_fp, sep = "\t", skip = 1, stringsAsFactors = FALSE)
tax  <- read.delim(covid_mex_tax_fp, sep = "\t", stringsAsFactors = FALSE)
phy  <- read.tree(covid_mex_phylotree_fp)

# Convert OTU table to matrix and set OTU IDs as row names
otu_mat <- as.matrix(otu[, -1])
rownames(otu_mat) <- otu$X.OTU.ID
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
# Create phyloseq OTU table 
tax_mat <- tax %>%
  select(-Confidence) %>%
  separate(Taxon, into = c("Domain","Phylum","Class","Order","Family","Genus","Species"),
           sep = ";", fill = "right", extra = "merge") %>%
  # Remove rank prefixes (e.g., g__, f__)
  mutate(across(everything(), ~ str_replace_all(.x, "^[a-z]__", ""))) %>%
  as.data.frame(stringsAsFactors = FALSE)
# Convert taxonomy to matrix and assign OTU IDs
tax_mat2 <- as.matrix(tax_mat[, -1])
rownames(tax_mat2) <- tax$Feature.ID
# Create phyloseq taxonomy table
TAX <- tax_table(tax_mat2)
# Prepare metadata for phyloseq
rownames(meta) <- meta$sample.id
META <- sample_data(meta)
# Combine OTU, taxonomy, metadata, and phylogeny into one object --> unrarefied phyloseq
ps <- phyloseq(OTU, TAX, META, phy)

# --- 4. Define sex-severity groups ---
# Replace blank or invalid sex entries with NA
sample_data(ps)$sex[sample_data(ps)$sex %in% c("", " ", "NA")] <- NA
# Initialize Severity variable
sample_data(ps)$Severity <- NA
# Assign disease severity based on Group field
sample_data(ps)$Severity[grepl("asymptomatic", sample_data(ps)$Group, ignore.case = TRUE)] <- "Control"
sample_data(ps)$Severity[grepl("ambulatory positive", sample_data(ps)$Group, ignore.case = TRUE)] <- "Mild"
sample_data(ps)$Severity[grepl("hospitalized positive", sample_data(ps)$Group, ignore.case = TRUE)] <- "Severe"
sample_data(ps)$Severity[grepl("Deceased hospitalized", sample_data(ps)$Group, ignore.case = TRUE)] <- "Deceased"
# Combine sex and severity into a single grouping variable
sample_data(ps)$Group_combined <- paste(sample_data(ps)$sex, sample_data(ps)$Severity, sep = "_")
# Define the desired ordering of groups
desired_levels <- c("male_Control","male_Mild","male_Severe","male_Deceased",
                    "female_Control","female_Mild","female_Severe","female_Deceased")
# Keep only group levels present in the data
present_levels <- intersect(desired_levels, unique(sample_data(ps)$Group_combined))

# Convert Group_combined to a factor with controlled ordering
sample_data(ps)$Group_combined <- factor(sample_data(ps)$Group_combined, levels = present_levels)
# Remove samples lacking valid group information
ps <- subset_samples(ps, !is.na(Group_combined))

# --- 5. Filter sparse taxa ---
# Minimum number of samples in which a taxon must appear
min_prevalence <- 3
# Minimum count threshold per sample
min_count <- 5
# Retain taxa that meet prevalence and total abundance criteria
ps_filt <- filter_taxa(ps, function(x) sum(x >= min_count) >= min_prevalence & sum(x) > 10, TRUE)

# --- 6. Collapse to genus level ---
# Aggregate OTUs to the genus rank
ps_genus <- tax_glom(ps_filt, taxrank = "Genus", NArm = TRUE)
# Extract genus labels and keep only classified genera
tax_df_genus <- as.data.frame(tax_table(ps_genus)) %>%
  rownames_to_column("OTU") %>%
  mutate(Label = ifelse(!is.na(Genus) & Genus != "", Genus, NA)) %>%
  filter(!is.na(Label)) %>%
  select(OTU, Label)

# Keep only OTUs with valid genus
ps_genus <- prune_taxa(tax_df_genus$OTU, ps_genus)
# Replace taxonomy table with genus-only labels
tax_table(ps_genus) <- tax_table(as.matrix(tax_df_genus %>% column_to_rownames("OTU")))

# --- 7. DESeq2 analysis ---
# Convert phyloseq object to DESeq2 format
dds <- phyloseq_to_deseq2(ps_genus, ~ Group_combined)
# Estimate size factors using a method robust to zeros
dds <- estimateSizeFactors(dds, type = "poscounts")
# Fit the DESeq2 model
dds <- DESeq(dds, fitType = "parametric")

# --- 8. Pairwise contrasts + LFC shrinkage ---
# Extract all group levels used in the model
groups <- levels(droplevels(sample_data(ps_genus)$Group_combined))
# Generate all pairwise group comparisons
comparisons <- combn(groups, 2, simplify = FALSE)
# Store results for each contrast
results_list <- list()

for (comp in comparisons) {
  # Define comparison groups
  a <- comp[1]; b <- comp[2]
  cname <- paste0(a, "_vs_", b)
  # Compute raw DESeq2 results
  res_raw <- tryCatch(results(dds, contrast = c("Group_combined", a, b)), error = function(e) NULL)
  if (is.null(res_raw)) next
  # Convert results to data frame
  res_raw_df <- as.data.frame(res_raw) %>% rownames_to_column("OTU")
  # Apply log2 fold-change shrinkage
  res_shrunk <- tryCatch(lfcShrink(dds, contrast = c("Group_combined", a, b), type = "normal"), error = function(e) NULL)
  # Merge shrunken LFCs if available
  if (!is.null(res_shrunk)) {
    res_shrunk_df <- as.data.frame(res_shrunk) %>% rownames_to_column("OTU") %>% select(OTU, log2FoldChange)
    res_combined <- res_raw_df %>%
      left_join(res_shrunk_df, by = "OTU", suffix = c("_raw", "")) %>%
      mutate(contrast = cname) %>%
      select(OTU, log2FoldChange, everything())
  } else {
    res_combined <- res_raw_df %>% mutate(log2FoldChange = log2FoldChange) %>% mutate(contrast = cname)
  }
  # Store results
  results_list[[cname]] <- res_combined
}
# Combine all contrasts into a single table
deseq_all <- bind_rows(results_list)

# --- 9. Annotate results with genus labels ---
deseq_all <- left_join(deseq_all, tax_df_genus, by = "OTU")

# --- 10. Filter significant genera  ---
# Use adjusted p-values if available; otherwise fall back to raw p-values
if(!"padj" %in% colnames(deseq_all)) {
  sig_deseq <- deseq_all %>% filter(!is.na(pvalue) & pvalue <= 0.05)
} else {
  sig_deseq <- deseq_all %>% filter(!is.na(padj) & padj <= 0.05)
}
# Remove entries lacking effect size or genus annotation
sig_deseq <- sig_deseq %>% filter(!is.na(log2FoldChange), !is.na(Label)) %>% mutate(Label = as.character(Label))
sig_deseq_plot <- sig_deseq

# Keep only male vs female contrasts within severity
severity_levels <- c("Control", "Mild", "Severe", "Deceased")

desired_contrasts <- c(
  paste0("male_", severity_levels, "_vs_female_", severity_levels),
  paste0("female_", severity_levels, "_vs_male_", severity_levels)
)

# Retain only contrasts present in the results
filtered_contrasts <- intersect(desired_contrasts, unique(sig_deseq_plot$contrast))
sig_deseq_plot <- sig_deseq_plot %>% filter(contrast %in% filtered_contrasts)
contrasts <- unique(sig_deseq_plot$contrast)


# --- 11a. Log2Fold Bar Plots ---
if (length(contrasts) == 0) {
  print("Skipping Bar Plots: No significant genera found for the Male vs Female same-severity comparisons.")
} else {
  for (c in contrasts) {
    df <- sig_deseq_plot %>% filter(contrast == c)
    
    # Extract groups and sex for dynamic labeling
    groups <- unlist(strsplit(c, "_vs_"))
    group1 <- groups[1]; group2 <- groups[2]
    sex1 <- strsplit(group1, "_")[[1]][1] 
    sex2 <- strsplit(group2, "_")[[1]][1]
    
    df <- df %>% mutate(Label = fct_reorder(Label, log2FoldChange),
                        direction = ifelse(log2FoldChange > 0, "Up", "Down"))
    
    df <- df %>% mutate(direction_label = case_when(
      direction == "Up" ~ paste("Higher in", sex1),
      direction == "Down" ~ paste("Higher in", sex2)
    ))
    # Customize colour
    color_mapping <- setNames(c("cornflowerblue", "plum2"),
                              c(paste("Higher in", sex1), paste("Higher in", sex2)))
    
    p <- ggplot(df, aes(x = Label, y = log2FoldChange, fill = direction_label)) +
      geom_col(show.legend = TRUE) +
      coord_flip() +
      theme_minimal() +
      labs(title = paste0("Significant Genera — ", c),
           x = "Genus",
           y = "log2 Fold Change",
           fill = "Abundance Difference") +
      scale_fill_manual(values = color_mapping) +
      theme(axis.text.y = element_text(size = 9),
            axis.text.x = element_text(size = 9),
            plot.title = element_text(hjust = 0.5))
    
    ggsave(filename = paste0("ISA_and_DESeq2_Plots/DESeq2_Genus_", c, ".png"),
           plot = p, width = 10, height = 6, dpi = 300)
  }
}




#### ---- Integrated Absolute abundance plot ---- ####

# --- 2. Get OTUs significant in ISA and DESeq2 ---
isa_otus <- significant_indicators$OTU 

deseq_sig_otus <- deseq_all %>%
  filter(!is.na(padj) & padj <= 0.05) %>%
  pull(OTU) %>%
  unique()

both_otus <- intersect(isa_otus, deseq_sig_otus)

# Prune phyloseq to only these significant OTUs
ps_both <- prune_taxa(both_otus, ps)

# --- 3. Clean metadata and define Group_combined ---
meta_df <- data.frame(sample_data(ps_both))

# Remove any samples with missing sex or Severity
meta_df <- meta_df %>% filter(!is.na(sex) & !is.na(Severity))

# Build Group_combined
meta_df <- meta_df %>% mutate(Group_combined = paste(sex, Severity, sep = "_"))

# Prune samples in phyloseq to match cleaned metadata
ps_both <- prune_samples(rownames(meta_df), ps_both)
sample_data(ps_both) <- sample_data(meta_df)

# --- 4. Clean taxonomy (Genus level) ---
tax_df <- as.data.frame(tax_table(ps_both)) %>%
  rownames_to_column("OTU") %>%
  # Keep only valid genera (remove unclassified)
  filter(!is.na(Genus) & Genus != "" & Genus != "g__") %>%
  mutate(Label = gsub("^g__", "", Genus))

# Keep only OTUs with valid genera
ps_both <- prune_taxa(tax_df$OTU, ps_both)

# --- 5. Melt phyloseq and join genus labels ---
ps_df <- psmelt(ps_both) %>%
  left_join(tax_df %>% select(OTU, Label), by = "OTU") %>%
  filter(!is.na(Label) & !is.na(Group_combined))

# --- 6. Summarize ABSOLUTE abundance per sample ---
plot_data_abs <- ps_df%>% # **** Using the filtered object ****
  group_by(Sample, Group_combined, Label) %>%
  summarize(sample_abundance = sum(Abundance), .groups = "drop") %>%
  ungroup() %>%
  # Summarize the mean ABSOLUTE count per group and genus
  group_by(Group_combined, Label) %>%
  summarize(mean_absolute_abundance = mean(sample_abundance), .groups = "drop") %>%
  # Add a pseudocount of 1 to handle zero values for the Log transformation
  mutate(log10_absolute_abundance = log10(mean_absolute_abundance + 1))

# --- 7. Set x-axis order ---
group_levels <- c(
  "male_Control", "female_Control",
  "male_Mild", "female_Mild",
  "male_Severe", "female_Severe",
  "male_Deceased", "female_Deceased"
)
plot_data_abs$Group_combined <- factor(plot_data_abs$Group_combined, levels = group_levels)

# --- 8. Color palette ---
unique_genera <- sort(unique(plot_data_abs$Label))
palette <- grDevices::colorRampPalette(c( "lightgreen", "#0072B8", "#009E73", "turquoise", "#56B4E8", "#3366CC", "lightblue"))(length(unique_genera))

# --- 9. Plot Log10 Stacked Bar Chart (Updated Title) ---
log_stacked_plot <- ggplot(plot_data_abs, aes(
  x = Group_combined,
  # Use the log10-transformed mean absolute abundance for the stack height
  y = log10_absolute_abundance,
  fill = Label
)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  labs(
    title = "Absolute Abundance of Significant (Log10 Stack)",
    x = "Sex–Severity Group",
    # Y-axis label reflects the transformation
    y = "Mean Absolute Abundance (Log10(Counts + 1))",
    fill = "Genus"
  ) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(face = "bold", size = 11),
    axis.title.y = element_text(face = "bold", size = 11),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9)
  )


print(log_stacked_plot)
ggsave("ISA_and_DESeq2_Plots/stacked_genus_log10_absolute_abundance.png", log_stacked_plot, width = 12, height = 9, dpi = 300)
