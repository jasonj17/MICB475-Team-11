library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(dplyr)
library(picante)
library(ggpubr)
library(rstatix)

#### Load Data ####
covid_mex_metadata_fp <- "COVID_Mexico_Final_Updated/covid_mex_metadata_edited.tsv"
covid_mex_metadata <- read_delim(covid_mex_metadata_fp, delim="\t")

covid_mex_otu_fp <- "COVID_Mexico_Final_Updated/covid_mexico_final_updated_export/table_export/feature-table.txt"
covid_mex_otu <- read_delim(covid_mex_otu_fp, delim="\t", skip=1)

covid_mex_tax_fp <- "COVID_Mexico_Final_Updated/covid_mexico_final_updated_export/taxonomy_export/taxonomy.tsv"
covid_mex_taxonomy <- read_delim(covid_mex_tax_fp, delim="\t")

covid_mex_phylotree_fp <- "COVID_Mexico_Final_Updated/covid_mexico_final_updated_export/rooted_tree_export/tree.nwk"
covid_mex_phylotree <- read.tree(covid_mex_phylotree_fp)

#### Format OTU table ####

# OTU tables should be a matrix with rownames and colnames as OTUs and SampleIDs respectively
# save everything except first column (OTU ID) into a matrix
covid_mex_otu_matrix <- as.matrix(covid_mex_otu[, -1])

# Make first column (#OTU ID) the rownames of the new matrix
rownames(covid_mex_otu_matrix) <- covid_mex_otu$`#OTU ID`

# Use the "otu_table" function to make an OTU table
OTU <- otu_table(covid_mex_otu_matrix, taxa_are_rows = TRUE)
class(OTU)

#### Format sample metadata ####

# Save everything except sample-id as new data frame
covid_mex_metadata_df <- as.data.frame(covid_mex_metadata[, -1])

# Make sample-id the rownames
rownames(covid_mex_metadata_df) <- covid_mex_metadata$'sample-id'

# Changing "Group" column name to "severity" and standardizing "Group" column into 5 categories
# which are: Asymptomatic, Hospitalized Negative, Hospitalized Positive, and Deceased
covid_mex_metadata_edited <- covid_mex_metadata_df %>%
  mutate(Group = case_when(
    grepl("asymptomatic", Group, ignore.case = TRUE) ~ "Asymptomatic",
    grepl("ambulatory negative", Group, ignore.case = TRUE) ~ "Mild negative",
    grepl("ambulatory positive", Group, ignore.case = TRUE) ~ "Mild positive",
    grepl("hospitalized positive", Group, ignore.case = TRUE) &
      !grepl("deceased", Group, ignore.case = TRUE) ~ "Severe positive",
    grepl("deceased", Group, ignore.case = TRUE) ~ "Deceased",
    TRUE ~ NA_character_     # <-- make unmatched truly missing
  )) %>%
  rename(severity = Group)  # rename "Group" to "severity"

# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(covid_mex_metadata_edited)
class(SAMP)

#### Formatting taxonomy ####

# Convert taxon strings to a table with separate taxa rank columns
covid_mex_taxonomy_matrix <- covid_mex_taxonomy %>% select(-Confidence)%>%
  separate(col=Taxon, sep=";"
           , into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  as.matrix() # Saving as a matrix

# Save everything except feature IDs
covid_mex_taxonomy_matrix <- covid_mex_taxonomy_matrix[, -1]

# Make sample-ids the rownames
rownames(covid_mex_taxonomy_matrix) <- covid_mex_taxonomy$`Feature ID`

# Make taxa table
TAX <- tax_table(covid_mex_taxonomy_matrix)
class(TAX)

#### Create phyloseq object ####

# Merge all into a phyloseq object
covid_mexico_phyloseq <- phyloseq(OTU, SAMP, TAX, covid_mex_phylotree)

#### Looking at phyloseq object ####

# View components of phyloseq object with the following commands
otu_table(covid_mexico_phyloseq)
sample_data(covid_mexico_phyloseq)
tax_table(covid_mexico_phyloseq)
phy_tree(covid_mexico_phyloseq)

# View sample names in our data set
sample_names(covid_mexico_phyloseq)

# Number of samples in our data set
nsamples(covid_mexico_phyloseq)

# Sum of reads in each sample
sample_sums(covid_mexico_phyloseq)

# 3 samples with the most reads from our data set
samps_most_reads <- names(sort(sample_sums(covid_mexico_phyloseq), decreasing = TRUE)[1:3])
get_taxa(covid_mexico_phyloseq, samps_most_reads)

# Names of taxa from our data set
taxa_names(covid_mexico_phyloseq)

# Number of taxa we have from our data set
ntaxa(covid_mexico_phyloseq)

# 3 most abundant taxa from our data set
taxa_most_abundant <- names(sort(taxa_sums(covid_mexico_phyloseq), decreasing = TRUE)[1:3])
get_sample(covid_mexico_phyloseq, taxa_most_abundant)

#### Filtering phyloseq object ####

# Remove any non-bacterial sequences, if any
covid_mexico_phyloseq_filt <- subset_taxa(
  covid_mexico_phyloseq,
  !(Family %in% c("Mitochondria") |
      Order %in% c("Chloroplast")))

# Remove ASVs that have less than 5 counts total
covid_mexico_phyloseq_filt_nolow <- filter_taxa(covid_mexico_phyloseq_filt,
                                                function(x) sum(x)>5,
                                                prune = TRUE)

# Remove samples with less than 100 reads
covid_mexico_phyloseq_filt_nolow_samps <- prune_samples(sample_sums(covid_mexico_phyloseq_filt_nolow)>100,
                                                        covid_mexico_phyloseq_filt_nolow)

# Remove samples where sex is NA, severity is NA, and filter out "Mild negative" group in severity column
covid_mexico_phyloseq_final <- subset_samples(
  covid_mexico_phyloseq_filt_nolow_samps,
  !is.na(sex) & !is.na(severity) & severity != "Mild negative")

#### Generating rarefaction curve ####
# rngseed sets a random number. To be able to reproduce this exact analysis each time, need to set rngseed the same number each time
# t transposes the table to use rarecurve function
# cex decreases font size
rarecurve(t(as.data.frame(otu_table(covid_mexico_phyloseq_final))), cex=0.1)

#### Rarefying phyloseq object ####
covid_mexico_phyloseq_rare <- rarefy_even_depth(covid_mexico_phyloseq_final, rngseed = 11, sample.size = 48792)

#### Saving final and rarefied phyloseq objects ####
save(covid_mexico_phyloseq_final, file="covid_mexico_phyloseq_final.RData")
save(covid_mexico_phyloseq_rare, file="covid_mexico_phyloseq_rare.RData")

#### Alpha diversity ####

### Shannon Diversity ###
## 1) Build the base Shannon diversity plot
base_alpha <- plot_richness(covid_mexico_phyloseq_rare, x = "sex", measures = "Shannon") +
  xlab("Sex") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~severity) +
  theme_minimal()

base_alpha$data <- subset(base_alpha$data, !is.na(severity) & severity != "NA" & !is.na(sex))

## 2) Compute Kruskal–Wallis p-values within each severity (Shannon ~ sex)
alpha_df <- estimate_richness(covid_mexico_phyloseq_rare, measures = "Shannon") %>%
  tibble::rownames_to_column("SampleID")

meta_df <- as.data.frame(sample_data(covid_mexico_phyloseq_rare)) %>%
  tibble::rownames_to_column("SampleID")

df <- left_join(alpha_df, meta_df, by = "SampleID") %>%
  mutate(
    severity = as.factor(severity),
    sex = as.factor(sex)
  ) %>%
  filter(!is.na(severity), severity != "NA", !is.na(sex))

kw_df <- base_alpha$data %>%
  dplyr::filter(!is.na(severity), severity != "NA", !is.na(sex)) %>%
  dplyr::group_by(severity) %>%
  dplyr::summarise(
    p = if (dplyr::n_distinct(sex) < 2) NA_real_
    else stats::kruskal.test(value ~ sex, data = dplyr::cur_data())$p.value,
    y = max(value, na.rm = TRUE) * 1.06,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    label = ifelse(is.na(p), "KW p = NA", paste0("KW p = ", signif(p, 3))),
    sex = "male"
  )

## 3) Add facet-wise labels + significance brackets

# Compute max y per severity for bracket placement
gg_richness_severity <- base_alpha +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("female", "male")),
    label = "p.signif",
    hide.ns = FALSE,
    bracket.size = 0.6,
    tip.length = 0.05,
    vjust = 1.2
  ) +
  geom_text(
    data = kw_df,
    aes(x = sex, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.15,
    vjust = 0,
    size = 3
  ) +
  coord_cartesian(clip = "off")

gg_richness_severity_final <- gg_richness_severity +
  theme_classic(base_size = 11) +
  theme(
    # Panel appearance
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.spacing = unit(1.2, "lines"),
    
    # Facet strip (severity labels)
    strip.background = element_blank(),
    strip.text = element_text(
      size = 11,
      face = "bold",
      margin = margin(b = 6)
    ),
    
    # Axis text & titles
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 8)),
    axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 8)),
    axis.text = element_text(size = 10),
    
    # Remove legend
    legend.position = "none",
    
    # Clean look
    axis.line = element_line(color = "black"),
    plot.margin = margin(10, 12, 10, 10)
  )

gg_richness_severity_final

ggsave(filename = "Alpha_and_Beta_Diversity_Plots/plot_richness_severity.png",
       gg_richness_severity_final, height=6, width=10)

estimate_richness(covid_mexico_phyloseq_rare)

## Phylogenetic diversity ##
# Calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(
  t(otu_table(covid_mexico_phyloseq_rare)),
  phy_tree(covid_mexico_phyloseq_rare),
  include.root = FALSE
)

# Add PD to metadata table
sample_data(covid_mexico_phyloseq_rare)$PD <- phylo_dist$PD

# Convert sample_data to a data.frame for ggplot + stats
# Force sample_data to a plain data.frame first (avoids <sample_data> pipe issues)
df_pd <- data.frame(sample_data(covid_mexico_phyloseq_rare), check.names = FALSE) %>%
  tibble::rownames_to_column("SampleID") %>%
  dplyr::mutate(
    severity = as.factor(severity),
    sex = as.factor(sex)
  ) %>%
  dplyr::filter(!is.na(severity), severity != "NA", !is.na(sex), !is.na(PD))

# 1) Kruskal–Wallis per severity (PD ~ sex)
kw_brackets <- df_pd %>%
  group_by(severity) %>%
  kruskal_test(PD ~ sex) %>%
  add_significance("p") %>%
  # Define the groups for the brackets
  mutate(
    group1 = "female",
    group2 = "male",
    # Create the label for the text: "KW p = 0.XXX"
    kw_label = paste0("KW p = ", format.pval(p, digits = 3, eps = 0.001))
  )

# 2) Calculate position for brackets and text per facet
pos_df <- df_pd %>%
  group_by(severity) %>%
  summarise(
    # INCREASED spacing from max PD value:
    # Set y.position for the bracket/ns label
    y_ns_signif = max(PD, na.rm = TRUE) + 0.15 * IQR(PD, na.rm = TRUE),
    # Set y.position for the KW p-value text (slightly higher to prevent overlap)
    y_kw_text   = max(PD, na.rm = TRUE) + 0.22 * IQR(PD, na.rm = TRUE),
    .groups = "drop"
  )

# Join the positions to the brackets data frame
kw_brackets <- kw_brackets %>%
  left_join(pos_df, by = "severity")

# 3) Create a data frame for the Kruskal-Wallis p-value text label
# This will be positioned specifically for the KW p-value.
kw_text_df <- kw_brackets %>%
  mutate(
    # Use the separate higher position
    y.position = y_kw_text,
    # Pin the text position to the 'male' boxplot column (position x=2)
    xmin = "male",
    xmax = "male",
    # Use the combined label for the text
    label = kw_label
  )

# Generate the final plot
plot_pd_severity_final <- ggplot(df_pd, aes(x = sex, y = PD)) +
  geom_jitter(
    width = 0.15,          # Maximum horizontal jitter. Keep it narrow.
    shape = 21,            # Circle shape (21) is often preferred for visibility
    size = 2,              # Point size
    fill = "gray70",       # Fill color for the points
    color = "black",       # Border color for the points
    alpha = 0.7            # Slightly transparent
  ) +
  geom_boxplot(width = 0.45, outlier.shape = NA, fill = NA, color = "black") +
  facet_wrap(~severity, scales = "free_y") +
  
  # Add significance brackets (using p.signif like ns, *, **)
  stat_pvalue_manual(
    kw_brackets,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y_ns_signif", # Use the lower position for the bracket/ns
    tip.length = 0.03,
    bracket.size = 0.5
  ) +
  
  # Add Kruskal-Wallis p-value text (KW p = 0.XXX)
  # Using stat_pvalue_manual for precise placement
  stat_pvalue_manual(
    kw_text_df,
    label = "label",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    tip.length = 0,
    bracket.size = 0,
    size = 3,
    hjust = 1.05
  ) +
  
  # Tweak scales for better appearance and spacing
  scale_x_discrete(expand = expansion(mult = c(0.6, 0.6))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  coord_cartesian(clip = "off") + 
  
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.spacing = unit(0.4, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold", margin = margin(b = 6)),
    legend.position = "none"
  )

plot_pd_severity_final

ggsave(filename = "Alpha_and_Beta_Diversity_Plots/plot_pd_severity.png",
       plot_pd_severity_final, height = 6, width = 10)

#### Beta diversity #####
# 1) Distance matrix
dm_unifrac <- UniFrac(covid_mexico_phyloseq_rare, weighted = TRUE)

# 2) Metadata table with explicit SampleID so we can subset the distance correctly
samp_dat_wdiv <- data.frame(sample_data(covid_mexico_phyloseq_rare),
                            estimate_richness(covid_mexico_phyloseq_rare)) %>%
  tibble::rownames_to_column("SampleID") %>%
  mutate(
    severity = as.factor(severity),
    sex      = as.factor(sex)
  )

# Helper: subset a 'dist' object to just a set of SampleIDs
subset_dist <- function(d, ids) {
  labs <- attr(d, "Labels")
  keep <- labs %in% ids
  as.dist(as.matrix(d)[keep, keep, drop = FALSE])
}

# 3) PERMANOVA: sex effect within each severity
permanova_by_severity <- samp_dat_wdiv %>%
  filter(!is.na(severity), !is.na(sex)) %>%
  group_by(severity) %>%
  group_modify(~{
    dat <- .x
    
    # Need at least 2 sex levels and enough samples
    if (n_distinct(dat$sex) < 2 || nrow(dat) < 4) {
      return(tibble(R2 = NA_real_, p = NA_real_, n = nrow(dat)))
    }
    
    d_sub <- subset_dist(dm_unifrac, dat$SampleID)
    
    # Ensure sample order matches dist labels
    dat2 <- dat %>%
      filter(SampleID %in% attr(d_sub, "Labels")) %>%
      arrange(match(SampleID, attr(d_sub, "Labels"))) %>%
      as.data.frame()
    
    a <- adonis2(d_sub ~ sex, data = dat2, permutations = 999)
    
    tibble(
      R2 = unname(a$R2[1]),
      p  = unname(a$`Pr(>F)`[1]),
      n  = nrow(dat2)
    )
  }) %>%
  ungroup()

# 4) Labels to add to each facet (one label per severity)
label_df <- permanova_by_severity %>%
  mutate(
    label = ifelse(
      is.na(p),
      "PERMANOVA (Weighted UniFrac)\nSex: NA",
      sprintf("PERMANOVA (Weighted UniFrac)\nSex: R² = %.3f, p = %.3g", R2, p)
    ),
    x = Inf, y = Inf
  ) %>%
  select(severity, x, y, label)

pcoa_wu <- ordinate(covid_mexico_phyloseq_rare, method="PCoA", distance=dm_unifrac)

base_beta <- plot_ordination(covid_mexico_phyloseq_rare, pcoa_wu, color="sex", shape="sex") +
  labs(shape = "Sex", colour = "Sex") +
  theme_minimal()

gg_pcoa <- base_beta +
  stat_ellipse(
    data = base_beta$data,
    aes(x = Axis.1, y = Axis.2, group = sex, color = sex),
    type = "norm",
    level = 0.68,
    linewidth = 0.7,
    alpha = 0.6,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ severity) +
  coord_cartesian(clip = "off")

gg_pcoa_final <- gg_pcoa +
  # Sex ellipse fill (within each severity facet)
  stat_ellipse(
    aes(group = sex, fill = sex),
    geom = "polygon",
    type = "norm",
    level = 0.68,
    linewidth = 0.0,
    alpha = 0.12,
    show.legend = FALSE
  ) +
  # Sex ellipse outline (within each severity facet)
  stat_ellipse(
    aes(group = sex, color = sex),
    type = "norm",
    level = 0.68,
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  # Make points slightly larger + readable
  geom_point(size = 2.4, alpha = 0.9) +
  # PERMANOVA labels per facet
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.15,
    size = 3.5
  ) +
  # Clean axes/labels
  labs(
    x = "PCoA1 [96.4%]",
    y = "PCoA2 [3.5%]",
    color = "Sex",
    fill  = "Sex",
    shape = "Sex"
  ) +
  # Make legends compact and non-distracting
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 3, alpha = 1)),
    shape = guide_legend(order = 2, override.aes = list(size = 3, alpha = 1)),
    fill  = "none"
  ) +
  # Change theme
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(3, "pt"),
    plot.margin = margin(10, 12, 10, 10),
    plot.title.position = "plot"
  ) +
  coord_cartesian(clip = "off")

gg_pcoa_final

ggsave(filename = "Alpha_and_Beta_Diversity_Plots/plot_pcoa.png",
       gg_pcoa_final, height=7, width=11)
