# PTB virome project
# config.R
# Central configuration for all analysis scripts.
# Run scripts from the repository root directory.

# -----------------------------
# Project directories
# -----------------------------
PROJECT_DIR <- "."
DATA_DIR <- file.path(PROJECT_DIR, "data")
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
FIGURE_DIR <- file.path(RESULTS_DIR, "figures")
TABLE_DIR <- file.path(RESULTS_DIR, "tables")
MAASLIN_DIR <- file.path(RESULTS_DIR, "maaslin")
ANNOTATION_DIR <- file.path(DATA_DIR, "annotations")

# -----------------------------
# Input data
# -----------------------------
METADATA_FILE <- file.path(DATA_DIR, "metadata.tsv")
VOTU_FILE <- file.path(DATA_DIR, "votu_abundance.tsv")
VIRUS_FAMILY_FILE <- file.path(DATA_DIR, "virus_family_abundance.tsv")
BACTERIAL_FILE <- file.path(DATA_DIR, "bacterial_species_abundance.tsv")
KO_ABUNDANCE_FILE <- file.path(DATA_DIR, "votu_ko_abundance.tsv")
VOTU_KO_MAPPING_FILE <- file.path(ANNOTATION_DIR, "votu_ko_mapping.tsv")

# Optional annotation/network inputs supplied separately.
VOTU_FAMILY_ANNOTATION <- file.path(ANNOTATION_DIR, "vOTU.family.tsv")
VOTU_HOST_ANNOTATION <- file.path(ANNOTATION_DIR, "vOTU.host.tax.genus.tsv")
HC_NODE_FILE <- file.path(ANNOTATION_DIR, "network", "hc_node.csv")
TB_NODE_FILE <- file.path(ANNOTATION_DIR, "network", "tb_node.csv")

# -----------------------------
# Common analysis parameters
# -----------------------------
GROUP_LEVELS <- c("Health", "Patient")
GROUP_COLORS <- c("Patient" = "#D44D5C", "Health" = "#33a02c")

MAASLIN_MIN_ABUNDANCE <- 0.01
MAASLIN_MIN_PREVALENCE <- 0.01
MAASLIN_NORMALIZATION <- "NONE"
MAASLIN_FIXED_EFFECTS <- c("Group", "Sex", "Age")
MAASLIN_REFERENCE <- c("Group,Patient")
MAASLIN_Q_THRESHOLD <- 0.25

OCCURRENCE_THRESHOLD <- 0.10
CORRELATION_THRESHOLD <- 0.60
FDR_THRESHOLD <- 0.05

BETADISPER_PERMUTATIONS <- 999

# -----------------------------
# Helper functions/scripts
# -----------------------------
# These files will be added to script/ when provided by the user:
#   calculate_alpha.R
#   plot_aplha.R
#   plot_beta.R
#   plot_family_barplot.R
#   plot_roc.R
#   plot_virf_host.R
#   calc_corr.py
#   rf.py

# Ensure output directories exist.
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MAASLIN_DIR, showWarnings = FALSE, recursive = TRUE)
