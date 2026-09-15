# Shared configuration for the PTB virome reproducibility package.
# The public metadata.tsv already uses Group = Health / Patient.

PROJECT_DIR <- "."
DATA_DIR <- file.path(PROJECT_DIR, "data")
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
FIGURE_DIR <- file.path(RESULTS_DIR, "figures")
TABLE_DIR <- file.path(RESULTS_DIR, "tables")
MAASLIN_DIR <- file.path(RESULTS_DIR, "maaslin")
ANNOTATION_DIR <- file.path(DATA_DIR, "annotations")

METADATA_FILE <- file.path(DATA_DIR, "metadata.tsv")
VOTU_FILE <- file.path(DATA_DIR, "votu_abundance.tsv")
VIRUS_FAMILY_FILE <- file.path(DATA_DIR, "virus_family_abundance.tsv")
BACTERIAL_FILE <- file.path(DATA_DIR, "bacterial_species_abundance.tsv")
KO_ABUNDANCE_FILE <- file.path(DATA_DIR, "votu_ko_abundance.tsv")

VOTU_KO_MAPPING_FILE <- file.path(ANNOTATION_DIR, "votu_ko_mapping.tsv")
VOTU_FAMILY_ANNOTATION <- file.path(ANNOTATION_DIR, "vOTU.family.tsv")
VOTU_HOST_ANNOTATION <- file.path(ANNOTATION_DIR, "vOTU.host.tax.genus.tsv")
HC_NODE_FILE <- file.path(ANNOTATION_DIR, "network", "hc_node.csv")
TB_NODE_FILE <- file.path(ANNOTATION_DIR, "network", "tb_node.csv")

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
