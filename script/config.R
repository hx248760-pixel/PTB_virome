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
KO_DESCRIPTION_FILE <- file.path(ANNOTATION_DIR,"ko_description.tsv")
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

# ------------------------------------------------------------------------------
# Random Forest Configuration
# ------------------------------------------------------------------------------
RF_DIR        <- file.path(RESULTS_DIR, "random_forest")
RF_SCRIPT     <- file.path("script", "rf.py")
RF_SHELL      <- file.path("script", "run_rf.sh")
RF_PYTHON     <- "python3"
RF_THREADS    <- 10
RF_N_EST      <- 1000
RF_K          <- 5

# 特征数网格
RF_TOP_N_VIR   <- c(seq(2, 100, 5), seq(101, 806, 50))
RF_TOP_N_BAC   <- seq(2, 81, 5)
RF_TOP_N_MERGE <- c(seq(2, 100, 5), seq(101, 887, 50))
RF_SEEDS       <- 1:10

RF_TOP_FEATURES <- 20
RF_COLORS <- c(
  "Bacteria only"           = "#869bc3",
  "Viruses only"            = "#9fcc56",
  "Viruses + Bacteria" = "#bc6363"
)

# ------------------------------------------------------------------------------
# External Validation Configuration
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# External Validation Configuration
# ------------------------------------------------------------------------------
# rf.py predict --fillna 会自动把训练模型中存在、外部数据缺失的特征补 0，
# 因此外部输入文件无需在 R 端做特征筛选，直接用原始丰度表即可。
RF_EXTERNAL_RESULT_DIR <- file.path(RF_DIR, "external")
RF_EXTERNAL_DATA_DIR   <- file.path(DATA_DIR, "external")

# 外部数据输入文件
RF_EXTERNAL_VOTU_FILE    <- file.path(RF_EXTERNAL_DATA_DIR, "votu.profile")
RF_EXTERNAL_BAC_FILE     <- file.path(RF_EXTERNAL_DATA_DIR, "species_relative_abundance.tsv")
RF_EXTERNAL_SAMPLE_INFO  <- file.path(RF_EXTERNAL_DATA_DIR, "sample_info.txt")

# 外部 ROC 曲线颜色
RF_EXTERNAL_COLORS <- c(
  Bacteria = "#869bc3",
  Virus    = "#9fcc56",
  Merged   = "#bc6363"
)