# config.R

project_dir <- normalizePath(".")

data_dir <- file.path(project_dir, "data")
script_dir <- file.path(project_dir, "scripts")
function_dir <- file.path(project_dir, "functions")
external_dir <- file.path(project_dir, "external")

result_dir <- file.path(project_dir, "results")
figure_dir <- file.path(result_dir, "figures")
table_dir <- file.path(result_dir, "tables")

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

metadata_file <- file.path(data_dir, "metadata.tsv")
votu_file <- file.path(data_dir, "votu_abundance.tsv")
virus_family_file <- file.path(data_dir, "virus_family_abundance.tsv")
bacterial_species_file <- file.path(data_dir, "bacterial_species_abundance.tsv")

source(file.path(function_dir, "calculate_alpha.R"))
source(file.path(function_dir, "plot_aplha.R"))
source(file.path(function_dir, "plot_beta.R"))
