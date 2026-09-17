# ==============================================================================
# PTB Virome Project
# 06_random_forest.R
# Random forest: input preparation, model training (via run_rf.sh),
# visualization (AUC curve, ROC curve, top-20 feature importance),
# and formal paired comparison (bootstrap + DeLong).
# ==============================================================================

rm(list = ls())
source("script/config.R")

library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)

source("script/plot_roc.R")

# ------------------------------------------------------------------------------
# 1. Prepare random forest input files
# ------------------------------------------------------------------------------
sig_bac <- read.delim(file.path(MAASLIN_DIR, "bac", "significant_results.tsv"),
                      check.names = FALSE, stringsAsFactors = FALSE) %>%
  filter(metadata == "Group")

sig_vir <- read.delim(file.path(MAASLIN_DIR, "vir", "significant_results.tsv"),
                      check.names = FALSE, stringsAsFactors = FALSE) %>%
  filter(metadata == "Group")

vir_spe <- read.delim(VOTU_FILE,      check.names = FALSE, row.names = 1)
bac_spe <- read.delim(BACTERIAL_FILE, check.names = FALSE, row.names = 1)

sig_vir_spe <- vir_spe[rownames(vir_spe) %in% sig_vir$feature, , drop = FALSE]
sig_bac_spe <- bac_spe[rownames(bac_spe) %in% sig_bac$feature, , drop = FALSE]

common_samples <- intersect(colnames(sig_vir_spe), colnames(sig_bac_spe))
sig_vir_spe <- sig_vir_spe[, common_samples, drop = FALSE]
sig_bac_spe <- sig_bac_spe[, common_samples, drop = FALSE]
merged_profile <- rbind(sig_vir_spe, sig_bac_spe)

metadata <- read.delim(METADATA_FILE, check.names = FALSE, stringsAsFactors = FALSE)
sample_info <- metadata %>% select(Sample, Group)

dir.create(RF_DIR, recursive = TRUE, showWarnings = FALSE)
write.table(sig_vir_spe,    file.path(RF_DIR, "sig_vir_spe.profile"),
            quote = FALSE, sep = "\t")
write.table(sig_bac_spe,    file.path(RF_DIR, "sig_bac_spe.profile"),
            quote = FALSE, sep = "\t")
write.table(merged_profile, file.path(RF_DIR, "merged.profile"),
            quote = FALSE, sep = "\t")
write.table(sample_info,    file.path(RF_DIR, "sample_info.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# ------------------------------------------------------------------------------
# 2. Run the random forest pipeline (shell script)
# ------------------------------------------------------------------------------
shell_script <- file.path("script", "run_rf.sh")
if (!file.exists(shell_script)) {
  stop("Shell script not found: ", shell_script)
}

message("Running random forest pipeline (this may take a while)...")
status <- system2("bash",
                  args   = c(shQuote(shell_script), shQuote(normalizePath("."))),
                  stdout = "",
                  stderr = "")
if (status != 0) {
  stop("Random forest pipeline failed with exit status ", status)
}

# ==============================================================================
# 3. AUC-vs-number-of-signatures curve
# ==============================================================================
read_kf <- function(fname) {
  read.delim(file.path(RF_DIR, fname),
             check.names = FALSE, stringsAsFactors = FALSE)
}

bac_kf   <- read_kf("bac_auc_kf.txt")
vir_kf   <- read_kf("vir_auc_kf.txt")
merge_kf <- read_kf("merge_auc_kf.txt")

calculate_auc <- function(data, nspecies_values) {
  sapply(nspecies_values, function(n) {
    sub <- data[data$nspecies == n, , drop = FALSE]
    if (nrow(sub) > 0 && length(unique(sub$Actual)) >= 2) {
      roc_obj <- pROC::roc(response  = sub$Actual,
                           predictor = sub$Patient,
                           quiet     = TRUE,
                           direction = "<")
      as.numeric(pROC::auc(roc_obj))
    } else {
      NA_real_
    }
  })
}

auc_df <- bind_rows(
  data.frame(nspecies = unique(merge_kf$nspecies),
             AUC      = calculate_auc(merge_kf, unique(merge_kf$nspecies)),
             type     = "Viruses + Bacteria"),
  data.frame(nspecies = unique(bac_kf$nspecies),
             AUC      = calculate_auc(bac_kf, unique(bac_kf$nspecies)),
             type     = "Bacteria only"),
  data.frame(nspecies = unique(vir_kf$nspecies),
             AUC      = calculate_auc(vir_kf, unique(vir_kf$nspecies)),
             type     = "Viruses only")
) %>% na.omit()

max_auc_df <- auc_df %>%
  group_by(type) %>%
  slice_max(AUC, n = 1, with_ties = FALSE) %>%
  ungroup()

p_auc <- ggplot(auc_df, aes(x = nspecies, y = AUC,
                            color = type, group = type)) +
  geom_line() +
  geom_point(size = 1.5) +
  geom_hline(data = max_auc_df,
             aes(yintercept = AUC), linetype = "dashed",
             color = "grey50") +
  geom_text(data = max_auc_df,
            aes(label = sprintf("%.3f", AUC),
                x = nspecies, y = AUC + 0.02),
            vjust = 0, size = 3, show.legend = FALSE) +
  scale_color_manual(values = RF_COLORS) +
  labs(x = "Number of signatures", y = "AUC", color = "") +
  theme_bw()

ggsave(file.path(FIGURE_DIR, "rf_auc_curve.pdf"),
       p_auc, width = 6, height = 4)

# ==============================================================================
# 4. Mean ROC curves across seeds  →  同时也拿到 best_nspecies
# ==============================================================================
files_vec <- c(
  Bacteria = file.path(RF_DIR, "bac_auc_kf.txt"),
  Virus    = file.path(RF_DIR, "vir_auc_kf.txt"),
  Combine  = file.path(RF_DIR, "merge_auc_kf.txt")
)

res <- plot_roc_flexible(
  files_vec     = files_vec,
  output_pdf    = file.path(FIGURE_DIR, "rf_roc.pdf"),
  predictor_col = "Patient",
  color_vec     = unname(RF_COLORS[c("Bacteria only",
                                     "Viruses only",
                                     "Viruses + Bacteria")]),
  nspecies_mode = "max"
)

best_ns <- res$best_nspecies
print(best_ns)

# ==============================================================================
# 5. Formal comparison of Virus vs Bacteria
#    - Primary result: 10-seed mean AUC
#    - Primary test:   paired bootstrap on seed-averaged predictions
#    - Sensitivity:    paired DeLong on seed-averaged predictions
# ==============================================================================
data_list <- lapply(files_vec, function(f) {
  read.table(f,
             header = TRUE, sep = "\t",
             stringsAsFactors = FALSE,
             check.names = FALSE)
})

sample_col    <- "Sample"
predictor_col <- "Patient"
actual_col    <- "Actual"
seed_col      <- "seed"
nspecies_col  <- "nspecies"

seeds <- c(75, 6, 86, 54, 92, 78, 43, 63, 100, 96)

# ----------------------------------------------------------
# 5.1 主结果：每个 seed 单独算 AUC，再对 10 个 seed 取平均
# ----------------------------------------------------------
seed_auc <- function(nm) {
  
  df <- data_list[[nm]] %>%
    filter(.data[[nspecies_col]] == best_ns[[nm]])
  
  do.call(rbind, lapply(seeds, function(s) {
    
    d <- df[df[[seed_col]] == s, , drop = FALSE]
    
    if (nrow(d) < 2 ||
        length(unique(d[[actual_col]])) < 2) {
      return(NULL)
    }
    
    d[[actual_col]] <- factor(d[[actual_col]],
                              levels = c("Health", "Patient"))
    
    r <- roc(d[[actual_col]], d[[predictor_col]],
             levels = c("Health", "Patient"),
             direction = "<", quiet = TRUE)
    
    data.frame(model = nm, seed = s, AUC = as.numeric(auc(r)))
  }))
}

seed_auc_all <- do.call(rbind, lapply(names(files_vec), seed_auc))

mean_auc_table <- seed_auc_all %>%
  group_by(model) %>%
  summarise(
    Mean_AUC = mean(AUC),
    SD_AUC   = sd(AUC),
    Min_AUC  = min(AUC),
    Max_AUC  = max(AUC),
    N_seeds  = n(),
    .groups  = "drop"
  )

cat("\n========== 10-seed mean AUC (main result) ==========\n")
print(mean_auc_table, digits = 4)
cat("====================================================\n\n")

# ----------------------------------------------------------
# 5.2 构造 seed-averaged OOF predictions
#     每个样本跨 10 个 seed 平均预测概率
# ----------------------------------------------------------
avg_pair <- function(nm) {
  data_list[[nm]] %>%
    filter(.data[[nspecies_col]] == best_ns[[nm]]) %>%
    group_by(.data[[sample_col]], .data[[actual_col]]) %>%
    summarise(pred = mean(.data[[predictor_col]], na.rm = TRUE),
              .groups = "drop") %>%
    as.data.frame()
}

vir_avg <- avg_pair("Virus")
bac_avg <- avg_pair("Bacteria")

# 按 Virus 的样本顺序对齐 Bacteria（paired 的前提）
bac_avg <- bac_avg[match(vir_avg[[sample_col]], bac_avg[[sample_col]]), ]
stopifnot(identical(bac_avg[[sample_col]], vir_avg[[sample_col]]))

vir_avg[[actual_col]] <- factor(vir_avg[[actual_col]],
                                levels = c("Health", "Patient"))
bac_avg[[actual_col]] <- factor(bac_avg[[actual_col]],
                                levels = c("Health", "Patient"))

get_auc <- function(df) {
  r <- roc(df[[actual_col]], df$pred,
           levels = c("Health", "Patient"),
           direction = "<", quiet = TRUE)
  as.numeric(auc(r))
}

auc_vir_avg <- get_auc(vir_avg)
auc_bac_avg <- get_auc(bac_avg)
delta_obs   <- auc_vir_avg - auc_bac_avg

cat(sprintf("\nSeed-averaged AUC: Virus = %.4f, Bacteria = %.4f, delta = %.4f\n\n",
            auc_vir_avg, auc_bac_avg, delta_obs))

# ----------------------------------------------------------
# 5.3 主比较：paired bootstrap（两个模型共用 idx）
# ----------------------------------------------------------
set.seed(123)
B <- 2000
n <- nrow(vir_avg)
delta_boot <- numeric(B)

for (i in seq_len(B)) {
  
  idx <- sample.int(n, n, replace = TRUE)
  
  d_bac <- bac_avg[idx, , drop = FALSE]
  d_vir <- vir_avg[idx, , drop = FALSE]
  
  if (length(unique(d_vir[[actual_col]])) < 2) {
    delta_boot[i] <- NA_real_
    next
  }
  
  delta_boot[i] <- get_auc(d_vir) - get_auc(d_bac)
}

delta_boot <- delta_boot[is.finite(delta_boot)]

ci_boot <- quantile(delta_boot, c(0.025, 0.975), names = FALSE)
p_boot  <- 2 * min(mean(delta_boot <= 0), mean(delta_boot >= 0))
p_boot  <- min(p_boot, 1)

# ----------------------------------------------------------
# 5.4 敏感性分析：paired DeLong
# ----------------------------------------------------------
roc_vir_avg <- roc(vir_avg[[actual_col]], vir_avg$pred,
                   levels = c("Health", "Patient"),
                   direction = "<", quiet = TRUE)
roc_bac_avg <- roc(bac_avg[[actual_col]], bac_avg$pred,
                   levels = c("Health", "Patient"),
                   direction = "<", quiet = TRUE)

delong <- roc.test(roc_vir_avg, roc_bac_avg,
                   paired = TRUE, method = "delong")

# ----------------------------------------------------------
# 5.5 汇总主表
# ----------------------------------------------------------
main_result <- data.frame(
  Comparison = "Virus vs Bacteria",
  
  Mean_AUC_Virus    = mean_auc_table$Mean_AUC[mean_auc_table$model == "Virus"],
  Mean_AUC_Bacteria = mean_auc_table$Mean_AUC[mean_auc_table$model == "Bacteria"],
  
  AUC_Virus_avgPred    = auc_vir_avg,
  AUC_Bacteria_avgPred = auc_bac_avg,
  delta_AUC            = delta_obs,
  
  bootstrap_SE = sd(delta_boot),
  CI_lower     = ci_boot[1],
  CI_upper     = ci_boot[2],
  p_bootstrap  = p_boot,
  
  DeLong_Z = as.numeric(delong$statistic),
  p_DeLong = delong$p.value
)

cat("\n========== Main comparison result ==========\n")
print(main_result, digits = 4)
cat("============================================\n\n")

# ----------------------------------------------------------
# 5.6 输出
# ----------------------------------------------------------
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

write.table(seed_auc_all,
            file.path(TABLE_DIR, "rf_10seed_auc_per_seed.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(mean_auc_table,
            file.path(TABLE_DIR, "rf_10seed_mean_auc.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(main_result,
            file.path(TABLE_DIR, "rf_virus_vs_bacteria_main_comparison.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ==============================================================================
# 6. Top-20 feature importance barplots
# ==============================================================================
process_sig <- function(filepath) {
  read.delim(filepath, check.names = FALSE, stringsAsFactors = FALSE) %>%
    filter(metadata == "Group") %>%
    mutate(Group = case_when(
      coef > 0 ~ "Health",
      coef < 0 ~ "Patient",
      TRUE     ~ "Not significant"
    )) %>%
    select(feature, Group)
}

sig_vir2 <- process_sig(file.path(MAASLIN_DIR, "vir", "significant_results.tsv"))
sig_bac2 <- process_sig(file.path(MAASLIN_DIR, "bac", "significant_results.tsv"))
all_sig  <- bind_rows(sig_bac2, sig_vir2)

process_importance <- function(filepath, sig_df, default_type = NULL) {
  read.delim(filepath, check.names = FALSE, stringsAsFactors = FALSE) %>%
    arrange(desc(Mean_Importance)) %>%
    slice(1:RF_TOP_FEATURES) %>%
    rename(feature = Feature) %>%
    left_join(sig_df, by = "feature") %>%
    mutate(
      type = if (is.null(default_type)) {
        if_else(startsWith(as.character(feature), "v"), "Vir", "Bac")
      } else {
        default_type
      },
      y_axis_color = case_when(
        Group == "Patient" ~ "#bb5757",
        Group == "Health"  ~ "#2d8152",
        TRUE               ~ "black"
      )
    )
}

top_20_v   <- process_importance(file.path(RF_DIR, "vir"),   sig_vir2, "Vir")
top_20_b   <- process_importance(file.path(RF_DIR, "bac"),   sig_bac2, "Bac")
top_20_all <- process_importance(file.path(RF_DIR, "merge"), all_sig)

plot_top20 <- function(df) {
  df <- df %>% mutate(feature = forcats::fct_reorder(feature, Mean_Importance))
  ordered_features <- levels(df$feature)
  color_order <- df$y_axis_color[match(ordered_features, df$feature)]
  
  ggplot(df, aes(x = feature, y = Mean_Importance, fill = type)) +
    geom_col(width = 0.8) +
    labs(x = "", y = "Mean Decrease Accuracy") +
    scale_fill_manual(
      name   = "Microbial Type",
      values = c("Bac" = RF_COLORS[["Bacteria"]],
                 "Vir" = RF_COLORS[["Viruses"]]),
      labels = c("Bac" = "Bacteria", "Vir" = "Virus")
    ) +
    coord_flip() +
    theme_bw() +
    theme(
      axis.text.y      = element_text(size = 10, color = color_order),
      panel.background = element_blank(),
      panel.grid.major = element_line(color = "gray90")
    )
}

ggsave(file.path(FIGURE_DIR, "rf_top20_virus.pdf"),
       plot_top20(top_20_v),   width = 6, height = 5)
ggsave(file.path(FIGURE_DIR, "rf_top20_bacteria.pdf"),
       plot_top20(top_20_b),   width = 6, height = 5)
ggsave(file.path(FIGURE_DIR, "rf_top20_merged.pdf"),
       plot_top20(top_20_all), width = 6, height = 6)

message("Random forest analysis completed.")