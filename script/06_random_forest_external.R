# ==============================================================================
# PTB Virome Project
# 06_random_forest_external.R
# External validation of random forest models:
#   1) call script/run_rf_external.sh to build models & predict on external cohort,
#   2) draw external ROC curves with AUC and 95% CI.
#
# Note: rf.py predict --fillna automatically fills missing features with 0,
#       so no feature selection is needed on the R side.
# ==============================================================================

rm(list = ls())
source("script/config.R")

library(pROC)

# ------------------------------------------------------------------------------
# 1. Run the external validation shell script
# ------------------------------------------------------------------------------
shell_script <- file.path("script", "run_rf_external.sh")
if (!file.exists(shell_script)) {
  stop("Shell script not found: ", shell_script)
}

message("Running external validation pipeline...")
status <- system2("bash",
                  args   = c(shQuote(shell_script), shQuote(normalizePath("."))),
                  stdout = "",
                  stderr = "")
if (status != 0) {
  stop("External validation pipeline failed with exit status ", status)
}

# ------------------------------------------------------------------------------
# 2. Read prediction results
# ------------------------------------------------------------------------------
files <- c(
  Bacteria = file.path(RF_EXTERNAL_RESULT_DIR, "predict_re_bac.txt"),
  Virus    = file.path(RF_EXTERNAL_RESULT_DIR, "predict_re_virus.txt"),
  Merged   = file.path(RF_EXTERNAL_RESULT_DIR, "predict_re_merge.txt")
)

data_list <- lapply(files, function(f) {
  read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
})

# ------------------------------------------------------------------------------
# 3. Compute ROC and AUC (with 95% CI)
# ------------------------------------------------------------------------------
roc_list    <- list()
auc_results <- data.frame(
  Model    = character(),
  AUC      = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  stringsAsFactors = FALSE
)

for (nm in names(data_list)) {
  df <- data_list[[nm]]
  
  # Health = 对照；Patient = 病例
  df$Group <- factor(df$Group, levels = c("Health", "Patient"))
  
  roc_obj <- pROC::roc(
    response  = df$Group,
    predictor = df$Patient,
    levels    = c("Health", "Patient"),
    direction = "<",
    percent   = TRUE,
    ci        = TRUE,
    quiet     = TRUE
  )
  
  roc_list[[nm]] <- roc_obj
  
  auc_value <- as.numeric(pROC::auc(roc_obj))
  ci_value  <- as.numeric(pROC::ci.auc(roc_obj))
  
  auc_results <- rbind(
    auc_results,
    data.frame(
      Model    = nm,
      AUC      = auc_value,
      CI_lower = ci_value[1],
      CI_upper = ci_value[3]
    )
  )
}

message("External validation AUC summary:")
print(auc_results)

write.table(auc_results,
            file.path(TABLE_DIR, "external_validation_auc.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Plot external ROC curves
# ------------------------------------------------------------------------------
pdf(file.path(FIGURE_DIR, "external_validation_ROC.pdf"),
    width = 5.5, height = 5.5)

curve_colors <- RF_EXTERNAL_COLORS[names(roc_list)]

plot(
  roc_list[[1]],
  percent      = TRUE,
  legacy.axes  = FALSE,
  col          = curve_colors[1],
  lwd          = 2,
  xlim         = c(100, 0),
  ylim         = c(0, 100),
  xlab         = "100 - Specificity (%)",
  ylab         = "Sensitivity (%)",
  main         = "External validation"
)

if (length(roc_list) > 1) {
  for (i in 2:length(roc_list)) {
    plot(
      roc_list[[i]],
      percent = TRUE,
      add     = TRUE,
      col     = curve_colors[i],
      lwd     = 2
    )
  }
}

abline(a = 0, b = 1, lty = 2, col = "grey70")

legend_text <- sprintf(
  "%s: AUC %.1f%% (95%% CI %.1f-%.1f%%)",
  auc_results$Model,
  auc_results$AUC,
  auc_results$CI_lower,
  auc_results$CI_upper
)

legend(
  "bottomright",
  legend = legend_text,
  col    = curve_colors,
  lwd    = 2,
  bty    = "n",
  cex    = 0.8
)

dev.off()

message("External validation ROC saved to ",
        file.path(FIGURE_DIR, "external_validation_ROC.pdf"))

