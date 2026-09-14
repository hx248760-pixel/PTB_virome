#' Plot ROC curves using mean AUC across independent random seeds
#'
#' AUC is calculated separately for each seed and feature number,
#' followed by averaging AUCs across seeds.
#'
#' @param files_vec
#'   Named character vector. Names are displayed as model labels,
#'   values are corresponding prediction files.
#'
#' @param output_pdf
#'   Output PDF file.
#'
#' @param predictor_col
#'   Prediction probability column. Default = "MHO".
#'
#' @param actual_col
#'   True class label column. Default = "Actual".
#'
#' @param seed_col
#'   Random seed column. Default = "seed".
#'
#' @param nspecies_col
#'   Number of selected features. Default = "nspecies".
#'
#' @param color_vec
#'   Curve colors.
#'
#' @param nspecies_mode
#'   "manual" or "max".
#'   "manual" is recommended for the final figure.
#'
#' @param nspecies_selected
#'   Selected feature number for each model.
#'
#' @param show_auc_range
#'   Whether to display 95% CI based on seed-specific AUCs.
#'
#' @param auc_text_position
#'   Position of AUC labels.
#'
#' @param n_grid
#'   Number of points used to interpolate ROC curves.
#'
#' @return
#'   A list containing AUC summary, selected feature numbers,
#'   and ROC data.
#'
#' @export
plot_roc_flexible <- function(
    files_vec,
    output_pdf,
    predictor_col = "MHO",
    actual_col = "Actual",
    seed_col = "seed",
    nspecies_col = "nspecies",
    color_vec = NULL,
    nspecies_mode = c("manual", "max"),
    nspecies_selected = NULL,
    show_auc_range = TRUE,
    auc_text_position = "topright",
    n_grid = 101,
    ...) {
  
  nspecies_mode <- match.arg(nspecies_mode)
  
  n_curves <- length(files_vec)
  model_names <- names(files_vec)
  
  # ------------------------------------------------------------
  # 1. Basic checks
  # ------------------------------------------------------------
  
  if (is.null(model_names) || any(model_names == "")) {
    stop(
      "files_vec must be a named vector, e.g. ",
      "c(Bacteria='bac.txt', Virus='vir.txt')"
    )
  }
  
  if (is.null(color_vec)) {
    
    default_palette <- c(
      "#49B192",
      "#D95F02",
      "#417CB4",
      "#83639F",
      "#E7298A"
    )
    
    if (n_curves > length(default_palette)) {
      stop("Not enough default colors for the number of curves.")
    }
    
    color_vec <- default_palette[seq_len(n_curves)]
    
  } else if (length(color_vec) < n_curves) {
    
    stop(
      "The number of colors in color_vec is smaller than ",
      "the number of input files."
    )
  }
  
  if (nspecies_mode == "manual") {
    
    if (is.null(nspecies_selected) ||
        length(nspecies_selected) != n_curves) {
      
      stop(
        "When nspecies_mode='manual', ",
        "nspecies_selected must have the same length as files_vec."
      )
    }
  }
  
  # ------------------------------------------------------------
  # 2. Read prediction files
  # ------------------------------------------------------------
  
  data_list <- list()
  
  for (nm in model_names) {
    
    file_path <- files_vec[[nm]]
    
    if (!file.exists(file_path)) {
      stop("File not found: ", file_path)
    }
    
    df <- read.table(
      file_path,
      header = TRUE,
      sep = "\t",
      fill = TRUE,
      stringsAsFactors = FALSE
    )
    
    required_cols <- c(
      actual_col,
      predictor_col,
      seed_col,
      nspecies_col
    )
    
    missing_cols <- setdiff(
      required_cols,
      colnames(df)
    )
    
    if (length(missing_cols) > 0) {
      
      stop(
        "Dataset [", nm, "] is missing: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    
    data_list[[nm]] <- df
  }
  
  # ------------------------------------------------------------
  # 3. Calculate seed-specific AUC
  #
  # IMPORTANT:
  # One seed = one AUC
  #
  # For each seed:
  #   77 OOF predictions -> 1 AUC
  #
  # Then:
  #   10 seed-specific AUCs -> mean AUC
  # ------------------------------------------------------------
  
  auc_summary_list <- list()
  seed_auc_list <- list()
  
  for (nm in model_names) {
    
    df <- data_list[[nm]]
    
    # Unique combinations of seed and feature number
    combinations <- unique(
      df[, c(seed_col, nspecies_col)]
    )
    
    seed_results <- list()
    
    for (i in seq_len(nrow(combinations))) {
      
      current_seed <- combinations[[seed_col]][i]
      current_n <- combinations[[nspecies_col]][i]
      
      sub_df <- df[
        df[[seed_col]] == current_seed &
          df[[nspecies_col]] == current_n,
        ,
        drop = FALSE
      ]
      
      # Need at least two classes
      if (
        nrow(sub_df) < 2 ||
        length(unique(sub_df[[actual_col]])) < 2
      ) {
        next
      }
      
      roc_obj <- pROC::roc(
        response = sub_df[[actual_col]],
        predictor = sub_df[[predictor_col]],
        percent = TRUE,
        quiet = TRUE,
        direction = "<"
      )
      
      auc_value <- as.numeric(
        pROC::auc(roc_obj)
      )
      
      seed_results[[length(seed_results) + 1]] <- data.frame(
        seed = current_seed,
        nspecies = current_n,
        auc = auc_value
      )
    }
    
    seed_auc <- do.call(
      rbind,
      seed_results
    )
    
    seed_auc_list[[nm]] <- seed_auc
    
    # ----------------------------------------------------------
    # Mean AUC across seeds
    # ----------------------------------------------------------
    
    auc_summary <- seed_auc %>%
      dplyr::group_by(nspecies) %>%
      dplyr::summarise(
        mean_auc = mean(auc, na.rm = TRUE),
        sd_auc = sd(auc, na.rm = TRUE),
        min_auc = min(auc, na.rm = TRUE),
        max_auc = max(auc, na.rm = TRUE),
        n_seeds = dplyr::n(),
        .groups = "drop"
      )
    
    auc_summary_list[[nm]] <- auc_summary
  }
  
  # ------------------------------------------------------------
  # 4. Select feature number
  # ------------------------------------------------------------
  
  if (nspecies_mode == "manual") {
    
    plot_ns <- nspecies_selected
    names(plot_ns) <- model_names
    
  } else {
    
    plot_ns <- sapply(
      auc_summary_list,
      function(df) {
        
        df$nspecies[
          which.max(df$mean_auc)
        ]
      }
    )
  }
  
  # ------------------------------------------------------------
  # 5. Generate mean ROC curve
  #
  # IMPORTANT:
  #
  # We DO NOT pool predictions from all seeds.
  #
  # Instead:
  #
  # seed 1 -> ROC
  # seed 2 -> ROC
  # ...
  # seed 10 -> ROC
  #
  # Then interpolate each ROC to a common FPR grid and
  # calculate the mean TPR.
  #
  # Therefore the displayed ROC represents the average
  # cross-validation ROC performance across seeds.
  # ------------------------------------------------------------
  
  roc_data_list <- list()
  
  fpr_grid <- seq(
    0,
    1,
    length.out = n_grid
  )
  
  for (nm in model_names) {
    
    df <- data_list[[nm]]
    
    ns <- plot_ns[[nm]]
    
    sub_all <- df[
      df[[nspecies_col]] == ns,
      ,
      drop = FALSE
    ]
    
    seed_values <- unique(
      sub_all[[seed_col]]
    )
    
    interpolated_rocs <- list()
    
    for (current_seed in seed_values) {
      
      sub_df <- sub_all[
        sub_all[[seed_col]] == current_seed,
        ,
        drop = FALSE
      ]
      
      if (
        nrow(sub_df) < 2 ||
        length(unique(sub_df[[actual_col]])) < 2
      ) {
        next
      }
      
      roc_obj <- pROC::roc(
        response = sub_df[[actual_col]],
        predictor = sub_df[[predictor_col]],
        percent = FALSE,
        quiet = TRUE,
        direction = "<"
      )
      
      # pROC ROC coordinates
      roc_coords <- pROC::coords(
        roc_obj,
        x = "all",
        ret = c("specificity", "sensitivity"),
        transpose = FALSE
      )
      
      specificity <- roc_coords$specificity
      sensitivity <- roc_coords$sensitivity
      
      # Convert to FPR
      fpr <- 1 - specificity
      tpr <- sensitivity
      
      # Remove duplicated FPR values
      tmp <- data.frame(
        fpr = fpr,
        tpr = tpr
      )
      
      tmp <- tmp[
        is.finite(tmp$fpr) &
          is.finite(tmp$tpr),
        ,
        drop = FALSE
      ]
      
      tmp <- tmp[
        order(tmp$fpr, tmp$tpr),
        ,
        drop = FALSE
      ]
      
      tmp <- tmp[
        !duplicated(tmp$fpr),
        ,
        drop = FALSE
      ]
      
      # Interpolate TPR on common FPR grid
      interp_tpr <- approx(
        x = tmp$fpr,
        y = tmp$tpr,
        xout = fpr_grid,
        method = "linear",
        rule = 2
      )$y
      
      interpolated_rocs[[length(interpolated_rocs) + 1]] <-
        interp_tpr
    }
    
    if (length(interpolated_rocs) == 0) {
      stop(
        "No valid ROC curves found for ",
        nm,
        " at nspecies = ",
        ns
      )
    }
    
    tpr_matrix <- do.call(
      cbind,
      interpolated_rocs
    )
    
    mean_tpr <- rowMeans(
      tpr_matrix,
      na.rm = TRUE
    )
    
    sd_tpr <- apply(
      tpr_matrix,
      1,
      sd,
      na.rm = TRUE
    )
    
    roc_data_list[[nm]] <- data.frame(
      fpr = fpr_grid,
      mean_tpr = mean_tpr,
      sd_tpr = sd_tpr
    )
  }
  
  # ------------------------------------------------------------
  # 6. Plot
  # ------------------------------------------------------------
  
  pdf(
    output_pdf,
    width = 5.5,
    height = 5.5
  )
  
  first_model <- model_names[1]
  
  plot(
    roc_data_list[[first_model]]$fpr * 100,
    roc_data_list[[first_model]]$mean_tpr * 100,
    type = "l",
    col = color_vec[1],
    lwd = 2,
    xlim = c(0, 100),
    ylim = c(0, 100),
    xlab = "100 - Specificity (%)",
    ylab = "Sensitivity (%)",
    axes = FALSE
  )
  
  axis(
    1,
    at = seq(0, 100, 10),
    labels = seq(0, 100, 10)
  )
  
  axis(
    2,
    at = seq(0, 100, 10),
    labels = seq(0, 100, 10)
  )
  
  box()
  
  abline(
    0,
    1,
    lty = 2,
    col = "grey70"
  )
  
  # Additional ROC curves
  if (n_curves > 1) {
    
    for (i in 2:n_curves) {
      
      nm <- model_names[i]
      
      lines(
        roc_data_list[[nm]]$fpr * 100,
        roc_data_list[[nm]]$mean_tpr * 100,
        col = color_vec[i],
        lwd = 2
      )
    }
  }
  
  # ------------------------------------------------------------
  # 7. Legend
  # ------------------------------------------------------------
  
  legend(
    "bottomright",
    legend = model_names,
    col = color_vec[seq_len(n_curves)],
    lwd = 2,
    bty = "n"
  )
  
  # ------------------------------------------------------------
  # 8. AUC labels
  #
  # AUC = mean of seed-specific AUCs
  #
  # CI = t-based 95% CI across seed-specific AUCs
  # ------------------------------------------------------------
  
  auc_labels <- character(n_curves)
  
  for (i in seq_len(n_curves)) {
    
    nm <- model_names[i]
    
    ns <- plot_ns[[nm]]
    
    row <- auc_summary_list[[nm]] %>%
      dplyr::filter(
        nspecies == ns
      )
    
    if (nrow(row) != 1) {
      stop(
        "Cannot find unique AUC summary for ",
        nm,
        " at nspecies = ",
        ns
      )
    }
    
    mean_auc <- row$mean_auc
    sd_auc <- row$sd_auc
    n_seeds <- row$n_seeds
    
    if (
      show_auc_range &&
      n_seeds >= 2 &&
      is.finite(sd_auc)
    ) {
      
      t_value <- qt(
        0.975,
        df = n_seeds - 1
      )
      
      se <- sd_auc / sqrt(n_seeds)
      
      ci_lower <- mean_auc - t_value * se
      ci_upper <- mean_auc + t_value * se
      
      # Bound to valid AUC range
      ci_lower <- max(
        0,
        ci_lower
      )
      
      ci_upper <- min(
        100,
        ci_upper
      )
      
      auc_labels[i] <- sprintf(
        "%s (N=%d): %.1f%% (95%% CI: %.1f%%–%.1f%%)",
        nm,
        ns,
        mean_auc,
        ci_lower,
        ci_upper
      )
      
    } else {
      
      auc_labels[i] <- sprintf(
        "%s (N=%d): %.1f%%",
        nm,
        ns,
        mean_auc
      )
    }
  }
  
  legend(
    auc_text_position,
    legend = auc_labels,
    bty = "n",
    text.col = color_vec[seq_len(n_curves)],
    cex = 0.85
  )
  
  dev.off()
  
  # ------------------------------------------------------------
  # 9. Return results
  # ------------------------------------------------------------
  
  invisible(
    list(
      auc_summary = auc_summary_list,
      seed_auc = seed_auc_list,
      best_nspecies = plot_ns,
      roc_data = roc_data_list
    )
  )
}