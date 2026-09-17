library(ggplot2)
library(dplyr)

plot_alpha_diversity <- function(data1, data2 = NULL, # 将data2设置为可选参数，默认为NULL
                                 index = c("Simpson", "Shannon", "Richness"),
                                 group_col = "Group",
                                 disease_group = "disease",
                                 control_group = "Control",
                                 cohort1_name = "Cohort1",
                                 cohort2_name = "Cohort2", # 仅在双队列模式下使用
                                 color_vec = c(
                                   "#FF6B6B",  # Cohort1 control
                                   "#4ECDC4",  # Cohort1 disease
                                   "#FFD166",  # Cohort2 control
                                   "#06D6A0"   # Cohort2 disease
                                 ),
                                 richness_col = "Richness",
                                 y_limits = NULL) {
  
  # 检查参数有效性
  index <- match.arg(index)
  required_cols <- c(group_col, index)
  if (index == "Richness") required_cols <- c(required_cols, richness_col)
  
  # 判断是否为单队列模式
  is_single_cohort <- is.null(data2)
  
  if (!all(required_cols %in% colnames(data1))) {
    stop("Cohort1 输入数据必须包含分组列和指定的多样性指数列")
  }
  if (!is_single_cohort && !all(required_cols %in% colnames(data2))) {
    stop("Cohort2 输入数据必须包含分组列和指定的多样性指数列")
  }
  
  # 1. 数据处理：添加队列标签并统一列名
  data1 <- data1 %>%
    mutate(Cohort = cohort1_name) %>%
    rename(Group = {{group_col}}) %>%
    {if(index == "Richness") rename(., Richness = all_of(richness_col)) else .}
  
  if (!is_single_cohort) {
    data2 <- data2 %>%
      mutate(Cohort = cohort2_name) %>%
      rename(Group = {{group_col}}) %>%
      {if(index == "Richness") rename(., Richness = all_of(richness_col)) else .}
    combined_data <- rbind(data1, data2)
  } else {
    combined_data <- data1
  }
  
  combined_data <- combined_data %>%
    mutate(
      Cohort_Group = paste(Cohort, Group, sep = "_"),
      Index_Value = pull(., {{index}})
    )
  
  # 2. Fisher方法合并p值函数
  fisher_merge_p <- function(p_values) {
    # 过滤掉NA值，因为单队列模式下会有NA
    p_values <- p_values[!is.na(p_values)]
    if (length(p_values) == 0) return(NA)
    if (length(p_values) == 1) return(p_values[1]) # 如果只有一个p值，直接返回
    
    chi_square <- -2 * sum(log(p_values))
    df <- 2 * length(p_values)
    p_combined <- pchisq(chi_square, df, lower.tail = FALSE)
    return(p_combined)
  }
  
  # 3. 计算各组Wilcoxon检验p值
  extract_index <- function(data, group_val) {
    data %>% filter(Group == group_val) %>% pull({{index}})
  }
  
  cohort1_disease <- extract_index(data1, disease_group)
  cohort1_con <- extract_index(data1, control_group)
  
  # 检查是否有足够的样本进行检验
  if (length(cohort1_disease) < 2 || length(cohort1_con) < 2) {
    cohort1_p <- NA
    warning(paste("Cohort1: Not enough samples in disease or control group to perform Wilcoxon test for", index, "."))
  } else {
    cohort1_p <- wilcox.test(cohort1_disease, cohort1_con, paired = FALSE)$p.value
  }
  
  cohort2_p <- NA
  combined_p <- NA
  adjusted_p <- NA
  
  if (!is_single_cohort) {
    cohort2_disease <- extract_index(data2, disease_group)
    cohort2_con <- extract_index(data2, control_group)
    
    if (length(cohort2_disease) < 2 || length(cohort2_con) < 2) {
      cohort2_p <- NA
      warning(paste("Cohort2: Not enough samples in disease or control group to perform Wilcoxon test for", index, "."))
    } else {
      cohort2_p <- wilcox.test(cohort2_disease, cohort2_con, paired = FALSE)$p.value
    }
    
    # 只有当两个队列的p值都有效时才合并
    valid_p_values <- c(cohort1_p, cohort2_p)
    valid_p_values <- valid_p_values[!is.na(valid_p_values)]
    
    if (length(valid_p_values) > 0) { # 至少有一个有效p值
        combined_p <- fisher_merge_p(valid_p_values)
        if (!is.na(combined_p)) {
            adjusted_p <- p.adjust(combined_p, method = "BH")
        }
    }
  }
  
  # 4. 输出统计结果
  if (!is_single_cohort) {
    results <- data.frame(
      Cohort = c(cohort1_name, cohort2_name, "Combined"),
      p_value = c(cohort1_p, cohort2_p, combined_p),
      adjusted_p = c(NA, NA, adjusted_p)
    )
  } else {
    results <- data.frame(
      Cohort = cohort1_name,
      p_value = cohort1_p,
      adjusted_p = NA
    )
  }
  cat(paste0(index, "指数的Wilcoxon检验和Fisher合并p值结果:\n"))
  print(results)
  
  # 5. 绘图数据准备
  # 确定y轴范围（用户指定或自动计算）
  if (is.null(y_limits)) {
    data_min <- min(combined_data$Index_Value, na.rm = TRUE)
    data_max <- max(combined_data$Index_Value, na.rm = TRUE)
    data_range <- data_max - data_min
    y_limits <- c(data_min - data_range * 0.1, data_max + data_range * 0.2)
  }
  
  # 计算p值标签位置（基于实际y轴范围）
  max_cohort1_val <- max(data1 %>% pull({{index}}), na.rm = TRUE)
  global_max <- max_cohort1_val
  if (!is_single_cohort) {
    max_cohort2_val <- max(data2 %>% pull({{index}}), na.rm = TRUE)
    global_max <- max(max_cohort1_val, max_cohort2_val, na.rm = TRUE)
  }
  
  y_range <- y_limits[2] - y_limits[1]
  offset1 <- y_range * 0.05  # 队列内p值偏移
  offset2 <- y_range * 0.1   # 合并p值偏移
  
  # 颜色映射
  if (!is_single_cohort) {
    color_names <- c(
      paste(cohort1_name, control_group, sep = "_"),
      paste(cohort1_name, disease_group, sep = "_"),
      paste(cohort2_name, control_group, sep = "_"),
      paste(cohort2_name, disease_group, sep = "_")
    )
    if (length(color_vec) < 4) stop("在双队列模式下，color_vec 必须至少包含4种颜色。")
    color_map_values <- color_vec[1:4]
  } else {
    color_names <- c(
      paste(cohort1_name, control_group, sep = "_"),
      paste(cohort1_name, disease_group, sep = "_")
    )
    if (length(color_vec) < 2) stop("在单队列模式下，color_vec 必须至少包含2种颜色。")
    color_map_values <- color_vec[1:2]
  }
  names(color_map_values) <- color_names
  
  # p值格式化
  format_p_value <- function(p) {
    if (is.na(p)) {
      return("N.A.")
    } else if (p < 0.001) {
      return("p < 0.001")
    } else {
      return(paste("p =", format(p, digits = 3)))
    }
  }
  
  # 6. 绘制箱线图（支持自定义y轴范围）
  y_label <- case_when(
    index == "Simpson" ~ "Simpson index",
    index == "Shannon" ~ "Shannon index",
    index == "Richness" ~ "Observed vOTUs number"
  )
  
  p <- ggplot(combined_data, aes(x = Cohort, y = Index_Value, fill = Cohort_Group)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(0.8)) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                 position = position_dodge(0.8), fill = "white") +
    labs(
      title = "",
      x = "",
      y = y_label,
      fill = "Group"
    ) +
    scale_fill_manual(
      values = color_map_values,
      labels = color_names
    ) +
    coord_cartesian(ylim = y_limits) +
    theme_minimal() +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(size = 10),
      legend.position = "right"
    ) +
    theme_bw()
  
  # 添加p值注释
  if (!is.na(cohort1_p)) {
    p <- p + annotate("text", x = 1, y = max_cohort1_val + offset1,
                      label = format_p_value(cohort1_p), size = 4)
  }
  
  if (!is_single_cohort && !is.na(cohort2_p)) {
    p <- p + annotate("text", x = 2, y = max_cohort2_val + offset1,
                      label = format_p_value(cohort2_p), size = 4)
    if (!is.na(adjusted_p)) {
      p <- p + annotate("text", x = 1.5, y = global_max + offset2,
                        label = paste("Combined p:", format_p_value(adjusted_p)),
                        size = 4, fontface = "bold")
    }
  }
  
  return(p)
}

# --- 使用示例 ---

# 假设你的数据结构如下：
# co1_a_group_res 和 co2_a_group_res 都是数据框，
# 包含 'Group' 列（例如 "Control", "disease"）和多样性指数列（例如 "Simpson", "Shannon", "Richness" 或 "richness"）

# 模拟数据 (请替换为你的实际数据)
#set.seed(123)
#co1_a_group_res <- data.frame(
#  Group = sample(c("Control", "disease"), 50, replace = TRUE),
#  Simpson = runif(50, 0.9, 0.99),
#  Shannon = runif(50, 2.5, 4.0),
#  Richness = sample(100:300, 50, replace = TRUE),
#  richness = sample(100:300, 50, replace = TRUE) # 兼容 richness_col 参数
#)

#co2_a_group_res <- data.frame(
#  Group = sample(c("Control", "disease"), 40, replace = TRUE),
#  Simpson = runif(40, 0.85, 0.98),
#  Shannon = runif(40, 2.0, 3.8),
#  Richness = sample(80:280, 40, replace = TRUE),
#  richness = sample(80:280, 40, replace = TRUE)
#)

# --- 1. 绘制两个队列的Alpha多样性图片 (与之前功能相同) ---
# 自动y轴范围
#p_simpson_dual_auto <- plot_alpha_diversity(
#  data1 = co1_a_group_res,
#  data2 = co2_a_group_res,
#  index = "Simpson",
#  cohort1_name = "Cohort A",
#  cohort2_name = "Cohort B"
#)
#print(p_simpson_dual_auto)

# 自定义y轴范围
#p_richness_dual_custom <- plot_alpha_diversity(
#  data1 = co1_a_group_res,
#  data2 = co2_a_group_res,
#  index = "Richness",
#  richness_col = "richness", # 如果你的丰富度列名是 'richness'
#  cohort1_name = "Cohort A",
#  cohort2_name = "Cohort B",
#  y_limits = c(50, 350)
#)
#print(p_richness_dual_custom)

# --- 2. 绘制单个队列的Alpha多样性图片 (data2 = NULL) ---
# 自动y轴范围
#p_shannon_single_auto <- plot_alpha_diversity(
#  data1 = co1_a_group_res,
#  index = "Shannon",
#  cohort1_name = "My Cohort" # 此时 cohort2_name 不起作用
#)
#print(p_shannon_single_auto)

# 自定义y轴范围
#p_simpson_single_custom <- plot_alpha_diversity(
#  data1 = co2_a_group_res,
#  index = "Simpson",
#  cohort1_name = "Another Cohort",
#  y_limits = c(0.8, 1.0)
#)
#print(p_simpson_single_custom)

# --- 3. 示例：当某组样本不足时 ---
#co_insufficient <- data.frame(
#  Group = c("Control", "Control", "disease"), # disease组只有1个样本
#  Simpson = c(0.95, 0.96, 0.90)
#)
#p_insufficient <- plot_alpha_diversity(
#  data1 = co_insufficient,
#  index = "Simpson",
#  cohort1_name = "Insufficient Samples"
#)
#print(p_insufficient) # 会有警告信息，p值显示为 N.A.
