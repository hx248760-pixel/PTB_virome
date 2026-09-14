calculate_alpha_diversity <- function(abundance_data, sample_col, group_col, group_data, result_name = NULL) {
  # 检查必要的包是否安装
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("需要安装vegan包，请运行 install.packages('vegan')")
  }
  
  abundance_matrix <- t(abundance_data)
  # 计算各种Alpha多样性指数
  shannon <- vegan::diversity(abundance_matrix, index = "shannon")
  simpson <- vegan::diversity(abundance_matrix, index = "simpson")
  richness <- vegan::specnumber(abundance_matrix)
  
  # 构建Alpha多样性数据框
  alpha_df <- data.frame(
    Sample = names(shannon),
    Simpson = simpson,
    Shannon = shannon,
    Richness = richness,
    check.names = FALSE
  )
  
  # 合并分组信息
  group_subset <- group_data[, c(sample_col, group_col)]
  colnames(group_subset)[colnames(group_subset) == sample_col] <- "Sample"
  result_df <- merge(alpha_df, group_subset, by = "Sample", all.x = TRUE)
  colnames(result_df)[colnames(result_df) == group_col] <- group_col
  
  # 将结果保存到当前环境（如果指定了变量名）
  if (!is.null(result_name)) {
    assign(result_name, result_df, envir = .GlobalEnv)
    message("结果已保存到当前环境，变量名为: ", result_name)
  }
  
  # 返回结果（即使保存到环境，仍返回数据框方便即时查看）
  return(result_df)
}

# 使用示例：
# calculate_alpha_diversity(
#   abundance_data = coQ_abd,
#   sample_col = "Sample",
#   group_col = "Group",
#   group_data = sample_coQ,
#   result_name = "coQ_a_group_res"  # 保存到环境中的变量名
# )

