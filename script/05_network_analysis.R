# ==============================================================================
# PTB Virome Project
# 05_network_analysis.R
# Network Analysis & Visualization (Co-occurrence, Degree Top20, Shared/Enriched)
# ==============================================================================

rm(list = ls())

source("script/config.R")

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(igraph)
library(ggraph)
library(ggtext)

# --------- 函数：拆分分类字符串生成表格 ----------
tax_table <- function(tax_vec) {
  tax_list <- strsplit(tax_vec, "\\|")
  max_level <- max(sapply(tax_list, length))

  # 补齐 NA
  tax_mat <- t(sapply(tax_list, function(x) {
    length(x) <- max_level
    x
  }))

  # 列名前缀映射
  rank_map <- c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order",
                f = "Family", g = "Genus", s = "Species", t = "Tribe")

  col_names <- sapply(1:ncol(tax_mat), function(i) {
    x <- tax_mat[which(!is.na(tax_mat[,i]))[1], i]
    if(is.na(x)) return(NA)
    prefix <- sub("__.*", "", x)
    rank_map[prefix]
  })

  tax_df <- as.data.frame(tax_mat, stringsAsFactors = FALSE)
  colnames(tax_df) <- col_names

  # 去掉前缀
  tax_df <- as.data.frame(lapply(tax_df, function(x) sub("^[a-z]__","", x)), stringsAsFactors = FALSE)

  return(tax_df)
}

# --------- 函数：读取节点信息 ----------
read_node <- function(node_path, enriched_df) {
  read.table(node_path, sep = ",", header = TRUE, check.names = FALSE) %>%
    select(name, Degree) %>%
    left_join(enriched_df %>% select(feature, Group), by = c("name" = "feature")) %>%
    mutate(type = ifelse(grepl("^v\\d+", name), "vir", "bac"))
}

# --------- 函数：处理 MaAsLin2 显著结果 ----------
read_and_process <- function(path) {
  read.table(path, header = TRUE) %>%
    filter(metadata == "Group") %>%
    mutate(Group = ifelse(coef > 0, "HC", "TB"))
}

# ------------------------------------------------------------------------------
# 1. 主流程：数据读取与整合
# ------------------------------------------------------------------------------
# 读取显著结果
enriched_imfor <- bind_rows(
  read_and_process(file.path("results", 'maaslin','vir', 'significant_results.tsv')),
  read_and_process(file.path("results", 'maaslin','bac', 'significant_results.tsv'))
)

# 读取节点信息 (直接读取，无需判断)
tb_node <- read_node(TB_NODE_FILE, enriched_imfor)
hc_node <- read_node(HC_NODE_FILE, enriched_imfor)

# 读取菌种分类表格
bac <- read.delim("data/metaphlan4.profile", header = TRUE)
tax_df <- tax_table(bac$name)

# 读取边信息
tb_con <- read.delim(file.path(TABLE_DIR, 'tb_sig.re'), header = TRUE)
hc_con <- read.delim(file.path(TABLE_DIR, 'hc_sig.re'), header = TRUE)

hc_node_bac <- hc_node %>% filter(type == "bac")
tb_node_bac <- tb_node %>% filter(type == "bac")

species_family <- tax_df %>%
  select(Species, Family) %>%
  filter(!is.na(Species))

hc_node <- hc_node %>% left_join(species_family, by = c("name" = "Species"))
tb_node <- tb_node %>% left_join(species_family, by = c("name" = "Species"))

# ------------------------------------------------------------------------------
# 2. 基础主题定义
# ------------------------------------------------------------------------------
base_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black"),
    strip.background = element_rect(fill = "grey85"),
    strip.text = element_text(size = 12)
  )

# ------------------------------------------------------------------------------
# 3. HC/TB Enrich 图数据处理与绘制
# ------------------------------------------------------------------------------
hc_con2 <- hc_con %>%
  left_join(hc_node %>% select(name, Family), by = c("B" = "name")) %>%
  left_join(hc_node %>% filter(type == "vir") %>% select(name, Group),
            by = c("A" = "name")) %>%
  rename(enrich = Group) %>%
  distinct(A, Family, enrich)

hc_family_count <- hc_con2 %>%
  group_by(Family, enrich) %>%
  summarise(vOTU_count = n(), .groups = "drop") %>%
  mutate(Group = "HC")

tb_con2 <- tb_con %>%
  left_join(tb_node %>% select(name, Family), by = c("B" = "name")) %>%
  left_join(tb_node %>% filter(type == "vir") %>% select(name, Group),
            by = c("A" = "name")) %>%
  rename(enrich = Group) %>%
  distinct(A, Family, enrich)

tb_family_count <- tb_con2 %>%
  group_by(Family, enrich) %>%
  summarise(vOTU_count = n(), .groups = "drop") %>%
  mutate(Group = "TB")

all_family_count <- bind_rows(hc_family_count, tb_family_count)

# 补全缺失值
complete_family <- expand.grid(
  Family = unique(all_family_count$Family),
  Group = c("HC", "TB"),
  enrich = c("HC", "TB"),
  stringsAsFactors = FALSE
)

all_family_count <- complete_family %>%
  left_join(all_family_count, by = c("Family", "Group", "enrich")) %>%
  mutate(vOTU_count = ifelse(is.na(vOTU_count), 0, vOTU_count))

# x 轴排序
family_order <- all_family_count %>%
  group_by(Family) %>%
  summarise(total = sum(vOTU_count), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(Family)

all_family_count$Family <- factor(all_family_count$Family, levels = family_order)

# 绘制上下拼接堆叠图
p_hc <- ggplot(all_family_count %>% filter(Group == "HC"),
               aes(x = Family, y = vOTU_count, fill = enrich)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(HC = "#2E8153", TB = "#BC5858")) +
  labs(x = NULL, y = "Number of vOTUs", title = "HC") +
  base_theme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top")

p_tb <- ggplot(all_family_count %>% filter(Group == "TB"),
               aes(x = Family, y = vOTU_count, fill = enrich)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(HC = "#2E8153", TB = "#BC5858")) +
  labs(x = "Bacterial Family", y = "Number of vOTUs", title = "TB") +
  base_theme

final_enrich_plot <- p_hc / p_tb

# ------------------------------------------------------------------------------
# 4. Shared/Special 图数据处理与绘制
# ------------------------------------------------------------------------------
hc_votu <- hc_con2 %>% distinct(A, Family)
tb_votu <- tb_con2 %>% distinct(A, Family)

shared_votu <- inner_join(hc_votu, tb_votu, by = c("A", "Family")) %>%
  mutate(enrich = "shared") %>%
  select(A, Family, enrich)

hc_only <- anti_join(hc_votu, shared_votu, by = c("A", "Family")) %>%
  mutate(enrich = "HC-special")
hc_plot_df <- bind_rows(hc_only, shared_votu)

tb_only <- anti_join(tb_votu, shared_votu, by = c("A", "Family")) %>%
  mutate(enrich = "TB-special")
tb_plot_df <- bind_rows(tb_only, shared_votu)

hc_family_count_s <- hc_plot_df %>%
  group_by(Family, enrich) %>%
  summarise(vOTU_count = n(), .groups = "drop") %>%
  mutate(Group = "HC")

tb_family_count_s <- tb_plot_df %>%
  group_by(Family, enrich) %>%
  summarise(vOTU_count = n(), .groups = "drop") %>%
  mutate(Group = "TB")

# 设置 factor 顺序：shared 在底
hc_family_count_s$enrich <- factor(hc_family_count_s$enrich,
                                   levels = c("shared","HC-special"))
tb_family_count_s$enrich <- factor(tb_family_count_s$enrich,
                                   levels = c("shared","TB-special"))

all_family_count_s <- bind_rows(hc_family_count_s, tb_family_count_s)

# 补全缺失值
complete_family_s <- expand.grid(
  Family = unique(all_family_count_s$Family),
  Group = c("HC","TB"),
  enrich = c("shared","HC-special","TB-special"),
  stringsAsFactors = FALSE
)

all_family_count_s <- complete_family_s %>%
  left_join(all_family_count_s, by = c("Family","Group","enrich")) %>%
  mutate(vOTU_count = ifelse(is.na(vOTU_count),0,vOTU_count))

# x 轴排序按照总数
family_order_s <- all_family_count_s %>%
  group_by(Family) %>%
  summarise(total=sum(vOTU_count), .groups="drop") %>%
  arrange(desc(total)) %>%
  pull(Family)

all_family_count_s$Family <- factor(all_family_count_s$Family, levels = family_order_s)

plot_family_bar <- function(df, title_text){
  
  # 关键：指定 fill 的因子顺序，shared 放最前（即堆叠在最上方）
  df$enrich <- factor(df$enrich, levels = c("shared", "HC-special", "TB-special"))
  
  ggplot(df, aes(x = Family, y = vOTU_count, fill = enrich)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(
      values = c("HC-special" = "#9db7a5",
                 "TB-special" = "#9db7a5",
                 "shared"     = "#fbb1a2"),
      breaks = c("shared", "HC-special", "TB-special")   # 图例顺序也一起改
    ) +
    labs(x = "Bacterial Family", y = "Number of vOTUs",
         fill = "Category", title = title_text) +
    base_theme
}

p_hc_shared <- plot_family_bar(all_family_count_s %>% filter(Group=="HC"), "HC Shared/HC-special") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="top")
p_tb_shared <- plot_family_bar(all_family_count_s %>% filter(Group=="TB"), "TB Shared/TB-special")

final_shared_plot <- p_hc_shared / p_tb_shared

# 保存柱状图结果
dir.create(file.path("results", "figures"), recursive = TRUE, showWarnings = FALSE)
ggsave(file.path("results", "figures", "share_plot.pdf"), final_shared_plot, width = 5, height = 6)
ggsave(file.path("results", "figures", "enrich_network.pdf"), final_enrich_plot, width = 5, height = 6)

# ------------------------------------------------------------------------------
# 5. 绘制网络圆圈图 (Circle Network Layout)
# --------------------------------------------------------------------

# --------- 通用网络图绘制函数 ----------
plot_circle_network <- function(node_file, edge_file, title_text) {
  
  # 1. 读取并处理节点属性
  node_attr <- read.table(node_file, sep = ",", header = TRUE, check.names = FALSE) %>%
    select(name, Degree) %>%
    left_join(enriched_imfor %>% select(feature, Group), by = c("name" = "feature")) %>%
    mutate(type = ifelse(grepl("^v\\d+", name), "vir", "bac"))
  
  # 排序与构造虚拟节点层级
  node_attr <- node_attr[order(-nchar(node_attr$name)), ]
  node_attr$order <- paste0(node_attr$name, "_order_")
  
  # 构造树结构 (origin -> order -> real node)
  d1 <- data.frame(from = 'origin', to = unique(node_attr$order))
  d2 <- data.frame(from = node_attr$order, to = node_attr$name)
  tree_df <- rbind(d1, d2)
  
  # 补充虚拟节点表
  temp_dt <- as.data.frame(matrix(NA,
                                  nrow = nrow(d1) + 1,
                                  ncol = ncol(node_attr),
                                  dimnames = list(NULL, colnames(node_attr))))
  temp_dt$name <- c("origin", d1$to)
  node_attr2 <- rbind(temp_dt, node_attr)
  rownames(node_attr2) <- node_attr2$name
  
  # 2. 构建 igraph 对象并生成 circular dendrogram 布局
  graph_data <- graph_from_data_frame(tree_df, vertices = node_attr2)
  
  layout_df <- create_layout(graph_data, layout = "dendrogram", circular = TRUE) %>%
    left_join(node_attr2 %>% select(name, Group, type), by = "name") %>%
    mutate(Group = Group.y, type = type.y) %>%
    select(-Group.y, -type.y) %>%
    mutate(is_real = !grepl("_order_", name) & name != "origin")
  
  # 3. 读取并生成边 Bundle 数据
  conn_df <- read.delim(edge_file, header = TRUE)[, c("A", "B", "corr")]
  from_idx <- match(conn_df$A, layout_df$name)
  to_idx   <- match(conn_df$B, layout_df$name)
  conn_df$edge_color <- ifelse(conn_df$corr > 0, "positive", "negative")
  
  edge_data <- get_con(
    from = from_idx,
    to = to_idx,
    corr = conn_df$corr,
    edge_color = conn_df$edge_color
  )
  
  # 4. 绘图
  p <- ggraph(layout_df) +
    geom_conn_bundle(
      data = edge_data,
      aes(edge_colour = edge_color),
      edge_width = 0.1,
      tension = 0.8,
      alpha = 0.3
    ) +
    geom_node_point(aes(fill = type, shape = type, size = type, filter = is_real)) +
    geom_node_text(
      aes(x = x * 1.04, y = y * 1.04, filter = is_real, label = name,
          colour = Group,
          angle = atan2(y, x) * 360 / (2 * pi),
          hjust = 'outward'),
      size = 1, alpha = 1
    ) +
    scale_color_manual(values = c(HC = "#2E8153", TB = "#BC5858")) +
    scale_edge_colour_manual(values = c(positive = "#BC5858", negative = "#39A2AE")) +
    scale_shape_manual(values = c(vir = 21, bac = 24)) +
    scale_fill_manual(values = c(vir = "#A0CC58", bac = "#8A9CC4")) +
    scale_size_manual(values = c(vir = 1, bac = 1)) +
    ggtitle(title_text) +
    coord_fixed() +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  return(p)
}

# --------- 分别绘制 HC 和 TB 网络图 ----------
p_net_hc <- plot_circle_network(
  node_file  = HC_NODE_FILE,
  edge_file  = file.path(TABLE_DIR, 'hc_sig.re'),
  title_text = "Health (HC) Network"
)

p_net_tb <- plot_circle_network(
  node_file  = TB_NODE_FILE,
  edge_file  = file.path(TABLE_DIR, 'tb_sig.re'),
  title_text = "Tuberculosis (TB) Network"
)

p_net_combined <- p_net_hc | p_net_tb
ggsave(file.path("results", "figures", "cor_combined.pdf"), p_net_combined, width = 20, height = 10)
# ------------------------------------------------------------------------------
# 6. 绘制 Importance Degree Top20 图
# ------------------------------------------------------------------------------
top_20_hc$formatted_label <- ifelse(
  top_20_hc$Group == "TB",
  paste0('<span style="color:#E64B35;">', top_20_hc$name, '</span>'),
  paste0('<span style="color:#33a02c;">', top_20_hc$name, '</span>')
)

top_20_tb$formatted_label <- ifelse(
  top_20_tb$Group == "TB",
  paste0('<span style="color:#E64B35;">', top_20_tb$name, '</span>'),
  paste0('<span style="color:#33a02c;">', top_20_tb$name, '</span>')
)

# =========================
# 2. 把富文本列设为 factor，并按 Degree 排序（保证 y 轴顺序）
# =========================
top_20_hc$name_md <- factor(
  top_20_hc$formatted_label,
  levels = top_20_hc$formatted_label[order(top_20_hc$Degree)]
)

top_20_tb$name_md <- factor(
  top_20_tb$formatted_label,
  levels = top_20_tb$formatted_label[order(top_20_tb$Degree)]
)

# =========================
# 3. 绘图函数（TB / HC 共用样式）
# =========================
make_deg_plot <- function(df, title) {
  ggplot(df, aes(x = Degree, y = name_md, fill = type)) +
    geom_col() +
    theme_minimal() +
    scale_fill_manual(values = c("vir" = "#9fcc56", "bac" = "#889bc3")) +
    labs(title = title, x = "Degree", y = "Name") +
    theme(
      axis.text.y      = ggtext::element_markdown(size = 10),
      axis.text.y.left = ggtext::element_markdown(size = 10),
      plot.title       = element_text(hjust = 0.5)
    )
}

# =========================
# 4. 分别出图
# =========================
p_deg_hc <- make_deg_plot(top_20_hc, "Top 20 Species by Degree HC")
p_deg_tb <- make_deg_plot(top_20_tb, "Top 20 Species by Degree TB")

# =========================
# 5. 拼图
# =========================
p_deg_final <- p_deg_hc | p_deg_tb
ggsave(file.path("results", "figures", "degree.pdf"), p_deg_final, width = 10, height = 4)
