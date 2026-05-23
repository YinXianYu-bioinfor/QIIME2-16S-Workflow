#!/usr/bin/env Rscript
#
# QIIME2 16S 扩增子分析 — 综合可视化与数据整理脚本
# ================================================================================
# 本脚本适用于标准 QIIME2 流程导出的结果文件，
# 生成出版级可视化图表，并准备 FAPROTAX 和 PICRUSt2 下游功能预测的输入文件。
#
# 使用方式：
#   1. 将本脚本放入 QIIME2 项目根目录 (与 metadata.txt 同级, export/ 由 pipeline 生成)
#   2. 修改下方 setwd() 路径后直接运行:
#      Rscript QIIME2_16S_visualization.R
#   3. 脚本自动读取 export/ 中的数据, 输出也在该目录下各子目录
#
# 运行后目录结构（项目根目录下）：
#   export/
#   ├── alpha/           -- α多样性箱线图
#   ├── beta/            -- β多样性 PCoA 图
#   ├── taxa/            -- 门水平物种组成图
#   ├── heatmap/         -- 属水平热图
#   ├── faprotax/        -- FAPROTAX 输入表
#   ├── picrust2/        -- PICRUSt2 输入文件
#   └── feature_tables/  -- 处理后的分类学与丰度表
# ================================================================================

# ============================================================
# 0. 配置区 — 请根据项目修改以下路径
# ============================================================

# 【重要】设置工作目录到 QIIME2 项目根目录
# 请将下方路径修改为您电脑上项目文件夹的实际路径
# 项目根目录应包含: metadata.txt、export/ (由 qiime2-16s-pipeline.sh 生成)
setwd("D:/Drivers/桌面/scripts/amplicon_analysis_script/QIIME2-16S-Workflow/examples")

# metadata.txt 路径（相对于工作目录或绝对路径）
metadata_file <- "metadata.txt"

# 分组列名（metadata.txt 中用于分组的列）
group_col <- "Group"

# 可视化参数
# 所有输入文件由 qiime2-16s-pipeline.sh 导出至 export/ 目录
alpha_diversity_file <- "export/evenness_vector.tsv"
feature_table_file <- "export/feature-table.tsv"
taxonomy_file <- "export/taxonomy.tsv"
dna_sequences_file <- "export/dna-sequences.fasta"
rarefied_biom_file <- "export/rarefied_table.biom"

# 输出子目录（统一在 export/ 下, 由 qiime2-16s-pipeline.sh 自动创建）
output_dirs <- c(
  "export/alpha", "export/beta", "export/taxa",
  "export/heatmap", "export/faprotax",
  "export/picrust2", "export/feature_tables"
)

# 热图显示 Top N 属
top_n_genera <- 30

# 堆叠柱状图显示 Top N 门（其余合并为 Others）
top_n_phyla_stacked <- 8

# 配色方案（ggsci 期刊配色）
# 可选: "npg" (Nature), "aaas" (Science), "nejm" (NEJM),
#       "lancet" (Lancet), "jco" (JCO), "jama" (JAMA), "uchicago"
color_palette <- "npg"

# ============================================================
# 1. 加载所需的 R 包
# ============================================================

cat("========================================\n")
cat("QIIME2 16S 可视化流程\n")
cat("========================================\n\n")

required_packages <- c(
  "ggplot2", "tidyr", "dplyr", "tibble", "readr",
  "vegan", "ape", "pheatmap", "RColorBrewer", "grDevices",
  "patchwork", "ggsignif", "ggsci"
)

check_packages <- function(pkgs) {
  missing <- c()
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    }
  }
  if (length(missing) > 0) {
    cat("正在安装缺失的包:", paste(missing, collapse = ", "), "\n")
    install.packages(missing, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  invisible(lapply(pkgs, library, character.only = TRUE, quietly = TRUE))
}

cat("[1/8] 加载 R 包...\n")
check_packages(required_packages)
cat("       所有包加载成功。\n\n")

# 配色辅助函数：根据分组数从 ggsci 期刊色板获取颜色
get_group_colors <- function(n, palette_name = color_palette) {
  base <- switch(palette_name,
    npg    = pal_npg()(10),
    aaas   = pal_aaas()(10),
    nejm   = pal_nejm()(8),
    lancet = pal_lancet()(9),
    jco    = pal_jco()(10),
    jama   = pal_jama()(7),
    uchicago = pal_uchicago()(9),
    pal_npg()(10)
  )
  if (n <= length(base)) base[1:n] else grDevices::colorRampPalette(base)(n)
}

# ============================================================
# 2. 创建输出目录
# ============================================================

cat("[2/8] 创建输出目录...\n")
for (dir_name in output_dirs) {
  if (!dir.exists(dir_name)) {
    dir.create(dir_name, recursive = TRUE)
    cat("       已创建:", dir_name, "\n")
  }
}
cat("       输出目录准备就绪。\n\n")

# ============================================================
# 3. 读取元数据
# ============================================================

cat("[3/8] 读取输入数据...\n")

if (!file.exists(metadata_file)) {
  stop("未找到元数据文件: ", metadata_file,
    "\n请检查上方配置区中的 metadata_file 路径。")
}

metadata <- read.table(
  metadata_file,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  strip.white = TRUE,
  stringsAsFactors = FALSE
)

# 清理列名（去掉可能的 # 前缀）
colnames(metadata) <- gsub("^#", "", colnames(metadata))

# 检查分组列是否存在
if (!group_col %in% colnames(metadata)) {
  stop("元数据中未找到分组列 '", group_col, "'。\n",
    "可用的列: ", paste(colnames(metadata), collapse = ", "))
}

# 将分组列转为因子
metadata[[group_col]] <- as.factor(metadata[[group_col]])
cat("       元数据加载成功:", nrow(metadata), "个样本,",
  ncol(metadata), "列\n")
cat("       分组情况:", paste(levels(metadata[[group_col]]), collapse = " / "), "\n")

# ============================================================
# 4. 读取特征表（Feature Table）
# ============================================================

if (!file.exists(feature_table_file)) {
  stop("未找到特征表: ", feature_table_file)
}

# 跳过第一行注释（# Constructed from biom file）
feature_table <- read.table(
  feature_table_file,
  sep = "\t",
  skip = 1,
  header = TRUE,
  check.names = FALSE,
  row.names = 1,
  comment.char = ""
)

# 确保数据为数值矩阵
feature_table <- as.matrix(feature_table)
mode(feature_table) <- "numeric"

# 只保留元数据中存在的样本
common_samples <- intersect(colnames(feature_table), rownames(metadata))
feature_table <- feature_table[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

cat("       特征表加载成功:", nrow(feature_table), "个ASV,",
  ncol(feature_table), "个样本\n")
cat("       总 reads 数:", sum(feature_table), "\n")

# ============================================================
# 5. 读取并处理分类学注释（Taxonomy）
# ============================================================

if (!file.exists(taxonomy_file)) {
  stop("未找到分类学文件: ", taxonomy_file)
}

taxonomy_raw <- read.table(
  taxonomy_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 解析分类学字符串为各层级
parse_taxonomy <- function(taxon_string) {
  result <- c("Unassigned", "Unassigned", "Unassigned",
              "Unassigned", "Unassigned", "Unassigned", "Unassigned")
  names(result) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  parts <- strsplit(trimws(taxon_string), ";")[[1]]
  for (part in parts) {
    if (nchar(part) < 3) next
    prefix <- substr(part, 1, 1)
    value <- substr(part, 4, nchar(part))
    if (value == "") value <- "Unassigned"

    level_name <- switch(prefix,
      "d" = "Kingdom",
      "p" = "Phylum",
      "c" = "Class",
      "o" = "Order",
      "f" = "Family",
      "g" = "Genus",
      "s" = "Species",
      NA
    )
    if (!is.na(level_name)) {
      result[level_name] <- value
    }
  }
  return(result)
}

# 解析所有分类学字符串
tax_parsed <- t(sapply(taxonomy_raw$Taxon, parse_taxonomy, USE.NAMES = FALSE))
taxonomy <- as.data.frame(tax_parsed, stringsAsFactors = FALSE)
taxonomy$FeatureID <- taxonomy_raw[["Feature ID"]]
taxonomy$Confidence <- taxonomy_raw$Confidence

# 构建完整分类学路径（用于 FAPROTAX）
taxonomy$FullPath <- apply(tax_parsed, 1, function(x) {
  paste(x, collapse = ";")
})

# 构建精简路径：去掉末端的 Unassigned
build_compact_path <- function(levels_vec) {
  non_na_idx <- which(levels_vec != "Unassigned")
  if (length(non_na_idx) == 0) return("Unassigned")
  last_idx <- max(non_na_idx)
  paste(levels_vec[1:last_idx], collapse = ";")
}
taxonomy$CompactPath <- apply(tax_parsed, 1, build_compact_path)

# 创建 Feature ID 到分类学信息的映射
tax_map <- setNames(taxonomy$FullPath, taxonomy$FeatureID)

cat("       分类学注释条目:", nrow(taxonomy), "\n")
cat("       已注释的门:", length(unique(taxonomy$Phylum[taxonomy$Phylum != "Unassigned"])), "\n")

# ============================================================
# 6. 保存处理后的数据表
# ============================================================

cat("\n[4/8] 保存处理后的数据表...\n")

# 6a. 处理后的分类学表（各层级拆分，Unassigned 填充）
taxonomy_out <- taxonomy[, c("FeatureID", "Kingdom", "Phylum", "Class",
                              "Order", "Family", "Genus", "Species",
                              "FullPath", "Confidence")]
colnames(taxonomy_out)[1] <- "Feature_ID"
write.table(
  taxonomy_out,
  file = "export/feature_tables/taxonomy_processed.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
cat("       [OK] export/feature_tables/taxonomy_processed.tsv\n")

# 6b. 带分类学注释的特征表
merged_data <- as.data.frame(feature_table)
merged_data$FeatureID <- rownames(merged_data)
merged_data <- merge(merged_data, taxonomy[, c("FeatureID", "FullPath", "CompactPath")],
                     by = "FeatureID", all.x = TRUE)

# 重排列序
merged_data <- merged_data[, c("FeatureID", "FullPath", "CompactPath",
                                setdiff(colnames(merged_data),
                                        c("FeatureID", "FullPath", "CompactPath")))]

write.table(
  merged_data,
  file = "export/feature_tables/feature_table_with_taxonomy.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
cat("       [OK] export/feature_tables/feature_table_with_taxonomy.tsv\n")

# ============================================================
# 7. α 多样性分析
# ============================================================

cat("\n[5/8] 生成可视化图表...\n")
cat("       正在处理 α 多样性...\n")

# 7a. 计算 α 多样性指标
# ------------------------------------------------------------------
# 从 feature_table 用 vegan 计算多个 α 多样性指标
ft_t <- t(feature_table)  # 转置：vegan 需要样本为行、ASV 为列

observed_features <- specnumber(ft_t)
shannon_div <- diversity(ft_t, index = "shannon")
simpson_div <- diversity(ft_t, index = "simpson")
est_rich <- estimateR(ft_t)
chao1_div <- est_rich["S.chao1", ]
ace_div   <- est_rich["S.ACE", ]

# Pielou evenness（从 alpha-diversity.tsv 读取，可选）
pielou_div <- setNames(rep(NA_real_, nrow(ft_t)), rownames(ft_t))
if (file.exists(alpha_diversity_file)) {
  alpha_data <- read.table(alpha_diversity_file, sep = "\t",
    header = TRUE, row.names = 1, check.names = FALSE)
  if ("pielou_evenness" %in% colnames(alpha_data)) {
    pielou_div <- alpha_data$pielou_evenness
  }
}

# 整合所有指标为数据框
alpha_metrics_df <- data.frame(
  SampleID = names(observed_features),
  observed_features = observed_features,
  shannon = shannon_div,
  simpson = simpson_div,
  chao1 = chao1_div,
  ACE = ace_div,
  pielou_evenness = pielou_div[names(observed_features)],
  row.names = NULL, stringsAsFactors = FALSE
)

# 合并元数据
alpha_plot_data <- merge(
  alpha_metrics_df,
  metadata[, group_col, drop = FALSE],
  by.x = "SampleID", by.y = "row.names"
)

# 保存所有 α 多样性指标到文件
write.table(alpha_plot_data,
  file = "export/feature_tables/alpha_diversity_metrics.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE)
cat("       [OK] export/feature_tables/alpha_diversity_metrics.tsv\n")

# 绘图指标（4 个主流指标）
plot_metrics <- c("observed_features", "shannon", "simpson", "chao1")
metric_labels <- c(
  observed_features = "Observed ASVs",
  shannon = "Shannon",
  simpson = "Simpson",
  chao1 = "Chao1"
)

# 转为长格式，过滤非有限值
alpha_long <- pivot_longer(
  alpha_plot_data,
  cols = all_of(plot_metrics),
  names_to = "Metric", values_to = "Value"
)
alpha_long[[group_col]] <- as.factor(alpha_long[[group_col]])
alpha_long <- alpha_long[is.finite(alpha_long$Value), ]

# 记录被过滤的样本
for (m in plot_metrics) {
  n_inf <- sum(!is.finite(alpha_plot_data[[m]]))
  if (n_inf > 0) {
    cat("        [注意]", m, "包含", n_inf, "个非有限值\n")
  }
}

# 7b. α 多样性统计检验
# ------------------------------------------------------------------
cat("       Alpha diversity statistics:\n")

for (metric in plot_metrics) {
  df_sub <- alpha_plot_data[is.finite(alpha_plot_data[[metric]]), ]

  # Kruskal-Wallis
  kw <- kruskal.test(as.formula(paste(metric, "~", group_col)), data = df_sub)
  cat("         ", metric, ": Kruskal-Wallis p =",
      format.pval(kw$p.value, digits = 4), "\n")

  # 成对 Wilcoxon + BH 校正
  pw <- suppressWarnings(
    pairwise.wilcox.test(df_sub[[metric]], df_sub[[group_col]],
                         p.adjust.method = "BH"))
  p_mat <- pw$p.value
  if (!is.null(p_mat)) {
    g_levels <- levels(df_sub[[group_col]])
    for (i in seq_len(length(g_levels) - 1)) {
      for (j in (i + 1):length(g_levels)) {
        pv <- tryCatch(
          if (g_levels[j] %in% rownames(p_mat) && g_levels[i] %in% colnames(p_mat))
            p_mat[g_levels[j], g_levels[i]]
          else if (g_levels[i] %in% rownames(p_mat) && g_levels[j] %in% colnames(p_mat))
            p_mat[g_levels[i], g_levels[j]]
          else NA_real_,
          error = function(e) NA_real_)
        if (!is.na(pv)) {
          star <- if (pv < 0.001) "***" else if (pv < 0.01) "**" else if (pv < 0.05) "*" else "ns"
          cat("             ", g_levels[i], "vs", g_levels[j],
              ": p =", format.pval(pv, digits = 4), star, "\n")
        }
      }
    }
  }
}

# 7c. 带显著性星号的 α 多样性箱线图（4 面板拼接）
# ------------------------------------------------------------------
plot_list <- list()

for (metric in plot_metrics) {
  df_sub <- alpha_long[alpha_long$Metric == metric, ]
  g_levels <- levels(df_sub[[group_col]])

  # 成对 Wilcoxon 检验
  pw <- suppressWarnings(
    pairwise.wilcox.test(df_sub$Value, df_sub[[group_col]],
                         p.adjust.method = "BH"))
  p_mat <- pw$p.value

  # 构建 ggsignif 参数
  comp_list <- list()
  annot_vec <- c()
  for (i in seq_len(length(g_levels) - 1)) {
    for (j in (i + 1):length(g_levels)) {
      pv <- tryCatch(
        if (g_levels[j] %in% rownames(p_mat) && g_levels[i] %in% colnames(p_mat))
          p_mat[g_levels[j], g_levels[i]]
        else if (g_levels[i] %in% rownames(p_mat) && g_levels[j] %in% colnames(p_mat))
          p_mat[g_levels[i], g_levels[j]]
        else NA_real_,
        error = function(e) NA_real_)
      if (!is.na(pv)) {
        comp_list <- c(comp_list, list(c(g_levels[i], g_levels[j])))
        star <- if (pv < 0.001) "***" else if (pv < 0.01) "**" else if (pv < 0.05) "*" else "ns"
        annot_vec <- c(annot_vec, star)
      }
    }
  }

  p <- ggplot(df_sub, aes(x = !!sym(group_col), y = Value, fill = !!sym(group_col))) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
    geom_jitter(aes(color = !!sym(group_col)), width = 0.15, size = 1.5, alpha = 0.8) +
    scale_fill_manual(values = get_group_colors(nlevels(df_sub[[group_col]]))) +
    scale_color_manual(values = get_group_colors(nlevels(df_sub[[group_col]]))) +
    labs(title = metric_labels[metric], x = NULL, y = "Value") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank())

  if (length(comp_list) > 0) {
    p <- p + geom_signif(comparisons = comp_list, annotation = annot_vec,
                         step_increase = 0.1, tip_length = 0.02, textsize = 3.5)
  }
  plot_list[[metric]] <- p
}

# 2x2 拼接
p_alpha_combined <- wrap_plots(plotlist = plot_list, ncol = 2, nrow = 2)
ggsave(filename = "export/alpha/alpha_diversity_boxplot.pdf",
       plot = p_alpha_combined, width = 10, height = 8, device = cairo_pdf)
cat("       [OK] export/alpha/alpha_diversity_boxplot.pdf\n")

# 7d. 稀释曲线（按组平均）
# ------------------------------------------------------------------
cat("       正在绘制稀疏曲线...\n")

# 手动 rarefaction：用 apply + rarefy 计算矩阵，避免 lapply/do.call 的兼容性问题
rarefy_step <- 100
min_depth <- min(rowSums(ft_t))
rarefy_depths <- seq(from = rarefy_step, to = min_depth, by = rarefy_step)

rare_mat <- suppressWarnings(apply(ft_t, 1, function(row) {
  rarefy(row, rarefy_depths)
}))
# rare_mat: 矩阵，行 = 测序深度，列 = 样本
rare_df <- data.frame(
  SampleID = rep(colnames(rare_mat), each = nrow(rare_mat)),
  Reads    = rep(rarefy_depths, ncol(rare_mat)),
  ASVs     = as.vector(rare_mat),
  stringsAsFactors = FALSE
)
# 按 sampleID 匹配分组信息
rare_df[[group_col]] <- metadata[as.character(rare_df$SampleID), group_col, drop = TRUE]

# 按 Reads 深度和组计算均值与标准误（base R aggregate）
rare_agg <- aggregate(
  x = rare_df$ASVs,
  by = list(Group = rare_df[[group_col]], Reads = rare_df$Reads),
  FUN = function(x) c(mean = mean(x, na.rm = TRUE), se = sd(x, na.rm = TRUE) / sqrt(length(x)))
)
rare_summary <- do.call(data.frame, rare_agg)
colnames(rare_summary) <- c(group_col, "Reads", "mean_ASVs", "se_ASVs")

p_rarefaction <- ggplot(rare_summary, aes(x = Reads, y = mean_ASVs,
                                           color = !!sym(group_col),
                                           fill = !!sym(group_col))) +
  geom_ribbon(aes(ymin = mean_ASVs - se_ASVs, ymax = mean_ASVs + se_ASVs),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = get_group_colors(nlevels(rare_summary[[group_col]]))) +
  scale_fill_manual(values = get_group_colors(nlevels(rare_summary[[group_col]]))) +
  labs(title = "Rarefaction Curves",
       x = "Number of Reads", y = "Number of ASVs",
       color = group_col, fill = group_col) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right", panel.grid.minor = element_blank())

ggsave(filename = "export/alpha/rarefaction_curves.pdf",
       plot = p_rarefaction, width = 8, height = 6, device = cairo_pdf)
cat("       [OK] export/alpha/rarefaction_curves.pdf\n")

# ============================================================
# 8. β 多样性 — PCoA 和 PERMANOVA
# ============================================================

cat("       正在处理 β 多样性...\n")

# 过滤掉在所有样本中均为 0 的 ASV
ft_sub <- feature_table[rowSums(feature_table) > 0, ]

# ---------- 8a. Bray-Curtis PCoA ----------
cat("       Bray-Curtis PCoA...\n")
bc_dist <- vegdist(t(ft_sub), method = "bray")

pcoa_result <- cmdscale(bc_dist, k = min(10, ncol(ft_sub) - 1), eig = TRUE)
var_exp <- round(pcoa_result$eig[1:4] / sum(pcoa_result$eig[pcoa_result$eig > 0]) * 100, 1)

pcoa_df <- as.data.frame(pcoa_result$points)
colnames(pcoa_df) <- paste0("PCo", seq_len(ncol(pcoa_df)))
pcoa_df$SampleID <- rownames(pcoa_df)
pcoa_df <- merge(pcoa_df, metadata[, group_col, drop = FALSE],
                 by.x = "SampleID", by.y = "row.names")

# PERMANOVA
adonis_formula <- as.formula(paste("bc_dist ~", group_col))
adonis_result <- adonis2(adonis_formula, data = metadata, permutations = 999)
print(adonis_result)
r2 <- round(adonis_result[1, "R2"], 3)
p_col <- grep("Pr", colnames(adonis_result))[1]
pval <- adonis_result[1, p_col]
pval_str <- ifelse(pval < 0.001, format.pval(pval, digits = 4),
                   paste0("= ", format.pval(pval, digits = 4)))

# PCoA 图
p_pcoa <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2, color = !!sym(group_col))) +
  geom_point(size = 3, alpha = 0.8)

min_group_size <- min(table(pcoa_df[[group_col]]))
if (min_group_size >= 4) {
  p_pcoa <- p_pcoa +
    stat_ellipse(aes(fill = !!sym(group_col)),
      geom = "polygon", alpha = 0.1, level = 0.95, show.legend = FALSE)
}

p_pcoa <- p_pcoa +
  scale_color_manual(values = get_group_colors(nlevels(pcoa_df[[group_col]]))) +
  scale_fill_manual(values = get_group_colors(nlevels(pcoa_df[[group_col]]))) +
  labs(title = "Beta Diversity - Bray-Curtis PCoA",
       x = paste0("PCo1 (", var_exp[1], "%)"),
       y = paste0("PCo2 (", var_exp[2], "%)"),
       color = group_col,
       caption = paste0("PERMANOVA: R² = ", r2, ", p ", pval_str)) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right", panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0.5, size = 10))

ggsave(filename = "export/beta/beta_diversity_pcoa_bray_curtis.pdf",
       plot = p_pcoa, width = 8, height = 6, device = cairo_pdf)
cat("       [OK] export/beta/beta_diversity_pcoa_bray_curtis.pdf\n")

# ---------- 8b. Jaccard PCoA ----------
cat("       Jaccard PCoA...\n")
jacc_dist <- vegdist(t(ft_sub), method = "jaccard", binary = TRUE)

pcoa_jacc <- cmdscale(jacc_dist, k = min(10, ncol(ft_sub) - 1), eig = TRUE)
var_exp_jacc <- round(pcoa_jacc$eig[1:4] / sum(pcoa_jacc$eig[pcoa_jacc$eig > 0]) * 100, 1)

pcoa_jacc_df <- as.data.frame(pcoa_jacc$points)
colnames(pcoa_jacc_df) <- paste0("PCo", seq_len(ncol(pcoa_jacc_df)))
pcoa_jacc_df$SampleID <- rownames(pcoa_jacc_df)
pcoa_jacc_df <- merge(pcoa_jacc_df, metadata[, group_col, drop = FALSE],
                      by.x = "SampleID", by.y = "row.names")

# PERMANOVA
adonis_jacc <- adonis2(as.formula(paste("jacc_dist ~", group_col)),
                       data = metadata, permutations = 999)
print(adonis_jacc)
r2_jacc <- round(adonis_jacc[1, "R2"], 3)
p_col_jacc <- grep("Pr", colnames(adonis_jacc))[1]
pval_jacc <- adonis_jacc[1, p_col_jacc]
pval_jacc_str <- ifelse(pval_jacc < 0.001, format.pval(pval_jacc, digits = 4),
                        paste0("= ", format.pval(pval_jacc, digits = 4)))

# Jaccard PCoA 图
p_jacc <- ggplot(pcoa_jacc_df, aes(x = PCo1, y = PCo2, color = !!sym(group_col))) +
  geom_point(size = 3, alpha = 0.8)

if (min_group_size >= 4) {
  p_jacc <- p_jacc +
    stat_ellipse(aes(fill = !!sym(group_col)),
      geom = "polygon", alpha = 0.1, level = 0.95, show.legend = FALSE)
}

p_jacc <- p_jacc +
  scale_color_manual(values = get_group_colors(nlevels(pcoa_jacc_df[[group_col]]))) +
  scale_fill_manual(values = get_group_colors(nlevels(pcoa_jacc_df[[group_col]]))) +
  labs(title = "Beta Diversity - Jaccard PCoA",
       x = paste0("PCo1 (", var_exp_jacc[1], "%)"),
       y = paste0("PCo2 (", var_exp_jacc[2], "%)"),
       color = group_col,
       caption = paste0("PERMANOVA: R² = ", r2_jacc, ", p ", pval_jacc_str)) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right", panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0.5, size = 10))

ggsave(filename = "export/beta/beta_diversity_pcoa_jaccard.pdf",
       plot = p_jacc, width = 8, height = 6, device = cairo_pdf)
cat("       [OK] export/beta/beta_diversity_pcoa_jaccard.pdf\n")
cat("\n")

# ============================================================
# 9. 门水平堆叠柱状图
# ============================================================

cat("       正在处理门水平物种组成...\n")

# 将每个 ASV 映射到其门分类
asv_phylum <- setNames(taxonomy$Phylum, taxonomy$FeatureID)

# 只保留特征表中存在的 ASV
common_asvs <- intersect(names(asv_phylum), rownames(feature_table))
asv_phylum <- asv_phylum[common_asvs]
ft_phylum <- feature_table[common_asvs, , drop = FALSE]

if (length(common_asvs) > 0) {
  # 按门汇总丰度
  phylum_levels <- unique(asv_phylum)
  phylum_counts <- matrix(0, nrow = length(phylum_levels), ncol = ncol(ft_phylum))
  rownames(phylum_counts) <- phylum_levels
  colnames(phylum_counts) <- colnames(ft_phylum)

  for (i in seq_along(asv_phylum)) {
    phylum_counts[asv_phylum[i], ] <- phylum_counts[asv_phylum[i], ] + ft_phylum[i, ]
  }

  # 转换为相对丰度 (%)
  phylum_rel <- sweep(phylum_counts, 2, colSums(phylum_counts), "/") * 100

  # 转为长格式
  phylum_df <- as.data.frame(phylum_rel)
  phylum_df$Phylum <- rownames(phylum_df)

  phylum_long <- pivot_longer(phylum_df, cols = -Phylum, names_to = "SampleID", values_to = "Abundance")
  phylum_long <- merge(phylum_long, metadata[, group_col, drop = FALSE],
                       by.x = "SampleID", by.y = "row.names")

  # 将低丰度门合并为 Others
  phylum_mean_abund <- aggregate(Abundance ~ Phylum, data = phylum_long, FUN = mean)
  phylum_mean_abund <- phylum_mean_abund[order(phylum_mean_abund$Abundance, decreasing = TRUE), ]
  top_phyla_stacked <- head(phylum_mean_abund$Phylum, top_n_phyla_stacked)
  phylum_long$Phylum_display <- ifelse(
    phylum_long$Phylum %in% top_phyla_stacked,
    phylum_long$Phylum,
    "Others"
  )
  phylum_long$Phylum_display <- factor(phylum_long$Phylum_display,
    levels = c(top_phyla_stacked, "Others"))

  # 配色：Top 门用 ggsci 色板，Others 用灰色
  n_top <- length(top_phyla_stacked)
  stk_colors <- get_group_colors(max(n_top, 3))[1:n_top]
  stk_colors <- setNames(c(stk_colors, "#BBBBBB"), c(top_phyla_stacked, "Others"))

  # 样本按分组排序
  sample_levels <- unique(as.character(phylum_long$SampleID))
  sample_levels <- sample_levels[order(phylum_long[[group_col]][match(sample_levels, phylum_long$SampleID)])]
  phylum_long$SampleID <- factor(phylum_long$SampleID, levels = sample_levels)

  p_phylum_bar <- ggplot(phylum_long, aes(x = SampleID, y = Abundance, fill = Phylum_display)) +
    geom_bar(stat = "identity", width = min(0.95, 18 / ncol(feature_table))) +
    scale_fill_manual(values = stk_colors) +
    labs(
      title = "Phylum Composition (Relative Abundance)",
      x = "Sample",
      y = "Relative Abundance (%)",
      fill = "Phylum"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
      legend.position = "bottom",
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.4, "cm"),
      panel.grid.major.x = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 4)) +
    facet_grid(as.formula(paste("~", group_col)), scales = "free_x", space = "free")

  ggsave(
    filename = "export/taxa/phylum_stacked_barplot.pdf",
    plot = p_phylum_bar,
    width = max(8, min(ncol(feature_table) * 0.35, 25)),
    height = 6,
    device = cairo_pdf
  )
  cat("       [OK] export/taxa/phylum_stacked_barplot.pdf\n")

  # ============================================================
  # 10. 门水平丰度条形图（聚合）
  # ============================================================

  # 计算每组的平均相对丰度
  phylum_group_mean <- phylum_long %>%
    group_by(!!sym(group_col), Phylum) %>%
    summarise(MeanAbundance = mean(Abundance), SD = sd(Abundance), .groups = "drop")

  # 计算总体平均
  phylum_overall <- phylum_long %>%
    group_by(Phylum) %>%
    summarise(MeanAbundance = mean(Abundance), SD = sd(Abundance), .groups = "drop") %>%
    arrange(desc(MeanAbundance))
  phylum_overall[[group_col]] <- "Overall"

  # 取 Top N，其余合并为 Others
  top_n_phyla <- 7
  top_phyla <- phylum_overall$Phylum[1:min(top_n_phyla, nrow(phylum_overall))]
  top_phyla <- top_phyla[top_phyla != "Unassigned"]
  top_phyla <- c(top_phyla, "Unassigned")

  phylum_group_mean$Phylum_grouped <- ifelse(
    phylum_group_mean$Phylum %in% top_phyla,
    phylum_group_mean$Phylum,
    "Others"
  )

  phylum_grouped <- phylum_group_mean %>%
    group_by(!!sym(group_col), Phylum_grouped) %>%
    summarise(MeanAbundance = sum(MeanAbundance), .groups = "drop")

  # 按总体相对丰度从高到低排序，Others 恒在底部
  phylum_grouped$Phylum_grouped <- factor(phylum_grouped$Phylum_grouped,
    levels = c(top_phyla, "Others"))

  phyla_group_levels <- levels(phylum_grouped$Phylum_grouped)
  phyla_group_levels <- phyla_group_levels[phyla_group_levels != "Others"]
  phy_colors <- pal_jco()(10)
  phy_colors <- grDevices::colorRampPalette(phy_colors)(length(phyla_group_levels))
  names(phy_colors) <- phyla_group_levels
  all_colors <- c(phy_colors, c("Others" = "#A9A9A9"))

  p_phylum_abundance <- ggplot(phylum_grouped,
                               aes(x = !!sym(group_col), y = MeanAbundance,
                                   fill = Phylum_grouped)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = all_colors) +
    labs(
      title = paste0("Top ", top_n_phyla, " Phyla (Mean Relative Abundance)"),
      x = group_col,
      y = "Mean Relative Abundance (%)",
      fill = "Phylum"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )

  ggsave(
    filename = "export/taxa/phylum_abundance_barchart.pdf",
    plot = p_phylum_abundance,
    width = 7,
    height = 6,
    device = cairo_pdf
  )
  cat("       [OK] export/taxa/phylum_abundance_barchart.pdf\n")

  # ============================================================
  # 11. 属水平热图
  # ============================================================

  cat("       正在生成属水平热图...\n")

  # 将 ASV 映射到属
  asv_genus <- setNames(taxonomy$Genus, taxonomy$FeatureID)
  asv_genus <- asv_genus[common_asvs]

  # 按属汇总
  genus_levels <- unique(asv_genus)
  genus_counts <- matrix(0, nrow = length(genus_levels), ncol = ncol(ft_phylum))
  rownames(genus_counts) <- genus_levels
  colnames(genus_counts) <- colnames(ft_phylum)

  for (i in seq_along(asv_genus)) {
    genus_counts[asv_genus[i], ] <- genus_counts[asv_genus[i], ] + ft_phylum[i, ]
  }

  # 移除 Unassigned
  if ("Unassigned" %in% rownames(genus_counts)) {
    genus_counts <- genus_counts[rownames(genus_counts) != "Unassigned", , drop = FALSE]
  }

  # 选取 Top N 属
  genus_total_abundance <- rowSums(genus_counts)
  genus_sorted <- sort(genus_total_abundance, decreasing = TRUE)
  top_genera <- names(genus_sorted)[1:min(top_n_genera, length(genus_sorted))]

  if (length(top_genera) >= 3) {
    genus_heatmap_data <- genus_counts[top_genera, , drop = FALSE]

    # Log10 转换（加伪计数避免 log(0)）
    genus_heatmap_log <- log10(genus_heatmap_data + 1)

    # 样本分组注释
    sample_groups <- metadata[colnames(genus_heatmap_log), group_col, drop = FALSE]
    annotation_colors <- list()
    annotation_colors[[group_col]] <- setNames(
      get_group_colors(nlevels(sample_groups[[group_col]])),
      levels(sample_groups[[group_col]])
    )

    cairo_pdf("export/heatmap/genus_heatmap.pdf",
        width = max(8, ncol(genus_heatmap_log) * 0.4),
        height = max(6, nrow(genus_heatmap_log) * 0.35))
    pheatmap(
      genus_heatmap_log,
      annotation_col = sample_groups,
      annotation_colors = annotation_colors,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 8,
      fontsize_col = 7,
      color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
      main = paste0("Top ", length(top_genera), " Genera (log10 Abundance)"),
      border_color = NA
    )
    dev.off()
    cat("       [OK] export/heatmap/genus_heatmap.pdf\n")
  } else {
    cat("       [跳过] 已注释的属太少，无法生成热图(< 3)\n")
  }

  # ============================================================
  # 12. 保存属水平丰度表
  # ============================================================

  genus_abundance_df <- as.data.frame(genus_counts)
  genus_abundance_df$Genus <- rownames(genus_abundance_df)
  genus_abundance_df <- genus_abundance_df[, c("Genus", setdiff(colnames(genus_abundance_df), "Genus"))]

  write.table(
    genus_abundance_df,
    file = "export/feature_tables/genus_abundance.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  cat("       [OK] export/feature_tables/genus_abundance.tsv\n")
} else {
  cat("       [跳过] 特征表与分类学注释之间无共同 ASV\n")
}

# ============================================================
# 13. FAPROTAX 输入文件
# ============================================================

cat("\n[6/8] 准备 FAPROTAX 输入文件...\n")

# FAPROTAX 需要稀疏特征表 (BIOM 格式) 和物种注释文件
# 直接复制 rarefied_table.biom 和 taxonomy.tsv 到 export/faprotax/ 目录

if (file.exists(rarefied_biom_file)) {
  file.copy(rarefied_biom_file, "export/faprotax/rarefied_table.biom", overwrite = TRUE)
  cat("       [OK] export/faprotax/rarefied_table.biom\n")
} else {
  cat("       [跳过] 未找到稀疏表 (", rarefied_biom_file, ")\n")
}

if (file.exists(taxonomy_file)) {
  tax_raw <- read.table(taxonomy_file, sep = "\t", header = TRUE,
    check.names = FALSE, stringsAsFactors = FALSE)
  tax_faprotax <- tax_raw[, c("Feature ID", "Taxon")]
  # 以二进制模式写入，避免 Windows CRLF 换行符
  con <- file("export/faprotax/taxonomy.tsv", "wb")
  write.table(tax_faprotax, file = con,
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  close(con)
  cat("       [OK] export/faprotax/taxonomy.tsv (已去除表头和 Confidence 列)\n")
} else {
  cat("       [跳过] 未找到物种注释文件 (", taxonomy_file, ")\n")
}

# ============================================================
# 14. PICRUSt2 输入文件准备
# ============================================================

cat("\n[7/8] 准备 PICRUSt2 输入文件...\n")

# 复制处理后的 feature-table.tsv（PICRUSt2 输入）
# 去除第一行 "# Constructed from biom file"，并将 "#OTU ID" 改为 "#OTUID"
if (file.exists(feature_table_file)) {
  ft_lines <- readLines(feature_table_file)
  ft_lines <- ft_lines[-1]
  ft_lines[1] <- sub("#OTU ID", "#OTUID", ft_lines[1])
  # 以二进制模式写入，避免 Windows CRLF 换行符
  con <- file("export/picrust2/feature-table.tsv", "wb")
  writeLines(ft_lines, con)
  close(con)
  cat("       [OK] export/picrust2/feature-table.tsv\n")
} else {
  cat("       [跳过] 未找到特征表\n")
}

# 复制代表序列
if (file.exists(dna_sequences_file)) {
  file.copy(dna_sequences_file, "export/picrust2/dna-sequences.fasta", overwrite = TRUE)
  cat("       [OK] export/picrust2/dna-sequences.fasta\n")
} else {
  cat("       [跳过] dna-sequences.fasta 未找到\n")
}


# ============================================================
# 15. 完成总结
# ============================================================

cat("\n[8/8] 流程执行完毕！\n")
cat("========================================\n")
cat("输出文件汇总:\n")
cat("========================================\n")
cat("  export/alpha/alpha_diversity_boxplot.pdf           — α 多样性箱线图（4 指标 × 星号标注）\n")
cat("  export/alpha/rarefaction_curves.pdf               — 稀释曲线\n")
cat("  export/beta/beta_diversity_pcoa_bray_curtis.pdf   — β 多样性 Bray-Curtis PCoA（PERMANOVA 标注）\n")
cat("  export/beta/beta_diversity_pcoa_jaccard.pdf       — β 多样性 Jaccard PCoA（PERMANOVA 标注）\n")
cat("  export/taxa/phylum_stacked_barplot.pdf              — 门水平堆叠柱状图\n")
cat("  export/taxa/phylum_abundance_barchart.pdf           — 门水平丰度条形图\n")
cat("  export/heatmap/genus_heatmap.pdf                    — 属水平热图\n")
cat("  export/faprotax/rarefied_table.biom                 — FAPROTAX 输入 (稀疏特征表)\n")
cat("  export/faprotax/taxonomy.tsv                        — FAPROTAX 输入 (物种注释)\n")
cat("  export/picrust2/feature-table.tsv                   — PICRUSt2 特征表\n")
cat("  export/picrust2/dna-sequences.fasta                 — PICRUSt2 代表序列\n")
cat("  export/feature_tables/taxonomy_processed.tsv        — 分类学层级拆分表\n")
cat("  export/feature_tables/genus_abundance.tsv           — 属水平丰度表\n")
cat("  export/feature_tables/alpha_diversity_metrics.tsv   — α 多样性指标表\n")
cat("  export/feature_tables/feature_table_with_taxonomy.tsv — 附分类学特征表\n")
cat("========================================\n")
cat("脚本运行成功！\n")
