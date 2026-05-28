#!/usr/bin/env Rscript
#
# FAPROTAX 功能注释分析 + 可视化工作流
# ==============================================================================
# 覆盖: 全部 FAPROTAX 功能（自动归类为生态循环）
#       全局概览 → 循环级差异 → Beta多样性 → 高CV热图 → OTU贡献追溯
# 用法: Rscript faprotax_visualization.R [选项]
# 数据: results/export/faprotax/ 下需有:
#         faprotax.txt          - 功能 × 样品丰度矩阵
#         faprotax_report.mat   - OTU × 功能 0/1 矩阵（用于 OTU 贡献追溯）
#         taxonomy.tsv          - OTU 分类学注释
# 输出: PDF + TSV → results/export/faprotax/faprotax_visualization/
# 帮助: Rscript faprotax_visualization.R --help / -h
# ==============================================================================

# ==============================================================================
# 0. 参数设置
# ==============================================================================

# ---- CLI 参数解析 ----
parse_cli_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  params <- list()

  short_map <- c(
    m = "metadata",
    f = "faprotax_dir",
    o = "output_dir",
    s = "sample_col",
    g = "group_col",
    "1" = "run_global_view",
    "2" = "run_diff_cycle",
    "3" = "run_beta",
    "4" = "run_heatmap",
    "5" = "run_otu_trace",
    t = "n_top_stack",
    T = "n_top_heatmap",
    p = "padj_cutoff",
    a = "min_abund_filter",
    u = "n_otu_func",
    x = "n_otu_tax",
    w = "fig_w_base",
    z = "fig_h_base",
    c = "group_palette"
  )

  i <- 1
  while (i <= length(args)) {
    if (args[i] %in% c("--help", "-h")) {
      cat("FAPROTAX 功能注释可视化脚本\n")
      cat("用法: Rscript faprotax_visualization.R [选项]\n")
      cat("所有参数均有默认值，仅需指定需覆盖的参数。\n\n")
      cat("路径参数:\n")
      cat("  -m, --metadata=<file>        元数据文件路径 (默认: metadata.txt)\n")
      cat("  -f, --faprotax-dir=<dir>     FAPROTAX 输入目录 (默认: results/export/faprotax)\n")
      cat("  -o, --output-dir=<dir>       输出目录 (默认: results/export/faprotax/faprotax_visualization)\n\n")
      cat("数据参数:\n")
      cat("  -s, --sample-col=<col>       样本ID列名 (默认: SampleID)\n")
      cat("  -g, --group-col=<col>        分组列名 (默认: Group)\n\n")
      cat("分析模块开关 (TRUE/FALSE):\n")
      cat("  -1, --run-global-view=<T/F>  全局功能组成概览 (默认: TRUE)\n")
      cat("  -2, --run-diff-cycle=<T/F>   生态循环级差异分析 (默认: TRUE)\n")
      cat("  -3, --run-beta=<T/F>         Beta 多样性 PCoA (默认: TRUE)\n")
      cat("  -4, --run-heatmap=<T/F>      高变异功能热图 (默认: TRUE)\n")
      cat("  -5, --run-otu-trace=<T/F>    OTU 贡献追溯 (默认: TRUE)\n\n")
      cat("可视化参数:\n")
      cat("  -t, --n-top-stack=<n>        全局堆叠图展示生态循环数 (默认: 12)\n")
      cat("  -T, --n-top-heatmap=<n>      热图展示高变异功能数 (默认: 40)\n")
      cat("  -p, --padj-cutoff=<n>        差异显著性阈值 (默认: 0.05)\n")
      cat("  -a, --min-abund-filter=<n>   功能相对丰度均值过滤阈值 (默认: 0.01)\n")
      cat("  -u, --n-otu-func=<n>         OTU 追溯展示功能数 (默认: 6)\n")
      cat("  -x, --n-otu-tax=<n>          每个功能展示 OTU 数 (默认: 15)\n")
      cat("  -w, --fig-w-base=<n>         图片宽度基准 (默认: 7)\n")
      cat("  -z, --fig-h-base=<n>         图片高度基准 (默认: 5)\n")
      cat("  -c, --group-palette=<name>   分组配色方案 (默认: npg)\n")
      quit(save = "no", status = 0)
    } else if (grepl("^--", args[i])) {
      kv <- sub("^--", "", args[i])
      if (grepl("=", kv, fixed = TRUE)) {
        parts <- strsplit(kv, "=", fixed = TRUE)[[1]]
        name <- gsub("-", "_", parts[1])
        params[[name]] <- parts[2]
      }
    } else if (grepl("^-", args[i])) {
      short <- sub("^-", "", args[i])
      flag <- substr(short, 1, 1)
      if (flag != "h") {
        if (grepl("=", short, fixed = TRUE)) {
          val <- sub("^[^=]+=", "", short)
          full <- short_map[flag]
          if (!is.na(full)) params[[full]] <- val
        } else if (nchar(short) == 1) {
          full <- short_map[flag]
          if (!is.na(full) && i < length(args)) {
            i <- i + 1
            params[[full]] <- args[i]
          }
        }
      }
    }
    i <- i + 1
  }
  return(params)
}

cli <- parse_cli_args()

# 参数辅助函数：取 CLI 值（若有）或默认值
use_param <- function(cli_val, default) {
  if (is.null(cli_val)) default else cli_val
}
as_flag <- function(x) {
  if (is.logical(x)) return(x)
  toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES")
}

# ---- 路径参数 ----
setwd(".")
metadata_file     <- use_param(cli[["metadata"]], "metadata.txt")
faprotax_dir      <- use_param(cli[["faprotax_dir"]], "results/export/faprotax")
output_dir        <- use_param(cli[["output_dir"]], "results/export/faprotax/faprotax_visualization")

# ---- 数据参数 ----
sample_col        <- use_param(cli[["sample_col"]], "SampleID")
group_col         <- use_param(cli[["group_col"]], "Group")

# ---- 分析模块开关 ----
run_global_view   <- as_flag(use_param(cli[["run_global_view"]], TRUE))
run_diff_cycle    <- as_flag(use_param(cli[["run_diff_cycle"]], TRUE))
run_beta          <- as_flag(use_param(cli[["run_beta"]], TRUE))
run_heatmap       <- as_flag(use_param(cli[["run_heatmap"]], TRUE))
run_otu_trace     <- as_flag(use_param(cli[["run_otu_trace"]], TRUE))

# ---- 可视化参数 ----
n_top_stack       <- as.numeric(use_param(cli[["n_top_stack"]], 12))
n_top_heatmap     <- as.numeric(use_param(cli[["n_top_heatmap"]], 40))
padj_cutoff       <- as.numeric(use_param(cli[["padj_cutoff"]], 0.05))
min_abund_filter  <- as.numeric(use_param(cli[["min_abund_filter"]], 0.01))

# ---- OTU 追溯参数 ----
n_otu_func        <- as.numeric(use_param(cli[["n_otu_func"]], 6))
n_otu_tax         <- as.numeric(use_param(cli[["n_otu_tax"]], 15))

# ---- 配色 ----
group_palette     <- use_param(cli[["group_palette"]], "npg")

# ---- 图片尺寸 ----
fig_w_base        <- as.numeric(use_param(cli[["fig_w_base"]], 7))
fig_h_base        <- as.numeric(use_param(cli[["fig_h_base"]], 5))

# ==============================================================================
# 1. 包加载与工具函数
# ==============================================================================
cat("[1/8] 加载 R 包...\n")

required_packages <- c(
  "ggplot2", "tidyr", "dplyr", "tibble", "readr",
  "vegan", "ape", "pheatmap", "RColorBrewer", "grDevices",
  "ggrepel", "ggsci", "viridis", "reshape2", "FSA", "grid",
  "scales"
)

check_packages <- function(pkgs) {
  missing <- c()
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    }
  }
  if (length(missing) > 0) {
    cat("  安装缺失的 CRAN 包:", paste(missing, collapse = ", "), "\n")
    install.packages(missing, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN",
                     quiet = TRUE, Ncpus = parallel::detectCores())
  }
  invisible(lapply(pkgs, library, character.only = TRUE, quietly = TRUE))
}

check_packages(required_packages)
cat("  所有包加载成功。\n\n")

# ---- 工具函数 ----

# 分组配色
get_group_colors <- function(n, palette = "npg") {
  if (palette == "npg")       return(rep(pal_npg()(10), length.out = n))
  if (palette == "aaas")      return(rep(pal_aaas()(10), length.out = n))
  if (palette == "nejm")      return(rep(pal_nejm()(8),  length.out = n))
  if (palette == "lancet")    return(rep(pal_lancet()(9), length.out = n))
  if (palette == "jco")       return(rep(pal_jco()(10),  length.out = n))
  if (palette == "viridis")   return(viridis(n))
  if (palette == "Set1")      return(brewer.pal(max(3, min(n, 9)), "Set1")[1:n])
  if (palette == "Set2")      return(brewer.pal(max(3, min(n, 8)), "Set2")[1:n])
  if (palette == "Set3")      return(brewer.pal(max(3, min(n, 12)), "Set3")[1:n])
  return(brewer.pal(max(3, min(n, 8)), "Set2")[1:n])
}

# 相对丰度转换（百分比）
normalize_abundance <- function(mat) {
  t(t(mat) / colSums(mat, na.rm = TRUE)) * 100
}

# Kruskal-Wallis + Dunn post-hoc
multi_group_test <- function(df, val_col, grp_col, p_adj_method = "bh") {
  df <- df[!is.na(df[[val_col]]) & !is.na(df[[grp_col]]), ]
  df[[grp_col]] <- as.factor(df[[grp_col]])

  kw <- kruskal.test(df[[val_col]] ~ df[[grp_col]])
  kw_p <- signif(kw$p.value, 4)

  n_grps <- length(unique(na.omit(df[[grp_col]])))
  if (n_grps >= 3) {
    dunn <- dunnTest(df[[val_col]] ~ df[[grp_col]], method = p_adj_method)
    dunn_df <- as.data.frame(dunn$res)
    dunn_df$Comparison <- gsub(" ", "", dunn_df$Comparison)
    dunn_df$sig_label <- ifelse(dunn_df$P.adj < 0.001, "***",
                         ifelse(dunn_df$P.adj < 0.01,  "**",
                         ifelse(dunn_df$P.adj < 0.05,  "*", "ns")))
    dunn_df$P.adj <- signif(dunn_df$P.adj, 4)
  } else {
    grp_vals <- levels(droplevels(df[[grp_col]]))
    dunn_df <- data.frame(
      Comparison = paste(grp_vals, collapse = " - "),
      P.adj = kw_p,
      sig_label = ifelse(kw_p < 0.05, ifelse(kw_p < 0.01, ifelse(kw_p < 0.001, "***", "**"), "*"), "ns"),
      stringsAsFactors = FALSE
    )
  }
  list(kw_p = kw_p, dunn = dunn_df)
}

# Wilcoxon 检验（两组）
wilcoxon_test <- function(df, val_col, grp_col) {
  df <- df[!is.na(df[[val_col]]) & !is.na(df[[grp_col]]), ]
  grps <- unique(as.character(df[[grp_col]]))
  if (length(grps) != 2) return(list(p_value = NA, stat = NA))
  x <- df[[val_col]][df[[grp_col]] == grps[1]]
  y <- df[[val_col]][df[[grp_col]] == grps[2]]
  if (length(x) < 2 || length(y) < 2 || sd(x) + sd(y) == 0) {
    return(list(p_value = NA, stat = NA))
  }
  wt <- wilcox.test(x, y, exact = FALSE)
  list(p_value = signif(wt$p.value, 4), stat = round(wt$statistic, 1))
}

# 导出 TSV
save_tsv_file <- function(data, file_path) {
  write.table(data, file_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# 合并低丰度项为 Others
merge_others <- function(rel_mat, n_keep, metadata, grp_col) {
  groups <- unique(as.character(metadata[[grp_col]]))
  grp_mean_list <- lapply(groups, function(g) {
    samps <- intersect(colnames(rel_mat), rownames(metadata[metadata[[grp_col]] == g, ]))
    if (length(samps) > 0) rowMeans(rel_mat[, samps, drop = FALSE]) else rep(0, nrow(rel_mat))
  })
  grp_means <- do.call(cbind, grp_mean_list)
  overall_rank <- order(rowMeans(grp_means), decreasing = TRUE)
  top_idx <- overall_rank[1:min(n_keep, length(overall_rank))]
  other_idx <- overall_rank[-(1:min(n_keep, length(overall_rank)))]
  if (length(other_idx) == 0) return(rel_mat[top_idx, , drop = FALSE])
  top_mat <- rel_mat[top_idx, , drop = FALSE]
  other_vec <- colSums(rel_mat[other_idx, , drop = FALSE])
  rbind(top_mat, Others = other_vec)
}

# 两两 PERMANOVA
pairwise_permanova <- function(dist_obj, group_vec, n_perm = 999) {
  grp_names <- unique(as.character(group_vec))
  results <- data.frame(Comparison = character(), R2 = numeric(),
                        p_value = numeric(), p_adj = numeric(), stringsAsFactors = FALSE)
  for (i in 1:(length(grp_names) - 1)) {
    for (j in (i + 1):length(grp_names)) {
      sel <- group_vec %in% c(grp_names[i], grp_names[j])
      sub_dist <- as.dist(as.matrix(dist_obj)[sel, sel])
      sub_meta <- data.frame(Group = factor(group_vec[sel]), row.names = names(group_vec)[sel])
      set.seed(123)
      perm <- adonis2(sub_dist ~ Group, data = sub_meta, permutations = n_perm)
      results <- rbind(results, data.frame(
        Comparison = paste0(grp_names[i], " vs ", grp_names[j]),
        R2 = round(perm$R2[1], 4),
        p_value = perm$`Pr(>F)`[1],
        p_adj = NA,
        stringsAsFactors = FALSE
      ))
    }
  }
  results$p_adj <- p.adjust(results$p_value, method = "BH")
  results$p_value <- signif(results$p_value, 3)
  results$p_adj <- signif(results$p_adj, 3)
  results
}

# ==============================================================================
# 1b. FAPROTAX 功能分类词典
# ==============================================================================

# 将 FAPROTAX 功能名称归类为生态循环
classify_function <- function(func_names) {
  # 分类规则: 关键词 + 白名单
  classification <- list(

    methane = c(
      "methanotrophy",
      "methanol_oxidation",
      "methylotrophy",
      "methanogenesis",
      "acetoclastic_methanogenesis",
      "methanogenesis_by_disproportionation_of_methyl_groups",
      "methanogenesis_using_formate",
      "methanogenesis_by_CO2_reduction_with_H2",
      "methanogenesis_by_reduction_of_methyl_compounds_with_H2",
      "hydrogenotrophic_methanogenesis"
    ),

    nitrogen = c(
      "nitrification",
      "aerobic_ammonia_oxidation",
      "aerobic_nitrite_oxidation",
      "denitrification",
      "nitrate_denitrification",
      "nitrite_denitrification",
      "nitrous_oxide_denitrification",
      "nitrogen_fixation",
      "nitrate_reduction",
      "nitrate_respiration",
      "nitrite_respiration",
      "nitrogen_respiration",
      "nitrate_ammonification",
      "nitrite_ammonification",
      "anammox"
    ),

    sulfur = c(
      "sulfate_respiration",
      "sulfur_respiration",
      "sulfite_respiration",
      "thiosulfate_respiration",
      "respiration_of_sulfur_compounds",
      "dark_sulfite_oxidation",
      "dark_sulfide_oxidation",
      "dark_sulfur_oxidation",
      "dark_thiosulfate_oxidation",
      "dark_oxidation_of_sulfur_compounds"
    ),

    iron_manganese = c(
      "iron_respiration",
      "dark_iron_oxidation",
      "manganese_oxidation",
      "manganese_respiration"
    ),

    carbon_degradation = c(
      "cellulolysis",
      "xylanolysis",
      "chitinolysis",
      "ligninolysis",
      "aromatic_compound_degradation",
      "aromatic_hydrocarbon_degradation",
      "aliphatic_non_methane_hydrocarbon_degradation",
      "hydrocarbon_degradation",
      "plastic_degradation"
    ),

    fermentation = c(
      "fermentation"
    ),

    heterotrophy = c(
      "aerobic_chemoheterotrophy",
      "chemoheterotrophy"
    ),

    hydrogen = c(
      "knallgas_bacteria",
      "dark_hydrogen_oxidation",
      "reductive_acetogenesis"
    ),

    phototrophy = c(
      "phototrophy",
      "photoautotrophy",
      "oxygenic_photoautotrophy",
      "anoxygenic_photoautotrophy",
      "anoxygenic_photoautotrophy_H2_oxidizing",
      "anoxygenic_photoautotrophy_S_oxidizing",
      "anoxygenic_photoautotrophy_Fe_oxidizing",
      "photoheterotrophy",
      "aerobic_anoxygenic_phototrophy",
      "photosynthetic_cyanobacteria",
      "nonphotosynthetic_cyanobacteria",
      "chloroplasts"
    ),

    host_associated = c(
      "human_pathogens_all",
      "human_pathogens_septicemia",
      "human_pathogens_pneumonia",
      "human_pathogens_nosocomia",
      "human_pathogens_meningitis",
      "human_pathogens_gastroenteritis",
      "human_pathogens_diarrhea",
      "human_associated",
      "human_gut",
      "mammal_gut",
      "animal_parasites_or_symbionts",
      "invertebrate_parasites",
      "fish_parasites",
      "plant_pathogen",
      "intracellular_parasites",
      "predatory_or_exoparasitic"
    ),

    other = c(
      "ureolysis",
      "fumarate_respiration",
      "chlorate_reducers",
      "arsenate_detoxification",
      "arsenate_respiration",
      "dissimilatory_arsenate_reduction",
      "arsenite_oxidation_detoxification",
      "arsenite_oxidation_energy_yielding",
      "dissimilatory_arsenite_oxidation",
      "oil_bioremediation"
    )
  )

  # 反向查找
  cycle_names <- c(
    "methane"            = "Methane / Methylotrophy",
    "nitrogen"           = "Nitrogen",
    "sulfur"             = "Sulfur",
    "iron_manganese"     = "Iron / Manganese",
    "carbon_degradation" = "Carbon Degradation",
    "fermentation"       = "Fermentation",
    "heterotrophy"       = "Heterotrophy",
    "hydrogen"           = "Hydrogen",
    "phototrophy"        = "Phototrophy",
    "host_associated"    = "Host-associated",
    "other"              = "Other"
  )

  # 构建 function → cycle 映射表
  func_to_cycle <- list()
  for (key in names(classification)) {
    for (f in classification[[key]]) {
      func_to_cycle[[f]] <- cycle_names[key]
    }
  }

  # 为每个输入名返回分类结果
  sapply(func_names, function(f) {
    if (f %in% names(func_to_cycle)) func_to_cycle[[f]] else "Unclassified"
  }, USE.NAMES = FALSE)
}

# 生态循环配色方案（跨模块复用）
get_cycle_colors <- function(n) {
  # 使用 Set3 作为基础，融合 npg 扩展
  base <- brewer.pal(min(12, max(3, n)), "Set3")
  if (n <= 12) return(base[1:n])
  c(base, viridis(n - 12))
}

# ==============================================================================
# 2. 数据读取与预处理
# ==============================================================================
cat("[2/8] 读取数据...\n")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 2a. 读取元数据 ----
metadata <- read.delim(metadata_file, sep = "\t", stringsAsFactors = FALSE)
stopifnot(sample_col %in% colnames(metadata), group_col %in% colnames(metadata))
metadata <- metadata[!is.na(metadata[[sample_col]]) & metadata[[sample_col]] != "", ]
rownames(metadata) <- metadata[[sample_col]]
metadata[[group_col]] <- factor(metadata[[group_col]])
groups  <- levels(metadata[[group_col]])
n_group <- length(groups)
cat(sprintf("  检测到 %d 个分组: %s\n", n_group, paste(groups, collapse = ", ")))

group_colors <- get_group_colors(n_group, group_palette)
names(group_colors) <- groups

sample_annot <- data.frame(Group = metadata[[group_col]], row.names = rownames(metadata))
ann_colors <- list(Group = group_colors)

# ---- 2b. 读取 FAPROTAX 功能丰度表 ----
faprotax_file <- file.path(faprotax_dir, "faprotax.txt")
if (!file.exists(faprotax_file)) stop("未找到 faprotax.txt: ", faprotax_file)

faprotax <- read.delim(faprotax_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
# 第一列名为 "group"，存储功能名
rownames(faprotax) <- faprotax$group
faprotax <- faprotax[, -1, drop = FALSE]  # 移除 group 列

# 与 metadata 对齐样品
common <- intersect(colnames(faprotax), rownames(metadata))
if (length(common) == 0) stop("样品名与 metadata 不匹配，请检查 SampleID 列")
faprotax <- faprotax[, common, drop = FALSE]
metadata <- metadata[common, , drop = FALSE]
cat(sprintf("  功能丰度表: %d 个功能 × %d 个样品\n", nrow(faprotax), ncol(faprotax)))

# 过滤所有行和均为 0 的功能
zero_func <- rowSums(faprotax) == 0
if (any(zero_func)) {
  cat(sprintf("  过滤 %d 个在所有样品中均为 0 的功能\n", sum(zero_func)))
  faprotax <- faprotax[!zero_func, , drop = FALSE]
}

# ---- 2c. 功能分类 ----
func_cycles <- classify_function(rownames(faprotax))
cat(sprintf("  功能归类为 %d 个生态循环\n", length(unique(func_cycles))))

# ---- 2d. 相对丰度标准化 ----
faprotax_rel <- normalize_abundance(faprotax)

# ---- 2e. 读取 OTU 贡献矩阵（可选） ----
mat_file <- file.path(faprotax_dir, "faprotax_report.mat")
otu_mat <- NULL
if (file.exists(mat_file)) {
  otu_mat <- read.delim(mat_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  cat(sprintf("  OTU 功能矩阵: %d 个 OTU × %d 个功能\n", nrow(otu_mat), ncol(otu_mat) - 1))
} else {
  cat("  未找到 faprotax_report.mat，OTU 追溯将不可用\n")
  if (run_otu_trace) {
    cat("  设置 run_otu_trace <- FALSE\n")
    run_otu_trace <- FALSE
  }
}

# ---- 2f. 读取 OTU 分类学注释（可选） ----
tax_file <- file.path(faprotax_dir, "taxonomy.tsv")
taxonomy <- NULL
if (file.exists(tax_file)) {
  taxonomy <- read.delim(tax_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(taxonomy) <- c("OTUID", "Taxonomy")
  cat(sprintf("  分类学注释: %d 个 OTU\n", nrow(taxonomy)))
} else {
  cat("  未找到 taxonomy.tsv，OTU 追溯将不可用\n")
  if (run_otu_trace) run_otu_trace <- FALSE
}

cat("  数据读取完成。\n\n")

# ==============================================================================
# 3. 全局功能组成概览
# ==============================================================================
cat("[3/8] 全局功能组成概览...\n")

if (run_global_view && nrow(faprotax_rel) > 0) {

  # ---- 3a. 按生态循环汇总 ----
  cycle_agg <- rowsum(faprotax_rel, group = func_cycles)
  cat(sprintf("  生态循环汇总: %d 个循环 × %d 个样品\n", nrow(cycle_agg), ncol(cycle_agg)))

  # 合并低丰度循环
  cycle_merged <- merge_others(cycle_agg, n_top_stack, metadata, group_col)

  # 转为长格式
  cycle_melt <- reshape2::melt(as.matrix(cycle_merged))
  colnames(cycle_melt) <- c("Cycle", "SampleID", "Abundance")
  cycle_melt <- merge(cycle_melt, metadata, by.x = "SampleID", by.y = sample_col)
  cycle_melt$SampleID <- factor(cycle_melt$SampleID, levels = colnames(cycle_merged))

  # 堆叠柱状图
  cycle_colors <- get_cycle_colors(nrow(cycle_merged))
  names(cycle_colors) <- rownames(cycle_merged)

  p <- ggplot(cycle_melt, aes(x = SampleID, y = Abundance, fill = Cycle)) +
    geom_bar(stat = "identity", position = "stack", width = 0.85) +
    facet_grid(. ~ get(group_col), scales = "free_x", space = "free") +
    scale_fill_manual(values = cycle_colors) +
    labs(title = "FAPROTAX Ecological Cycle Composition",
         x = "", y = "Relative Abundance (%)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
  ggsave(file.path(output_dir, "01_Ecological_Cycle_Stackplot.pdf"), p,
         width = max(10, ncol(faprotax) * 0.35), height = 6)
  save_tsv_file(as.data.frame(cycle_agg),
                file.path(output_dir, "01_Cycle_Abundance.txt"))
  cat("  01_Ecological_Cycle_Stackplot.pdf 已生成\n")

  # ---- 3b. 功能丰度排序柱状图（Top N） ----
  func_means <- rowMeans(faprotax_rel)
  top_func <- names(sort(func_means, decreasing = TRUE))[1:min(n_top_heatmap, length(func_means))]
  top_df <- data.frame(
    Function = factor(top_func, levels = rev(top_func)),
    MeanRelAbund = func_means[top_func],
    Cycle = func_cycles[top_func]
  )

  p <- ggplot(top_df, aes(x = MeanRelAbund, y = Function, fill = Cycle)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = cycle_colors[intersect(names(cycle_colors), unique(top_df$Cycle))]) +
    labs(title = paste0("Top ", length(top_func), " Functions by Mean Relative Abundance"),
         x = "Mean Relative Abundance (%)", y = "") +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 7)
    )
  ggsave(file.path(output_dir, "02_Top_Functions_Barplot.pdf"), p, width = 8, height = max(6, length(top_func) * 0.2))
  save_tsv_file(top_df, file.path(output_dir, "02_Top_Functions_Abundance.txt"))
  cat("  02_Top_Functions_Barplot.pdf 已生成\n")

} else {
  cat("  跳过：全局概览已禁用或无数据\n")
}

# ==============================================================================
# 4. 生态循环级差异分析
# ==============================================================================
cat("[4/8] 生态循环差异分析...\n")

if (run_diff_cycle && nrow(faprotax_rel) > 0) {

  # ---- 准备绘图数据（长格式 + 分类 + 分组） ----
  long_df <- reshape2::melt(as.matrix(faprotax_rel))
  colnames(long_df) <- c("Function", "SampleID", "Abundance")
  long_df$Cycle <- func_cycles[long_df$Function]
  long_df <- merge(long_df, metadata[, c(sample_col, group_col), drop = FALSE],
                   by.x = "SampleID", by.y = sample_col)

  # 收集所有检验结果
  all_test_results <- data.frame()

  # 获取有数据的循环
  active_cycles <- unique(long_df$Cycle)
  active_cycles <- intersect(active_cycles, names(which(table(func_cycles) > 0)))

  # ---- 逐循环绘图 ----
  for (cyc in active_cycles) {
    cyc_df <- long_df[long_df$Cycle == cyc, ]
    # 替换零值为极小值（log10 变换需要）
    cyc_df$Abundance <- ifelse(cyc_df$Abundance == 0, 0.001, cyc_df$Abundance)
    cyc_funcs <- unique(cyc_df$Function)
    n_func <- length(cyc_funcs)

    # 跳过只有 0-1 个功能的循环
    if (n_func < 1) next

    # ---- 4a. 逐功能统计检验 ----
    if (n_group == 2) {
      # 两组用 Wilcoxon
      test_list <- lapply(cyc_funcs, function(f) {
        sub <- cyc_df[cyc_df$Function == f, ]
        wt <- wilcoxon_test(sub, "Abundance", group_col)
        data.frame(Function = f, p_value = wt$p_value, Method = "Wilcoxon", stringsAsFactors = FALSE)
      })
    } else {
      # 多组用 Kruskal-Wallis
      test_list <- lapply(cyc_funcs, function(f) {
        sub <- cyc_df[cyc_df$Function == f, ]
        kw <- kruskal.test(sub$Abundance ~ as.factor(sub[[group_col]]))
        data.frame(Function = f, p_value = signif(kw$p.value, 4), Method = "Kruskal-Wallis", stringsAsFactors = FALSE)
      })
    }
    cyc_tests <- do.call(rbind, test_list)
    cyc_tests$p_adj <- p.adjust(cyc_tests$p_value, method = "BH")
    cyc_tests$sig_label <- ifelse(is.na(cyc_tests$p_adj), "ns",
      ifelse(cyc_tests$p_adj < 0.001, "***",
      ifelse(cyc_tests$p_adj < 0.01,  "**",
      ifelse(cyc_tests$p_adj < 0.05,  "*", "ns"))))
    cyc_tests$Cycle <- cyc
    all_test_results <- rbind(all_test_results, cyc_tests)

    # 记录显著性统计
    n_sig <- sum(cyc_tests$sig_label %in% c("*", "**", "***"), na.rm = TRUE)
    cat(sprintf("  %s: %d 功能, %d 显著 (BH-p < 0.05)\n", cyc, n_func, n_sig))

    # 无显著差异时跳过箱线图（避免生成无信息量的长图）
    if (n_sig == 0) {
      cat("    无显著差异，跳过箱线图\n")
      next
    }

    # 过滤掉非显著功能，仅保留有差异的功能绘制箱线图
    sig_funcs <- cyc_tests$Function[cyc_tests$sig_label %in% c("*", "**", "***")]
    cyc_df <- cyc_df[cyc_df$Function %in% sig_funcs, ]
    cyc_funcs <- unique(cyc_df$Function)
    n_func <- length(cyc_funcs)
    cat(sprintf("    过滤后保留 %d 个显著功能\n", n_func))

    # ---- 4b. 箱线图 ----
    # 长功能名换行
    cyc_df$Function_label <- gsub("_", "\n", cyc_df$Function)

    # 确定图的高度
    h <- max(4, n_func * 1.8)

    # 分组均值排序（便于阅读）
    func_order <- cyc_df %>%
      group_by(Function) %>%
      summarise(MeanAbund = mean(Abundance), .groups = "drop") %>%
      arrange(MeanAbund)
    cyc_df$Function <- factor(cyc_df$Function, levels = func_order$Function)

    # 构建显著性标签位置
    sig_pos <- cyc_df %>%
      group_by(Function) %>%
      summarise(max_y = max(Abundance, na.rm = TRUE) * 1.1, .groups = "drop")
    sig_pos <- merge(sig_pos, cyc_tests[, c("Function", "sig_label")], by = "Function", all.x = TRUE)
    sig_pos <- sig_pos[!is.na(sig_pos$sig_label) & sig_pos$sig_label != "ns", , drop = FALSE]

    p <- ggplot(cyc_df, aes(x = Function, y = Abundance, fill = .data[[group_col]])) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6, position = position_dodge(0.8)) +
      geom_jitter(aes(color = .data[[group_col]]),
                  position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                  size = 1.2, alpha = 0.6) +
      geom_text(data = sig_pos, aes(x = Function, y = max_y, label = sig_label),
                inherit.aes = FALSE, size = 3.5, vjust = 0.3) +
      scale_fill_manual(values = group_colors) +
      scale_color_manual(values = group_colors) +
      scale_y_log10(limits = c(0.001, NA)) +
      labs(title = paste0(cyc, " (", n_func, " functions)"),
           x = "", y = "Relative Abundance (%) [log10]") +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 8),
        legend.position = "top",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid.minor = element_blank()
      )

    safe_cyc <- gsub("[^a-zA-Z0-9_-]", "_", cyc)
    ggsave(file.path(output_dir, paste0("03_", safe_cyc, "_Boxplot.pdf")), p,
           width = max(6, n_func * 0.8), height = h)
  }

  # ---- 导出全部检验结果 ----
  save_tsv_file(all_test_results, file.path(output_dir, "03_All_Functions_Test_Results.txt"))
  cat("  03_*_Boxplot.pdf 已生成\n")

  # ---- 4c. 差异功能概览火山图 ----
  # 对所有功能计算组间均值差（仅两组时）
  if (n_group == 2) {
    cat("  生成差异功能概览图...\n")
    grp1 <- groups[1]
    grp2 <- groups[2]
    grp1_samps <- rownames(metadata[metadata[[group_col]] == grp1, ])
    grp2_samps <- rownames(metadata[metadata[[group_col]] == grp2, ])

    func_mean_diff <- data.frame(
      Function = rownames(faprotax_rel),
      Mean_grp1 = rowMeans(faprotax_rel[, intersect(colnames(faprotax_rel), grp1_samps), drop = FALSE]),
      Mean_grp2 = rowMeans(faprotax_rel[, intersect(colnames(faprotax_rel), grp2_samps), drop = FALSE]),
      stringsAsFactors = FALSE
    )
    func_mean_diff$log2FC <- log2((func_mean_diff$Mean_grp2 + 0.01) / (func_mean_diff$Mean_grp1 + 0.01))
    func_mean_diff <- merge(func_mean_diff, all_test_results[, c("Function", "p_adj", "sig_label")],
                            by = "Function", all.x = TRUE)
    func_mean_diff$sig_label[is.na(func_mean_diff$sig_label)] <- "ns"
    func_mean_diff$sig <- ifelse(func_mean_diff$sig_label %in% c("*", "**", "***"), "Significant", "Not Significant")

    # 标注 top 差异功能
    sig_funcs <- func_mean_diff[func_mean_diff$sig == "Significant", ]
    n_label <- min(10, nrow(sig_funcs))
    if (n_label > 0) {
      sig_funcs <- sig_funcs[order(abs(sig_funcs$log2FC), decreasing = TRUE), ]
      label_df <- sig_funcs[1:n_label, ]
    } else {
      label_df <- data.frame()
    }

    p <- ggplot(func_mean_diff, aes(x = log2FC, y = -log10(p_adj), color = sig)) +
      geom_point(alpha = 0.7, size = 1.5) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", linewidth = 0.3) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", linewidth = 0.3) +
      scale_color_manual(values = c("Significant" = "firebrick3", "Not Significant" = "#A9A9A9")) +
      labs(title = paste0(grp2, " vs ", grp1, "  |  ", sum(func_mean_diff$sig == "Significant"), " significant"),
           x = "log2(Fold Change)", y = expression(-log[10](adjusted~p))) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 10))

    if (n_label > 0) {
      p <- p + geom_text_repel(data = label_df, aes(label = Function),
                                size = 2.5, box.padding = 0.5, max.overlaps = 15, show.legend = FALSE)
    }

    ggsave(file.path(output_dir, "04_Function_Volcano.pdf"), p, width = 8, height = 6)
    save_tsv_file(func_mean_diff, file.path(output_dir, "04_Function_Diff_Results.txt"))
    cat("  04_Function_Volcano.pdf 已生成\n")
  }

} else {
  cat("  跳过：循环差异分析已禁用或无数据\n")
}

# ==============================================================================
# 5. Beta 多样性（PCoA + PERMANOVA）
# ==============================================================================
cat("[5/8] Beta 多样性分析...\n")

if (run_beta && nrow(faprotax_rel) > 0) {

  # ---- 5a. Bray-Curtis PCoA ----
  bray_dist <- vegdist(t(faprotax_rel), method = "bray")

  pcoa <- cmdscale(bray_dist, k = 2, eig = TRUE)
  pcoa_df <- data.frame(
    SampleID = rownames(pcoa$points),
    PCoA1 = pcoa$points[, 1],
    PCoA2 = pcoa$points[, 2],
    stringsAsFactors = FALSE
  )
  pcoa_df <- merge(pcoa_df, metadata, by.x = "SampleID", by.y = sample_col)

  set.seed(123)
  permanova <- adonis2(bray_dist ~ metadata[[group_col]], permutations = 999)
  perma_r2 <- round(permanova$R2[1], 3)
  perma_p <- permanova$`Pr(>F)`[1]

  # betadisper
  bd <- betadisper(bray_dist, metadata[[group_col]])
  bd_anova <- anova(bd)
  bd_p <- bd_anova$`Pr(>F)`[1]

  # pairwise PERMANOVA
  grp_vec <- setNames(metadata[[group_col]], rownames(metadata))
  pw_perma <- pairwise_permanova(bray_dist, grp_vec)

  eig1 <- round(pcoa$eig[1] / sum(abs(pcoa$eig)) * 100, 1)
  eig2 <- round(pcoa$eig[2] / sum(abs(pcoa$eig)) * 100, 1)

  # 形状向量
  shapes_avail <- c(16, 17, 15, 18, 8, 3, 4, 11, 1, 2)
  shapes <- shapes_avail[1:n_group]
  names(shapes) <- groups

  p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = .data[[group_col]], shape = .data[[group_col]])) +
    geom_point(size = 2.5, alpha = 0.8) +
    stat_ellipse(level = 0.95, linewidth = 0.6) +
    scale_color_manual(values = group_colors) +
    scale_shape_manual(values = shapes) +
    labs(title = paste0("FAPROTAX PCoA (Bray-Curtis)\nPERMANOVA R² = ", perma_r2,
                        ", p = ", signif(perma_p, 3),
                        "  |  betadisper p = ", signif(bd_p, 3)),
         x = paste0("PCoA1 (", eig1, "%)"),
         y = paste0("PCoA2 (", eig2, "%)")) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 10))
  ggsave(file.path(output_dir, "05_FAPROTAX_PCoA.pdf"), p, width = 7, height = 6)
  save_tsv_file(pcoa_df, file.path(output_dir, "05_PCoA_Results.txt"))
  save_tsv_file(as.data.frame(permanova), file.path(output_dir, "05_PERMANOVA.txt"))
  save_tsv_file(pw_perma, file.path(output_dir, "05_Pairwise_PERMANOVA.txt"))
  cat("  05_FAPROTAX_PCoA.pdf 已生成\n")

  # ---- 5b. 按生态循环分别做 PCoA ----
  cat("  逐循环 PCoA...\n")
  active_cycles_pcoa <- unique(func_cycles)
  for (cyc in active_cycles_pcoa) {
    funcs_in_cycle <- names(func_cycles[func_cycles == cyc])
    funcs_in_cycle <- intersect(funcs_in_cycle, rownames(faprotax_rel))

    if (length(funcs_in_cycle) < 3) next  # 至少需要 3 个功能

    cyc_mat <- faprotax_rel[funcs_in_cycle, , drop = FALSE]
    cyc_dist <- vegdist(t(cyc_mat), method = "bray")

    cyc_pcoa <- cmdscale(cyc_dist, k = 2, eig = TRUE)
    cyc_pcoa_df <- data.frame(
      SampleID = rownames(cyc_pcoa$points),
      PCoA1 = cyc_pcoa$points[, 1],
      PCoA2 = cyc_pcoa$points[, 2],
      stringsAsFactors = FALSE
    )
    cyc_pcoa_df <- merge(cyc_pcoa_df, metadata, by.x = "SampleID", by.y = sample_col)

    set.seed(123)
    cyc_perm <- adonis2(cyc_dist ~ metadata[[group_col]], permutations = 999)
    cyc_r2 <- round(cyc_perm$R2[1], 3)
    cyc_p <- cyc_perm$`Pr(>F)`[1]

    ceig1 <- round(cyc_pcoa$eig[1] / sum(abs(cyc_pcoa$eig)) * 100, 1)
    ceig2 <- round(cyc_pcoa$eig[2] / sum(abs(cyc_pcoa$eig)) * 100, 1)

    p <- ggplot(cyc_pcoa_df, aes(x = PCoA1, y = PCoA2, color = .data[[group_col]], shape = .data[[group_col]])) +
      geom_point(size = 2, alpha = 0.8) +
      stat_ellipse(level = 0.95, linewidth = 0.5) +
      scale_color_manual(values = group_colors) +
      scale_shape_manual(values = shapes) +
      labs(title = paste0(cyc, "\nPERMANOVA R² = ", cyc_r2, ", p = ", signif(cyc_p, 3)),
           x = paste0("PCoA1 (", ceig1, "%)"),
           y = paste0("PCoA2 (", ceig2, "%)")) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 9))

    safe_cyc <- gsub("[^a-zA-Z0-9_-]", "_", cyc)
    ggsave(file.path(output_dir, paste0("05_", safe_cyc, "_PCoA.pdf")), p, width = 6, height = 5.5)
  }
  cat("  逐循环 PCoA 完成\n")

} else {
  cat("  跳过：Beta 多样性已禁用或无数据\n")
}

# ==============================================================================
# 6. 高变异功能热图
# ==============================================================================
cat("[6/8] 高变异功能热图...\n")

if (run_heatmap && nrow(faprotax_rel) > 1) {

  # ---- 计算 CV ----
  func_cv <- apply(faprotax_rel, 1, function(x) sd(x) / mean(x))
  func_cv <- func_cv[is.finite(func_cv)]

  if (length(func_cv) > 0) {
    top_cv <- names(sort(func_cv, decreasing = TRUE))[1:min(n_top_heatmap, length(func_cv))]
    heatmap_mat <- faprotax_rel[top_cv, , drop = FALSE]

    # 行缩放（Z-score）
    heatmap_scaled <- t(scale(t(heatmap_mat)))

    # 行为功能添加分类注释
    heatmap_annot_row <- data.frame(
      Cycle = func_cycles[top_cv],
      row.names = top_cv
    )

    # 循环颜色
    uni_cycles <- unique(func_cycles[top_cv])
    cycle_cmap <- get_cycle_colors(length(uni_cycles))
    names(cycle_cmap) <- uni_cycles

    ann_colors_hm <- list(
      Group = group_colors,
      Cycle = cycle_cmap
    )

    # 确定尺寸
    hm_w <- max(7, ncol(heatmap_mat) * 0.35)
    hm_h <- max(5, nrow(heatmap_mat) * 0.25)

    ph <- pheatmap(heatmap_scaled,
                   annotation_col = sample_annot,
                   annotation_row = heatmap_annot_row,
                   annotation_colors = ann_colors_hm,
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                   scale = "none",
                   clustering_distance_rows = "correlation",
                   clustering_distance_cols = vegdist(t(heatmap_mat), method = "bray"),
                   show_rownames = TRUE, show_colnames = TRUE,
                   main = paste0("Top ", length(top_cv), " Variable FAPROTAX Functions"),
                   treeheight_row = 15, treeheight_col = 15,
                   fontsize_row = 7, fontsize_col = 6,
                   cutree_rows = max(1, min(6, length(top_cv) %/% 5)),
                   silent = TRUE)

    pdf(file.path(output_dir, "06_Function_Heatmap.pdf"),
        width = hm_w, height = hm_h)
    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
    dev.off()

    save_tsv_file(heatmap_mat, file.path(output_dir, "06_TopCV_Function_Abundance.txt"))
    cat("  06_Function_Heatmap.pdf 已生成\n")
  } else {
    cat("  跳过：CV 计算失败（数据不足）\n")
  }

} else {
  cat("  跳过：热图已禁用或数据不足\n")
}

# ==============================================================================
# 7. OTU 贡献追溯
# ==============================================================================
cat("[7/8] OTU 贡献追溯...\n")

if (run_otu_trace && !is.null(otu_mat) && !is.null(taxonomy)) {

  # ---- 7a. 确定要追溯的功能 ----
  # 选择丰度最高的功能中，在 .mat 中也有记录的功能
  func_means <- rowMeans(faprotax_rel)
  top_funcs <- names(sort(func_means, decreasing = TRUE))

  # 与 otu_mat 列名求交集（跳过 OTUID 列）
  available_funcs <- intersect(top_funcs, colnames(otu_mat)[-1])
  trace_funcs <- head(available_funcs, n_otu_func)

  if (length(trace_funcs) == 0) {
    cat("  跳过：无可追溯的功能（功能名与 .mat 列名不匹配）\n")
  } else {
    cat(sprintf("  追溯 %d 个功能: %s\n", length(trace_funcs),
                paste(trace_funcs, collapse = ", ")))

    # 解析分类学字符串为层级列
    parse_taxonomy <- function(tax_strings) {
      levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      # 单字母缩写 → 全名映射（d__ → Domain, p__ → Phylum, etc.）
      abbrev <- setNames(levels, c("d", "p", "c", "o", "f", "g", "s"))
      result <- matrix("", nrow = length(tax_strings), ncol = length(levels))
      colnames(result) <- levels

      for (i in seq_along(tax_strings)) {
        parts <- strsplit(tax_strings[i], ";")[[1]]
        for (j in seq_along(parts)) {
          kv <- strsplit(trimws(parts[j]), "__")[[1]]
          if (length(kv) == 2) {
            level_name <- abbrev[tolower(kv[1])]
            if (!is.na(level_name)) {
              result[i, level_name] <- kv[2]
            }
          }
        }
      }
      as.data.frame(result, stringsAsFactors = FALSE)
    }

    # 解析分类学
    tax_parsed <- parse_taxonomy(taxonomy$Taxonomy)
    taxonomy <- cbind(taxonomy, tax_parsed)

    # 遍历每个功能
    all_contrib_list <- list()

    for (func_name in trace_funcs) {
      # 在 .mat 中找到对应该功能的 OTU
      otu_col <- otu_mat[[func_name]]
      if (is.null(otu_col)) next

      contrib_otus <- otu_mat$OTUID[otu_col == 1]
      if (length(contrib_otus) == 0) {
        cat(sprintf("    %s: 无 OTU 记录\n", func_name))
        next
      }

      # 关联分类学信息
      contrib_tax <- taxonomy[taxonomy$OTUID %in% contrib_otus, ]
      cat(sprintf("    %s: %d 个 OTU\n", func_name, nrow(contrib_tax)))

      # 门水平汇总
      phyla_summary <- as.data.frame(table(contrib_tax$Phylum), stringsAsFactors = FALSE)
      colnames(phyla_summary) <- c("Phylum", "Count")
      phyla_summary <- phyla_summary[phyla_summary$Phylum != "", , drop = FALSE]
      phyla_summary <- phyla_summary[order(phyla_summary$Count, decreasing = TRUE), ]
      phyla_summary$Proportion <- round(phyla_summary$Count / sum(phyla_summary$Count) * 100, 1)
      phyla_summary$Function <- func_name

      # 属水平汇总（如果有）
      genus_summary <- as.data.frame(table(contrib_tax$Genus), stringsAsFactors = FALSE)
      colnames(genus_summary) <- c("Genus", "Count")
      genus_summary <- genus_summary[genus_summary$Genus != "", , drop = FALSE]
      genus_summary <- genus_summary[order(genus_summary$Count, decreasing = TRUE), ]
      genus_summary$Proportion <- round(genus_summary$Count / sum(genus_summary$Count) * 100, 1)
      genus_summary$Function <- func_name

      # Top OTU 列表（直接展示 OTU ID + 分类）
      top_n <- min(n_otu_tax, nrow(contrib_tax))
      top_otus <- contrib_tax[1:top_n, c("OTUID", "Phylum", "Class", "Order", "Family", "Genus")]
      top_otus$Function <- func_name

      # 保存
      all_contrib_list[[func_name]] <- list(
        phyla = phyla_summary,
        genus = genus_summary,
        top_otus = top_otus,
        n_otu = nrow(contrib_tax)
      )

      # 绘制门水平条形图
      n_phyla <- nrow(phyla_summary)
      if (n_phyla > 0) {
        phyla_summary$Phylum <- factor(phyla_summary$Phylum, levels = rev(phyla_summary$Phylum))
        p <- ggplot(phyla_summary, aes(x = Count, y = Phylum)) +
          geom_col(fill = "#3C5488", width = 0.7) +
          labs(title = paste0(func_name, " (", nrow(contrib_tax), " OTUs)"),
               x = "Number of OTUs", y = "") +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
                panel.grid.major.y = element_blank())

        ggsave(file.path(output_dir,
                         paste0("07_OTU_", gsub("[^a-zA-Z0-9_-]", "_", func_name), ".pdf")),
               p, width = 7, height = max(3, n_phyla * 0.35))
      }
    }

    # ---- 7b. OTU 丰富度总览（每个功能有多少 OTU）----
    func_otu_richness <- data.frame(
      Function = names(all_contrib_list),
      N_OTU = sapply(all_contrib_list, function(x) x$n_otu),
      stringsAsFactors = FALSE
    )
    # 补充总丰度
    func_otu_richness$MeanRelAbund <- round(func_means[func_otu_richness$Function], 2)
    func_otu_richness <- func_otu_richness[order(func_otu_richness$N_OTU, decreasing = TRUE), ]

    save_tsv_file(func_otu_richness, file.path(output_dir, "07_OTU_Richness_per_Function.txt"))

    if (nrow(func_otu_richness) > 1) {
      p <- ggplot(func_otu_richness, aes(x = N_OTU, y = reorder(Function, N_OTU))) +
        geom_col(aes(fill = MeanRelAbund), width = 0.7) +
        scale_fill_gradient(low = "steelblue", high = "firebrick") +
        labs(title = "OTU Richness per Function",
             x = "Number of OTUs", y = "", fill = "Mean Rel.\nAbund. (%)") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              panel.grid.major.y = element_blank())

      ggsave(file.path(output_dir, "07_OTU_Richness_Overview.pdf"), p, width = 7, height = max(4, nrow(func_otu_richness) * 0.4))
    }

    cat("  07_OTU_*.pdf 已生成\n")
  }

} else {
  cat("  跳过：OTU 追溯已禁用或数据不足\n")
}

# ==============================================================================
# 8. 结果汇总
# ==============================================================================
cat("\n======================================================================\n")
cat("  FAPROTAX 分析可视化完成!\n")
cat("  输出目录:", output_dir, "\n")
cat("  分组:", paste(groups, collapse = ", "), "\n")
cat("  样本数:", nrow(metadata), "\n")
cat("  功能数:", nrow(faprotax), "\n")
cat("----------------------------------------------------------------------\n")

out_files <- list.files(output_dir)
for (f in out_files) {
  cat("    -", f, "\n")
}
cat("======================================================================\n")
