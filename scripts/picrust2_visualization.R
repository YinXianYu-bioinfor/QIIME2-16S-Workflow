#!/usr/bin/env Rscript
#
# PICRUSt2 功能宏基因组多组可视化（覆盖全部功能层级）
# ==============================================================================
# 覆盖层级: NSTI, KEGG L1/L2, KEGG Pathway, KO, EC
# 用法: Rscript picrust2_visualization.R [选项]
#       数据读取自 results/export/picrust2/out/，输出至 results/export/picrust2/
# 输出: PDF图片 + TSV表格，存放在 results/export/picrust2/picrust2_visualization/ 下
# 帮助: Rscript picrust2_visualization.R --help / -h
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
    p = "picrust2_dir",
    o = "output_dir",
    s = "sample_col",
    g = "group_col",
    "1" = "run_pathway",
    "2" = "run_ko",
    "3" = "run_ec",
    t = "n_top_stack",
    T = "n_top_heatmap",
    d = "n_top_diff",
    P = "padj_cutoff",
    l = "lfc_cutoff",
    w = "fig_w_base",
    z = "fig_h_base",
    c = "group_palette"
  )

  i <- 1
  while (i <= length(args)) {
    if (args[i] %in% c("--help", "-h")) {
      cat("PICRUSt2 功能预测可视化脚本\n")
      cat("用法: Rscript picrust2_visualization.R [选项]\n")
      cat("所有参数均有默认值，仅需指定需覆盖的参数。\n\n")
      cat("路径参数:\n")
      cat("  -m, --metadata=<file>        元数据文件路径 (默认: metadata.txt)\n")
      cat("  -p, --picrust2-dir=<dir>     PICRUSt2 输出目录 (默认: results/export/picrust2/out)\n")
      cat("  -o, --output-dir=<dir>       输出目录 (默认: results/export/picrust2/picrust2_visualization)\n\n")
      cat("数据参数:\n")
      cat("  -s, --sample-col=<col>       样本ID列名 (默认: SampleID)\n")
      cat("  -g, --group-col=<col>        分组列名 (默认: Group)\n\n")
      cat("分析模块开关 (TRUE/FALSE):\n")
      cat("  -1, --run-pathway=<T/F>      KEGG Pathway + L1/L2 分析 (默认: TRUE)\n")
      cat("  -2, --run-ko=<T/F>           KO 级分析 (默认: TRUE)\n")
      cat("  -3, --run-ec=<T/F>           EC 级分析 (默认: TRUE)\n\n")
      cat("可视化参数:\n")
      cat("  -t, --n-top-stack=<n>        堆叠图展示功能数 (默认: 15)\n")
      cat("  -T, --n-top-heatmap=<n>      热图展示高变异功能数 (默认: 40)\n")
      cat("  -d, --n-top-diff=<n>         差异分析展示 top n (默认: 25)\n")
      cat("  -P, --padj-cutoff=<n>        差异显著性阈值 (默认: 0.05)\n")
      cat("  -l, --lfc-cutoff=<n>         火山图 log2FC 截断 (默认: 1.0)\n")
      cat("  -w, --fig-w-base=<n>         图片宽度基准 (默认: 6)\n")
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
metadata_file  <- use_param(cli[["metadata"]], "metadata.txt")
picrust2_dir   <- use_param(cli[["picrust2_dir"]], "results/export/picrust2/out")
output_dir     <- use_param(cli[["output_dir"]], "results/export/picrust2/picrust2_visualization")

# ---- 数据参数 ----
sample_col     <- use_param(cli[["sample_col"]], "SampleID")
group_col      <- use_param(cli[["group_col"]], "Group")

# ---- 功能层级开关 ----
run_pathway    <- as_flag(use_param(cli[["run_pathway"]], TRUE))
run_ko         <- as_flag(use_param(cli[["run_ko"]], TRUE))
run_ec         <- as_flag(use_param(cli[["run_ec"]], TRUE))

# ---- 可视化参数 ----
n_top_stack    <- as.numeric(use_param(cli[["n_top_stack"]], 15))
n_top_heatmap  <- as.numeric(use_param(cli[["n_top_heatmap"]], 40))
n_top_diff     <- as.numeric(use_param(cli[["n_top_diff"]], 25))
padj_cutoff    <- as.numeric(use_param(cli[["padj_cutoff"]], 0.05))
lfc_cutoff     <- as.numeric(use_param(cli[["lfc_cutoff"]], 1.0))

# ---- 配色 ----
group_palette  <- use_param(cli[["group_palette"]], "npg")

# ---- 图片尺寸 ----
fig_w_base     <- as.numeric(use_param(cli[["fig_w_base"]], 6))
fig_h_base     <- as.numeric(use_param(cli[["fig_h_base"]], 5))

# ==============================================================================
# 1. 包加载与工具函数
# ==============================================================================
cat("[1/8] 加载 R 包...\n")

required_packages <- c(
  "ggplot2", "tidyr", "dplyr", "tibble", "readr",
  "vegan", "ape", "pheatmap", "RColorBrewer", "grDevices",
  "ggrepel", "ggsci", "viridis", "reshape2", "FSA"
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
    install.packages(missing, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN", quiet = TRUE, Ncpus = parallel::detectCores())
  }
  invisible(lapply(pkgs, library, character.only = TRUE, quietly = TRUE))
}

check_bioc_packages <- function(pkgs) {
  missing <- c()
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    }
  }
  if (length(missing) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
    }
    cat("  安装缺失的 Bioconductor 包:", paste(missing, collapse = ", "), "\n")
    BiocManager::install(missing, update = FALSE, ask = FALSE, quiet = TRUE)
  }
  invisible(lapply(pkgs, library, character.only = TRUE, quietly = TRUE))
}

check_packages(required_packages)
check_bioc_packages("DESeq2")
cat("  所有包加载成功。\n\n")

# ---- 工具函数 ----

# 分组配色（默认 npg 期刊色板）
get_group_colors <- function(n, palette = "npg") {
  if (palette == "npg") {
    return(rep(pal_npg()(10), length.out = n))
  } else if (palette == "aaas") {
    return(rep(pal_aaas()(10), length.out = n))
  } else if (palette == "nejm") {
    return(rep(pal_nejm()(8), length.out = n))
  } else if (palette == "lancet") {
    return(rep(pal_lancet()(9), length.out = n))
  } else if (palette == "jco") {
    return(rep(pal_jco()(10), length.out = n))
  } else if (palette == "viridis") {
    return(viridis(n))
  } else if (palette == "Set1") {
    return(brewer.pal(max(3, min(n, 9)), "Set1")[1:n])
  } else if (palette == "Set2") {
    return(brewer.pal(max(3, min(n, 8)), "Set2")[1:n])
  } else if (palette == "Set3") {
    return(brewer.pal(max(3, min(n, 12)), "Set3")[1:n])
  }
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

  grp_lev <- unique(as.character(df[[grp_col]]))
  if (length(grp_lev) <= 2) {
    # 只有两组时 K-W 已足够，dunnTest 在两组时内部会报维度错误
    dunn_df <- data.frame(Comparison = paste(grp_lev, collapse = " - "),
                          Z = NA, P.unadj = NA, P.adj = kw_p,
                          sig_label = ifelse(kw_p < 0.001, "***",
                                      ifelse(kw_p < 0.01,  "**",
                                      ifelse(kw_p < 0.05,  "*", "ns"))),
                          stringsAsFactors = FALSE)
  } else {
    dunn <- dunnTest(df[[val_col]] ~ df[[grp_col]], method = p_adj_method)
    dunn_df <- as.data.frame(dunn$res)
    dunn_df$Comparison <- gsub(" ", "", dunn_df$Comparison)
    dunn_df$sig_label <- ifelse(dunn_df$P.adj < 0.001, "***",
                         ifelse(dunn_df$P.adj < 0.01,  "**",
                         ifelse(dunn_df$P.adj < 0.05,  "*", "ns")))
    dunn_df$P.adj <- signif(dunn_df$P.adj, 4)
  }
  list(kw_p = kw_p, dunn = dunn_df)
}

# 导出 TSV
save_tsv_file <- function(data, file_path) {
  write.table(data, file_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# 读取 gzipped PICRUSt2 表格
read_picrust_tsv <- function(file_path, ...) {
  read.delim(gzfile(file_path), sep = "\t", check.names = FALSE, ...)
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
  results <- data.frame(Comparison = character(), R2 = numeric(), p_value = numeric(), p_adj = numeric(), stringsAsFactors = FALSE)
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

# DESeq2 LRT 差异分析通用函数
run_deseq2_lrt <- function(count_mat, metadata, grp_col, lfc_cutoff = 1.0, padj_cutoff = 0.05) {
  common <- intersect(colnames(count_mat), rownames(metadata))
  count_mat <- count_mat[, common, drop = FALSE]
  meta <- metadata[common, , drop = FALSE]
  meta[[grp_col]] <- factor(meta[[grp_col]])

  int_mat <- round(count_mat)
  int_mat <- int_mat[rowSums(int_mat) > 0, ]

  dds <- DESeqDataSetFromMatrix(countData = int_mat, colData = meta, design = as.formula(paste0("~ ", grp_col)))
  dds <- DESeq(dds, test = "LRT", reduced = ~ 1, quiet = TRUE)
  res <- results(dds, tidy = TRUE)
  colnames(res)[1] <- "Feature"
  res <- res[order(res$padj), ]
  return(list(dds = dds, result = res))
}

# 绘制 DESeq2 LRT 柱状图
plot_lrt_bar <- function(res_df, n_top, padj_cutoff, title, file_path, width = 9, height = 7) {
  top_df <- res_df[1:min(n_top, nrow(res_df)), ]
  top_df$Feature <- factor(top_df$Feature, levels = rev(unique(top_df$Feature)))
  p <- ggplot(top_df, aes(x = -log10(padj), y = Feature)) +
    geom_bar(stat = "identity", aes(fill = -log10(padj)), width = 0.7) +
    geom_vline(xintercept = -log10(padj_cutoff), linetype = "dashed", color = "red", linewidth = 0.4) +
    scale_fill_gradient(low = "steelblue", high = "firebrick") +
    labs(title = title, x = expression(-log[10](adjusted~p)), y = "") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 10))
  ggsave(file_path, p, width = width, height = height)
}

# 绘制火山图
plot_volcano <- function(dds, grp_col, cmp, ref, padj_cutoff, lfc_cutoff, top_n_label, title, file_path, width = 7, height = 6) {
  res_pair <- results(dds, contrast = c(grp_col, cmp, ref), tidy = TRUE)
  colnames(res_pair)[1] <- "Feature"
  res_pair$sig <- "Not Significant"
  res_pair$sig[res_pair$padj < padj_cutoff & res_pair$log2FoldChange > lfc_cutoff & !is.na(res_pair$padj)] <- "Up"
  res_pair$sig[res_pair$padj < padj_cutoff & res_pair$log2FoldChange < -lfc_cutoff & !is.na(res_pair$padj)] <- "Down"
  res_pair$sig <- factor(res_pair$sig, levels = c("Up", "Down", "Not Significant"))

  n_up <- sum(res_pair$sig == "Up")
  n_down <- sum(res_pair$sig == "Down")

  top_up <- head(res_pair[res_pair$sig == "Up", ], top_n_label)
  top_down <- head(res_pair[res_pair$sig == "Down", ], top_n_label)
  label_df <- rbind(top_up, top_down)

  p <- ggplot(res_pair, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.7, size = 1.0) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "gray50", linewidth = 0.3) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "gray50", linewidth = 0.3) +
    scale_color_manual(values = c("Up" = "firebrick3", "Down" = "navy", "Not Significant" = "#A9A9A9")) +
    geom_text_repel(data = label_df, aes(label = Feature), size = 2.5, box.padding = 0.5, max.overlaps = 20, show.legend = FALSE) +
    labs(title = title, x = "log2(Fold Change)", y = expression(-log[10](adjusted~p)), color = "Regulation") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 10))
  ggsave(file_path, p, width = width, height = height)

  return(res_pair)
}

# 差异热图（pheatmap + vst 标准化）
plot_diff_heatmap <- function(dds, feature_names, sample_annot, ann_colors, title, file_path, width = 9, height = 8) {
  feat <- feature_names[feature_names %in% rownames(dds)]
  if (length(feat) < 2) {
    warning("差异热图: 显著功能数不足 2 个，跳过。")
    return(NULL)
  }
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  vsd_mat <- assay(vsd)[feat, , drop = FALSE]
  vsd_mat <- vsd_mat[rowSums(is.na(vsd_mat)) == 0, , drop = FALSE]
  vsd_mat <- t(scale(t(vsd_mat)))
  vsd_mat[is.nan(vsd_mat)] <- 0

  ph <- pheatmap(vsd_mat, annotation_col = sample_annot, annotation_colors = ann_colors,
                 color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                 show_rownames = TRUE, show_colnames = TRUE, main = title,
                 treeheight_row = 15, treeheight_col = 15,
                 fontsize_row = 7, fontsize_col = 6,
                 clustering_distance_rows = "correlation",
                 clustering_distance_cols = "euclidean", silent = TRUE)
  pdf(file_path, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()
}

# ==============================================================================
# 2. 数据读取与预处理
# ==============================================================================
cat("[2/8] 读取数据...\n")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 读取元数据 ----
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

# ---- 读取通路丰度 (pathways_out) ----
pathway_file <- file.path(picrust2_dir, "pathways_out/path_abun_unstrat.tsv.gz")
if (file.exists(pathway_file)) {
  pathway <- read_picrust_tsv(pathway_file)
  rownames(pathway) <- pathway$pathway
  pathway <- pathway[, -1, drop = FALSE]
  # 去除可能的重复行名
  if (any(duplicated(rownames(pathway)))) {
    pathway <- as.data.frame(rowsum(pathway, rownames(pathway)))
  }
  common <- intersect(colnames(pathway), rownames(metadata))
  pathway <- pathway[, common, drop = FALSE]
  metadata <- metadata[common, , drop = FALSE]
  cat(sprintf("  通路丰度表: %d × %d\n", nrow(pathway), ncol(pathway)))
} else {
  cat("  ! 未找到通路丰度表，跳过 Pathway 分析:", pathway_file, "\n")
  run_pathway <- FALSE
}

# ---- 读取 NSTI ----
nsti_file <- file.path(picrust2_dir, "KO_metagenome_out/weighted_nsti.tsv.gz")
nsti_meta <- NULL
if (file.exists(nsti_file)) {
  nsti <- read_picrust_tsv(nsti_file)
  colnames(nsti) <- c("SampleID", "weighted_NSTI")
  nsti <- nsti[nsti$SampleID %in% rownames(metadata), , drop = FALSE]
  nsti_meta <- merge(nsti, metadata, by.x = "SampleID", by.y = sample_col)
  cat(sprintf("  NSTI 数据: %d 样本\n", nrow(nsti)))
}

# ---- 读取 KEGG L1/L2 ----
l1_file <- file.path(picrust2_dir, "KEGG.PathwayL1.raw.txt")
l2_file <- file.path(picrust2_dir, "KEGG.PathwayL2.raw.txt")
kegg_l1 <- kegg_l2 <- NULL
if (run_pathway && file.exists(l1_file)) {
  kegg_l1 <- read.delim(l1_file, sep = "\t", check.names = FALSE)
  rownames(kegg_l1) <- trimws(kegg_l1$PathwayL1)
  kegg_l1 <- kegg_l1[, -1, drop = FALSE]
  common <- intersect(colnames(kegg_l1), rownames(metadata))
  kegg_l1 <- kegg_l1[, common, drop = FALSE]
  cat(sprintf("  KEGG L1: %d × %d\n", nrow(kegg_l1), ncol(kegg_l1)))
}
if (run_pathway && file.exists(l2_file)) {
  kegg_l2 <- read.delim(l2_file, sep = "\t", check.names = FALSE)
  rownames(kegg_l2) <- trimws(kegg_l2$PathwayL2)
  kegg_l2 <- kegg_l2[, -1, drop = FALSE]
  # 合并重复行（KEGG 层级中断行）
  if (any(duplicated(rownames(kegg_l2)))) {
    kegg_l2 <- as.data.frame(rowsum(kegg_l2, rownames(kegg_l2)))
  }
  common <- intersect(colnames(kegg_l2), rownames(metadata))
  kegg_l2 <- kegg_l2[, common, drop = FALSE]
  cat(sprintf("  KEGG L2: %d × %d\n", nrow(kegg_l2), ncol(kegg_l2)))
}

# ---- 读取 KO 数据 ----
ko_data <- ko_count <- NULL
ko_file <- file.path(picrust2_dir, "KO.tsv")
ko_meta_file <- file.path(picrust2_dir, "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")
if (run_ko && file.exists(ko_file) && file.exists(ko_meta_file)) {
  ko_data <- read.delim(ko_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  rownames(ko_data) <- ko_data[["function"]]
  ko_desc <- ko_data[, "description", drop = FALSE]
  ko_count <- read_picrust_tsv(ko_meta_file)
  rownames(ko_count) <- ko_count[[1]]
  ko_count <- ko_count[, -1, drop = FALSE]
  common <- intersect(colnames(ko_count), rownames(metadata))
  ko_count <- ko_count[, common, drop = FALSE]
  ko_desc <- ko_desc[intersect(rownames(ko_desc), rownames(ko_count)), , drop = FALSE]
  cat(sprintf("  KO 丰度: %d × %d\n", nrow(ko_count), ncol(ko_count)))
}

# ---- 读取 EC 数据 ----
ec_data <- ec_count <- NULL
ec_file <- file.path(picrust2_dir, "EC.tsv")
ec_meta_file <- file.path(picrust2_dir, "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz")
if (run_ec && file.exists(ec_file) && file.exists(ec_meta_file)) {
  ec_data <- read.delim(ec_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  rownames(ec_data) <- ec_data[["function"]]
  ec_desc <- ec_data[, "description", drop = FALSE]
  ec_count <- read_picrust_tsv(ec_meta_file)
  rownames(ec_count) <- ec_count[[1]]
  ec_count <- ec_count[, -1, drop = FALSE]
  common <- intersect(colnames(ec_count), rownames(metadata))
  ec_count <- ec_count[, common, drop = FALSE]
  ec_desc <- ec_desc[intersect(rownames(ec_desc), rownames(ec_count)), , drop = FALSE]
  cat(sprintf("  EC 丰度: %d × %d\n", nrow(ec_count), ncol(ec_count)))
}

cat("  数据读取完成。\n\n")

# ==============================================================================
# 3. NSTI 质量评估
# ==============================================================================
cat("[3/8] NSTI 质量评估...\n")

if (!is.null(nsti_meta) && nrow(nsti_meta) > 0) {
  nsti_kw <- multi_group_test(nsti_meta, "weighted_NSTI", group_col)

  p <- ggplot(nsti_meta, aes_string(x = group_col, y = "weighted_NSTI", fill = group_col)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    scale_fill_manual(values = group_colors) +
    labs(title = bquote(atop("Weighted NSTI Distribution", "Kruskal-Wallis p = " ~ .(nsti_kw$kw_p))),
         x = "", y = "Weighted NSTI") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                       axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(output_dir, "01_NSTI_Boxplot.pdf"), p, width = fig_w_base + n_group * 0.3, height = fig_h_base)
  save_tsv_file(nsti_meta, file.path(output_dir, "01_NSTI_Values.txt"))
  save_tsv_file(nsti_kw$dunn, file.path(output_dir, "01_NSTI_Dunn_posthoc.txt"))
  cat("  01_NSTI_Boxplot.pdf 已生成\n")
} else {
  cat("  跳过：NSTI 数据不可用\n")
}

# ==============================================================================
# 4. KEGG 层级全局概览
# ==============================================================================
cat("[4/8] KEGG 层级概览...\n")

if (run_pathway && !is.null(kegg_l1) && nrow(kegg_l1) > 0) {
  # ---- 4a. L1 大类堆叠柱状图 ----
  l1_rel <- normalize_abundance(kegg_l1)
  l1_melt <- reshape2::melt(as.matrix(l1_rel))
  colnames(l1_melt) <- c("Category", "SampleID", "Abundance")
  l1_melt <- merge(l1_melt, metadata, by.x = "SampleID", by.y = sample_col)
  l1_melt$SampleID <- factor(l1_melt$SampleID, levels = colnames(l1_rel))

  p <- ggplot(l1_melt, aes(x = SampleID, y = Abundance, fill = Category)) +
    geom_bar(stat = "identity", position = "stack", width = 0.85) +
    facet_grid(. ~ get(group_col), scales = "free_x", space = "free") +
    scale_fill_manual(values = rep(pal_npg()(10), length.out = nrow(l1_rel))) +
    labs(title = "KEGG L1 Category Composition", x = "", y = "Relative Abundance (%)") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
                       plot.title = element_text(hjust = 0.5), strip.text = element_text(size = 8))
  ggsave(file.path(output_dir, "02_KEGG_L1_Stackplot.pdf"), p,
         width = max(10, ncol(kegg_l1) * 0.4), height = 6)
  cat("  02_KEGG_L1_Stackplot.pdf 已生成\n")

  # ---- 4b. L2 高变异分类热图 ----
  if (!is.null(kegg_l2) && nrow(kegg_l2) > 1) {
    l2_rel <- normalize_abundance(kegg_l2)
    l2_cv <- apply(l2_rel, 1, function(x) sd(x) / mean(x))
    l2_cv <- l2_cv[is.finite(l2_cv)]
    top_cv <- names(sort(l2_cv, decreasing = TRUE))[1:min(n_top_heatmap, length(l2_cv))]
    l2_mat <- l2_rel[top_cv, , drop = FALSE]

    ph <- pheatmap(l2_mat, annotation_col = sample_annot, annotation_colors = ann_colors,
                   scale = "row", clustering_distance_rows = "correlation",
                   clustering_distance_cols = vegdist(t(l2_mat), method = "bray"),
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                   show_rownames = TRUE, show_colnames = TRUE,
                   main = paste0("Top ", length(top_cv), " Variable KEGG L2 Categories"),
                   treeheight_row = 15, treeheight_col = 15,
                   fontsize_row = 7, fontsize_col = 6, silent = TRUE)
    pdf(file.path(output_dir, "03_KEGG_L2_Heatmap.pdf"),
        width = max(8, ncol(kegg_l2) * 0.3), height = max(6, nrow(l2_mat) * 0.25))
    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
    dev.off()
    cat("  03_KEGG_L2_Heatmap.pdf 已生成\n")
  } else {
    cat("  跳过 L2 热图: L2 数据不足\n")
  }
} else {
  cat("  跳过：KEGG 层级数据不可用\n")
}

# ==============================================================================
# 5. KEGG Pathway 级分析
# ==============================================================================
cat("[5/8] Pathway 级分析...\n")

if (run_pathway && !is.null(pathway) && nrow(pathway) > 0) {

  pathway_rel <- normalize_abundance(pathway)

  # ---- 5a. Top 通路分面水平柱状图（每组 top 10） ----
  top_pathway_list <- lapply(groups, function(g) {
    samps <- rownames(metadata[metadata[[group_col]] == g, ])
    samps <- intersect(colnames(pathway_rel), samps)
    if (length(samps) == 0) return(NULL)
    means <- rowMeans(pathway_rel[, samps, drop = FALSE])
    sorted <- sort(means, decreasing = TRUE)
    n <- min(10, length(sorted))
    data.frame(Pathway = names(sorted)[1:n],
               Abundance = sorted[1:n],
               Group = g,
               stringsAsFactors = FALSE)
  })
  top_df <- do.call(rbind, top_pathway_list)
  top_df$Pathway <- factor(top_df$Pathway, levels = rev(unique(top_df$Pathway)))

  p <- ggplot(top_df, aes(x = Abundance, y = Pathway)) +
    geom_col(fill = "#3C5488", width = 0.7) +
    facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
    labs(title = "Top 10 Pathways by Group",
         x = "Relative Abundance (%)", y = "") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          strip.background = element_rect(fill = "#F0F0F0"),
          strip.text = element_text(face = "bold"),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 7))
  ggsave(file.path(output_dir, "04_Pathway_Top10_Barplot.pdf"), p, width = 8, height = 7)
  cat("  04_Pathway_Top10_Barplot.pdf 已生成\n")



  # ---- 5b. 功能 beta 多样性（Bray-Curtis PCoA） ----
  bray_dist <- vegdist(t(pathway), method = "bray")
  pcoa <- cmdscale(bray_dist, k = 2, eig = TRUE)
  pcoa_df <- data.frame(SampleID = rownames(pcoa$points),
                        PCoA1 = pcoa$points[, 1],
                        PCoA2 = pcoa$points[, 2],
                        stringsAsFactors = FALSE)
  pcoa_df <- merge(pcoa_df, metadata, by.x = "SampleID", by.y = sample_col)

  set.seed(123)
  permanova <- adonis2(bray_dist ~ metadata[[group_col]], permutations = 999)
  perma_r2 <- round(permanova$R2[1], 3)
  perma_p <- permanova$`Pr(>F)`[1]

  # betadisper 检验
  bd <- betadisper(bray_dist, metadata[[group_col]])
  bd_anova <- anova(bd)
  bd_p <- bd_anova$`Pr(>F)`[1]

  # pairwise PERMANOVA
  grp_vec <- setNames(metadata[[group_col]], rownames(metadata))
  pw_perma <- pairwise_permanova(bray_dist, grp_vec[colnames(pathway)])

  eig1 <- round(pcoa$eig[1] / sum(abs(pcoa$eig)) * 100, 1)
  eig2 <- round(pcoa$eig[2] / sum(abs(pcoa$eig)) * 100, 1)

  p <- ggplot(pcoa_df, aes_string(x = "PCoA1", y = "PCoA2", color = group_col, shape = group_col)) +
    geom_point(size = 2.5, alpha = 0.8) +
    stat_ellipse(level = 0.95, linewidth = 0.6) +
    scale_color_manual(values = group_colors) +
    scale_shape_manual(values = (1:n_group) %% 25 + 1) +
    labs(title = paste0("PCoA (Bray-Curtis)\nPERMANOVA R² = ", perma_r2,
                        ", p = ", signif(perma_p, 3),
                        "  |  betadisper p = ", signif(bd_p, 3)),
         x = paste0("PCoA1 (", eig1, "%)"),
         y = paste0("PCoA2 (", eig2, "%)")) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(output_dir, "07_Pathway_PCoA.pdf"), p, width = 7, height = 6)
  save_tsv_file(pcoa_df, file.path(output_dir, "07_PCoA_Results.txt"))
  save_tsv_file(as.data.frame(permanova), file.path(output_dir, "07_PERMANOVA.txt"))
  save_tsv_file(pw_perma, file.path(output_dir, "07_Pairwise_PERMANOVA.txt"))
  cat("  07_Pathway_PCoA.pdf 已生成\n")

  # ---- 5c. DESeq2 LRT 多组差异分析 ----
  cat("  运行 DESeq2 LRT 差异分析...\n")
  de_res <- run_deseq2_lrt(pathway, metadata, group_col, lfc_cutoff, padj_cutoff)
  dds <- de_res$dds
  lrt_res <- de_res$result

  sig_res <- lrt_res[lrt_res$padj < padj_cutoff & !is.na(lrt_res$padj), ]
  n_sig <- nrow(sig_res)
  cat(sprintf("  LRT 检测到 %d 个显著差异通路 (padj < %.3f)\n", n_sig, padj_cutoff))

  top_diff <- if (n_sig > 0) head(sig_res$Feature, n_top_diff) else head(lrt_res$Feature, n_top_diff)

  # LRT 柱状图
  plot_lrt_bar(lrt_res, min(n_top_diff, nrow(lrt_res)), padj_cutoff,
               paste0("DESeq2 LRT: Top ", n_top_diff, " Differential Pathways"),
               file.path(output_dir, "08_DESeq2_LRT_Barplot.pdf"))
  save_tsv_file(lrt_res, file.path(output_dir, "08_DESeq2_LRT_Full_Results.txt"))
  cat("  08_DESeq2_LRT_Barplot.pdf 已生成\n")


  # ---- 5d. 火山图（两两比较） ----
  ref_grp <- names(sort(table(metadata[[group_col]]), decreasing = TRUE))[1]
  cmp_grps <- setdiff(groups, ref_grp)
  for (cmp in cmp_grps) {
    safe_cmp <- gsub("[^a-zA-Z0-9_]", "", cmp)
    safe_ref <- gsub("[^a-zA-Z0-9_]", "", ref_grp)
    volcano_title <- paste0(cmp, " vs ", ref_grp)
    volcano_file <- file.path(output_dir, paste0("11_Volcano_", safe_cmp, "_vs_", safe_ref, ".pdf"))
    pair_res <- plot_volcano(dds, group_col, cmp, ref_grp, padj_cutoff, lfc_cutoff,
                             3, volcano_title, volcano_file)
    save_tsv_file(pair_res, file.path(output_dir, paste0("11_Diff_", safe_cmp, "_vs_", safe_ref, ".txt")))
    cat(sprintf("  11_Volcano_%s_vs_%s.pdf 已生成\n", safe_cmp, safe_ref))
  }
} else {
  cat("  跳过：Pathway 数据不可用\n")
}

# ==============================================================================
# 6. KO 级分析
# ==============================================================================
cat("[6/8] KO 级分析...\n")

if (run_ko && !is.null(ko_count) && nrow(ko_count) > 0) {

  # ---- 6a. KO PCA ----
  # KO 维度高，用 top 500 高变异 KO 做 PCA
  ko_rel <- normalize_abundance(ko_count)
  ko_cv <- apply(ko_rel, 1, function(x) sd(x) / mean(x))
  ko_cv <- ko_cv[is.finite(ko_cv)]
  top_ko <- names(sort(ko_cv, decreasing = TRUE))[1:min(500, length(ko_cv))]
  ko_pca_mat <- log2(ko_rel[top_ko, , drop = FALSE] + 1)

  set.seed(123)
  pca <- prcomp(t(ko_pca_mat), center = TRUE, scale. = TRUE)
  pca_var <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)
  pca_df <- as.data.frame(pca$x)
  pca_df$SampleID <- rownames(pca_df)
  pca_df <- merge(pca_df, metadata, by.x = "SampleID", by.y = sample_col)

  p <- ggplot(pca_df, aes_string(x = "PC1", y = "PC2", color = group_col, shape = group_col)) +
    geom_point(size = 2.5, alpha = 0.8) +
    stat_ellipse(level = 0.95, linewidth = 0.6) +
    scale_color_manual(values = group_colors) +
    scale_shape_manual(values = (1:n_group) %% 25 + 1) +
    labs(title = paste0("KO PCA (top ", length(top_ko), " variable KO by CV)"),
         x = paste0("PC1 (", pca_var[1], "%)"),
         y = paste0("PC2 (", pca_var[2], "%)")) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(output_dir, "12_KO_PCA.pdf"), p, width = 7, height = 6)
  save_tsv_file(pca_df, file.path(output_dir, "12_KO_PCA_Results.txt"))
  cat("  12_KO_PCA.pdf 已生成\n")

  # ---- 6b. KO DESeq2 LRT ----
  cat("  运行 KO DESeq2 LRT 差异分析...\n")
  de_res <- run_deseq2_lrt(ko_count, metadata, group_col, lfc_cutoff, padj_cutoff)
  dds <- de_res$dds
  lrt_res <- de_res$result

  # 添加 KO 描述
  if (!is.null(ko_desc)) {
    lrt_res$Description <- ko_desc[lrt_res$Feature, "description"]
    lrt_res$Description[is.na(lrt_res$Description)] <- ""
  }

  sig_res <- lrt_res[lrt_res$padj < padj_cutoff & !is.na(lrt_res$padj), ]
  n_sig <- nrow(sig_res)
  cat(sprintf("  KO LRT 检测到 %d 个显著差异 KO (padj < %.3f)\n", n_sig, padj_cutoff))

  top_diff <- if (n_sig > 0) head(sig_res$Feature, n_top_diff) else head(lrt_res$Feature, n_top_diff)

  # LRT 柱状图（使用带描述的标签）
  bar_df <- lrt_res[1:min(n_top_diff, nrow(lrt_res)), ]
  if (!is.null(ko_desc) && "Description" %in% colnames(bar_df)) {
    bar_df$label <- ifelse(bar_df$Description != "" & !is.na(bar_df$Description),
                           paste0(bar_df$Feature, " (", bar_df$Description, ")"),
                           bar_df$Feature)
  } else {
    bar_df$label <- bar_df$Feature
  }
  bar_df$label <- factor(bar_df$label, levels = rev(unique(bar_df$label)))

  p <- ggplot(bar_df, aes(x = -log10(padj), y = label)) +
    geom_bar(stat = "identity", aes(fill = -log10(padj)), width = 0.7) +
    geom_vline(xintercept = -log10(padj_cutoff), linetype = "dashed", color = "red", linewidth = 0.4) +
    scale_fill_gradient(low = "steelblue", high = "firebrick") +
    labs(title = paste0("KO DESeq2 LRT: Top ", n_top_diff, " Differential KO"),
         x = expression(-log[10](adjusted~p)), y = "") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 10))
  ggsave(file.path(output_dir, "13_KO_LRT_Barplot.pdf"), p, width = 10, height = max(6, n_top_diff * 0.3))
  save_tsv_file(lrt_res, file.path(output_dir, "13_KO_LRT_Full_Results.txt"))
  cat("  13_KO_LRT_Barplot.pdf 已生成\n")


} else {
  cat("  跳过：KO 数据不可用\n")
}

# ==============================================================================
# 7. EC 级分析
# ==============================================================================
cat("[7/8] EC 级分析...\n")

if (run_ec && !is.null(ec_count) && nrow(ec_count) > 0) {

  # ---- 7a. EC 大类堆叠柱状图 ----
  # EC:1.x.x.x → Class 1 (Oxidoreductases), etc.
  ec_class_num <- function(ec_id) {
    as.integer(sub("EC:(\\d+).*", "\\1", ec_id))
  }
  ec_rel <- normalize_abundance(ec_count)
  ec_classes <- ec_class_num(rownames(ec_rel))
  ec_class_names <- c("1" = "Oxidoreductases", "2" = "Transferases",
                      "3" = "Hydrolases", "4" = "Lyases",
                      "5" = "Isomerases", "6" = "Ligases",
                      "7" = "Translocases")
  ec_class_labels <- ec_class_names[as.character(ec_classes)]
  ec_class_labels[is.na(ec_class_labels)] <- "Other"

  ec_class_agg <- rowsum(ec_rel, group = ec_class_labels)
  ec_class_melt <- reshape2::melt(as.matrix(ec_class_agg))
  colnames(ec_class_melt) <- c("EC_Class", "SampleID", "Abundance")
  ec_class_melt <- merge(ec_class_melt, metadata, by.x = "SampleID", by.y = sample_col)
  ec_class_melt$SampleID <- factor(ec_class_melt$SampleID, levels = colnames(ec_rel))

  p <- ggplot(ec_class_melt, aes(x = SampleID, y = Abundance, fill = EC_Class)) +
    geom_bar(stat = "identity", position = "stack", width = 0.85) +
    facet_grid(. ~ get(group_col), scales = "free_x", space = "free") +
    scale_fill_manual(values = rep(pal_npg()(10), length.out = nrow(ec_class_agg))) +
    labs(title = "EC Class Composition", x = "", y = "Relative Abundance (%)") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
                       plot.title = element_text(hjust = 0.5), strip.text = element_text(size = 8))
  ggsave(file.path(output_dir, "17_EC_Class_Stackplot.pdf"), p,
         width = max(10, ncol(ec_count) * 0.4), height = 6)
  save_tsv_file(ec_class_agg, file.path(output_dir, "17_EC_Class_Abundance.txt"))
  cat("  17_EC_Class_Stackplot.pdf 已生成\n")

  # ---- 7b. EC DESeq2 LRT ----
  cat("  运行 EC DESeq2 LRT 差异分析...\n")
  de_res <- run_deseq2_lrt(ec_count, metadata, group_col, lfc_cutoff, padj_cutoff)
  dds <- de_res$dds
  lrt_res <- de_res$result

  # 添加 EC 描述
  if (!is.null(ec_desc)) {
    lrt_res$Description <- ec_desc[lrt_res$Feature, "description"]
    lrt_res$Description[is.na(lrt_res$Description)] <- ""
  }

  sig_res <- lrt_res[lrt_res$padj < padj_cutoff & !is.na(lrt_res$padj), ]
  n_sig <- nrow(sig_res)
  cat(sprintf("  EC LRT 检测到 %d 个显著差异 EC (padj < %.3f)\n", n_sig, padj_cutoff))

  top_diff <- if (n_sig > 0) head(sig_res$Feature, n_top_diff) else head(lrt_res$Feature, n_top_diff)

  # LRT 柱状图（带描述）
  bar_df <- lrt_res[1:min(n_top_diff, nrow(lrt_res)), ]
  if (!is.null(ec_desc) && "Description" %in% colnames(bar_df)) {
    bar_df$label <- ifelse(bar_df$Description != "" & !is.na(bar_df$Description),
                           paste0(bar_df$Feature, " (", bar_df$Description, ")"),
                           bar_df$Feature)
  } else {
    bar_df$label <- bar_df$Feature
  }
  bar_df$label <- factor(bar_df$label, levels = rev(unique(bar_df$label)))

  p <- ggplot(bar_df, aes(x = -log10(padj), y = label)) +
    geom_bar(stat = "identity", aes(fill = -log10(padj)), width = 0.7) +
    geom_vline(xintercept = -log10(padj_cutoff), linetype = "dashed", color = "red", linewidth = 0.4) +
    scale_fill_gradient(low = "steelblue", high = "firebrick") +
    labs(title = paste0("EC DESeq2 LRT: Top ", n_top_diff, " Differential EC"),
         x = expression(-log[10](adjusted~p)), y = "") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 10))
  ggsave(file.path(output_dir, "18_EC_LRT_Barplot.pdf"), p, width = 10, height = max(6, n_top_diff * 0.3))
  save_tsv_file(lrt_res, file.path(output_dir, "18_EC_LRT_Full_Results.txt"))
  cat("  18_EC_LRT_Barplot.pdf 已生成\n")


} else {
  cat("  跳过：EC 数据不可用\n")
}

# ==============================================================================
# 8. 结果汇总
# ==============================================================================
cat("\n======================================================================\n")
cat("  PICRUSt2 多组可视化完成!\n")
cat("  输出目录:", output_dir, "\n")
cat("  分组:", paste(groups, collapse = ", "), "\n")
cat("  样本数:", nrow(metadata), "\n")
cat("----------------------------------------------------------------------\n")

out_files <- list.files(output_dir)
for (f in out_files) {
  cat("    -", f, "\n")
}
cat("======================================================================\n")
