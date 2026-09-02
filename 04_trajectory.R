#!/usr/bin/env Rscript
# ============================================================================
# 04_trajectory.R — 骨系分化轨迹与漂移连续性分析（步骤 3）
# 目的：验证 H-状态1"漂移是连续谱"——DriftIndex 沿分化轴的连续变化
# 方法：slingshot（monocle3 因 Bioc 版本不兼容未装，slingshot 为等位替代）
#       用 Harmony 嵌入做轨迹推断，检验 DriftIndex/pseudotime 单调性
# 输入：data/processed/seurat_gse169396_drift.rds
# 输出：
#   results/figures/04_trajectory_umap.pdf      — 轨迹 + DriftIndex 着色
#   results/figures/04_drift_vs_pseudotime.pdf  — DriftIndex ~ pseudotime
#   results/figures/04_genes_along_trajectory.pdf — 关键基因沿轨迹
#   results/tables/04_trajectory_summary.csv
# 用法：Rscript scripts/04_trajectory.R
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(slingshot)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

base_dir <- "D:/WorkBuddySpace/output/Lu-2026-mesenchymal-drift"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "results/figures")
tab_dir  <- file.path(base_dir, "results/tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

cat("===== 步骤 3：分化轨迹与漂移连续性（slingshot）=====\n")

seu <- readRDS(file.path(proc_dir, "seurat_gse169396_drift.rds"))
cat(sprintf("骨系对象: %d 细胞\n", ncol(seu)))

# ---- 1. 准备降维坐标（slingshot 用 harmony 嵌入）---------------------------
# 确保 harmony 嵌入存在（步骤 1 已算）
emb <- seu[["harmony"]]@cell.embeddings[, 1:15]
cat(sprintf("Harmony 嵌入: %d × %d\n", nrow(emb), ncol(emb)))

# 聚类标签（slingshot 需要 cluster 作为轨迹骨架）
clusters <- as.character(Idents(seu))
cat("聚类数:", length(unique(clusters)), "\n")

# ---- 2. slingshot 轨迹推断 --------------------------------------------------
set.seed(42)
sds <- slingshot(emb, clusterLabels = clusters,
                 start.clus = NULL,  # 自动选择起点（MSC 富集簇）
                 approx_points = 200)
cat("检测到轨迹曲线:", ncol(slingCurves(sds)), "条\n")

# 提取 pseudotime（第一曲线）
pt <- slingPseudotime(sds, na = FALSE)[, 1]
names(pt) <- rownames(emb)
seu$pseudotime <- pt[colnames(seu)]
cat(sprintf("pseudotime 范围: [%.2f, %.2f]\n", min(pt, na.rm=TRUE), max(pt, na.rm=TRUE)))

# ---- 3. 关键分析：DriftIndex 沿 pseudotime 的连续性 -------------------------
# 3.1 Spearman 相关（细胞水平）
cor_res <- cor.test(seu$DriftIndex, seu$pseudotime, method = "spearman")
cat(sprintf("\nDriftIndex ~ pseudotime: rho = %.3f, p = %.2e\n",
            cor_res$estimate, cor_res$p.value))

# 3.2 分段趋势（按 pseudotime 十分位看均值——检验单调性而非线性）
pt_bin <- cut(seu$pseudotime, breaks = 10, labels = FALSE)
trend_df <- data.frame(
  pseudotime_bin = pt_bin,
  DriftIndex     = seu$DriftIndex,
  AD_score       = seu$AD_drift_ucell,
  OS_score       = seu$OS_identity_ucell,
  MD_score       = seu$MD_score
) %>%
  group_by(pseudotime_bin) %>%
  summarise(across(everything(), mean), .groups = "drop")
write.csv(trend_df, file.path(tab_dir, "04_trajectory_summary.csv"), row.names = FALSE)

# ---- 4. 绘图 ---------------------------------------------------------------
theme_set(theme_classic(base_size = 12))

# 图 1：轨迹 + DriftIndex 着色
plot_df <- data.frame(
  umap1 = seu[["umap.bone"]]@cell.embeddings[,1],
  umap2 = seu[["umap.bone"]]@cell.embeddings[,2],
  DriftIndex = seu$DriftIndex,
  pseudotime = seu$pseudotime,
  group = seu$group
)
# 去除 NA pseudotime
plot_df <- plot_df[!is.na(plot_df$pseudotime), ]

p1 <- ggplot(plot_df, aes(umap1, umap2, color = DriftIndex)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_gradient2(low = "#377EB8", mid = "#FFFFFF", high = "#E41A1C",
                        midpoint = median(plot_df$DriftIndex)) +
  labs(title = "Trajectory with cell-level DriftIndex") +
  theme(legend.position = "right")

p2 <- ggplot(plot_df, aes(umap1, umap2, color = pseudotime)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_viridis_c(option = "C") +
  labs(title = "slingshot pseudotime") +
  theme(legend.position = "right")

pdf(file.path(fig_dir, "04_trajectory_umap.pdf"), width = 12, height = 6)
print(p1 | p2)
dev.off()

# 图 2：DriftIndex 沿 pseudotime（散点 + 平滑 + 分段均值）
pdf(file.path(fig_dir, "04_drift_vs_pseudotime.pdf"), width = 8, height = 6)
ggplot(plot_df, aes(pseudotime, DriftIndex)) +
  geom_point(size = 0.3, alpha = 0.3, color = "grey60") +
  geom_smooth(method = "loess", se = TRUE, color = "#E41A1C", linewidth = 1) +
  labs(title = "DriftIndex along pseudotime",
       x = "Pseudotime (osteoblast trajectory)",
       y = "Cell-level DriftIndex") +
  annotate("text", x = min(plot_df$pseudotime, na.rm=T) + 0.1,
           y = max(plot_df$DriftIndex, na.rm=T),
           label = sprintf("Spearman rho = %.2f\np = %.1e",
                           cor_res$estimate, cor_res$p.value),
           hjust = 0, size = 4)
dev.off()

# 图 3：关键基因沿轨迹（成骨 vs 成脂程序反转点）
genes_show <- c("RUNX2","SP7","ALPL","BGLAP","SOST","PPARG","CEBPA","ADIPOQ","FABP4","VIM")
genes_show <- intersect(genes_show, rownames(seu))
if (length(genes_show) > 0) {
  expr <- GetAssayData(seu, assay = "SCT", layer = "data")[genes_show, ]
  # 按 pseudotime 排序后滑动平均
  ord <- order(seu$pseudotime, na.last = NA)
  pt_ord <- seu$pseudotime[ord]
  expr_ord <- as.matrix(expr)[, ord]
  # 10 个窗口均值
  n_win <- 20
  win <- floor(seq(1, length(ord), length.out = n_win + 1))
  plot_list <- list()
  for (g in genes_show) {
    means <- sapply(1:n_win, function(i) mean(expr_ord[g, win[i]:win[i+1]]))
    plot_list[[g]] <- ggplot(data.frame(x = 1:n_win, y = means), aes(x, y)) +
      geom_line(color = "#185FA5", linewidth = 0.8) +
      labs(title = g, x = "Pseudotime window", y = "Expr") +
      theme_minimal(base_size = 9)
  }
  pdf(file.path(fig_dir, "04_genes_along_trajectory.pdf"), width = 14, height = 10)
  print(wrap_plots(plot_list, ncol = 4))
  dev.off()
}

# ---- 5. 结论判定 -----------------------------------------------------------
cat("\n=== H-状态1 判定（漂移是否沿分化轴连续变化） ===\n")
if (cor_res$p.value < 0.05) {
  cat(sprintf("  ✅ DriftIndex 与 pseudotime 显著相关 (rho=%.3f, p<0.05)\n",
              cor_res$estimate))
} else {
  cat(sprintf("  ⚠️ DriftIndex 与 pseudotime 相关性不显著 (rho=%.3f, p=%.2f)\n",
              cor_res$estimate, cor_res$p.value))
}
# 单调性检查：分段均值是否递增/递减趋势
first_3 <- mean(trend_df$DriftIndex[1:3])
last_3  <- mean(trend_df$DriftIndex[8:10])
cat(sprintf("  早期(前3段)均值=%.3f vs 晚期(后3段)均值=%.3f → %s\n",
            first_3, last_3,
            ifelse(abs(last_3 - first_3) > 0.01, "存在趋势", "无明确趋势")))

saveRDS(seu, file.path(proc_dir, "seurat_gse169396_traj.rds"))
cat("\n===== 步骤 3 完成 =====\n")
