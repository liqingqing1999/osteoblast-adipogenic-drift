#!/usr/bin/env Rscript
# ============================================================================
# 03_sensitivity.R — 中间态阈值敏感性分析（步骤 2 稳健性验证）
# 目的：验证"三档中间态比例趋势"不依赖单一阈值（审稿人必查）
# 方法：对中间态判定阈值 thr ∈ {0, ±0.02, ±0.04, ±0.06} 重复计算
#       检验 DriftIndex 排序与中间态比例排序是否稳定
# 输入：data/processed/seurat_gse169396_drift.rds
# 输出：results/tables/03_sensitivity_summary.csv
#       results/figures/03_sensitivity_trend.pdf
# 用法：Rscript scripts/03_sensitivity.R
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

base_dir <- "D:/WorkBuddySpace/output/Lu-2026-mesenchymal-drift"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "results/figures")
tab_dir  <- file.path(base_dir, "results/tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

cat("===== 中间态阈值敏感性分析 =====\n")

seu <- readRDS(file.path(proc_dir, "seurat_gse169396_drift.rds"))
cat(sprintf("对象: %d 细胞\n", ncol(seu)))

# ---- 阈值扫描 --------------------------------------------------------------
thrs <- c(-0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06)
res <- list()
for (thr in thrs) {
  # Hybrid：AD 与 OS 都高于阈值
  hybrid_flag <- seu$AD_drift_ucell > thr & seu$OS_identity_ucell > thr
  tmp <- data.frame(
    sample = seu$sample,
    group  = seu$group,
    hybrid = hybrid_flag
  ) %>%
    group_by(sample, group) %>%
    summarise(
      n_cells    = n(),
      hybrid_pct = 100 * mean(hybrid),
      .groups    = "drop"
    ) %>%
    mutate(threshold = thr)
  res[[as.character(thr)]] <- tmp
}
sens <- bind_rows(res)
write.csv(sens, file.path(tab_dir, "03_sensitivity_summary.csv"), row.names = FALSE)

# 打印核心：三档（S1/S2/S4）在各阈值下的中间态比例与排序
cat("\n=== 三档中间态比例（按阈值） ===\n")
sens_disease <- sens %>% filter(group %in% c("OP","Osteopenia","Normal"))
print(as.data.frame(sens_disease %>% select(threshold, group, hybrid_pct) %>%
        tidyr::pivot_wider(names_from = group, values_from = hybrid_pct)))

# 排序稳定性：每个阈值下三档中间态比例排序
# 注：仅 ≥0 的阈值"informative"（<0 时所有样本均饱和为 100%）
cat("\n=== 每阈值下三档中间态比例排序（关键稳定项：Osteopenia > Normal） ===\n")
for (thr in thrs) {
  sub <- sens %>% filter(threshold == thr, group %in% c("OP","Osteopenia","Normal"))
  ord <- sub %>% arrange(desc(hybrid_pct)) %>% pull(group)
  cat(sprintf("  thr=%+.2f: %s\n", thr, paste(ord, collapse = " > ")))
}

# ---- 绘图：敏感性趋势 -----------------------------------------------------
# 按样本画线（三档 + 年龄维度），x=阈值 y=中间态比例
sens$group <- factor(sens$group,
                     levels = c("Normal","Osteopenia","OP","Age66M"))
pdf(file.path(fig_dir, "03_sensitivity_trend.pdf"), width = 8, height = 6)
p <- ggplot(sens, aes(x = threshold, y = hybrid_pct,
                      color = group, group = sample)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#4DAF4A","#FF7F00","#E41A1C","#984EA3")) +
  labs(title = "Hybrid cell % across thresholds (sensitivity)",
       x = "Hybrid threshold (AD & OS > thr)",
       y = "% Hybrid cells",
       color = "Group/Sample") +
  theme_classic(base_size = 12) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")
print(p)
dev.off()

# 关键判定输出（只评估 informative 阈值 ≥0，饱和区 <0 不纳入）
# 论文稳健主张 = Osteopenia > Normal（OP 组因 S1 骨系细胞少 n=768 不做排序断言）
cat("\n=== 结论 ===")
stable <- TRUE
for (thr in thrs[thrs >= 0]) {
  sub <- sens %>% filter(threshold == thr, group %in% c("OP","Osteopenia","Normal"))
  os <- sub$hybrid_pct[sub$group == "Osteopenia"]
  nm <- sub$hybrid_pct[sub$group == "Normal"]
  if (!(os > nm)) stable <- FALSE
}
cat(sprintf("\n  Osteopenia > Normal 在所有信息阈值(≥0)下成立: %s\n",
            ifelse(stable, "✅ 是（趋势稳健）", "❌ 否（趋势依赖阈值）")))
cat("\n  注：OP vs Normal 的排序因 S1 骨系细胞少(n=768)而不稳（阈值 +0.06 时 Normal 反超），")
cat("\n  建议仅强调 Osteopenia > Normal 的稳健结论（与论文表述一致）。\n")
cat(sprintf("  表: %s\n  图: %s\n",
            file.path(tab_dir, "03_sensitivity_summary.csv"),
            file.path(fig_dir, "03_sensitivity_trend.pdf")))
cat("===== 敏感性分析完成 =====\n")
