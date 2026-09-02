#!/usr/bin/env Rscript
# ============================================================================
# 02_drift_analysis.R — GSE169396 骨系细胞漂移分析（步骤 2，核心分析）
# 项目：方向② 单细胞景观验证"成脂漂移是连续谱"
# 输入：data/processed/seurat_gse169396_bone.rds（步骤 1 产出）
# 输出：
#   results/figures/02_*_*.pdf        — 全部结果图
#   results/tables/02_drift_summary.csv — 每样本漂移统计表
#   results/tables/02_hybrid_table.csv  — 中间态检测表
# 内容：
#   1. UCell 计算细胞级 DriftIndex（AD_drift − OS_identity）
#   2. 中间态检测：AD_drift>0 & OS_identity>0 双高细胞
#   3. 三档趋势：OP vs Osteopenia vs Normal（S1/S2/S4）
#   4. 年龄维度：S3(66M) vs S4(31M)
# 用法：Rscript scripts/02_drift_analysis.R
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ---- 0. 路径 ---------------------------------------------------------------
base_dir <- "D:/WorkBuddySpace/output/Lu-2026-mesenchymal-drift"
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "results/figures")
tab_dir  <- file.path(base_dir, "results/tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

cat("===== 步骤 2：骨系漂移分析（DriftIndex / 中间态 / 三档趋势）=====\n")

# ---- 1. 读取骨系对象 -------------------------------------------------------
seu <- readRDS(file.path(proc_dir, "seurat_gse169396_bone.rds"))
cat(sprintf("骨系对象: %d 细胞 × %d 基因\n", ncol(seu), nrow(seu)))
print(table(seu$sample, seu$group))

# ---- 2. 基因集定义（与 bulk 协议一致，5 组核心基因集）-----------------------
gene_sets <- list(
  OS_identity = c("RUNX2","SP7","ALPL","BGLAP","SPP1","IBSP","DMP1",
                  "COL1A1","SOST","PHEX","BMP2","BMP4"),
  AD_drift    = c("PPARG","CEBPA","CEBPB","ADIPOQ","LEP","FABP4",
                  "PLIN1","LPL","GPD1","RETN","ACSL1"),
  MD_sig      = c("VIM","FN1","COL1A1","SPARC","ACTA2","TAGLN",
                  "POSTN","SNAI2","TWIST1","DKK3"),
  SASP        = c("IL6","IL1B","CXCL8","MMP3","MMP9","SERPINE1",
                  "CCL2","IGFBP3","CXCL1","CXCL2"),
  Osteo_lineage = c("RUNX2","SP7","ALPL","BGLAP","IBSP","SPP1",
                    "COL1A1","SOST","DMP1","PHEX")
)

# 只保留在数据中存在的基因
gene_sets_avail <- lapply(gene_sets, function(g) intersect(g, rownames(seu)))
for (nm in names(gene_sets_avail)) {
  cat(sprintf("  基因集 %s: %d/%d 基因可用\n", nm,
              length(gene_sets_avail[[nm]]), length(gene_sets[[nm]])))
}

# ---- 3. UCell 打分（单细胞级模块分数）--------------------------------------
cat("\n[3] UCell 打分...\n")
seu <- AddModuleScore_UCell(seu, features = gene_sets_avail, name = "_ucell")

# UCell 输出列名：OS_identity_ucell, AD_drift_ucell, ...
score_cols <- paste0(names(gene_sets_avail), "_ucell")
print(head(seu[[score_cols]]))

# 细胞级 DriftIndex = AD_drift − OS_identity
seu$DriftIndex  <- seu$AD_drift_ucell - seu$OS_identity_ucell
seu$MD_score    <- seu$MD_sig_ucell
seu$SASP_score  <- seu$SASP_ucell

# 细胞级"身份模糊度"（MD 特征）与漂移方向的关系
cat(sprintf("DriftIndex 范围: [%.3f, %.3f]\n",
            min(seu$DriftIndex), max(seu$DriftIndex)))

# ---- 4. 中间态检测 ----------------------------------------------------------
# 定义：AD_drift > 0 且 OS_identity > 0（双程序活跃 = 身份模糊/过渡态）
# 注：UCell 分数为秩次标准化，0 为全基因集的中间水平；双高 = 两个程序都高于中位
hybrid_thr <- 0
seu$cell_state <- ifelse(seu$AD_drift_ucell > hybrid_thr & seu$OS_identity_ucell > hybrid_thr,
                         "Hybrid",
                  ifelse(seu$OS_identity_ucell > hybrid_thr, "Osteo-biased",
                  ifelse(seu$AD_drift_ucell > hybrid_thr, "Adipo-biased", "Uncommitted")))
seu$cell_state <- factor(seu$cell_state,
                         levels = c("Uncommitted","Osteo-biased","Adipo-biased","Hybrid"))
cat("\n[4] 细胞状态分布:\n")
print(table(seu$cell_state))

# ---- 5. 三档趋势分析（S1=OP, S2=Osteopenia, S4=Normal；S3 单独年龄维度）-----
cat("\n[5] 每样本漂移统计:\n")
summary_by_sample <- seu@meta.data %>%
  group_by(sample, group) %>%
  summarise(
    n_cells       = n(),
    DriftIndex_mean = mean(DriftIndex),
    DriftIndex_sd   = sd(DriftIndex),
    MD_score_mean   = mean(MD_score),
    SASP_score_mean = mean(SASP_score),
    n_hybrid      = sum(cell_state == "Hybrid"),
    hybrid_pct    = 100 * n_hybrid / n(),
    n_osteo       = sum(cell_state == "Osteo-biased"),
    osteo_pct     = 100 * n_osteo / n(),
    n_adipo       = sum(cell_state == "Adipo-biased"),
    adipo_pct     = 100 * n_adipo / n(),
    .groups = "drop"
  )
print(as.data.frame(summary_by_sample))
write.csv(summary_by_sample, file.path(tab_dir, "02_drift_summary.csv"), row.names = FALSE)

# 三档（疾病梯度）：OP > Osteopenia > Normal
disease_samples <- c("S1", "S2", "S4")
seu_disease <- subset(seu, sample %in% disease_samples)
seu_disease$group <- factor(seu_disease$group,
                            levels = c("Normal","Osteopenia","OP"))

# 年龄维度：S3(66M) vs S4(31M)
seu_age <- subset(seu, sample %in% c("S3","S4"))
seu_age$age_label <- ifelse(seu_age$sample == "S3", "66M", "31M")

# ---- 6. 绘图 ---------------------------------------------------------------
theme_set(theme_classic(base_size = 12))

# 图 A1：UMAP 按 DriftIndex 着色（三档 + 年龄）
pdf(file.path(fig_dir, "02_UMAP_DriftIndex.pdf"), width = 16, height = 8)
p1 <- FeaturePlot(seu, features = "DriftIndex", reduction = "umap.bone",
                  cols = c("#377EB8", "#FFFFFF", "#E41A1C")) +
        ggtitle("Cell-level DriftIndex (AD - OS)")
p2 <- DimPlot(seu, reduction = "umap.bone", group.by = "group",
              cols = c("#4DAF4A", "#FF7F00", "#E41A1C", "#984EA3")) +
        ggtitle("Disease group")
p3 <- DimPlot(seu, reduction = "umap.bone", group.by = "sample") +
        ggtitle("Sample")
print(p1 | p2 | p3)
dev.off()

# 图 A2：细胞状态（Hybrid/成骨/成脂/未定型）在 UMAP 上
pdf(file.path(fig_dir, "02_UMAP_cell_state.pdf"), width = 12, height = 6)
DimPlot(seu, reduction = "umap.bone", group.by = "cell_state",
        cols = c("#BDBDBD","#377EB8","#984EA3","#E41A1C")) +
  ggtitle("Cell state: Hybrid = AD+ & OS+")
dev.off()

# 图 B1：三档 DriftIndex 分布（小提琴 + 箱线）
pdf(file.path(fig_dir, "02_DriftIndex_by_group.pdf"), width = 8, height = 6)
ggplot(seu_disease@meta.data, aes(x = group, y = DriftIndex, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3) +
  scale_fill_manual(values = c("#4DAF4A","#FF7F00","#E41A1C")) +
  labs(title = "DriftIndex across bone-mass gradient",
       x = "Group (by lumbar T-score)", y = "Cell-level DriftIndex") +
  theme(legend.position = "none")
dev.off()

# 图 B2：中间态比例按三档（堆叠条形）
hybrid_summary <- seu_disease@meta.data %>%
  group_by(group) %>%
  summarise(Hybrid = 100*sum(cell_state=="Hybrid")/n(),
            Osteo  = 100*sum(cell_state=="Osteo-biased")/n(),
            Adipo  = 100*sum(cell_state=="Adipo-biased")/n(),
            Uncommitted = 100*sum(cell_state=="Uncommitted")/n(),
            .groups="drop") %>%
  tidyr::pivot_longer(-group, names_to = "state", values_to = "pct")
hybrid_summary$state <- factor(hybrid_summary$state,
                               levels = c("Uncommitted","Osteo","Adipo","Hybrid"))

pdf(file.path(fig_dir, "02_cell_state_proportion.pdf"), width = 8, height = 6)
ggplot(hybrid_summary, aes(x = group, y = pct, fill = state)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("#BDBDBD","#377EB8","#984EA3","#E41A1C")) +
  labs(title = "Cell state proportions across bone-mass gradient",
       x = "Group", y = "% of bone-lineage cells") +
  theme(legend.position = "right")
dev.off()

# 图 C1：年龄维度 DriftIndex（S3 66M vs S4 31M）
pdf(file.path(fig_dir, "02_DriftIndex_age_dimension.pdf"), width = 7, height = 6)
ggplot(seu_age@meta.data, aes(x = age_label, y = DriftIndex, fill = age_label)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3) +
  scale_fill_manual(values = c("#66C2A5","#FC8D62")) +
  labs(title = "Age dimension: 66M vs 31M",
       x = "Sample", y = "Cell-level DriftIndex") +
  theme(legend.position = "none")
dev.off()

# 图 C2：年龄维度中间态比例
age_hybrid <- seu_age@meta.data %>%
  group_by(age_label) %>%
  summarise(Hybrid = 100*sum(cell_state=="Hybrid")/n(), .groups="drop")
pdf(file.path(fig_dir, "02_hybrid_pct_age.pdf"), width = 6, height = 5)
ggplot(age_hybrid, aes(x = age_label, y = Hybrid, fill = age_label)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("#66C2A5","#FC8D62")) +
  labs(title = "Hybrid cell % by age", x = "Sample", y = "% Hybrid cells") +
  theme(legend.position = "none")
dev.off()

# 图 D：中间态特征热图（Hybrid vs 其他状态的标志基因）
pdf(file.path(fig_dir, "02_hybrid_signature_heatmap.pdf"), width = 10, height = 8)
genes_heat <- unique(unlist(gene_sets_avail[c("OS_identity","AD_drift","MD_sig","SASP")]))
genes_heat <- intersect(genes_heat, rownames(seu))
if (length(genes_heat) > 5) {
  # 随机抽样加速（每状态最多 300 细胞）
  set.seed(42)
  cells_sub <- unlist(lapply(split(rownames(seu@meta.data), seu$cell_state), function(x) {
    sample(x, size = min(300, length(x)))
  }))
  print(DoHeatmap(seu, cells = cells_sub, features = genes_heat,
                  group.by = "cell_state", size = 4) +
          ggtitle("Hybrid state signature"))
}
dev.off()

# ---- 7. 保存 ---------------------------------------------------------------
saveRDS(seu, file.path(proc_dir, "seurat_gse169396_drift.rds"))
cat("\n===== 步骤 2 完成 =====")
cat(sprintf("\n  漂移对象: %s/seurat_gse169396_drift.rds", proc_dir))
cat(sprintf("\n  统计表:   %s/02_drift_summary.csv", tab_dir))
cat(sprintf("\n  图表:     %s/02_*.pdf（%d 张）\n", fig_dir,
            length(list.files(fig_dir, pattern = "^02_"))))
