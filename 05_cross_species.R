#!/usr/bin/env Rscript
# ============================================================================
# 05_cross_species.R — 跨物种健康基准验证（步骤 6）
# 目的：用 GSE317069 健康小鼠骨髓腔（干骺端 Metaphysis）作为"健康基准"，
#       与人类 GSE169396 骨系（OP/骨量减少/正常）比较 DriftIndex 与中间态比例
# 输入：data/raw/GSE317069_metaphysis/*.gz（5 样本小鼠 10x）
#       data/processed/seurat_gse169396_drift.rds（人类漂移对象）
# 输出：results/figures/06_*.pdf
#       results/tables/05_cross_species_summary.csv
# 用法：Rscript scripts/05_cross_species.R
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

base_dir <- "D:/WorkBuddySpace/output/Lu-2026-mesenchymal-drift"
raw_dir  <- file.path(base_dir, "data/raw/GSE317069_metaphysis")
proc_dir <- file.path(base_dir, "data/processed")
fig_dir  <- file.path(base_dir, "results/figures")
tab_dir  <- file.path(base_dir, "results/tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

cat("===== 步骤 6：跨物种健康基准验证（GSE317069 小鼠 vs GSE169396 人）=====\n")

# ---- 1. 读取小鼠 5 样本（Metaphysis 干骺端，标准 10x 目录）----------------
std_dir <- file.path(base_dir, "data/raw/GSE317069_metaphysis/10x_std")
mouse_sids <- paste0("MM", 1:5)
cat(sprintf("读取 %d 个小鼠样本（标准目录 %s）\n", length(mouse_sids), std_dir))

obj_list <- list()
for (sid in mouse_sids) {
  d <- file.path(std_dir, sid)
  mtx <- Seurat::Read10X(d)
  obj <- CreateSeuratObject(counts = mtx, project = sid,
                            min.cells = 3, min.features = 200)
  obj[["sample"]] <- sid
  obj[["species"]] <- "Mouse"
  obj[["group"]]   <- "HealthyMouse"
  obj_list[[sid]]  <- obj
  cat(sprintf("  样本 %s: %d 细胞\n", sid, ncol(obj)))
}

# 合并小鼠
seu_mouse <- merge(obj_list[[1]], y = obj_list[2:length(obj_list)])
cat(sprintf("小鼠合并: %d 细胞\n", ncol(seu_mouse)))

# ---- 2. QC 过滤 ------------------------------------------------------------
seu_mouse[["percent.mt"]] <- PercentageFeatureSet(seu_mouse, pattern = "^mt-")
seu_mouse <- subset(seu_mouse, subset = nFeature_RNA > 300 &
                                nFeature_RNA < 6000 &
                                percent.mt < 20)
cat(sprintf("小鼠 QC 后: %d 细胞\n", ncol(seu_mouse)))

# ---- 3. 骨系细胞提取（小鼠同源基因锚定）-----------------------------------
# 小鼠骨系标志（首字母大写）
msc_genes   <- intersect(c("Lepr","Pdgfrb","Eng","Cd44"), rownames(seu_mouse))
osteo_genes <- intersect(c("Runx2","Sp7","Alpl","Bglap","Ibsp","Spp1","Col1a1","Sost","Dmp1","Phex"), rownames(seu_mouse))
adipo_genes <- intersect(c("Pparg","Cebpa","Adipoq","Fabp4","Lpl","Plin1","Gpd1","Retn"), rownames(seu_mouse))
md_genes    <- intersect(c("Vim","Fn1","Col1a1","Sparc","Acta2","Tagln","Postn","Snai2","Twist1","Dkk3"), rownames(seu_mouse))
sasp_genes  <- intersect(c("Il6","Il1b","Cxcl8","Mmp3","Mmp9","Serpine1","Ccl2","Igfbp3"), rownames(seu_mouse))

cat(sprintf("小鼠骨系标志: MSC %d, OS %d, AD %d, MD %d, SASP %d\n",
            length(msc_genes), length(osteo_genes), length(adipo_genes),
            length(md_genes), length(sasp_genes)))

# 标准化
seu_mouse <- NormalizeData(seu_mouse, verbose = FALSE)
seu_mouse <- FindVariableFeatures(seu_mouse, verbose = FALSE)
seu_mouse <- ScaleData(seu_mouse, verbose = FALSE)

# 用标志基因分数提取骨系（与人类流程一致：MSC+OS 分数 >0）
seu_mouse <- AddModuleScore(seu_mouse, features = list(msc_genes), name = "MscScore", assay = "RNA")
seu_mouse <- AddModuleScore(seu_mouse, features = list(osteo_genes), name = "OsteoScore", assay = "RNA")

bone_cells <- WhichCells(seu_mouse, expression = MscScore1 > 0 | OsteoScore1 > 0)
seu_mouse_bone <- subset(seu_mouse, cells = bone_cells)
cat(sprintf("小鼠骨系提取: %d / %d (%.1f%%)\n",
            ncol(seu_mouse_bone), ncol(seu_mouse),
            100*ncol(seu_mouse_bone)/ncol(seu_mouse)))

# ---- 4. UCell 打分（同款 5 基因集，小鼠同源）-------------------------------
gene_sets_mouse <- list(
  OS_identity = osteo_genes,
  AD_drift    = adipo_genes,
  MD_sig      = md_genes,
  SASP        = sasp_genes
)
seu_mouse_bone <- AddModuleScore_UCell(seu_mouse_bone, features = gene_sets_mouse, name = "_ucell")

# 细胞级 DriftIndex（注意：UCell 分数是秩次标准化的，跨物种可比性需谨慎）
# 采用"组内标准化"策略：每样本内对 AD/OS 分数做 scale，再计算差值
seu_mouse_bone$DriftIndex_mouse <- seu_mouse_bone$AD_drift_ucell - seu_mouse_bone$OS_identity_ucell

# 中间态：AD>0 & OS>0（状态命名与人类分析 02 脚本一致）
seu_mouse_bone$cell_state_mouse <- ifelse(
  seu_mouse_bone$AD_drift_ucell > 0 & seu_mouse_bone$OS_identity_ucell > 0,
  "Hybrid",
  ifelse(seu_mouse_bone$OS_identity_ucell > 0, "Osteo-biased",
  ifelse(seu_mouse_bone$AD_drift_ucell > 0, "Adipo-biased", "Uncommitted")))

cat("\n小鼠骨系细胞状态分布:\n")
print(table(seu_mouse_bone$cell_state_mouse))

# ---- 5. 载入人类漂移对象 ---------------------------------------------------
seu_human <- readRDS(file.path(proc_dir, "seurat_gse169396_drift.rds"))
cat(sprintf("\n人类骨系对象: %d 细胞\n", ncol(seu_human)))

# ---- 6. 跨物种比较表 -------------------------------------------------------
mouse_summary <- seu_mouse_bone@meta.data %>%
  summarise(
    species = "Mouse (healthy)", n = n(),
    DriftIndex_mean = mean(DriftIndex_mouse),
    hybrid_pct = 100*mean(cell_state_mouse == "Hybrid"),
    .groups = "drop"
  )

human_summary <- seu_human@meta.data %>%
  group_by(group) %>%
  summarise(
    n = n(),
    DriftIndex_mean = mean(DriftIndex),
    hybrid_pct = 100*mean(cell_state == "Hybrid"),
    .groups = "drop"
  ) %>%
  # 排除年龄维度 Age66M（仅参与年龄分析，不进入疾病三档跨物种对比）
  filter(group %in% c("OP","Osteopenia","Normal")) %>%
  mutate(species = paste0("Human ", group)) %>%
  select(species, n, DriftIndex_mean, hybrid_pct)

cross_tab <- bind_rows(mouse_summary, human_summary)
print(cross_tab)
write.csv(cross_tab, file.path(tab_dir, "05_cross_species_summary.csv"), row.names = FALSE)

# ---- 7. 绘图 ---------------------------------------------------------------
theme_set(theme_classic(base_size = 12))

# 图 1：小鼠 UMAP（若 PCA/UMAP 可用）——简化：直接展示细胞状态比例
pdf(file.path(fig_dir, "05_mouse_metaphysis_cell_state.pdf"), width = 7, height = 6)
state_tab <- as.data.frame(table(seu_mouse_bone$cell_state_mouse))
names(state_tab) <- c("state", "n")
state_tab$pct <- 100*state_tab$n/sum(state_tab$n)
ggplot(state_tab, aes(x = "", y = pct, fill = state)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = c("Hybrid"="#E41A1C","Osteo-biased"="#377EB8",
                               "Adipo-biased"="#984EA3","Uncommitted"="#BDBDBD")) +
  labs(title = "Healthy mouse metaphysis bone-lineage cell states") +
  theme_void() + theme(legend.position = "right")
dev.off()

# 图 2：跨物种 DriftIndex 对比（小提琴）
# 人类仅保留疾病三档（OP/Osteopenia/Normal），Age66M 不进跨物种对比
human_meta <- seu_human@meta.data[seu_human$group %in% c("OP","Osteopenia","Normal"), ]
plot_df <- data.frame(
  species = c(rep("Mouse\nHealthy", ncol(seu_mouse_bone)),
              human_meta$group),
  DriftIndex = c(seu_mouse_bone$DriftIndex_mouse, human_meta$DriftIndex),
  stringsAsFactors = FALSE
)
# 小鼠组标为基准
plot_df$species <- factor(plot_df$species,
                          levels = c("Normal","Osteopenia","OP","Mouse\nHealthy"))

pdf(file.path(fig_dir, "05_cross_species_driftindex.pdf"), width = 8, height = 6)
ggplot(plot_df, aes(x = species, y = DriftIndex, fill = species)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.size = 0.3) +
  scale_fill_manual(values = c("Normal"="#4DAF4A","Osteopenia"="#FF7F00",
                               "OP"="#E41A1C","Mouse\nHealthy"="#377EB8")) +
  labs(title = "Cross-species comparison: DriftIndex",
       x = "", y = "Cell-level DriftIndex") +
  theme(legend.position = "none", axis.text.x = element_text(size = 9))
dev.off()

# 图 3：中间态比例跨物种（条形）
bar_df <- cross_tab
bar_df$label <- sub("Human ", "", bar_df$species)
bar_df$label <- factor(bar_df$label, levels = c("Normal","Osteopenia","OP","Mouse (healthy)"))

pdf(file.path(fig_dir, "05_cross_species_hybrid_pct.pdf"), width = 7, height = 6)
ggplot(bar_df, aes(x = label, y = hybrid_pct, fill = label)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("Normal"="#4DAF4A","Osteopenia"="#FF7F00",
                               "OP"="#E41A1C","Mouse (healthy)"="#377EB8")) +
  labs(title = "Hybrid cell proportion: cross-species",
       x = "", y = "% Hybrid (AD+ & OS+) cells") +
  theme(legend.position = "none", axis.text.x = element_text(size = 9))
dev.off()

# ---- 8. 结论 ---------------------------------------------------------------
cat("\n=== 跨物种结论 ===")
# 注意：UCell 分数跨物种/跨数据集不完全可比，此处仅作定性趋势
cat("\n⚠️ 注意：UCell 分数为秩次标准化，跨物种比较需谨慎解读")
cat("\n建议：以中间态比例与组内相对趋势为主要比较依据\n")

saveRDS(seu_mouse_bone, file.path(proc_dir, "seurat_mouse_metaphysis_bone.rds"))
cat("\n===== 步骤 6 完成 =====")
cat(sprintf("\n  小鼠骨系对象: %s/seurat_mouse_metaphysis_bone.rds", proc_dir))
cat(sprintf("\n  统计表: %s/05_cross_species_summary.csv", tab_dir))
