#!/usr/bin/env Rscript
# ============================================================================
# 01_qc_cluster.R — GSE169396 人股骨头 scRNA-seq
# 步骤 1：QC 聚类 + Harmony 整合 + 骨系提取
# 项目：方向② 单细胞景观验证"成脂漂移是连续谱"
# 数据：GSE169396（4 样本，GSM5201883-S1 ~ GSM5201886-S4，10x v3）
# 输出：data/processed/seurat_gse169396_qc.rds
#        results/figures/QC_*.pdf, UMAP_*.pdf
# 用法：Rscript scripts/01_qc_cluster.R
# ============================================================================

# ---- 0. 环境与路径 ---------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

base_dir  <- "D:/WorkBuddySpace/output/Lu-2026-mesenchymal-drift"
raw_dir   <- file.path(base_dir, "data/raw")
proc_dir  <- file.path(base_dir, "data/processed")
fig_dir   <- file.path(base_dir, "results/figures")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)

cat("===== GSE169396 步骤 1：QC + 整合 + 骨系提取 =====\n")

# ---- 1. 读取 4 个样本的 10x 矩阵 -------------------------------------------
# GSM5201883=S1, GSM5201884=S2, GSM5201885=S3, GSM5201886=S4
sample_ids <- c("S1", "S2", "S3", "S4")
gsm_ids    <- c("GSM5201883", "GSM5201884", "GSM5201885", "GSM5201886")

# 临床信息（2026-09-02 定稿：腰椎T值三档 + S3 年龄维度保留）
# S1 = GSM5201883: 61y 女，腰椎T=-3.0 → 骨质疏松(OP)
# S2 = GSM5201884: 45y 女，腰椎T=-1.3 → 骨量减少(Osteopenia)
# S3 = GSM5201885: 66y 男，T值NA → 不进疾病三档，作年龄维度样本（66M vs 31M）
# S4 = GSM5201886: 31y 男，腰椎T=+0.6 → 正常(Normal)
age_gender <- c(S1 = "61F", S2 = "45F", S3 = "66M", S4 = "31M")
t_lumbar   <- c(S1 = -3.0, S2 = -1.3, S3 = NA,  S4 = 0.6)
t_hip      <- c(S1 = -1.9, S2 = -1.2, S3 = NA,  S4 = -1.1)
group_map  <- c(S1 = "OP", S2 = "Osteopenia", S3 = "Age66M", S4 = "Normal")
cat("\n[样本信息] 腰椎T值三档(S1/S2/S4) + S3年龄维度:\n")
print(data.frame(AgeSex = age_gender, Lumbar_T = t_lumbar, Hip_T = t_hip, Group = group_map))

# Seurat 5.x 的 Read10X 不再支持 prefix 参数：
# 为每个样本建子目录并重命名为标准 10x 文件名（barcodes/features/matrix）
prep_dir <- file.path(raw_dir, "10x_standard")
dir.create(prep_dir, recursive = TRUE, showWarnings = FALSE)

obj_list <- list()
for (i in seq_along(sample_ids)) {
  sid  <- sample_ids[i]
  gsm  <- gsm_ids[i]

  # 建子目录：data/raw/10x_standard/S1/
  sample_dir <- file.path(prep_dir, sid)
  dir.create(sample_dir, showWarnings = FALSE)

  # 重命名/复制三件套为标准 10x 文件名
  for (ext in c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")) {
    src <- file.path(raw_dir, paste0(gsm, "_", sid, "_", ext))
    dst <- file.path(sample_dir, ext)
    if (!file.exists(dst)) {
      file.copy(src, dst, overwrite = TRUE)
    }
  }

  mtx  <- Seurat::Read10X(data.dir = sample_dir)
  obj  <- CreateSeuratObject(counts = mtx,
                             project = sid,
                             min.cells = 3,      # 至少在 3 个细胞中表达
                             min.features = 200) # 至少 200 个基因
  obj[["sample"]] <- sid
  obj[["group"]]  <- unname(group_map[sid])  # 用 [[ 或 unname 去掉名字，避免 Seurat5 metadata 匹配报错
  obj_list[[sid]] <- obj
  cat(sprintf("  样本 %s (%s): %d 细胞 × %d 基因\n",
              sid, gsm, ncol(obj), nrow(obj)))
}

# 合并
seu <- merge(obj_list[[1]], y = obj_list[2:length(obj_list)])
cat(sprintf("\n合并后: %d 细胞 × %d 基因\n", ncol(seu), nrow(seu)))

# ---- 2. QC 指标与过滤 ------------------------------------------------------
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
# 人骨髓/骨组织线粒体基因高表达常见，阈值适当放宽（<20%）

pdf(file.path(fig_dir, "01_QC_violin_prefilter.pdf"), width = 12, height = 5)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0, group.by = "sample") &
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# 过滤：基因数 500-6000，线粒体 <20%
seu <- subset(seu, subset = nFeature_RNA > 500 &
                        nFeature_RNA < 6000 &
                        percent.mt < 20)
cat(sprintf("QC 过滤后: %d 细胞\n", ncol(seu)))

pdf(file.path(fig_dir, "01_QC_violin_postfilter.pdf"), width = 12, height = 5)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0, group.by = "sample") &
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# ---- 3. 标准化 + 高变基因 + 整合 ------------------------------------------
# 3.1 SCTransform（v2 更快）——按样本独立标准化
seu <- SCTransform(seu, vars.to.regress = "percent.mt",
                   verbose = FALSE)

# 3.2 Harmony 整合（去除 4 样本技术批次）
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- RunHarmony(seu, group.by.vars = "sample",
                  reduction.save = "harmony",
                  verbose = FALSE)
cat("Harmony 整合完成\n")

# 3.3 UMAP + 聚类（基于 harmony 嵌入）
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30,
               reduction.name = "umap.harmony", verbose = FALSE)
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
cat(sprintf("聚类完成: %d 个 cluster\n", length(unique(Idents(seu)))))

# ---- 4. 可视化（整合前后对比）---------------------------------------------
pdf(file.path(fig_dir, "01_UMAP_integration_compare.pdf"), width = 14, height = 7)
p1 <- DimPlot(seu, reduction = "umap.harmony", group.by = "sample",
              label = FALSE) + ggtitle("After Harmony integration - by sample")
p2 <- DimPlot(seu, reduction = "umap.harmony", group.by = "group",
              label = FALSE) + ggtitle("After Harmony integration - by group")
p3 <- DimPlot(seu, reduction = "umap.harmony", label = TRUE) +
        ggtitle("After Harmony integration - clusters")
print(p1 + p2 + p3)
dev.off()

# 经典标志基因点图（辅助注释）
markers_check <- c("PTPRC","CD3D","CD79A","LYZ","HBB",      # 免疫/血液
                   "VWF","PECAM1",                            # 内皮
                   "COL1A1","LEPR","RUNX2","SP7","ALPL",      # 骨系：MSC/成骨
                   "BGLAP","IBSP","SPP1","SOST","DMP1","PHEX",# 骨系：成骨/骨细胞
                   "PPARG","ADIPOQ","FABP4","LPL","PLIN1",    # 成脂
                   "COL2A1","SOX9")                           # 软骨
markers_check <- intersect(markers_check, rownames(seu))
pdf(file.path(fig_dir, "01_DotPlot_lineage_markers.pdf"), width = 14, height = 8)
DotPlot(seu, features = markers_check, group.by = "seurat_clusters") +
  RotatedAxis() + ggtitle("Canonical lineage marker DotPlot")
dev.off()

# ---- 5. 骨系细胞提取 ------------------------------------------------------
# 5.1 用标志基因分数锚定骨系（LEPR/RUNX2/SP7/BGLAP/SOST 等）
osteo_genes <- intersect(c("RUNX2","SP7","ALPL","BGLAP","IBSP","SPP1","COL1A1",
                           "SOST","DMP1","PHEX"), rownames(seu))
adipo_genes <- intersect(c("PPARG","ADIPOQ","FABP4","LPL","PLIN1","CEBPA"), rownames(seu))
msc_genes   <- intersect(c("LEPR","PDGFRB","NT5E","ENG"), rownames(seu))

seu <- AddModuleScore(seu, features = list(osteo_genes),
                      name = "OsteoScore", assay = "SCT")
seu <- AddModuleScore(seu, features = list(adipo_genes),
                      name = "AdipoScore", assay = "SCT")
seu <- AddModuleScore(seu, features = list(msc_genes),
                      name = "MscScore", assay = "SCT")

# 5.2 骨系候选：MSC 分数 > 0 或 成骨分数 > 0（排除免疫/红细胞/内皮/软骨）
bone_cells <- WhichCells(seu, expression = MscScore1 > 0 | OsteoScore1 > 0)
seu_bone <- subset(seu, cells = bone_cells)
cat(sprintf("\n骨系细胞提取: %d / %d (%.1f%%)\n",
            ncol(seu_bone), ncol(seu), 100 * ncol(seu_bone) / ncol(seu)))

# 5.3 骨系亚群再聚类（分辨率更高）
seu_bone <- RunUMAP(seu_bone, reduction = "harmony", dims = 1:30,
                    reduction.name = "umap.bone", verbose = FALSE)
seu_bone <- FindNeighbors(seu_bone, reduction = "harmony", dims = 1:30, verbose = FALSE)
seu_bone <- FindClusters(seu_bone, resolution = 0.8, verbose = FALSE)

pdf(file.path(fig_dir, "01_UMAP_bone_lineage.pdf"), width = 14, height = 7)
p1 <- DimPlot(seu_bone, reduction = "umap.bone", group.by = "seurat_clusters",
              label = TRUE) + ggtitle("Bone-lineage cell clusters")
p2 <- FeaturePlot(seu_bone, features = c("RUNX2","SP7","BGLAP","SOST","LEPR","PPARG"),
                  reduction = "umap.bone", ncol = 3) &
        theme(plot.title = element_text(size = 10))
print(p1 / p2)
dev.off()

# 5.4 保存
saveRDS(seu,       file.path(proc_dir, "seurat_gse169396_all.rds"))
saveRDS(seu_bone,  file.path(proc_dir, "seurat_gse169396_bone.rds"))
cat("\n===== 步骤 1 完成 =====")
cat(sprintf("\n  全细胞对象: %s/seurat_gse169396_all.rds", proc_dir))
cat(sprintf("\n  骨系对象:   %s/seurat_gse169396_bone.rds", proc_dir))
cat(sprintf("\n  图表输出:   %s\n", fig_dir))
