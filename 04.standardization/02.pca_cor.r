#!/usr/bin/env Rscript
library(optparse)
library(pheatmap)

#--------------------- 1. 解析命令行参数
option_list <- list(
  make_option(c("-g", "--group"), type="character", default="group",
              help="Grouping column in sample_info: group, genotype, week [default %default]",
              metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
group_by <- opt$group

cat("Grouping by column:", group_by, "\n")

#--------------------- 2. 加载数据
outdir <- "/disk2/cai113/data/stateTrans/04.standardization"
vst_file <- file.path(outdir, "03.vst_matrix.txt")
cat("Loading VST matrix from ", vst_file, "...\n")
vst_mat <- as.matrix(read.delim(vst_file, row.names=1))
cat("VST matrix loaded: ", nrow(vst_mat), " genes, ", ncol(vst_mat), " samples\n")

clinic_info <- read.delim("/disk2/cai113/data/stateTrans/00.clinic.txt", check.names=FALSE)
sample_info <- clinic_info[match(colnames(vst_mat), clinic_info$sample), ]
rownames(sample_info) <- sample_info$sample
sample_info$group <- paste(sample_info$genotype, sample_info$week, sep="_")
sample_info$group <- factor(sample_info$group)
sample_info$genotype <- factor(sample_info$genotype)
sample_info$week <- factor(sample_info$week)

# 检查 group_by 是否存在
if (!group_by %in% colnames(sample_info)) {
  stop(paste("Column", group_by, "not found in sample_info"))
}

cat("Sample table:\n")
print(sample_info)

#--------------------- 3. PCA
pca <- prcomp(t(vst_mat))
percentVar <- pca$sdev^2 / sum(pca$sdev^2)

pdf(file.path(outdir, paste0("PCA_plot_", group_by, ".pdf")), width=7, height=6)
plot(pca$x[,1], pca$x[,2],
     xlab=paste0("PC1 (", round(percentVar[1]*100,1),"%)"),
     ylab=paste0("PC2 (", round(percentVar[2]*100,1),"%)"),
     pch=19, col=sample_info[[group_by]])
text(pca$x[,1], pca$x[,2], labels=rownames(sample_info), pos=3, cex=0.7)
legend("topright", legend=levels(sample_info[[group_by]]),
       col=1:length(levels(sample_info[[group_by]])), pch=19)
dev.off()

#--------------------- 4. 样本相关性热图
cor_mat <- cor(vst_mat)
pdf(file.path(outdir, paste0("correlation_heatmap_", group_by, ".pdf")), width=7, height=6)
pheatmap(cor_mat, main=paste("Sample Correlation (VST) by", group_by),
         annotation_col=sample_info[group_by])
dev.off()

cat("PCA and heatmap generated using '", group_by, "' for sample grouping.\n", sep="")
