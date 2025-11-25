library(tidyverse)
library(ggfortify)

# 加载数据
basedir <- "/disk2/cai113/data/stateTrans"
vst_file <- file.path(basedir, "04.standardization/03.vst_matrix.txt")
cat("Loading VST matrix from ", vst_file, "...\n")
vst_mat <- as.matrix(read.delim(vst_file, row.names=1, check.names = FALSE))
cat("VST matrix loaded: ", nrow(vst_mat), " genes, ", ncol(vst_mat), " samples\n")

clinic_info <- read.delim("/disk2/cai113/data/stateTrans/00.clinic.txt", check.names=FALSE)
sample_info <- clinic_info[match(colnames(vst_mat), clinic_info$sample), ]
rownames(sample_info) <- sample_info$sample

# 转置：样本为行，基因为列（prcomp 要求）
X <- t(vst_mat)  # n x p, n=样本数, p=基因数

# 中心化（prcomp 默认 center=TRUE，但显式写出更清晰）
X_centered <- scale(X, center = TRUE, scale = FALSE)

# PCA via SVD (prcomp uses SVD under the hood)
pca_result <- prcomp(X_centered, center = FALSE, scale. = FALSE)

# 提取主成分得分（样本在PC空间的坐标）
pcs <- pca_result$x  # n x k

# 查看方差解释比例
var_explained <- summary(pca_result)$importance[2, ]
cum_var <- cumsum(var_explained)

# 打印前10个PC的累计方差
print(data.frame(PC = 1:10, CumVar = cum_var[1:10]))

# 添加 PCA 得分到 sample_info
sample_info_pcs <- cbind(sample_info, pcs)

# 计算每个PC在KO组的时间趋势（线性回归斜率 & p值）
pc_candidates <- colnames(pcs)
ko_samples <- sample_info_pcs$genotype == "Ncf2"
wt_samples <- sample_info_pcs$genotype == "WT"

candidate_scores <- data.frame(
  PC = character(),
  ko_pval = numeric(),
  ko_slope = numeric(),
  wt_pval = numeric(),
  wt_slope = numeric(),
  stringsAsFactors = FALSE
)

for (pc in pc_candidates[1:which(cum_var >= 0.8)[1]]) {  # 只考虑前80%方差内的PC
  # KO 组：lung_weight 对 PC 回归
  fit_ko <- lm(pcs[ko_samples, pc] ~ sample_info_pcs$lung_weight[ko_samples])
  ko_p <- summary(fit_ko)$coefficients[2, 4]
  ko_slope <- coef(fit_ko)[2]
  
  # WT 组
  fit_wt <- lm(pcs[wt_samples, pc] ~ sample_info_pcs$lung_weight[wt_samples])
  wt_p <- summary(fit_wt)$coefficients[2, 4]
  wt_slope <- coef(fit_wt)[2]
  
  candidate_scores <- rbind(candidate_scores, data.frame(
    PC = pc,
    ko_pval = ko_p,
    ko_slope = ko_slope,
    wt_pval = wt_p,
    wt_slope = wt_slope
  ))
}

# 排序：优先选 KO 显著（p < 0.05）、WT 不显著（p > 0.05）、|slope| 大的
candidate_scores <- candidate_scores %>%
  mutate(ko_sig = ko_pval < 0.05,
         wt_not_sig = wt_pval > 0.05,
         abs_slope = abs(ko_slope)) %>%
  arrange(desc(ko_sig), desc(wt_not_sig), desc(abs_slope))

print(candidate_scores)

# 原 PC1 与 lung_weight 正相关 → 越大越病
# 现在取反：x_cgd 越小 = 越病
x_cgd <- -pcs[, "PC1"]  # 关键修改！

sample_info_pcs$x_cgd <- x_cgd

# 验证：KO 中 x_cgd 应与 lung_weight 负相关
cor.test(x_cgd[ko_samples], sample_info_pcs$lung_weight[ko_samples])
# 得到显著负相关（p < 0.05, r < 0）

save(x_cgd,sample_info_pcs,file = file.path(basedir, "05.model/sample_info_with_PCs.rdata"))
