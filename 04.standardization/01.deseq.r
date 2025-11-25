library(DESeq2)
library(tidyverse)

#--------------------- 1. 路径设置
    base_dir <- "/disk2/cai113/data/stateTrans"
    count_file <- file.path(base_dir, "03.quantify", "gene_counts.txt")
    clinic_file <- file.path(base_dir, "00.clinic.txt")
    outdir <- file.path(base_dir, "04.standardization")
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
#--------------------- 2. 读取数据，构建 count_mat
    clinic_info <- read.delim(clinic_file, check.names = FALSE) # load clinic info
    cat("Loading featureCounts matrix...\n")
    raw <- read.delim(count_file, comment.char="#",check.names = FALSE)
    count_mat <- raw[, 7:ncol(raw)]# featureCounts 前6列是注释
    # 处理列名
        extract_sample <- function(full_path) {
            # 取路径中最后一级目录
            dir <- basename(dirname(full_path))
            # 从目录中提取最后的 "-xxx" 样本名
            # 匹配 "-" 后面的部分，如 w4_1, n4_2, left_w5_1 等
            m <- regexpr("[-]([A-Za-z0-9_]+)$", dir)
            if (m == -1) return(NA)
            substring(dir, m + 1, nchar(dir))
        }# 1) 从 count_mat 列名中提取 sample 名称
        raw_colnames <- colnames(count_mat)
        sample_names <- sapply(raw_colnames, extract_sample)
        if (any(is.na(sample_names))) {
            warning("有列名无法从路径中提取 sample 名称：")
            print(raw_colnames[is.na(sample_names)])
            stop("请检查目录结构或正则表达式。")
        }# 2) 检查是否有无法提取的样本名
        missing <- setdiff(sample_names, clinic_info$sample)
        if (length(missing) > 0) {
            warning("以下 sample 无法在 clinic_info 中找到，请检查命名是否一致：")
            print(missing)
            stop("存在无法匹配的 sample，已停止列名替换。")
        }# 3) 检查提取的 sample 是否全部能在 clinic_info 中找到
        colnames(count_mat) <- sample_names# 4) 没问题则替换列名
    rownames(count_mat) <- raw$Geneid
    count_mat <- count_mat[!duplicated(rownames(count_mat)), ]# 去掉可能重复行
    cat("Matrix loaded: ", nrow(count_mat), " genes, ", ncol(count_mat), " samples\n")
#--------------------- 3. 构建样本信息表
    sample_info <- clinic_info[match(colnames(count_mat), clinic_info$sample), ]# 1) 按 count_mat 的列名顺序匹配 clinic_info
    if (any(is.na(sample_info$sample))) {
        missing <- colnames(count_mat)[is.na(sample_info$sample)]
        stop("以下样本在 clinic_info 中未找到，请检查 clinic.txt：\n",
            paste(missing, collapse = ", "))
    }# 2) 检查是否有无法匹配的样本
    rownames(sample_info) <- sample_info$sample# 3) 设置 rownames，方便 DESeq2 使用
    sample_info$group <- paste(sample_info$genotype, sample_info$week, sep = "_")# 4) 添加复合分组：genotype_week
    cat("Final sample_info:\n")
    print(sample_info)
#--------------------- 4. 构建 DESeq2 数据集
    dds <- DESeqDataSetFromMatrix(
        countData = count_mat,
        colData = sample_info,
        design = ~ genotype + week #不建议直接用group作为 DESeq2 的 design，因为它会产生太多比较，用这个分析NCF2 deficiency 在不同时间点的动态变化足够
    )
#--------------------- 5. 过滤低表达基因
    keep <- rowSums(counts(dds) >= 10) >= 3 # （至少 3 个样本 count >= 10）
    dds <- dds[keep, ]
    cat("Genes kept after filtering: ", nrow(dds), "\n")
#--------------------- 6. DESeq2 normalization（estimateSizeFactors）
    dds <- estimateSizeFactors(dds)
    norm_counts <- counts(dds, normalized = TRUE)
#--------------------- 7. vst 转换（用于 PCA / heatmap）
    vst_mat <- vst(dds, blind = TRUE) %>% assay()
#--------------------- 8. 保存结果
    write.table(raw,
                file.path(outdir, "00.gene_counts.txt"),
                sep="\t", quote=FALSE)#featureCounts 原始结果
    write.table(count_mat,
                file.path(outdir, "01.counts_raw.txt"),
                sep="\t", quote=FALSE)# 过滤后的原始 counts 矩阵,可用于 DESeq2 / edgeR 差异分析的输入
    write.table(norm_counts,
                file.path(outdir, "02.counts_norm.txt"),
                sep="\t", quote=FALSE)# DESeq2 归一化后的 counts 矩阵，可用于画校正后表达量（画单基因图）
    write.table(vst_mat,
                file.path(outdir, "03.vst_matrix.txt"),
                sep="\t", quote=FALSE)# VST 转换后的矩阵，可用于 PCA / heatmap，样本相关性分析，大多数下游可视化
    cat("Output saved to ", outdir, "\n")
