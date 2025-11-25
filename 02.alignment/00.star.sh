#!/bin/bash
# STAR alignment for mouse RNA-seq (GRCm39)
# Input: clean FASTQ from 01.cleandata/
# Output: BAM + gene counts in 02.alignment/<sample>/

set -euo pipefail

# ================== 配置区 ==================
CLEAN_DIR="/disk2/cai113/data/stateTrans/01.cleandata"
ALIGN_DIR="/disk2/cai113/data/stateTrans/02.alignment"
INDEX_DIR="/disk2/cai113/data/stateTrans/02.alignment/STAR_ref/mouse_GRCm39/STAR_index"
CONDA_ENV="star"
THREADS=16
# ===========================================

echo "[$(date)] Starting STAR alignment..."
echo "Clean data dir: $CLEAN_DIR"
echo "Output dir:     $ALIGN_DIR"
echo "STAR index:     $INDEX_DIR"
echo "Conda env:      $CONDA_ENV"
echo "Threads:        $THREADS"
echo

# 初始化 conda（兼容 nohup）
__conda_setup="$('/disk2/cai113/software/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/disk2/cai113/software/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/disk2/cai113/software/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/disk2/cai113/software/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup

# 激活环境
conda activate "$CONDA_ENV"

# 验证 STAR 可用
echo "[$(date)] Using STAR version: $(STAR --version)"

# 获取所有 R1 文件
R1_FILES=("$CLEAN_DIR"/*.R1.clean.fastq.gz)

if [ ${#R1_FILES[@]} -eq 0 ]; then
    echo "错误：未在 $CLEAN_DIR 中找到 .R1.clean.fastq.gz 文件！"
    exit 1
fi

echo "共检测到 ${#R1_FILES[@]} 个样本"

for r1 in "${R1_FILES[@]}"; do
    # 提取样本名（如 L1MKI0600921-w4_1）
    sample_name=$(basename "$r1" .R1.clean.fastq.gz)
    r2="${CLEAN_DIR}/${sample_name}.R2.clean.fastq.gz"
    out_dir="${ALIGN_DIR}/${sample_name}"

    # 检查 R2 是否存在
    if [ ! -f "$r2" ]; then
        echo "警告：R2 文件缺失，跳过 $sample_name"
        continue
    fi

    # 如果已完成，跳过
    if [ -f "$out_dir/Aligned.sortedByCoord.out.bam" ]; then
        echo "已完成: $sample_name"
        continue
    fi

    echo "正在比对: $sample_name"
    mkdir -p "$out_dir"

    STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$INDEX_DIR" \
        --readFilesIn "$r1" "$r2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$out_dir/" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outSAMattributes All \
        --limitBAMsortRAM 31000000000 \
        > "$out_dir/star.log" 2>&1

    echo "完成: $sample_name"
done

echo
echo "所有样本比对完成！"
echo "BAM 文件: $ALIGN_DIR/*/Aligned.sortedByCoord.out.bam"
echo "基因计数: $ALIGN_DIR/*/ReadsPerGene.out.tab"
echo "结束时间: $(date)"