#!/bin/bash
# 02.featurecounts.sh
# 使用 featureCounts 对所有 STAR BAM 文件进行基因定量

set -euo pipefail

# --- 配置区 ---
ALIGN_DIR="/disk2/cai113/data/stateTrans/02.alignment"
OUT_DIR="/disk2/cai113/data/stateTrans/03.quantify"
GTF="/disk2/cai113/data/stateTrans/02.alignment/STAR_ref/mouse_GRCm39/Mus_musculus.GRCm39.110.gtf"
THREADS=12    # CPU 核心数调整
# -------------
if [ ! -f "$GTF" ]; then # 检查 GTF 是否存在
    echo "错误：GTF 文件不存在！"
    echo "请检查路径是否正确：$GTF"
    exit 1
fi
echo "使用的注释文件: $GTF"
# -------------

mkdir -p "$OUT_DIR"

# 收集所有样本 BAM 文件
bam_files=()
sample_names=()

for sample_dir in "$ALIGN_DIR"/L1MKI*; do
    [ ! -d "$sample_dir" ] && continue
    bam="$sample_dir/Aligned.sortedByCoord.out.bam"
    if [ -f "$bam" ]; then
        bam_files+=("$bam")
        sample_names+=("$(basename "$sample_dir")")
    fi
done

if [ ${#bam_files[@]} -eq 0 ]; then
    echo "未找到任何 BAM 文件！检查: $ALIGN_DIR"
    exit 1
fi

echo "共检测到 ${#bam_files[@]} 个样本："
printf '  %s\n' "${sample_names[@]}"

# 运行 featureCounts
featureCounts \
  -T "$THREADS" \
  -p -B -C \
  -a "$GTF" \
  -o "$OUT_DIR/gene_counts.txt" \
  "${bam_files[@]}"

echo ""
echo "featureCounts 运行完成！"
echo "主输出: $OUT_DIR/gene_counts.txt"
echo "摘要:   $OUT_DIR/gene_counts.txt.summary"