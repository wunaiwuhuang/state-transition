#!/bin/bash
# 01.samtools.sh
# 对 STAR 生成的 BAM 文件进行 QC 统计

set -euo pipefail

# 输入目录（STAR 输出）
ALIGN_DIR="/disk2/cai113/data/stateTrans/02.alignment"

# 输出目录
OUT_DIR="/disk2/cai113/data/stateTrans/03.quantify"

# 创建输出目录
mkdir -p "$OUT_DIR"

# 遍历所有样本文件夹（排除非样本项如 STAR_ref, 00.*）
for sample_dir in "$ALIGN_DIR"/L1MKI*; do
    # 跳过非目录项（如文件）
    [ ! -d "$sample_dir" ] && continue

    sample_name=$(basename "$sample_dir")
    bam_file="$sample_dir/Aligned.sortedByCoord.out.bam"

    # 检查 BAM 是否存在
    if [ ! -f "$bam_file" ]; then
        echo "BAM not found for $sample_name, skipping..."
        continue
    fi

    echo "Processing: $sample_name"

    # 运行 samtools flagstat
    samtools flagstat "$bam_file" > "$OUT_DIR/${sample_name}.flagstat.txt"

    # 运行 samtools idxstats
    samtools idxstats "$bam_file" > "$OUT_DIR/${sample_name}.idxstats.txt"
done

echo "samtools QC 完成！结果保存在: $OUT_DIR"