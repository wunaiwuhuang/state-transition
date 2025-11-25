#!/bin/bash
# 00.fastp.sh
# 功能：对原始RNA-seq数据进行质控，生成clean data
# 输入：rawdata目录下的 *.R1.raw.fastq.gz 和 *.R2.raw.fastq.gz
# 输出：01.cleandata目录下的 clean fastq + QC报告

set -euxo pipefail

# =============== 配置区域 ===============
RAW_DIR="/disk2/cai113/data/stateTrans/00.rawdata/MJ20250904414-ZX-R-250918200-董阁-真核转录组纯测序-32个样本/rawdata"
OUT_DIR="/disk2/cai113/data/stateTrans/01.cleandata"

# 检查输出目录是否存在，不存在则报错退出
if [ ! -d "$OUT_DIR" ]; then
    echo "错误：输出目录不存在！请先手动创建："
    echo "  mkdir -p \"$OUT_DIR\""
    exit 1
fi

# 接头序列（来自测序公司报告）
ADAPTER_R1="AGATCGGAAGAGCACACGTC"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGT"

# 线程数
THREADS=8
# ======================================

# 获取所有R1文件
R1_FILES=("$RAW_DIR"/*.R1.raw.fastq.gz)

# 检查是否有R1文件
if [ ${#R1_FILES[@]} -eq 0 ]; then
    echo "错误：在 $RAW_DIR 中未找到 *.R1.raw.fastq.gz 文件！"
    exit 1
fi

echo "共找到 ${#R1_FILES[@]} 个样本，开始质控..."

# 遍历每个样本
for r1 in "${R1_FILES[@]}"; do
    # 提取样本名（去掉路径和后缀）
    base_name=$(basename "$r1" .R1.raw.fastq.gz)
    
    # 对应的R2文件
    r2="${RAW_DIR}/${base_name}.R2.raw.fastq.gz"
    
    # 检查R2是否存在
    if [ ! -f "$r2" ]; then
        echo "警告：R2文件不存在: $r2，跳过样本 $base_name"
        continue
    fi
    
    echo "处理样本: $base_name"
    
    # 输出文件路径
    out_r1="${OUT_DIR}/${base_name}.R1.clean.fastq.gz"
    out_r2="${OUT_DIR}/${base_name}.R2.clean.fastq.gz"
    html_report="${OUT_DIR}/${base_name}_fastp.html"
    json_report="${OUT_DIR}/${base_name}_fastp.json"
    
    # 运行fastp
    fastp \
        -i "$r1" \
        -I "$r2" \
        -o "$out_r1" \
        -O "$out_r2" \
        --adapter_sequence "$ADAPTER_R1" \
        --adapter_sequence_r2 "$ADAPTER_R2" \
        --detect_adapter_for_pe \
        --trim_front1 0 \
        --trim_front2 0 \
        --cut_front \
        --cut_tail \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 50 \
        --length_required 30 \
        --n_base_limit 5 \
        --detect_adapter_for_pe \
        --thread "$THREADS" \
        --html "$html_report" \
        --json "$json_report" \
        --report_title "Fastp QC Report - $base_name"
    
    echo "完成样本: $base_name"
done

echo "=== Fastp质控完成！ ==="
echo "Clean数据位置: $OUT_DIR/*.clean.fastq.gz"
echo "QC报告位置: $OUT_DIR/*_fastp.html"

# =============== 新增：使用 multiqc 环境生成汇总报告 ===============
if command -v conda &> /dev/null && conda env list | grep -q 'multiqc'; then
    echo "检测到 multiqc 环境，正在生成汇总报告..."
    
    # 使用 conda run 在指定环境中执行命令（无需手动激活/退出）
    conda run -n multiqc multiqc "$OUT_DIR" -o "$OUT_DIR"/multiqc_report
    
    echo "✅ MultiQC汇总报告已保存至: $OUT_DIR/multiqc_report/"
    echo "   报告文件: $OUT_DIR/multiqc_report/multiqc_report.html"
else
    echo "⚠️ 未找到 multiqc 环境，跳过汇总报告生成"
    echo "   如需生成，请先运行："
    echo "   conda create -n multiqc -c conda-forge -c bioconda python=3.10 multiqc"
fi