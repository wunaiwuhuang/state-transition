#!/bin/bash
# Build STAR index for mouse GRCm39 genome
# check fastq files to determine the appropriate --sjdbOverhang value
# sjdbOverhang = read length - 1
# I have checked all samples: 'sequencing': 'paired end (151 cycles + 151 cycles)', 'before_filtering': {'total_reads': 41258850, 'total_bases': 6230086350, 'q20_bases': 6166628905, 'q30_bases': 5989245780, 'q20_rate': 0.989814, 'q30_rate': 0.961342, 'read1_mean_length': 151, 'read2_mean_length': 151, 'gc_content': 0.490562}, 'after_filtering': {'total_reads': 41106358, 'total_bases': 6182147744, 'q20_bases': 6133344466, 'q30_bases': 5966255282, 'q20_rate': 0.992106, 'q30_rate': 0.965078, 'read1_mean_length': 150, 'read2_mean_length': 150, 'gc_content': 0.490123}}	
# sjdbOverhang = 150 - 1 = 149 is appropriate for all samples
# Uses conda environment "star"

set -euo pipefail
echo "[$(date)] Starting STAR index build..."

# Change to the STAR index directory
cd /disk2/cai113/data/stateTrans/02.alignment/STAR_ref/mouse_GRCm39/STAR_index

# Initialize conda in non-interactive shell (critical for nohup)
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

# Activate the 'star' environment

conda activate star
echo "[$(date)] Using STAR version: $(STAR --version)"

# Build index
STAR \
    --runMode genomeGenerate \
    --genomeDir ./ \
    --genomeFastaFiles ../Mus_musculus.GRCm39.dna.primary_assembly.fa \
    --sjdbGTFfile ../Mus_musculus.GRCm39.110.gtf \
    --sjdbOverhang 149 \
    --runThreadN 16

echo "[$(date)] ✅ STAR index build completed successfully!"