#!/bin/bash

# 检查输入参数的个数
if [ $# -ne 3 ]; then
    echo "Usage: $0 <combined_gvcf_file> <output_dir> <reference_genome>"
    exit 1
fi

# 读取输入参数
combined_gvcf=$1
output_dir=$2
reference_genome=$3

# 检查 GATK 是否安装
if ! command -v gatk &> /dev/null; then
    echo "Error: GATK is not installed."
    exit 1
fi

# 检查参考基因组文件是否存在
if [ ! -f "$reference_genome" ]; then
    echo "Error: Reference genome not found at $reference_genome"
    exit 1
fi

# 检查合并的 GVCF 文件是否存在
if [ ! -f "$combined_gvcf" ]; then
    echo "Error: Combined GVCF file not found at $combined_gvcf"
    exit 1
fi

# 创建输出目录（如果不存在的话）
mkdir -p "$output_dir"

# 生成输出 VCF 文件路径
output_vcf="${output_dir}/all_samples.vcf"

# 运行 GATK GenotypeGVCFs
echo "Running GATK GenotypeGVCFs..."
gatk GenotypeGVCFs \
    -R "$reference_genome" \
    -V "$combined_gvcf" \
    -O "$output_vcf"

# 检查 GATK 是否成功运行
if [ $? -ne 0 ]; then
    echo "Error: GATK GenotypeGVCFs failed."
    exit 1
fi

echo "GenotypeGVCFs completed successfully. Output VCF saved at $output_vcf."

