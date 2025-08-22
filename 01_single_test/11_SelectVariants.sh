#!/bin/bash

# 检查输入参数的个数
if [ $# -ne 3 ]; then
    echo "Usage: $0 <input_vcf> <output_dir> <output_name>"
    exit 1
fi

# 读取输入参数
input_vcf=$1
output_dir=$2
output_name=$3

# 检查 GATK 是否安装
if ! command -v gatk &> /dev/null; then
    echo "Error: GATK is not installed."
    exit 1
fi

# 创建输出目录（如果不存在的话）
mkdir -p "$output_dir"

# 定义 SNP 和 INDEL 输出文件路径
output_snp_vcf="${output_dir}/${output_name}.snp.vcf.gz"
output_indel_vcf="${output_dir}/${output_name}.indel.vcf.gz"

# 使用 GATK SelectVariants 选出 SNP
echo "Selecting SNP variants..."
gatk SelectVariants \
    -select-type SNP \
    -V "$input_vcf" \
    -O "$output_snp_vcf"
if [ $? -ne 0 ]; then
    echo "Error: SelectVariants for SNPs failed."
    exit 1
fi
echo "** SNP selection done **"

# 使用 GATK SelectVariants 选出 INDEL
echo "Selecting INDEL variants..."
gatk SelectVariants \
    -select-type INDEL \
    -V "$input_vcf" \
    -O "$output_indel_vcf"
if [ $? -ne 0 ]; then
    echo "Error: SelectVariants for INDELs failed."
    exit 1
fi
echo "** INDEL selection done **"

