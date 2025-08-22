#!/bin/bash

# 检查输入参数的个数
if [ $# -ne 4 ]; then
    echo "Usage: $0 <sample_list> <input_dir> <output_file> <reference_genome>"
    exit 1
fi

# 读取输入参数
sample_list=$1
input_dir=$2
output_dir=$3
reference_genome=$4

# 检查 GATK 是否安装
if ! command -v gatk &> /dev/null; then
    echo "Error: GATK is not installed."
    exit 1
fi

# 创建输出目录（如果不存在的话）
mkdir -p $output_dir

# 检查参考基因组文件是否存在
if [ ! -f "$reference_genome" ]; then
    echo "Error: Reference genome not found at $reference_genome"
    exit 1
fi

# 构建样本 GVCF 文件的参数列表
gvcf_files=()
while IFS= read -r sample_id; do
    gvcf_file="${input_dir}/${sample_id}.gvcf"
    
    # 检查 GVCF 文件是否存在
    if [ ! -f "$gvcf_file" ]; then
        echo "Warning: GVCF file for sample $sample_id is missing ($gvcf_file)"
        continue
    fi
    
    # 将 GVCF 文件添加到参数列表中
    gvcf_files+=("-V" "$gvcf_file")
done < "$sample_list"

# 检查是否有 GVCF 文件
if [ ${#gvcf_files[@]} -eq 0 ]; then
    echo "Error: No GVCF files found to combine."
    exit 1
fi

# 运行 GATK CombineGVCFs
echo "Running GATK CombineGVCFs..."
gatk CombineGVCFs -R $reference_genome "${gvcf_files[@]}" -O ${output_dir}/combined.g.vcf

# 检查 GATK 是否成功运行
if [ $? -ne 0 ]; then
    echo "Error: GATK CombineGVCFs failed."
    exit 1
fi

echo "CombineGVCFs completed successfully. Combined GVCF saved at $output_dir."

