#!/bin/bash

# 检查输入参数的个数
if [ $# -ne 4 ]; then
    echo "Usage: $0 <sample_list> <input_dir> <output_dir> <reference_genome>"
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

# 遍历 sample_list 中的每一行
while IFS= read -r sample_id; do
    # 获取对应的 BAM 文件
    input_bam="${input_dir}/${sample_id}.markdup.BQSR.bam"
    
    # 检查 BAM 文件是否存在
    if [ ! -f "$input_bam" ]; then
        echo "Warning: BAM file for sample $sample_id is missing ($input_bam)"
        continue
    fi

    # 运行 GATK HaplotypeCaller
    echo "Running GATK HaplotypeCaller for sample: $sample_id"
    gatk HaplotypeCaller \
        -R $reference_genome \
        --emit-ref-confidence GVCF \
        -I $input_bam \
        -O ${output_dir}/${sample_id}.gvcf

    # 检查 GATK 是否成功运行
    if [ $? -ne 0 ]; then
        echo "Error: GATK HaplotypeCaller failed for sample $sample_id"
        continue
    fi

done < "$sample_list"

echo "HaplotypeCaller GVCF generation completed for all samples. Results saved in $output_dir."

