#!/bin/bash

# 检查输入参数的个数
if [ $# -ne 3]; then
    echo "Usage: $0 <sample_list> <input_dir> <output_dir> <reference_genome>"
    exit 1
fi

# 读取输入参数
sample_list=$1
input_dir=$2
output_dir=$3

# 检查 gatk 是否安装
if ! command -v gatk &> /dev/null; then
    echo "Error: GATK is not installed."
    exit 1
fi

# 创建输出目录（如果不存在的话）
mkdir -p $output_dir

# 遍历 sample_list 中的每一行
while IFS= read -r sample_id; do
    # 使用 GATK MarkDuplicates 标记重复
    echo "Running GATK MarkDuplicates for sample: $sample_id"
    gatk MarkDuplicates \
        -I ${input_dir}/${sample_id}.bam \
        -O ${output_dir}/${sample_id}_markdup.bam \
        -M ${output_dir}/${sample_id}_metrics.csv \
	--CREATE_INDEX true \
	--REMOVE_DUPLICATES true

    # 检查 GATK 是否成功运行
    if [ $? -ne 0 ]; then
        echo "Error: GATK MarkDuplicates failed for sample $sample_id"
        continue
    fi

done < "$sample_list"

echo "BWA alignment and GATK MarkDuplicates completed. Results saved in $output_dir."

