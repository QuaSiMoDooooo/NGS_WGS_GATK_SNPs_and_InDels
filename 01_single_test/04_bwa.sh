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

# 检查 bwa 和 samtools 是否安装
if ! command -v bwa &> /dev/null; then
    echo "Error: bwa is not installed."
    exit 1
fi

if ! command -v samtools &> /dev/null; then
    echo "Error: samtools is not installed."
    exit 1
fi

index_files=(
    "${reference_genome}.amb"
    "${reference_genome}.ann"
    "${reference_genome}.bwt"
    "${reference_genome}.pac"
    "${reference_genome}.sa"
)
# 检查索引文件是否存在
missing_index=0
for index in "${index_files[@]}"; do
    if [ ! -f "$index" ]; then
        echo "Index file $index not found. Building bwa index..."
        missing_index=1
    fi
done
# 如果缺少索引文件，则使用 bwa index 创建索引
if [ $missing_index -eq 1 ]; then
    echo "Building bwa index for reference genome: $reference_genome"
    bwa index $reference_genome
    if [ $? -ne 0 ]; then
        echo "Error: Failed to build bwa index for reference genome $reference_genome"
        exit 1
    fi
fi

# 创建输出目录（如果不存在的话）
mkdir -p $output_dir

# 遍历 sample_list 中的每一行
while IFS= read -r sample_id; do
    # 生成双端测序文件名
    fq1="${input_dir}/${sample_id}_1_val_1.fq.gz"
    fq2="${input_dir}/${sample_id}_2_val_2.fq.gz"

    # 检查输入文件是否存在
    if [ ! -f "$fq1" ] || [ ! -f "$fq2" ]; then
        echo "Warning: One or both files for sample $sample_id are missing ($fq1 or $fq2)"
        continue
    fi

    read_group="@RG\tID:${sample_id}_id\tPL:illumina\tSM:${sample_id}"
    # 运行 bwa mem 进行比对
    echo "Running bwa mem for sample: $sample_id"
    bwa mem -t 30 -R "$read_group" $reference_genome $fq1 $fq2 | samtools sort -@ 30 -o ${output_dir}/${sample_id}.bam
    samtools flagstat -@ 30 ${output_dir}/${sample_id}.bam > ${output_dir}/${sample_id}.stats

done < "$sample_list"

echo "BWA mem alignment completed. Results saved in $output_dir."

