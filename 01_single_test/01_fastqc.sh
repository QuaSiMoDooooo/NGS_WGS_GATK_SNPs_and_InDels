#!/bin/bash

# 检查输入参数的个数
if [ $# -ne 3 ]; then
    echo "Usage: $0 <sample_list> <input_dir> <output_dir>"
    exit 1
fi

# 读取输入参数
sample_list=$1
input_dir=$2
output_dir=$3

# 检查 FastQC 是否安装
if ! command -v fastqc &> /dev/null; then
    echo "Error: FastQC is not installed."
    exit 1
fi

# 创建输出目录（如果不存在的话）
mkdir -p $output_dir

# 遍历 sample_list 中的每一行
while IFS= read -r sample_id; do
    # 生成双端测序文件名
    fq1="${input_dir}/${sample_id}_1.fastq.gz"
    fq2="${input_dir}/${sample_id}_2.fastq.gz"

    # 检查输入文件是否存在
    if [ ! -f "$fq1" ] || [ ! -f "$fq2" ]; then
        echo "Warning: One or both files for sample $sample_id are missing ($fq1 or $fq2)"
        continue
    fi

    # 使用 FastQC 对双端测序文件进行分析
    echo "Running FastQC for sample: $sample_id"
    fastqc -t 12 -o $output_dir $fq1
    fastqc -t 12 -o $output_dir $fq2
done < "$sample_list"

echo "FastQC analysis completed. Results saved in $output_dir."

