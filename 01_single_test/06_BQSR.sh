#!/bin/bash

# 检查输入参数的个数
if [ $# -ne 5 ]; then
    echo "Usage: $0 <sample_list> <input_dir> <output_dir> <reference_genome> <known_sites_dir>"
    exit 1
fi

# 读取输入参数
sample_list=$1
input_dir=$2
output_dir=$3
reference_genome=$4
known_sites_dir=$5

# 检查 GATK 是否安装
if ! command -v gatk &> /dev/null; then
    echo "Error: GATK is not installed."
    exit 1
fi

# 检查输入文件路径和目录
if [ ! -d "$known_sites_dir" ]; then
    echo "Error: Known sites directory does not exist: $known_sites_dir"
    exit 1
fi

# 创建输出目录（如果不存在的话）
mkdir -p $output_dir

# 确定已知变异位点文件的路径
dbsnp="${known_sites_dir}/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
mills_and_1000G="${known_sites_dir}/hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf"
phase1_1000G="${known_sites_dir}/hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf"

# 检查已知变异位点文件是否存在
if [ ! -f "$dbsnp" ] || [ ! -f "$mills_and_1000G" ] || [ ! -f "$phase1_1000G" ]; then
    echo "Error: One or more known sites VCF files are missing."
    exit 1
fi

# 遍历 sample_list 中的每一行
while IFS= read -r sample_id; do
    # 生成标记重复的 BAM 文件路径
    markdup_bam="${input_dir}/${sample_id}_markdup.bam"

    # 检查 BAM 文件是否存在
    if [ ! -f "$markdup_bam" ]; then
        echo "Warning: BAM file for sample $sample_id is missing ($markdup_bam)"
        continue
    fi

    # 运行 GATK BaseRecalibrator
    echo "Running GATK BaseRecalibrator for sample: $sample_id"
    gatk BaseRecalibrator \
        -R $reference_genome \
        -I $markdup_bam \
        --known-sites $dbsnp \
        --known-sites $mills_and_1000G \
        --known-sites $phase1_1000G \
        -O ${output_dir}/${sample_id}.table

    # 检查 GATK 是否成功运行
    if [ $? -ne 0 ]; then
        echo "Error: GATK BaseRecalibrator failed for sample $sample_id"
        continue
    fi

    # 运行 GATK ApplyBQSR
    echo "Running GATK ApplyBQSR for sample: $sample_id"
    gatk ApplyBQSR \
        -R $reference_genome \
        -I ${input_dir}/${sample_id}_markdup.bam \
        --bqsr-recal-file ${output_dir}/${sample_id}.table \
        -O ${output_dir}/${sample_id}.markdup.BQSR.bam

    # 检查 ApplyBQSR 是否成功运行
    if [ $? -ne 0 ]; then
        echo "Error: GATK ApplyBQSR failed for sample $sample_id"
        continue
    fi

done < "$sample_list"

echo "BQSR and ApplyBQSR completed for all samples. Results saved in $output_dir."
