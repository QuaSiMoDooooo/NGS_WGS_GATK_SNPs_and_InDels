#!/bin/bash

# 检查输入参数的个数
if [ $# -ne 5 ]; then
    echo "Usage: $0 <input_vcf> <output_dir> <reference_genome> <gatk_bundle_dir> <output_name>"
    exit 1
fi

# 读取输入参数
input_vcf=$1
output_dir=$2
reference_genome=$3
gatk_bundle_dir=$4
output_name=$5

# 检查 GATK 是否安装
if ! command -v gatk &> /dev/null; then
    echo "Error: GATK is not installed."
    exit 1
fi

# 创建输出目录（如果不存在的话）
mkdir -p "$output_dir"

# 定义输出路径
output_vcf_snps="${output_dir}/${output_name}.HC.snps.VQSR.vcf.gz"
output_vcf_final="${output_dir}/${output_name}.HC.VQSR.vcf.gz"

# SNP 模式：VariantRecalibrator 和 ApplyVQSR
echo "Running VQSR for SNPs..."
gatk VariantRecalibrator \
    -R "$reference_genome" \
    -V "$input_vcf" \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 "${gatk_bundle_dir}/hg38_v0_hapmap_3.3.hg38.vcf" \
    -resource:omini,known=false,training=true,truth=false,prior=12.0 "${gatk_bundle_dir}/hg38_v0_1000G_omni2.5.hg38.vcf" \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 "${gatk_bundle_dir}/hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf" \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "${gatk_bundle_dir}/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf" \
    -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
    -mode SNP \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
    --rscript-file "${output_dir}/${output_name}.HC.snps.plots.R" \
    --tranches-file "${output_dir}/${output_name}.HC.snps.tranches" \
    -O "${output_dir}/${output_name}.HC.snps.recal"

gatk ApplyVQSR \
    -R "$reference_genome" \
    -V "$input_vcf" \
    --ts-filter-level 99.0 \
    --tranches-file "${output_dir}/${output_name}.HC.snps.tranches" \
    --recal-file "${output_dir}/${output_name}.HC.snps.recal" \
    -mode SNP \
    -O "$output_vcf_snps"
echo "** SNPs VQSR done **"

# INDEL 模式：VariantRecalibrator 和 ApplyVQSR
echo "Running VQSR for INDELs..."
gatk VariantRecalibrator \
    -R "$reference_genome" \
    -V "$output_vcf_snps" \
    -resource:mills,known=true,training=true,truth=true,prior=12.0 "${gatk_bundle_dir}/hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf" \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "${gatk_bundle_dir}/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf" \
    -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
    -mode INDEL \
    --max-gaussians 6 \
    --rscript-file "${output_dir}/${output_name}.HC.snps.indels.plots.R" \
    --tranches-file "${output_dir}/${output_name}.HC.snps.indels.tranches" \
    -O "${output_dir}/${output_name}.HC.snps.indels.recal"

gatk ApplyVQSR \
    -R "$reference_genome" \
    -V "$output_vcf_snps" \
    --ts-filter-level 99.0 \
    --tranches-file "${output_dir}/${output_name}.HC.snps.indels.tranches" \
    --recal-file "${output_dir}/${output_name}.HC.snps.indels.recal" \
    -mode INDEL \
    -O "$output_vcf_final"
echo "** SNPs and Indels VQSR (${output_name}.HC.VQSR.vcf.gz) done **"

