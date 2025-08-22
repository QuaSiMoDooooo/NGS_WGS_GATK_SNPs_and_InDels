## 项目结构

```
NGS_WGS_GATK_SNPs_and_InDels/
├── 01_single_test/          # 逐步脚本版本
│   ├── 01_fastqc.sh        # 原始数据质控
│   ├── 02_trim_galore.sh   # 数据修剪
│   ├── 03_fastqc.sh        # 修剪后质控
│   ├── 04_bwa.sh          # BWA比对
│   ├── 05_bwa_MarkDuplicates.sh  # 去重
│   ├── 06_BQSR.sh         # 基础质量重校准
│   ├── 07_haplotypecaller_gvcf.sh  # 变异检测
│   ├── 08_CombineGVCFs.sh # GVCF合并
│   ├── 09_GenotypeGVCFs.sh # 联合基因分型
│   ├── 10_VQSR.sh         # 变异质量校准
│   ├── 11_SelectVariants.sh # SNP/INDEL分离
│   └── sample_list.txt    # 样本列表
└── 02_batch_pipeline/     # Snakemake自动化流程
    ├── Snakefile          # 主流程文件
    ├── config.yaml        # 配置文件
    └── sample_list.txt    # 样本列表
```

## 功能特点

- **完整GATK最佳实践**：严格遵循GATK官方推荐流程
- **双模式支持**：提供逐步脚本和自动化流程两种选择
- **高质量结果**：包含质控、校准、过滤等完整步骤
- **灵活配置**：支持自定义参考基因组和参数调整
- **并行处理**：支持多线程加速计算

## 依赖要求

### 必需软件
- **BWA** : 序列比对工具
- **GATK** : 基因组分析工具包
- **Samtools** : SAM/BAM文件处理
- **FastQC** : 质量控制工具
- **Trim Galore** : 数据修剪工具
- **Snakemake** : 工作流管理（仅自动化流程需要）

### 参考数据
- **GRCh38参考基因组**：`GRCh38.NCBI.no_alt.no_contigs.toUCSC.fasta`
- **GATK资源包**：包含以下文件
  - `hg38_v0_hapmap_3.3.hg38.vcf`
  - `hg38_v0_1000G_omni2.5.hg38.vcf`
  - `hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf`
  - `hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf`
  - `hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf`

## 使用方法

### 方法一：逐步脚本（适合学习和调试）

1. **准备数据**
   ```bash
   # 编辑样本列表
   echo "sample1" > sample_list.txt
   echo "sample2" >> sample_list.txt
   ```

2. **运行流程**
   ```bash
   cd 01_single_test/
   
   # 1. 原始数据质控
   bash 01_fastqc.sh sample_list.txt /path/to/fastq ./01_fastqc_result
   
   # 2. 数据修剪
   bash 02_trim_galore.sh sample_list.txt /path/to/fastq ./02_trim_galore_result
   
   # 3. 修剪后质控
   bash 03_fastqc.sh sample_list.txt ./02_trim_galore_result ./03_fastqc_result
   
   # 4. BWA比对
   bash 04_bwa.sh sample_list.txt ./02_trim_galore_result ./04_bwa_result /path/to/reference
   
   # 5. 去重
   bash 05_bwa_MarkDuplicates.sh sample_list.txt ./04_bwa_result ./05_bwa_MarkDuplicates_result
   
   # 6. 基础质量重校准
   bash 06_BQSR.sh sample_list.txt ./05_bwa_MarkDuplicates_result ./06_BQSR_result /path/to/reference /path/to/gatk_bundle
   
   # 7. 变异检测
   bash 07_haplotypecaller_gvcf.sh sample_list.txt ./06_BQSR_result ./07_haplotypecaller_gvcf_result /path/to/reference
   
   # 8. GVCF合并
   bash 08_CombineGVCFs.sh sample_list.txt ./07_haplotypecaller_gvcf_result ./08_CombineGVCFs_result /path/to/reference
   
   # 9. 联合基因分型
   bash 09_GenotypeGVCFs.sh ./08_CombineGVCFs_result/combined.g.vcf.gz ./09_GenotypeGVCFs_result /path/to/reference
   
   # 10. 变异质量校准
   bash 10_VQSR.sh ./09_GenotypeGVCFs_result/all_samples.vcf ./10_VQSR_result /path/to/reference /path/to/gatk_bundle all_sample
   
   # 11. SNP/INDEL分离
   bash 11_SelectVariants.sh ./10_VQSR_result/all_sample.HC.VQSR.vcf.gz ./11_SelectVariants_result /path/to/reference
   ```

### 方法二：Snakemake自动化流程（推荐）

1. **配置参数**
   ```bash
   cd 02_batch_pipeline/
   
   # 编辑配置文件
   vim config.yaml
   ```

2. **修改配置文件**
   ```yaml
   # 路径配置
   input_dir: "/path/to/your/fastq/files"
   sample_list: "sample_list.txt"
   output_dir: "results"
   reference: "/path/to/GRCh38.fasta"
   gatk_bundle: "/path/to/gatk/resource"
   
   # 线程数
   threads: 10
   ```

3. **运行流程**
   ```bash
   # 试运行（检查流程）
   snakemake -n
   
   # 多线程运行
   snakemake -j 10
   ```

## 输出文件说明

### 主要结果文件
- `*.snp.vcf.gz`: 最终的高质量SNP文件
- `*.indel.vcf.gz`: 最终的高质量INDEL文件
- `*.g.vcf.gz`: 各样本的基因组VCF文件

### 质控文件
- `*_fastqc.html`: FastQC质控报告
- `*_trimming_report.txt`: 数据修剪报告
- `*.stats`: 比对统计信息
- `*_metrics.txt`: 重复序列统计

### 校准文件
- `*.table`: BQSR校准表
- `*.recal`: VQSR校准文件
- `*.tranches`: VQSR质量分层文件
- `*.plots.R.pdf`: VQSR可视化报告

## 注意事项

1. **数据准备**：确保FASTQ文件命名格式为`{sample}_1.fastq.gz`和`{sample}_2.fastq.gz`


## DEBUG日志

### 常见问题
1. **内存不足**：减少线程数或增加系统内存
2. **文件权限**：确保所有目录有读写权限
3. **依赖缺失**：使用conda环境管理依赖
4. **路径错误**：检查所有路径是否为绝对路径

### 日志查看
```bash
# 查看Snakemake日志
tail -f .snakemake/log/*.log

# 查看具体步骤日志
cat results/benchmarks/*.txt
```

        