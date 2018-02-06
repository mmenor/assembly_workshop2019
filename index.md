* [Server Information](#server-information)
* [Dataset Information](#dataset-information)
* [Raw QC](#raw-quality-control-(qc))
* [Quality Filtering and Trimming](#quality-filtering-and-trimming)
* [Reference Genome Index](#reference-genome-index)
* [Align Reads](#align-reads)
* [Alignment Metrics](#alignment-metrics)
* [Determine Genotypes](#determine-genotypes)
* [Variant Annotation-(#variant-annotation)


# Server Information

# Dataset Information

# Linux Command-line

# Raw Quality Control (QC)

```bash
mkdir -p raw_qc
fastqc -o raw_qc reads/SRR097849_1.fastq reads/SRR097849_2.fastq
```
# Quality Filtering and Trimming

```bash
mkdir -p filtered_reads
cutadapt -m 20 -q 20,20 --pair-filter=any -o filtered_reads/trim_SRR097849_1.fastq -p filtered_reads/trim_SRR097849_2.fastq reads/SRR097849_1.fastq reads/SRR097849_2.fastq
fastqc -o raw_qc filtered_reads/trim_SRR097849_1.fastq filtered_reads/trim_SRR097849_2.fastq
```

# Reference Genome Index

```bash
bowtie2-build reference/hg19chr21.fa reference/hg19chr21.fa
```

# Align Reads

```bash
mkdir -p alignment
bowtie2 -p 4 -x reference/hg19chr21.fa -1 filtered_reads/trim_SRR097849_1.fastq -2 filtered_reads/trim_SRR097849_2.fastq -S alignment/SRR097849.sam
#samtools view -bS alignment/SRR097849.sam | samtools sort - alignment/SRR097849
samtools view -bS alignment/SRR097849.sam -o alignment/temp.bam
samtools sort alignment/temp.bam alignment/SRR097849
samtools index alignment/SRR097849.bam
```

# Alignment Metrics

```bash
picard CollectAlignmentSummaryMetrics INPUT=alignment/SRR097849.bam OUTPUT=alignment/SRR097849_metrics.txt REFERENCE_SEQUENCE=reference/hg19chr21.fa
```

# Determine Genotypes

```bash
mkdir -p variants
freebayes -f reference/hg19chr21.fa alignment/SRR097849.bam > variants/SRR097849.vcf
bcftools norm -f reference/hg19chr21.fa variants/SRR097849.vcf > variants/SRR097849_norm.vcf
```

# Variant Annotation

```bash
snpEff -Xmx50g ann hg19 -s variants/stats.html variants/SRR097849_norm.vcf > variants/SRR097849_ann.vcf
bgzip -c variants/SRR097849_ann.vcf > variants/SRR097849_ann.vcf.gz
tabix -p vcf variants/SRR097849_ann.vcf.gz
```
