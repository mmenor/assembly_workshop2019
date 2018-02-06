* [Server Information](#server-information)
* [Dataset Information](#dataset-information)
* [Linux Command Line](#linux-command-line)
* [Raw QC](#raw-quality-control-(qc))
* [Quality Filtering and Trimming](#quality-filtering-and-trimming)
* [Reference Genome Index](#reference-genome-index)
* [Align Reads](#align-reads)
* [Alignment Metrics](#alignment-metrics)
* [Determine Genotypes](#determine-genotypes)
* [Variant Annotation](#variant-annotation)
* [Analysis in R](#analysis-in-r)


# Server Information

# Dataset Information

# Linux Command Line

A complete tutorial on the Linux command line is a full workshop on its own. Please refer to other resources, e.g. [Ryan's Tutorials](https://ryanstutorials.net/linuxtutorial/). For our purposes, we just need to know some basic navigation.

You are signed on as a guest on the tutorial server, so when you first log in, your current working directory is `/home/guest`. Everyone in the workshop should work in a different directory, so let's make a folder named after yourself.

```bash
mkdir write_your_UH_username_here
```

For example, `mkdir mmenor` if your UH username is _mmenor_. Let us confirm you created a directory by listing the contents of your working directory,

```bash
ls
```

You should see a folder with your UH username, among the other workshop participants. Let's now change our working directory to your newly created folder using the change directory command,

```bash
cd write_your_UH_username_here
```

For example, `cd mmenor` in my case. You can run the `ls` command again to confirm you've changed folders and your current working directory is empty. Now we're ready to start the data analysis.

# Raw Quality Control (QC)

Garbage in is garbage out, so let's see check if we have quality issues with our raw reads. We'll be using [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for this purpose.

```bash
mkdir raw_qc
fastqc -o raw_qc /home/bqhs/reads/SRR097849_1.fastq /home/bqhs/reads/SRR097849_2.fastq
```
The first command creates a new folder that will store our results `raw_qc`. The `fastqc` command takes in three arguments. `-o raw_qc` specifies where the output will be saved, in this case our newly created `raw_qc` folder. Since our reads our paired-end, the next two arguments specify the FASTQ files of the pair.

The output of FastQC is an HTML page with plots that unfortunately we cannot view with an SSH remote terminal. Let's download the results using (FileZilla)[https://filezilla-project.org/].

The accuracy of a base's identifiy is measured using [Phred quality scores](https://en.wikipedia.org/wiki/Phred_quality_score). The higher the score, the better. Roughly, a Phred score of 20 or above (99+% accuracy) is great.

# Quality Filtering and Trimming

As you saw, base quality degrades at the ends of the read, particularly at the 5' end. Let's trim low quality bases of the ends our reads using [cutadapt](https://cutadapt.readthedocs.io/en/stable/). Cutadapt can also trim off adapters, but as we didn't identify any with FastQC, it is not required today.

```bash
mkdir filtered_reads
cutadapt -m 20 -q 20,20 --pair-filter=any -o filtered_reads/trim_SRR097849_1.fastq -p filtered_reads/trim_SRR097849_2.fastq /home/bqhs/reads/SRR097849_1.fastq /home/bqhs/reads/SRR097849_2.fastq
fastqc -o raw_qc filtered_reads/trim_SRR097849_1.fastq filtered_reads/trim_SRR097849_2.fastq
```

As with the previously step, we first create a new folder to store our filtered and trimmed reads, `filtered_reads`. Next we call `cutadapt` with the following arguments. `-m 20` tells cutadapt to discard reads, after trimming, shorter than length 20 bases. `-q 20,20` specifies the Phred qualty score cutoffs for the low-quality trimming of the 5' and 3' ends of the read. In this case, we trim a base off the 5' or 3' end if its Phred score is below 20.

Since we are working with paired-end reads, we need to specify the pair filtering mode, `--pair-filter=any`. I selected the `any` mode that will discard the pair if one of the ends meets the filtering criteria. In this case, the only filter we specified was `-m 20`, so if any end of the pair is less the 20 bases long, the whole pair will be discarded.

The next two arguments, `-o` and `-p` specify the names of the cutadapt output FASTQ files, the filtered and trimmed reads for the first and second ends of the paired-end reads, respectively. The final two arguments specify the input raw FASTQ files.

The final `fastqc` command is optional in practice and is done here for this workshop to illustrate that cutadapt did in fact perform some trimming.

# Reference Genome Index

The next command is for your future reference. We've already built an index for the chromosome 21 of the human genome for this workshop, so you don't need to run it yourself. If you do attempt the command, it will fail as you don't have permission to alter this file.

The next step will be aligning the filtered and trimmed reads to a reference human genome (hg19). Aligners for short DNA reads require an index of the reference genome in order to align efficiently. We've chosen [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) as our aligner. Therefore, we need an index compatible with Bowtie2 and that can be created using the following command,

```bash
bowtie2-build /home/bqhs/reference/hg19chr21.fa /home/bqhs/reference/hg19chr21.fa
```

While the arguments look identical, they have different purposes. The first argument specifies a FASTA file of the sequence(s) we want to index. In this case contains only chromosome 21 of the human genome. The second argument specifies the name of the index, which I chose to name the same as the input.

# Align Reads

There are many progams to pick for aligning DNA reads to a reference. We've chosen [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) as our aligner. The output of Bowtie2 aligning the reads (FASTQ files) is a single file containing both the aligned and unaligned reads called a SAM file. Typically a compressed version of the SAM file, called a BAM file, is used in downstream analysis. We will use [Samtools](http://www.htslib.org/) for this conversion. Samtools also lets us sort the read by alignment position and create an index for the BAM file, both of which are usually required in downstream analysis.

```bash
mkdir alignment
bowtie2 -x /home/bqhs/reference/hg19chr21.fa -1 filtered_reads/trim_SRR097849_1.fastq -2 filtered_reads/trim_SRR097849_2.fastq -S alignment/SRR097849.sam
samtools view -bS alignment/SRR097849.sam -o alignment/temp.bam
samtools sort alignment/temp.bam alignment/SRR097849
samtools index alignment/SRR097849.bam
```

We again create a new folder, `alignment`, to save our results (BAM and SAM files). The next command, `bowtie2`, executes the alignment using Bowtie2. The first argument, `-x`, specifies the reference genome. In this case, it is chromosome 21 of the human genome. `-1` and `-2` specifies the input paired-end reads. We're using our filtered and trimmed reads here. Finally `-S` specifies the name of our output SAM file.

As mentioned previously, it is usually required in downstread analysis to convert the SAM file to a BAM file, sort it, and index it. The first command, `samtools view`, is used to convert the SAM file to a BAM file. The `-bS` flag specifies that we want our output to be a BAM file (`b`) and the input is a SAM file (`S`). The next argument is the name of our input file, which is the SAM file that Bowtie2 created for us. Lastly the `-o` argument specifies the name of the output BAM file.

To sort the read in the BAM file by alignment position, we use the `samtools sort` command. The first argument is the input BAM file and the second argument is the name for the output, sorted BAM file. Lastly, we created a index for our BAM file using the `samtools index` command that takes in a single argument, the name of the input BAM file.

# Alignment Metrics

```bash
picard CollectAlignmentSummaryMetrics INPUT=alignment/SRR097849.bam OUTPUT=alignment/SRR097849_metrics.txt REFERENCE_SEQUENCE=reference/hg19chr21.fa
```

# Determine Genotypes

```bash
mkdir variants
freebayes -f reference/hg19chr21.fa alignment/SRR097849.bam > variants/SRR097849.vcf
bcftools norm -f reference/hg19chr21.fa variants/SRR097849.vcf > variants/SRR097849_norm.vcf
```

# Variant Annotation

```bash
snpEff -Xmx1g ann hg19 -s variants/stats.html variants/SRR097849_norm.vcf > variants/SRR097849_ann.vcf
bgzip -c variants/SRR097849_ann.vcf > variants/SRR097849_ann.vcf.gz
tabix -p vcf variants/SRR097849_ann.vcf.gz
```

# Analysis in R
