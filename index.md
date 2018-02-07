* [Server Information](#server-information)
* [Dataset Information](#dataset-information)
* [Linux Command Line](#linux-command-line)
* [Raw QC](#raw-quality-control)
* [Quality Filtering and Trimming](#quality-filtering-and-trimming)
* [Reference Genome Index](#reference-genome-index)
* [Align Reads](#align-reads)
* [Alignment Metrics](#alignment-metrics)
* [Variant Detection](#variant-detection)
* [Variant Annotation](#variant-annotation)

# Server Information

We will be using a Linux server (Ubuntu) for this workshop. You will need to use [PuTTY](https://www.chiark.greenend.org.uk/~sgtatham/putty/) on Windows or `ssh` on the Terminal in macOS.

# Dataset Information

[MCF7 breast cancer cell line](https://www.ncbi.nlm.nih.gov/sra/SRX040531). On the server we have a 10% subsample of the data for tutorial purposes saved at `/home/bqhs/reads/SRR097849_1.fastq` and `/home/bqhs/reads/SRR097849_2.fastq`. This is a paired-end experiment and thus we have two FASTQ files, one per end. FASTQ is the file format for raw reads and it contains information on the sequence and the quality (accurracy) of each base call.

Alignment of raw reads to a full reference genome is resource consuming, so for this tutorial we will only align to chromosome 21 of the human genome (hg19). There reference can be viewed here, `/home/bqhs/reference/hg19chr21.fa`.

At the end of this tutorial, we would like to see which variants we've detected in the MCF7 dataset correspond to known variants in the dbSnp database. The database can be viewed here, `/home/bqhs/dbsnp.hg19.vcf`.

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

# Raw Quality Control

Garbage in is garbage out, so let's see check if we have quality issues with our raw reads. We'll be using [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for this purpose.

```bash
mkdir raw_qc
fastqc -o raw_qc /home/bqhs/reads/SRR097849_1.fastq /home/bqhs/reads/SRR097849_2.fastq
```
The first command creates a new folder that will store our results `raw_qc`. The `fastqc` command takes in three arguments. `-o raw_qc` specifies where the output will be saved, in this case our newly created `raw_qc` folder. Since our reads our paired-end, the next two arguments specify the FASTQ files of the pair.

The output of FastQC is an HTML page with plots that unfortunately we cannot view with an SSH remote terminal. Let's download the results using [FileZilla](https://filezilla-project.org/). Connect using the same guest credentials you used on SSH. Navigate to `/home/guest/write_your_UH_username_here/raw_qc` to find your FASTQC results.

The accuracy of a base's identifiy is measured using [Phred quality scores](https://en.wikipedia.org/wiki/Phred_quality_score). The higher the score, the better. Roughly, a Phred score of 20 or above (99+% accuracy) is great.

# Quality Filtering and Trimming

As you saw, base quality degrades at the ends of the read, particularly at the 5' end. Let's trim low quality bases of the ends our reads using [cutadapt](https://cutadapt.readthedocs.io/en/stable/). Cutadapt can also trim off adapters, but as we didn't identify any with FastQC, it is not required today.

```bash
mkdir filtered_reads
cutadapt -m 20 -q 20,20 --pair-filter=any -o filtered_reads/trim_SRR097849_1.fastq -p filtered_reads/trim_SRR097849_2.fastq /home/bqhs/reads/SRR097849_1.fastq /home/bqhs/reads/SRR097849_2.fastq
```

As with the previously step, we first create a new folder to store our filtered and trimmed reads, `filtered_reads`. Next we call `cutadapt` with the following arguments. `-m 20` tells cutadapt to discard reads, after trimming, shorter than length 20 bases. `-q 20,20` specifies the Phred qualty score cutoffs for the low-quality trimming of the 5' and 3' ends of the read. In this case, we trim a base off the 5' or 3' end if its Phred score is below 20.

Since we are working with paired-end reads, we need to specify the pair filtering mode, `--pair-filter=any`. I selected the `any` mode that will discard the pair if one of the ends meets the filtering criteria. In this case, the only filter we specified was `-m 20`, so if any end of the pair is less the 20 bases long, the whole pair will be discarded.

The next two arguments, `-o` and `-p` specify the names of the cutadapt output FASTQ files, the filtered and trimmed reads for the first and second ends of the paired-end reads, respectively. The final two arguments specify the input raw FASTQ files.

To illustrate that cutadapt did in fact perform some trimming, we can run FastQC on our new data. In practice this is optional.

```bash
fastqc -o raw_qc filtered_reads/trim_SRR097849_1.fastq filtered_reads/trim_SRR097849_2.fastq
```

# Reference Genome Index

The next command is for your future reference. We've already built an index for the chromosome 21 of the human genome for this workshop, so you don't need to run it yourself. If you do attempt the command, it will fail as you don't have permission to alter this file.

The next step will be aligning the filtered and trimmed reads to a reference human genome (hg19). Aligners for short DNA reads require an index of the reference genome in order to align efficiently. We've chosen [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) as our aligner. Therefore, we need an index compatible with Bowtie2 and that can be created using the following command,

```bash
bowtie2-build /home/bqhs/reference/hg19chr21.fa /home/bqhs/reference/hg19chr21.fa
```

While the arguments look identical, they have different purposes. The first argument specifies a FASTA file of the sequence(s) we want to index. In this case contains only chromosome 21 of the human genome. The second argument specifies the name of the index, which I chose to name the same as the input.

# Align Reads

There are many progams to pick for aligning DNA reads to a reference. We've chosen [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) as our aligner. The output of Bowtie2 aligning the reads (FASTQ files) is a single file containing both the aligned and unaligned reads called a SAM file.
```bash
mkdir alignment
bowtie2 -x /home/bqhs/reference/hg19chr21.fa -1 filtered_reads/trim_SRR097849_1.fastq -2 filtered_reads/trim_SRR097849_2.fastq -S alignment/SRR097849.sam
```

We again create a new folder, `alignment`, to save our results (BAM and SAM files). The next command, `bowtie2`, executes the alignment using Bowtie2. The first argument, `-x`, specifies the reference genome. In this case, it is chromosome 21 of the human genome. `-1` and `-2` specifies the input paired-end reads. We're using our filtered and trimmed reads here. Finally `-S` specifies the name of our output SAM file.

Typically a compressed version of the SAM file, called a BAM file, is used in downstream analysis. We will use [Samtools](http://www.htslib.org/) for this conversion. Samtools also lets us sort the read by alignment position and create an index for the BAM file, both of which are usually required in downstream analysis.

```bash
samtools view -bS alignment/SRR097849.sam -o alignment/temp.bam
samtools sort alignment/temp.bam alignment/SRR097849
samtools index alignment/SRR097849.bam
```
The first command, `samtools view`, is used to convert the SAM file to a BAM file. The `-bS` flag specifies that we want our output to be a BAM file (`b`) and the input is a SAM file (`S`). The next argument is the name of our input file, which is the SAM file that Bowtie2 created for us. Lastly the `-o` argument specifies the name of the output BAM file.

To sort the read in the BAM file by alignment position, we use the `samtools sort` command. The first argument is the input BAM file and the second argument is the name for the output, sorted BAM file. Lastly, we created a index for our BAM file using the `samtools index` command that takes in a single argument, the name of the input BAM file.

# Alignment Metrics

Bowtie2 does provide some metrics on how well the alignment went, e.g. % of reads aligned, but it is useful to save a table with such metrics. This is easy to do using [Picard](https://broadinstitute.github.io/picard/index.html).

```bash
picard CollectAlignmentSummaryMetrics INPUT=alignment/SRR097849.bam OUTPUT=alignment/SRR097849_metrics.txt REFERENCE_SEQUENCE=/home/bqhs/reference/hg19chr21.fa
```

The `picard CollectAlignmentSummaryMetrics` takes three arguments. `INPUT` specifies the input BAM file, as you'd expect. `OUTPUT` specifies your preferred name of the output table. Lastly `REFERECE_SEQUENCE` should specify the same FASTA file you used as the reference for the alignment.

# Variant Detection

We will use [FreeBayes](https://github.com/ekg/freebayes) for variant detection. From the input BAM file, FreeBayes will detect variants and report the genotype(s) of the sample(s) in a VCF file.

```bash
mkdir variants
freebayes -f reference/hg19chr21.fa alignment/SRR097849.bam > variants/SRR097849.vcf
```

We first create a new folder, `variants`, to store our variant data. We then call `freebayes`. The `-f` argument specifies the FASTA file for the reference genome. This should be the same file used for alignment. Next we specify the input BAM file. Normally FreeBayes will print the output to your terminal's screen. We can redirect it and save it in an output file instead, using `> variants/SRR097849.vcf` at the end of the command.

# Variant Annotation

The raw VCF isn't too interesting, as it is missing annotation of what that variant could be (missense mutation, etc.). We can add that information using [SnpEff](http://snpeff.sourceforge.net/). 

```bash
snpEff -Xmx2g ann hg19 -s variants/stats.html variants/SRR097849.vcf > variants/SRR097849_snpEff.vcf
```

SnpEff is a Java program and usually you'll need to specify the max amount of RAM it can use. Here I specified 2 GB of RAM with `-Xmx2g`, as in my experience that is sufficient for this dataset. The next argument `ann` specifies that we want to annote our input VCF file. `hg19` specifies the reference genome database to annotate with, which is human genome build hg19 in this case. `-s` specifies the name of the output for some basic stats on our dataset. Next we specify `variants/SRR097849.vcf` as our input VCF file. We save our output using the same redirection trick we did with FreeBayes, `> variants/SRR097849_snpEff.vcf`.

It would also be useful know if any of the SNPs are known and documented. We can compare to [dbSnp](ftp://ftp.ncbi.nih.gov/snp/organisms/) using [SnpSift](http://snpeff.sourceforge.net/SnpSift.html).

```bash
SnpSift annotate /home/bqhs/dbsnp.hg19.vcf variants/SRR097849_snpEff.vcf > variants/SRR097849_dbSnp.vcf
```
Here we call SnpSift to annotate using the dbSnp VCF we downloaded from the link above for hg19. We again use redirection to save the resulting VCF file.

Typically we would need to filter the VCF file by chromosome or quality using `SnpSift filter`, but our tutorial VCF is small enough that we can work it directly. So let's turn the VCF file into a TXT file for ease of importing to Excel. SnpSift can do this for us too.

```bash
SnpSift extractFields variants/SRR097849_dbSnp.vcf ID CHROM POS REF ALT QUAL DP RO AO TYPE ANN[0].GENE ANN[0].EFFECT ANN[0].IMPACT ANN[0].BIOTYPE ANN[0].HGVS_C ANN[0].HGVS_P GEN[0].GT > variants/SRR097849.txt
```

The `SnpSift extractFields`command is rather complicated. You first specify the input VCF file. And in the end you redirect the output to a TXT file. The arguments in between specificy which fields of the VCF you want. These are what I selected.

* ID: Will contain dbSnp IDs if available
* CHROM: Chromosome number. All variants will be on chromosome 21 in this tutorial.
* POS: Position number of variant
* REF: Reference allele
* ALT: Alternate allele
* QUAL: For FreeBayes, this is the probability of locus not being polymorphic in Phred scale. E.g. If QUAL > 20, then the probability of polymorphism > 99%.
* DP: Read depth. Number of reads covering locus.
* RO: Number of reads supporting reference allele
* AO: Number of reads supporting alternative allele
* TYPE: SNP, complex, etc.
* ANN[0].GENE: SnpEff can attach multiple annotations (ANN). Here I only use the first one, ANN[0]. GENE lists the gene symbol.
* ANN[0].EFFECT: Predicted effect, e.g. missense variant
* ANN[0].IMPACT: Predicted impact of mutation
* ANN[0].BIOTYPE: E.g. protein coding, pseudogene
* ANN[0].HGVS_C: Variant using HGVS notation (DNA level)
* ANN[0].HGVS_P: Variant using HGVS notation (protein level) if applicable
* GEN[0].GT: Typically you'd analyze multiple samples, and GEN[0], GEN[1], etc. would denote them. We only have one sample in this tutorial, so GEN[0].GT refers to this sample's genotype.
