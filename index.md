<center>
  <h2>Friday February 15 @ 9:00 AM - 12:00 PM</h2>
  <h2>JABSOM MEB, Computer Lab 107D</h2>
  <h3>Presenters: Dr. Youping Deng, Dr. Mark Menor, and Dr. Vedbar Khadka</h3>
  Bioinformatics and Biostatistics Cores<br>
  Complementary Integrative Medicine Department,<br>
  JABSOM, University of Hawaii<br>
  UH Cancer Center Genomics & Bioinformatics Shared Resource<br>
  <img src="images/inbre-iii.png" alt="INBRE-III" style="height: 100px;"/>&nbsp;
  <img src="images/pceidr.png" alt="PCEIDR" style="height: 100px;"/>&nbsp;
  <img src="images/rmatrix.png" alt="RMATRIX-II" style="height: 100px;"/>&nbsp;
  <img src="images/ola_hawaii.png" alt="Ola Hawaii" style="height: 100px;"/>
</center>
<br>

* [Server Information](#server-information)
* [Dataset Information](#dataset-information)
* [Linux Command Line](#linux-command-line)
* [Downloading Data](#downloading-data)
* [Quality Control](#quality-control)
* [Quality Filtering and Trimming](#quality-filtering-and-trimming)
* [Genome Assembly](#genome-assembly)
* [Assembly Evaluation](#assembly-evaluation)
* [Assembly Visualization](#assembly-visualization)
* [Genome Annotation](#genome-annotation)

# Server Information

We will be using a Linux server (Ubuntu) for this workshop. You will need to use [PuTTY](https://www.chiark.greenend.org.uk/~sgtatham/putty/) on Windows or `ssh` on the Terminal in macOS.

Windows software:
* This computer lab requires 32-bit versions
* [FileZilla](https://filezilla-project.org/download.php?type=client)
* [PuTTY](https://the.earth.li/~sgtatham/putty/latest/w32/putty-0.70-installer.msi)

With PuTTY and FileZilla you can connect to tutorial server:
* Address: 168.105.161.70
* Port: 22
* Log-in information to be provided at workshop


# Dataset Information

[Zaire ebolavirus genome sequencing from 2014 outbreak in Sierra Leone](https://www.ncbi.nlm.nih.gov/sra/SRR1553425). This is a paired-end experiment and thus we have two FASTQ files, one per end. FASTQ is the file format for raw reads and it contains information on the sequence and the quality (accurracy) of each base call. Check out some details of the FASTQ format on the [Wikipedia entry](https://en.wikipedia.org/wiki/FASTQ_format)

# Linux Command Line

A complete tutorial on the Linux command line is a full workshop on its own. Please refer to other resources, e.g. [Ryan's Tutorials](https://ryanstutorials.net/linuxtutorial/). For our purposes, we just need to know some basic navigation.

You are signed on as a guest on the tutorial server, so when you first log in, your current working directory is `/home/guest2019`. Everyone in the workshop should work in a different directory, so let's make a folder named after yourself.

```bash
mkdir write_your_UH_username_here
```

For example, `mkdir mmenor` if your UH username is _mmenor_. The command is short for "make directory." Let us confirm you created a directory by listing the contents of your working directory,

```bash
ls
```

`ls` is short for "list". You should see a folder with your UH username, among the other workshop participants. Let's now change our working directory to your newly created folder using the change directory command,

```bash
cd write_your_UH_username_here
```

For example, `cd mmenor` in my case. `cd` is short for "change directory." You can run the `ls` command again to confirm you've changed folders and your current working directory is empty. Now we're ready to start the data analysis.

# Downloading Data

In the interest of time, we'll only a 100,000 read subsample of the dataset (SRR1553425) from the Sequence Read Archive (SRA). NCBI provides tools in Linux to download and extract datasets called SRA Tools.

```bash
fastq-dump --split-files -X 100000 SRR1553425
```

The `fastq-dump` command downloads the first 100,000 reads from SRR1553425 and saves into a pair of FASTQ files. The `-X` parameter is used to specify the number of reads (default is all reads). The `--split-files` parameter specifies that the data should be saved into two separate FASTQ files, which is typical for paired-end sequencing. Each end gets its own file. Use the `ls` command to verify the files `SRR1553425_1.fastq` and `SRR1553425_2.fastq` have been downloaded.

You can browse the FASTQ files using the `less` command, e.g. for the first FASTQ file:

```bash
less SRR1553425_1.fastq
```

The arrow keys scroll the file up and down and the q key exits the viewer. You can view the second FASTQ file by changing the file name in the command.

# Quality Control

Garbage in is garbage out, so let's see check if we have quality issues with our raw reads. We'll be using [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for this purpose.

```bash
fastqc SRR1553425_1.fastq SRR1553425_2.fastq
```
The `fastqc` command simply requires you to list all the FASTQ files to analyze.

The output of FastQC is an HTML page with plots that unfortunately we cannot view with an SSH remote terminal. Let's download the results using [FileZilla](https://filezilla-project.org/). Connect using the same guest credentials you used on SSH. Navigate to `/home/guest2019/write_your_UH_username_here` to find your FASTQC results. The report files names should be SRR1553425_1_fastqc.html and SRR1553425_2_fastqc.html.

The accuracy of a base's identifiy is measured using [Phred quality scores](https://en.wikipedia.org/wiki/Phred_quality_score). The higher the score, the better. Roughly, a Phred score of 20 or above (99+% accuracy) is great.

Check out the [FastQC manual](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf) for more information on each plot.

# Quality Filtering and Trimming

Typically base quality degrades at the ends of the read, particularly at the 5' end. In the case of this dataset, the quality is already great, but for the sake of this exercise, we'll trim low quality bases of the ends our reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) anyway. Trimmomatic can also trim off adapters, which also has been detected by FASTQC.

First let's download a FASTA file of common adapters.

```bash
curl -OL https://raw.githubusercontent.com/BioInfoTools/BBMap/master/resources/adapters.fa
```

With that, we can now do quality and adapter trimming,

```bash
trimmomatic PE SRR1553425_1.fastq SRR1553425_2.fastq trimmed_1.fastq unpaired_1.fastq trimmed_2.fastq unpaired_2.fastq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:20 TRAILING:20 AVGQUAL:20 MINLEN:20
```

There is a lot going on in the `trimmomatic` command. First we specify `PE` to say we're working with paired-end data. Next we specify our pair of input FASTQ files (`SRR1553425_1.fastq` and `SRR1553425_2.fastq`). Trimmomatic outputs two FASTQ files per input FASTQ, which gives us a total of four output filenames to specify. You can name them whatever you want, but the names I've given them indicate what the file is.

During the trimming process, entire reads may be filtered out. Thus the output is split into two files: one file that contains reads that are still paired-end (both reads survived filtering) and one file that contains reads that are now single-end because their mate was filtered out. For the first mate file `SRR1553425_1.fastq`, the file `trimmed_1.fastq` contains the trimmed paired-end reads, while `unpaired_1.fastq` contains the now single-ends reads. 

The remaining parameters specify how the trimming should be executed. First we specify trimming of adapters with `ILLUMINACLIP:adapters.fa:2:30:1`. This specifies using the FASTA file of adapters we downloaded earlier. The three numbers control how strict the adapter matching will be.

The parameters `LEADING:20 TRAILING:20` specifies the Phred qualty score cutoffs for the low-quality trimming of the 5' and 3' ends of the read. In this case, we trim a base off the 5' or 3' end if its Phred score is below 20. On the other hand, `AVGQUAL:20` specifies the overall average quality of the read needs to be at least Phred score 20. If not, the read will be removed from the dataset completely. Lastly `MINLEN:20` specifies that after trimming, if the read is less the 20 bases long, we remove it from the dataset. Very short reads aren't very useful in genome assembly.


# Genome Assembly

There are a lot of genome assemblers out there for different purposes. [SPAdes](http://cab.spbu.ru/software/spades/) is designed for single- and multi-cell bacterial datasets. SPAdes is easier to use than Velvet, one of the more popular genome assemblers, since you can specify multiple k-mer lengths and it will handle combining the results for you.

```bash
cat unpaired_1.fastq unpaired_2.fastq > unpaired.fastq
spades.py -k 21,33,55,77,89,99 --careful -o spades_output -1 trimmed_1.fastq -2 trimmed_2.fastq -s unpaired.fastq
```

First the `cat` command is used to "concatenate" our two unnpaired FASTQ files into a single FASTQ file.

Next the `spades.py` command you specify a list of k-mer lengths with the `-k` parameter. The values must be positive odd integers. This parameter can be omitted and default k-mer lengths used, but I found the results subpar. Next `--careful` tells SPAdes to be conservative in its matching and reduce the number of mismatches and indels. This is a computationally taxing operation so it is only recommend for small genomes such as this virus. The `-o` for "output" is used to name the folder where all the SPAdes results will be saved. In this case a folder named "spades_output" will be created. SPAdes generates a lot of files, so it's good to put its results in a separate folder. Finally you specify your paired-end reads with the `-1` and `-2` parameters and you unpaired reads with `-s`.

Download the `spades_output` folder using FileZilla. Will be using its contents soon.

The important SPAdes results:
* `assembly_graph.fastg`: Useful for visualing the assembly graph. Some visualization tools let's you manually edit the assembly graph to resolve ambiguities that the assembler couldn't solve. In this case, we will have one very large contig dominating the visualization that's clearly the whole genome and no manual work is necessary.
* `contigs.fasta`: A fasta file of each assembled contig
* `scaffolds.fasta`: A fasta file of each assembled scaffold


# Assembly Evaluation

We will use [Quality Assement Tool for Genome Assemblies (QUAST)](http://quast.sourceforge.net/quast) to evaluate our assembly. This is done to compare different assemblies (results from different k-mer lengths, assembler programs, etc.). Typically in _de novo_ sequencing you would not have a reference genome to compare to, but we do in this case. If you provide QUAST a reference genome it can caculate additional stats like the number of mistmatches and indels compared to the reference.

```bash
quast -R /home/bqhs/workshop2019/KJ660346.fasta spades_output/scaffolds.fasta
```

The `quast` command is simple. You first specify the reference genome to compare to with `-R`. Then you specify the assembly you want to evaluate, in this case the final scaffolds from SPAdes. The results by default will save into a folder called `quast_results`. Download that folder using FileZilla. In the folder, there is a file called `report.html` that will list all the assembly metrics. You can hover your mouse cursor over a statistic's name to learn more about it.


# Assembly Visualization

We will use [Bandage](https://rrwick.github.io/Bandage/) to visualize the assembly graph. Download and install this on your computer. When it is open, hit File->Load graph in the main menu. Here pick your final SPAdes assembly graph that you downloaded, `spades_output/assembly_graph.fastg`. Then hit the "Draw graph" button on the left side. You'll see we got clean assembly results.

To see a messier assembly result, you can load the graph `spades_output/K21/assembly_graph.fastg`. This is the assembly graph using 21-mers and it didn't turn out as well.


# Genome Annotation

Lastly, let's annotate possible genes our assembled genome. First we'll find open reading frames (ORFs) using NCBI's [ORFfinder](https://www.ncbi.nlm.nih.gov/orffinder/). Copy and paste _only_ the first contig from your `spades_output/scaffolds.fasta` file. Submit with the default settings.

To download all ORF results into a FASTA file, on the bottom right box, press "Mark subset..." and pick "All ORFs." The press "Download marked set." By default, this will save the protein predictions. You can also download CDS predictions by clicking on the drop down and changing "Protein FASTA" to "CDS FASTA".

With the list of ORFs, you can now use [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to try to predict the function of the ORF. Use Nucleotide BLAST for the CDS regions and Protein BLAST for the protein file.
