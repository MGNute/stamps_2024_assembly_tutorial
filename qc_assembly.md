# Metagenome QC & Assembly Tutorial
This tutorial will walk us through several of the steps required to quality control a metagenomic dataset (including adapter trimming) and assemble the reads into contigs. Much of this tutorial is pulled from [Hadrien Gourle's tutorial site](https://www.hadriengourle.com/tutorials/meta_assembly/) including the particular dataset we will test on. Some of the commands have been changed to reflect the tools we have discussed, but many of the commands are unchanged from that version.

### Getting the Data
The tutorial linked above uses a dataset that is ideal for this lesson because it is small enough that a metagenome assembly can be run on it to reasonable completeness in a short time. It contains simulated Illumina HiSeq reads from 20 genomes discovered as part of the Tara Oceans Project. Let's get the data:

```
mkdir -p ~/data
cd ~/data
curl -O -J -L https://osf.io/th9z6/download
curl -O -J -L https://osf.io/k6vme/download
chmod -w tara_reads_R*
```

##### Basic Data Integrity Checking

Let's start with a couple of commands to just check that the data is what we expect it to be.

First let's check how many reads we have in each file. As with most raw read data, these files are stored as `gzip` archives, so we first `gunzip` the file to the standard output at the command line and then pipe that result to another command line function. Here we pipe a `gunzip -c` command (where the `-c` flag indicates the result should go to standard out) to a `sed -n '1~4p'` command, which prints every 4th line of the result. Then that output is further piped to `wc -l`, which counts the number of lines in the output:

```
gunzip -c tara_reads_R1.fastq.gz | sed -n '1~4p' | wc -l
gunzip -c tara_reads_R2.fastq.gz | sed -n '1~4p' | wc -l
```

First we note that both files have about the same number of reads, consistent with what we'd expect for paired-end data. Also, we have 1.5 million reads in total which is a mangeable number for MEGAHIT in a tutorial. 

Now let's double check the length of the reads in each file. Here we are using a bash loop that will iterate over a particular set of lines selected from the reads. Specifically, it will pipe the `gunzip -c` output to the command `head -n 20` which chooses the first 20 lines of the fastq file (i.e., the first 5 reads), and then to the command `sed -n '2~4p'` which will take every 4th line starting with the _2nd_ (i.e., the line containing the read sequence). Each of those lines is then piped to `wc -c` which will count the number of characters in the line and print it (i.e., the read length):

```
for i in $(gunzip -c tara_reads_R1.fastq.gz | head -n 20 | sed -n '2~4p'); do echo $i | wc -c; done;
for i in $(gunzip -c tara_reads_R2.fastq.gz | head -n 20 | sed -n '2~4p'); do echo $i | wc -c; done;
```

From this we can tell that we have paired-end reads where each read is 127 bases long. We can do more detailed QC checking in the following section using FastQC.

### Quality Control
Now let's move on to true sequencing QC with FastQC:

```
mkdir -p ~/results
cd ~/results
ln -s ~/data/tara_reads_* .
fastqc tara_reads_*.fastq.gz
```

Open the html file that was created in the results folder. What observations occur to you? Did any of the QC checks fail? 

### Adapter Trimming
##### Installing BBMap
We are going to attempt to install the BBMap tool suite before this tutorial, however if that does not work installation is simple and only requires downloading and extracting the source file from the JGI website. The following commands will do the trick. These commands do the following, in order: 1) Download the latest version from the JGI sourceforge site (which can be found [here](https://sourceforge.net/projects/bbmap/), which is linked to in the [download instructions](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/)). 2) Extract the tar file. 3) Create a local environment variable noting the directory we have extracted this file to. (Note: installing a software to a project-specific folder like this is a bit sloppy and not a best practice, but it will do for the moment.)

```
wget https://sourceforge.net/projects/bbmap/files/BBMap_39.08.tar.gz
tar -xzf BBMap_39.08.tar.gz
bbmappath=$(pwd)/bbmap
```

##### Running BBduk
Let's run BBduk on the reads using the default adapters file provided with the software:

```
$bbmappath/bbduk.sh in1=~/data/tara_reads_R1.fastq.gz in2=~/data/tara_reads_R2.fastq.gz out1=~/data/tara_reads_R1_trimmed.fastq.gz out2=~/data/tara_reads_R2_trimmed.fastq.gz ref=$bbmappath/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
```

What does the output tell us about how many reads were trimmed and/or removed entirely?

### Assembly



