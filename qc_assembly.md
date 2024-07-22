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

Let's do a quick check of how many reads we have in each file. As with most raw read data, these files are stored as `gzip` archives, so we first `gunzip` the file to the standard output at the command line and then pipe that result to another command line function. Here we pipe a `gunzip -c` command (where the `-c` flag indicates the result should go to standard out) to a `sed -n '1~4p'` command, which prints every 4th line of the result. Then that output is further piped to `wc -l`, which counts the number of lines in the output:

```
gunzip -c tara_reads_R1.fastq.gz | sed -n '1~4p' | wc -l
gunzip -c tara_reads_R2.fastq.gz | sed -n '1~4p' | wc -l
```

This should show us that there are roughly 1.5 million reads in each file, which is a mangeable number for MEGAHIT in a tutorial. Now let's double check the length of the reads in each file. Here we are using a bash loop that will iterate over a particular set of lines selected from the reads. Specifically, it will pipe the `gunzip -c` output to the command `head -n 20` which chooses the first 20 lines of the fastq file (i.e., the first 5 reads), and then to the command `sed -n '2~4p'` which will take every 4th line starting with the _2nd_ (i.e., the line containing the read sequence). Each of those lines is then piped to `wc -c` which will count the number of characters in the line and print it (i.e., the read length):

```
for i in $(gunzip -c tara_reads_R1.fastq.gz | head -n 20 | sed -n '2~4p'); do echo $i | wc -c; done;
for i in $(gunzip -c tara_reads_R2.fastq.gz | head -n 20 | sed -n '2~4p'); do echo $i | wc -c; done;
```

From this we can tell that we have paired-end reads where each read 127 bases long. We can do more detailed QC checking in the following section using FastQC.

### Quality Control
