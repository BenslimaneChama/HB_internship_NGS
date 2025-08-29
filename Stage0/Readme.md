# HackBio – Stage 0
This repository contains my solutions for Stage 0 of the HackBio Internship, under the NGS Fundamentals. 
The stage included two main projects:
## Project 1 – Linux & Bash Basics for Bioinformatics
In this project, I practiced the basics of bash scripting by applying them to biological data files (`.fna` and `.gbk`):

- File and directory management: creating folders (chama, biocomputing), moving and deleting files.

- Downloading datasets: using wget to retrieve wildtype.fna and wildtype.gbk.

- File manipulation:

    * hecking whether the sequence in the .fna file was mutant or wild type using grep.
    * Extracting matching lines into a new file.
    * Counting the number of sequence lines in .gbk excluding the header.
    * Extracting metadata such as sequence length (from the LOCUS tag), source organism (from the SOURCE tag), and a list of all genes (/gene=).

- Command history and listing: using history and ls to review what was done.

_This project showed me how even simple Linux commands can help answer biological questions, like checking for mutations or extracting features from GenBank files._

## Project 2 – Conda & Bioinformatics Tools

In this project, I worked with conda to create and manage environments and install core bioinformatics tools:

- Created and activated a conda environment named funtools.

- Installed widely used bioinformatics tools from the bioconda channel:

   * *bwa* – sequence alignment
   * *blast* – sequence similarity search
   * *samtools* – manipulation of SAM/BAM/CRAM files
   * *bedtools* – comparing and analyzing genomic regions
   * *spades* – genome assembly
   * *bcftools* – working with VCF/BCF variant files
   * *fastp* – FASTQ preprocessing and quality control
   * *multiqc* – aggregating QC and analysis reports

- Verified each tool was properly installed and accessible inside the funtools environment.

_This project gave me hands-on practice with setting up a reproducible bioinformatics environment, a key skill for any genomics project._
