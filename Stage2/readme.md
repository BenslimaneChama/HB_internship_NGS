# Unmasking Molecular Signatures of Bipolar II Disorder
## Parrt 1 : Linux (Pretreatement and Alignement)
In this part we worked on linux, in order to prepare and analyse sequencing dataset. The goal was to from raw FastQ, that we downloaded using SRA explorer files, to actually check their quality, to clean the sequences (trimming) to align them on the human genome reference whcih is (GRCh38à, and then to generate a gene count file that will be useful for further analyses using R.
### Folders organisation
We have created a clear tree structure to separate each stage of the pipeline:
```
mkdir -p project_human/{Dataset,qc,trim,reference,genomeIndex,mapped,counts,IGV}
cd project_human
```
`Dataset/` → downloaded FASTQ files

`qc/` → quality reports (FastQC/MultiQC)

`trim/` → cleaned reads (fastp)

`reference/` → genome and annotations (FASTA, GTF)

`genomeIndex/` → STAR index

`mapped/` → BAM alignments sorted by coordinates

`counts/` → count table (featureCounts)

`IGV/` → alignments + index for visualization

### Dataset downloading

The dataset is provided from the project `PRJNA1253883` (iPSC-derived microglia, healthy vs bipolar disorder). And it is available in ther EBI-ENA server.
we looked for them using SRA explorer, that gave us the code to download it.
Example of downloading using `curl` :
```
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/064/SRR33243164/SRR33243164_1.fastq.gz -o SRR33243164_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/064/SRR33243164/SRR33243164_2.fastq.gz -o SRR33243164_2.fastq.gz
```
We've got 8 datasetm each of them got the R1 and R2.

### Quality control

Before any treatement, the quality was verified using `FastQC`
```
fastqc Dataset/*.fastq.gz -o qc/
multiqc qc/ -o qc/
```
`FastQC`: provides individual reports per file (quality, adapters, enriched sequences).
`MultiQC`: aggregates all reports into a single interactive page → easier to compare.

### Trimming and filtering 

In order to improve the quality, we use `fastp`.
```
fastp \
    -i "$r1" -I "$r2" \
    -o "trim/${base}_1.trim.fastq.gz" \
    -O "trim/${base}_2.trim.fastq.gz" \
    -h "trim/${base}_report.html" \
    -j "trim/${base}_report.json" \
    --detect_adapter_for_pe
```
And then we make another fastp to visualize the changes. 

### Reference Genome

We used the Homo sapiens genome `GRCh38` (release 109 Ensembl).
```
wget -O reference/homo_sapiens.fa.gz  ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -O reference/homo_sapiens.gtf.gz ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
gunzip reference/homo_sapiens.fa.gz reference/homo_sapiens.gtf.gz
```
The `FASTA` file contains the genome sequence, and the `GTF` file describes the annotation : genes and exons.

### STAR indexation

for a fast alignement, STAR buildds the genome index.

```
STAR --runMode genomeGenerate \
     --genomeDir genomeIndex \
     --genomeFastaFiles reference/homo_sapiens.fa \
     --sjdbGTFfile reference/homo_sapiens.gtf \
     --sjdbOverhang 149 \
     --runThreadN 8
```

### Reads alignement 

each R1/R2 paires are aligned with STAR : 

```
STAR --runThreadN 8 \
       --genomeDir genomeIndex \
       --readFilesIn "$r1" "$r2" \
       --readFilesCommand zcat \
       --outFileNamePrefix mapped/${sample}_ \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMattributes All
```
we then can check the verification using
`samtools flagstat mapped/*bam | head
`
### BAM indexation and IGV

in order to visualise the alignements in IGV we needed to move all `bam` files into another folder called IGV, and then index all the 8 samples using `samtools`.
```
for bam in IGV/*bam; do
  samtools index "$bam"
done
```
### Gene counting 

the counting was realised using `featureCounts` (Subread) :
```
featureCounts -T 8 -p -B -C \
  -a reference/homo_sapiens.gtf \
  -o counts/gene_counts.txt \
  mapped/*bam
```
`-p `→ paired-end mode
`-B` → correctly paired reads
`-C` → exclusion of chimeric fragments

At this point, the data is ready for differential analysis in R.
