# Microbial Project
This project was done in order to be able to processes and analyzes Illumina paired-end sequencing data from 100 outbreak samples. The workflow includes raw data quality control, trimming, assembly, quality assessment, organism identification, and antimicrobial resistance profiling.

## Methodes

### Data Download
The first step was to actually download the script that would help us to download the 100 dataset, we used wget.
```
wget https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/SA_Polony_100_download.sh
```
And then changed it's mode so we can execute it easily. The output was then downloaded in another directory so we can manipulate it easily.
### Sample Preparation
Once the data is downloaded, we had to be sure that we have all the dataset and nothing is missing, for each sample we got the R1 and R2, we sorted them, and then created a loop to find if there is any sample with R1 or R2 missing
```
ls SA_Polony_data/*_1.fastq.gz \
  | sed -E 's#.*/##; s/_1\.fastq\.gz$//' \
  | sort > samples.txt
```
the loop 
```
while read S; do
  r1="SA_Polony_data/${S}_1.fastq.gz"
  r2="SA_Polony_data/${S}_2.fastq.gz"
  if [[ ! -f "$r1" ]] || [[ ! -f "$r2" ]]; then
    echo "INCOMPLETE: $S"
  fi
done < samples.txt
```
### Quality Control
Now the dataset is ready and we didn't left anything behind, we checked raw reads with FastQC and summarized with MultiQC.
```
mkdir -p qc/raw_fastqc qc/multiqc_raw
find SA_Polony_data -name "*.fastq.gz" -print0 \
 | xargs -0 -n 8 -P 8 fastqc -q -o qc/raw_fastqc
multiqc -q qc/raw_fastqc -o qc/multiqc_raw
```
### Trimming and Filtering
In order to assure a better quality, we had to trim, so we can remove the adapter and low-quality base trimming using the command fastp, then in order to compare we used another multiqc
```
while read -r S; do
  fastp \
    -i SA_Polony_data/${S}_1.fastq.gz \
    -I SA_Polony_data/${S}_2.fastq.gz \
    -o trimming/${S}_1.trim.fastq.gz \
    -O trimming/${S}_2.trim.fastq.gz \
    --detect_adapter_for_pe \
    --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
    --qualified_quality_phred 20 --length_required 50 \
    --thread 8 \
    -h trimming/report/${S}_fastp.html \
    -j trimming/report/${S}_fastp.json
done < samples.txt
```
### De Novo assembly
Now we are ready to actually assemblies it using SPAdes.
```
while read -r S; do
  spades.py \
    -1 trimming/${S}_1.trim.fastq.gz \
    -2 trimming/${S}_2.trim.fastq.gz \
    -o assembly/${S} \
    -t 8 -m 12
done < samples.txt
```
### Assembly Collection & Visualization
To make our work easier, we runnded a loop so we can assemble all the contigs and the fastg into another directory for further analysis and visualization (e.g., Bandage).
```
while read -r S; do
  [[ -s assembly/${S}/contigs.fasta ]] && cp assembly/${S}/contigs.fasta "contigs_ready/${S}.fasta"
  [[ -s assembly/${S}/assembly_graph.fastg ]] && cp assembly/${S}/assembly_graph.fastg "bandage_fastg/${S}.fastg"
done < samples.txt
```
### Assembly Quality Assessment
We then ran QUAST to evaluate genome assembly quality across all contigs.
```
quast.py contigs_ready/*.fasta -o quast_all --threads 8 -q
```
### Organism Identification
Representative contigs were analyzed with BLASTn against NCBI nt (remote).
Due to genome size, some contigs were split before blasting.

### Antimicrobial Resistance Profiling
We screened assemblies for AMR genes using ABRicate (ResFinder DB).


## Results
### First MultiQC
![First MultiQC](Images/mutliqc_1.png)
![First MultiQC](Images/multiqc1_1.png)

This plot shows the total number of reads per sample. Blue bars represent unique reads, while black bars indicate duplicates. Most samples have sufficient coverage for downstream assembly, although variability between samples is visible

![First MultiQC](Images/multiqc1_2.png)

Average read quality is generally high (>30, green zone), indicating reliable sequencing. A slight drop in quality is observed toward the end of the reads, which is typical for Illumina data.

![First MultiQC](Images/multiqc1_5.png)

The GC distribution peaks around ~38–40%, consistent with the genome of Listeria monocytogenes. Secondary peaks may suggest contamination or sequencing bias in some samples.

![First MultiQC](Images/multiqc1_3.png)

Most reads exhibit high Phred scores (>30). Very few reads fall below the acceptable threshold, confirming the dataset is overall high quality.

![First MultiQC](Images/multiqc1_4.png)

The proportions of A, T, G, and C stabilize across the read length after the initial bases. The small imbalance at the beginning is typical of Illumina priming and does not indicate a major problem. 

![First MultiQC](Images/multiqc1_6.png)

Shows the percentage of undetermined bases ("N") across all read positions. The low values confirm that base calling was accurate with minimal ambiguity.

![First MultiQC](Images/multiqc1_9.png)

Highlights sequences that appear more frequently than expected. Only a very small fraction (<1% of reads) are overrepresented, suggesting minimal contamination or technical artifacts.

![First MultiQC](Images/multiqc1_10.png)

Displays the presence of sequencing adapter sequences. The increasing signal at the end of the reads indicates some adapter contamination, which trimming tools (like fastp) are designed to remove.

![First MultiQC](Images/multiqc1_7.png)

Plots the read length across all samples. The sharp peak at ~300 bp reflects uniform sequencing library preparation, as expected for Illumina paired-end reads.

![First MultiQC](Images/multiqc1_8.png)

Indicates the percentage of duplicated reads per library. Most sequences have low duplication, but some samples show higher duplication levels, which can result from high coverage or PCR bias.
#### MultiQC summary

![First MultiQC](Images/multiqc1_11.png)

This heatmap provides an overview of quality metrics across all 200 samples. Each column corresponds to a QC module (e.g., per-base quality, GC content, duplication levels), and each row represents one sample. Green indicates passing QC, yellow means warnings, and red marks failures. Overall, most samples passed the key metrics, with warnings or failures mainly in GC content and duplication, which are common for bacterial genomes.

### Second MultiQC

![Second MultiQC](Images/multiqc2_.png)

![Second MultiQC](Images/multiqc2_1.png)

This bar plot shows the proportion of reads retained after filtering. The majority of reads passed quality filters (blue), with only small fractions removed due to short length, excess Ns, or low complexity. This indicates good sequencing quality overall.

![Second MultiQC](Images/multiqc2_2.png)

This plot represents the distribution of insert sizes across all libraries. The majority of insert sizes are within a consistent range, with a clear peak around the expected library size. A few variations are observed, but the global profile confirms correct library construction.

![Second MultiQC](Images/multiqc2_3.png) ![Second MultiQC](Images/multiqc2_4.png)

Both Read 1 and Read 2 maintain high-quality scores across most of the read length, with values above Q30. A minor decline appears toward the ends of reads, which is typical of Illumina sequencing, but overall quality is strong.

![Second MultiQC](Images/multiqc2_5.png) ![Second MultiQC](Images/multiqc2_6.png)

The GC content remains stable across reads, averaging around 40–45%, which is consistent with Listeria monocytogenes. This uniformity indicates no major contamination or technical bias.

![Second MultiQC](Images/multiqc2_7.png) ![Second MultiQC](Images/multiqc2_8.png) 

The proportion of ambiguous bases (N) is near zero across read positions for both read pairs. This confirms high base-calling accuracy and effective filtering.

### Example of assembly graph visualization in Bandage

![Bandage](Images/BANDAGE.png)

This graph illustrates contig connectivity from one of the assembled samples. Each color represents a different contig, and the nodes highlight branching points in the assembly. The visualization helps detect assembly quality issues such as repeats, breaks, or unresolved structures.

### Quast 

![Quast](Images/quast_3.png)

This Nx plot shows the contig length distribution across assemblies. Each line represents an isolate, and the stepwise decline illustrates the number and size of contigs needed to cover a percentage of the genome. Higher, flatter curves indicate assemblies with longer contigs and better continuity.

![Quast](Images/quast1_1.png)

The GC content distribution of all assemblies is plotted here. The sharp peak around ~38–40% corresponds to the expected GC content of Listeria monocytogenes. The consistent curves across isolates confirm homogeneity in genomic composition and absence of major contamination.

![Quast](Images/quast2.png)

This cumulative plot shows how total genome length is recovered as contigs are added from longest to shortest. Most assemblies reach ~3 Mbp, which is in line with the known genome size of Listeria monocytogenes. Curves with earlier plateaus represent assemblies with fewer, larger contigs (better assembly quality).

### Genome Identification 

![Genome Identification](Images/NCBI.png)

BLAST alignment of assembled contigs against the NCBI nt database confirmed Listeria monocytogenes as the organism of origin. Multiple reference strains (e.g., OB080183, 19-02390, QIO055, IZ-SAM, BL91/023) showed 100% query coverage and 100% identity, with negligible E-values (0.0). This provides strong evidence that all isolates analyzed belong to Listeria monocytogenes, consistent with the outbreak origin.

### AMR 

![AMR](Images/AMR.png)

```
     10 fosX_2
```
Example output of AMR screening for 10 Listeria monocytogenes outbreak isolates using ResFinder. All 10 samples carry the fosX_2 gene, which confers resistance to fosfomycin. The alignments show high sequence identity (~92.56%) and nearly complete coverage (~99.75%), confirming the consistent presence of this resistance determinant across isolate

## Public Health Discussion 

### Listeria monocytogenes as a Pathogen

![Listeria monocytogenes](Images/Listeria_monocytogenes.jpg)


Listeria monocytogenes is a Gram-positive, facultative intracellular bacterium and the causative agent of listeriosis, a serious foodborne disease. It is particularly dangerous for vulnerable populations, including pregnant women, newborns, the elderly, and immunocompromised individuals. Ingestion of contaminated food products such as dairy, ready-to-eat meats, or vegetables can lead to systemic infection, causing septicemia, meningitis, or neonatal infections. Due to its ability to survive and grow under refrigeration and in harsh environments, L. monocytogenes represents a major challenge for food safety and public health systems worldwide.

### The fosX Gene and Fosfomycin Resistance

![fosX Gene and Fosfomycin Resistance](Images/fosx_resistance.jpg)

*Mattioni Marchetti, V., Hrabak, J., & Bitar, I. (2023). Fosfomycin resistance mechanisms in Enterobacterales: an increasing threat. Frontiers in cellular and infection microbiology, 13, 1178547.*

Our analysis identified the fosX_2 gene in all tested isolates. The fosX gene encodes a glutathione transferase enzyme that inactivates fosfomycin, an antibiotic that targets bacterial cell wall synthesis. Although fosfomycin is not the primary therapeutic option against Listeria, the presence of this resistance marker indicates that the bacterium has the genetic capacity to withstand treatment with this molecule. Importantly, the fosX gene is considered intrinsic to Listeria monocytogenes, meaning that most isolates naturally carry it. Thus, fosfomycin is generally not considered effective in clinical management of listeriosis.

### Treatement with antibiotics

According to the StatPearls entry on Listeria monocytogenes (NCBI Bookshelf, 2023), the antibiotic agents of choice for invasive listeriosis are intravenous ampicillin or penicillin G. In patients with a penicillin allergy, trimethoprim-sulfamethoxazole is the preferred alternative. Listeria exhibits *intrinsic resistance to cephalosporins, and gentamicin is sometimes added (especially in neonatal or severe cases) for broader coverage or synergy.

https://www.ncbi.nlm.nih.gov/books/NBK534838/
