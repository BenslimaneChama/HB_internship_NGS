#Create a project folder and go inside
$mkdir -p ~/chama/project_microbes
cd ~/chama/project_microbes

#Download the dataset script
$wget https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/SA_Polony_100_download.sh

#Make it executable
$chmod +x SA_Polony_100_download.sh

#Choose a directory where files will be saved, and run the script from inside it
$mkdir -p SA_Polony_data
$cd SA_Polony_data && bash ../SA_Polony_100_download.sh

#Run this to generate a clean list of sample names
$ls SA_Polony_data/*_1.fastq.gz \
  | sed -E 's#.*/##; s/_1\.fastq\.gz$//' \
  | sort > samples.txt

#Count how many FASTQ files we have
ls SA_Polony_data/*.fastq.gz | wc -l

#Check if each sample has both R1 and R2 files
$while read S; do
  r1="SA_Polony_data/${S}_1.fastq.gz"
  r2="SA_Polony_data/${S}_2.fastq.gz"
  if [[ ! -f "$r1" ]] || [[ ! -f "$r2" ]]; then
    echo "INCOMPLETE: $S"
  fi
done < samples.txt


#QC of raw reads (first FastQC and MultiQC)
$mkdir -p qc/raw_fastqc qc/multiqc_raw
$find SA_Polony_data -name "*.fastq.gz" -print0 \
 | xargs -0 -n 8 -P 8 fastqc -q -o qc/raw_fastqc
$multiqc -q qc/raw_fastqc -o qc/multiqc_raw

#Trimming (fastp) + QC after trimming
$mkdir -p trimming/report trimming/qc
$while read -r S; do
  fastp \
    -i  SA_Polony_data/${S}_1.fastq.gz \
    -I  SA_Polony_data/${S}_2.fastq.gz \
    -o  trimming/${S}_1.trim.fastq.gz \
    -O  trimming/${S}_2.trim.fastq.gz \
    --detect_adapter_for_pe \
    --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
    --qualified_quality_phred 20 --length_required 50 \
    --thread 8 \
    -h trimming/report/${S}_fastp.html \
    -j trimming/report/${S}_fastp.json
done < samples.txt

$mkdir -p qc/trim_fastqc qc/multiqc_trim
$find trimming -name "*.trim.fastq.gz" -print0 \
 | xargs -0 -n 8 -P 8 fastqc -q -o qc/trim_fastqc
$multiqc -q qc/trim_fastqc -o qc/multiqc_trim

#De novo assembly (While using SPAdes)
$mkdir -p assembly
$while read -r S; do
  spades.py \
    -1 trimming/${S}_1.trim.fastq.gz \
    -2 trimming/${S}_2.trim.fastq.gz \
    -o assembly/${S} \
    -t 8 -m 12
done < samples.txt

#Collect contigs + assembly graphs.fastg in order to visulise them for badage application
$mkdir -p contigs_ready bandage_fastg
$while read -r S; do
  [[ -s assembly/${S}/contigs.fasta ]] && cp assembly/${S}/contigs.fasta "contigs_ready/${S}.fasta"
  [[ -s assembly/${S}/assembly_graph.fastg ]] && cp assembly/${S}/assembly_graph.fastg "bandage_fastg/${S}.fastg"
done < samples.txt

#Assembly quality (QUAST; multi-sample report)
$mkdir -p quast_all
$quast.py contigs_ready/*.fasta -o quast_all --threads 8 -q

#Organism identification (quick BLAST option) 
#Before using BLAST either by code or by web, we had to split, because the genome size was way too large (3Mb)
#We used SkyBLAST and we got Listeria monocytogenes with 100% for Per. Identity, another way to actually run blast with code is by this
$q="contigs_ready/$(head -n1 samples.txt).fasta"
$blastn -task megablast -query "$q" -db nt -remote \
       -entrez_query "Bacteria[Organism]" \
       -max_target_seqs 25 -evalue 1e-20 \
       -outfmt "6 qseqid sacc staxids ssciname pident length evalue bitscore stitle" \
       -out blast_results/$(basename "$q" .fasta).tsv

#AMR detection (ABRicate / ResFinder) + summaries
$mkdir -p amr
$abricate --db resfinder contigs_ready/*.fasta > amr/all_amr.tab
$abricate --summary amr/all_amr.tab > amr/all_amr_summary.tsv
$cut -f6 amr/all_amr.tab | grep -v '^GENE' | sort | uniq -c | sort -nr > amr/all_genes_count.txt #in order to list genes and their counts across isolates


