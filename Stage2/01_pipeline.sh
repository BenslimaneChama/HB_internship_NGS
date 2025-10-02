#Folder Organisation

mkdir project_human && cd project_human #main folder
mkdir Dataset qc trim genomeIndex mapped counts reference IGV

#Dataset downloading from SRA explorer

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/064/SRR33243164/SRR33243164_1.fastq.gz -o SRR33243164_RNA_seq_of_iPSC_derived_microglia_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/064/SRR33243164/SRR33243164_2.fastq.gz -o SRR33243164_RNA_seq_of_iPSC_derived_microglia_2.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/065/SRR33243165/SRR33243165_1.fastq.gz -o SRR33243165_RNA_seq_of_iPSC_derived_microglia_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/065/SRR33243165/SRR33243165_2.fastq.gz -o SRR33243165_RNA_seq_of_iPSC_derived_microglia_2.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/066/SRR33243166/SRR33243166_1.fastq.gz -o SRR33243166_RNA_seq_of_iPSC_derived_microglia_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/066/SRR33243166/SRR33243166_2.fastq.gz -o SRR33243166_RNA_seq_of_iPSC_derived_microglia_2.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/067/SRR33243167/SRR33243167_1.fastq.gz -o SRR33243167_RNA_seq_of_iPSC_derived_microglia_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/067/SRR33243167/SRR33243167_2.fastq.gz -o SRR33243167_RNA_seq_of_iPSC_derived_microglia_2.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/068/SRR33243168/SRR33243168_1.fastq.gz -o SRR33243168_RNA_seq_of_iPSC_derived_microglia_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/068/SRR33243168/SRR33243168_2.fastq.gz -o SRR33243168_RNA_seq_of_iPSC_derived_microglia_2.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/069/SRR33243169/SRR33243169_1.fastq.gz -o SRR33243169_RNA_seq_of_iPSC_derived_microglia_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/069/SRR33243169/SRR33243169_2.fastq.gz -o SRR33243169_RNA_seq_of_iPSC_derived_microglia_2.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/070/SRR33243170/SRR33243170_1.fastq.gz -o SRR33243170_RNA_seq_of_iPSC_derived_microglia_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/070/SRR33243170/SRR33243170_2.fastq.gz -o SRR33243170_RNA_seq_of_iPSC_derived_microglia_2.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/071/SRR33243171/SRR33243171_1.fastq.gz -o SRR33243171_RNA_seq_of_iPSC_derived_microglia_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/071/SRR33243171/SRR33243171_2.fastq.gz -o SRR33243171_RNA_seq_of_iPSC_derived_microglia_2.fastq.gz

cd ..
#fastQC and MultiQC

fastqc Dataset/*.fastq.gz -o qc/ 
multiqc qc/ -o qc/ #for better visualisation

#Trimming using fastp

for r1 in Dataset/*_1.fastq.gz; do
    r2=${r1/_1.fastq.gz/_2.fastq.gz}
    base=$(basename "$r1" _1.fastq.gz)

    fastp \
      -i "$r1" -I "$r2" \
      -o "trim/${base}_1.trim.fastq.gz" \
      -O "trim/${base}_2.trim.fastq.gz" \
      -h "trim/${base}_report.html" \
      -j "trim/${base}_report.json" \
      --detect_adapter_for_pe
done

multiqc trim/ -o qc/ #for summary and better visualisation

#Now we need the genome reference : Human genome

wget -O reference/homo_sapiens.fa.gz ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz #fasta format
gunzip reference/homo_sapiens.fa.gz #unzipping

wget -O reference/homo_sapiens.gtf.gz ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz #GTF annotation
gunzip reference/homo_sapiens.gtf.gz  #unzipping

#Generate STAR index

STAR --runMode genomeGenerate \
     --genomeDir genomeIndex \
     --genomeFastaFiles reference/homo_sapiens.fa \
     --sjdbGTFfile reference/homo_sapiens.gtf \
     --sjdbOverhang 149 \
     --runThreadN 8

#Begin mapping with STAR

for r1 in trim/*_1.trim.fastq.gz; do
    r2=${r1/_1.trim.fastq.gz/_2.trim.fastq.gz}
    sample=$(basename "$r1" _1.trim.fastq.gz)

    STAR --runThreadN 8 \
         --genomeDir genomeIndex \
         --readFilesIn "$r1" "$r2" \
         --readFilesCommand zcat \
         --outFileNamePrefix mapped/${sample}_ \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes All
done

#Now for indexing BAM For IGV

cp mapped/*_Aligned.sortedByCoord.out.bam IGV/

for bam in IGV/*_Aligned.sortedByCoord.out.bam; do
    samtools index "$bam"
done

#We can then visualise with IGV : https://igv.org/app/

#Counting genes with featureCounts
featureCounts -T 8 -p -B -C \
  -a reference/homo_sapiens.gtf \
  -o counts/gene_counts.txt \
  mapped/*_Aligned.sortedByCoord.out.bam

#summary
column -t counts/gene_counts.txt.summary | less -S

#Next step is made with R after exporting the counts.txt file to our computer 
