# Pipeline: *Saccharomyces cerevisiae* Differential Gene Expression Analysis

## 0. Environment Setup and Data Aquisition
```bash
# Create metadata.csv - manually
nano metadata.csv
stage,sample_id,sra_accession
early,IL20,SRR10551665
early,IL21,SRR10551664
early,IL22,SRR10551663
thin,IL23,SRR10551662
thin,IL24,SRR10551661
thin,IL25,SRR10551660
mature,IL29,SRR10551659
mature,IL30,SRR10551658
mature,IL31,SRR10551657

# Install entrez-direct to get SRA list
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
export PATH=${HOME}/edirect:${PATH}

# Get SRR list for project PRJNA592304 from NCBI
```bash
esearch -db sra -query "PRJNA592304" | \
efetch -format runinfo | \
cut -d "," -f1 | \
grep SRR > srr_list.txt

# Download fastq files corresponding to SRA numbers 
conda activate sra-tools
prefetch --option-file srr_list.txt --output-directory data/raw/

# Convert .sra format to .fastq and compress
for srr in $(cat srr_list.txt); do
echo "Processing $srr:"
fasterq-dump --split-3 --threads 10 --outdir data/raw/ data/raw/$srr/$srr.sra
gzip data/raw/${srr}*.fastq
rm -rf data/raw/$srr
done
```
## 1. Quality Control
```bash
# Quality check raw reads with FastQC and generate consolidated report with MultiQC
conda activate fastqc 
fastqc ./data/raw/*.fastq.gz -o ./qc_results/raw_reads

conda activate multiqc
multiqc ./qc_results/raw_reads -o ./qc_results/raw_reads
```
## 2. Read Alignment
## 3. Count Gene Expression
## 4. Statistical Analysis
## 5. Visualize Differentially Expressed Genes
## 6. Functional Enrichment Analysis
