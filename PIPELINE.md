# Pipeline: *Saccharomyces cerevisiae* Differential Gene Expression Analysis

## 0. Environment Setup and Data Aquisition
The following Conda environments were created and used for the analysis:
- **sra-tools**: for obtaining SRA list from project number, and converting NCBI sra files to fastq files with sra-tools
- **fastqc**: for quality checking raw and filtered reads with FastQC 
- **multiqc**: for generating combined reports with read statistics for all samples before and after quality control with MultiQC
- **fastp**: for filtering reads with low quality with fastp
- **salmon**: to quantify gene expression with Salmon

```bash
# Create metadata.csv - manually
nano ./data/metadata.csv
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

# Download reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/046/745/GCA_003046745.1_ASM304674v1/GCA_003046745.1_ASM304674v1_genomic.fna.gz -O ./data/references/scerev_ref_genome.fna.gz

# Download reference transcriptome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/046/745/GCA_003046745.1_ASM304674v1/GCA_003046745.1_ASM304674v1_rna_from_genomic.fna.gz -O ./data/references/scerev_ref_transcripts.fna.gz

# Install entrez-direct to get SRA list
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
export PATH=${HOME}/edirect:${PATH}

# Get SRR list for project PRJNA592304 from NCBI
```bash
esearch -db sra -query "PRJNA592304" | \
efetch -format runinfo | \
cut -d "," -f1 | \
grep SRR > ./data/metadata/srr_list.txt

# Download fastq files corresponding to SRA numbers 
conda activate sra-tools
prefetch --option-file ./data/metadata/srr_list.txt --output-directory data/raw/

# Convert .sra format to .fastq and compress
for srr in $(cat ./data/metadata/srr_list.txt); do
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
fastqc ./data/raw/*.fastq.gz -o ./results/qc_results/raw_reads

conda activate multiqc
multiqc ./results/qc_results/raw_reads -o ./results/qc_results/raw_reads

# Filter raw reads based on quality 30 using Fastp
conda create -n fastp
conda activate fastp
conda install -c bioconda fastp

# Install script to run filtering for all files at once 
wget https://raw.githubusercontent.com/OpenGene/fastp/master/parallel.py -O ./tools/

python3 ./tools/parallel.py -i ./data/raw/ -o ./data/filtered -r ./results/qc_results/filtered_reads/fastp/ --args="--disable_adapter_trimming --disable_length_filtering --qualified_quality_phred 30 --unqualified_percent_limit 20"

# Quality check after filtering with MultiQC 
multiqc ./ -o ./results/qc_results/filtered_reads
```

## 2. Quantify Gene Expression
```bash
conda create -n salmon
conda activate salmon
conda install salmon

# Build salmon index 
## Prepare metadata
### Extract genome targets (decoys)
grep "^>" <(gunzip -c ./data/references/scerev_ref_genome.fna.gz) | cut -d " " -f 1 > ./data/salmon_index/decoys.txt
sed -i.bak -e 's/>//g' ./data/salmon_index/decoys.txt

### Concatenate transcriptome and genome reference file
cat ./data/references/scerev_ref_transcripts.fna.gz ./data/references/scerev_ref_genome.fna.gz > ./data/salmon_index/scerev_gentrome.fa.gz

## Index 
cd data/salmon_index
salmon index -t scerev_gentrome.fa.gz -d decoys.txt -p 12 -i index

# Quantify gene expression for each sample with Salmon
brew install parallel






```
## 3. Statistical Analysis
## 4. Visualize Differentially Expressed Genes
## 5. Functional Enrichment Analysis
