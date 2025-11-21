# inst/scripts/generate_extdata.sh
# This pseudo-code describes how to generate the full dataset (only a subset is included in inst/extdata).
# Sources and licenses:
# - Ly 2024 HeLa cell cycle dataset: NCBI SRA (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA957808)
# - GENCODE v47 annotation and transcript FASTA
# - hg38 reference genome: UCSC Genome Browser
# Tools: SRA Toolkit, STAR, cutadapt, featureCounts, RIBOSS, DOTSeq orf_to_gtf.py script
# Assume these tools are installed.


# Step 1. Align Ribo-seq and RNA-seq reads

# 1. Clone DOTSeq repository
git clone https://github.com/compgenom/DOTSeq.git

# Assume FASTQ files are downloaded via SRA Toolkit and stored in:
# DOTSeq/inst/extdata/ly_2024
# STAR index will be generated in:
# DOTSeq/inst/extdata/hg38_star_index

# 2. Download GENCODE annotation and transcript FASTA files
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz \
  -O DOTSeq/inst/extdata/gencode.v47.annotation.gtf.gz

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.pc_transcripts.fa.gz \
  -O DOTSeq/inst/extdata/gencode.v47.pc_transcripts.fa.gz

# Assume reference genome FASTA is downloaded
# 3. Generate STAR index
STAR --runMode genomeGenerate \
  --runThreadN 32 \
  --genomeFastaFiles DOTSeq/inst/extdata/hg38.fa \
  --sjdbGTFfile DOTSeq/inst/extdata/gencode.v47.annotation.gtf.gz \
  --genomeDir DOTSeq/inst/extdata/hg38_star_index

# 4. Trim and align reads
for i in DOTSeq/inst/extdata/ly_2024/*/*/*.fastq.gz; do
  cutadapt -j 16 -m 15 -u 8 -e 0.1 --match-read-wildcards \
    -a TCGTATGCCGTCTTCTGCTTG -O 1 \
    -o $(dirname "$i")/$(basename "$i" .fastq.gz).trimmed.fasta.gz "$i"

  STAR --runMode alignReads \
    --runThreadN 32 \
    --outFilterType BySJout --outFilterMismatchNmax 2 \
    --genomeDir DOTSeq/inst/extdata/hg38_star_index \
    --readFilesIn $(dirname "$i")/$(basename "$i" .fastq.gz).trimmed.fasta.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix $(dirname "$i")/$(basename "$i" .fastq.gz) \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFilterMultimapNmax 1 --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd --outSAMattributes All
done


# Step 2. Prepare ORF-level GTF and BED files

# Recommended: See the vignette for a step-by-step guide to preparing ORF-level annotation natively in R.
# Input: GTF files from GENCODE, Ensembl, or Araport.

# Alternative: The flattened annotation files in inst/extdata were generated using the RIBOSS engine (https://github.com/lcscs12345/riboss).

# 1. Install Miniforge3 on Linux.
wget https://github.com/conda-forge/miniforge/releases/download/24.7.1-2/Miniforge3-24.7.1-2-Linux-x86_64.sh
bash Miniforge3-24.7.1-2-Linux-x86_64.sh -b -p $HOME/miniforge3
eval "$(/$HOME/miniforge3/bin/conda shell.bash hook)" # your terminal prompt will show (base) bash-5.1$

# 2. Create a conda environment and install RIBOSS
git clone https://github.com/lcscs12345/riboss.git
cd riboss
conda config --set channel_priority flexible # required if channel_priority was set to strict. See https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-channels.html
conda env create -f environment.yml

conda activate riboss # your terminal prompt will show (riboss) bash-5.1$
DIRNAME=`which python | xargs dirname`
cp bin/riboprof $DIRNAME
chmod +x $DIRNAME/riboprof
which python | awk 'sub(/python/,"pip3") {print $1, "install -e ."}' | sh # editable mode

# 3. Activate the conda environment for next time
eval "$(/$HOME/miniforge3/bin/conda shell.bash hook)"
conda activate riboss

# 4. Generate ORF-level GTF using DOTSeq's Python script
python DOTSeq/inst/scripts/orf_to_gtf.py \
  --gtf DOTSeq/inst/extdata/gencode.v47.annotation.gtf.gz \
  --transcripts DOTSeq/inst/extdata/gencode.v47.pc_transcripts.fa.gz \
  --output DOTSeq/inst/extdata/gencode.v47.annotation.orf_flattened


# Step 3: Count reads

featureCounts -f -O -s 0 -T 16 \
-F GTF -a DOTSeq/inst/extdata/gencode.v47.annotation.orf_flattened.gtf \
-o DOTSeq/inst/extdata/featureCounts.dotseq.out DOTSeq/inst/extdata/ly_2024/*/*/*Aligned.sortedByCoord.out.bam

