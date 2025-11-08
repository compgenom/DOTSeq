# DOTSeq

`DOTSeq` is an R package for identifying 
**differentially translated open reading frames (ORFs)** from ribosome 
profiling (Ribo-seq) and matched RNA-seq datasets. Unlike traditional 
gene-level approaches, `DOTSeq` performs analysis at the **ORF level**, 
enabling detection of:

- **Differential ORF Usage (DOU)** — changes in ORF usage within the same 
gene across conditions.
- **Differential Translation Efficiency (DTE)** — changes in ribosome 
loading relative to RNA level across conditions.

`DOTSeq` models Ribo-seq and RNA-seq read counts using a 
**beta-binomial generalised linear model (GLM)** implemented via 
[`glmmTMB`](https://CRAN.R-project.org/package=glmmTMB). 
It supports experimental designs with multiple conditions, and uses an 
interaction term (`condition:strategy`) to isolate translation-specific 
effects.

Post hoc contrasts are computed using 
[`emmeans`](https://CRAN.R-project.org/package=emmeans), 
and empirical Bayes shrinkage is applied via 
[`ashr`](https://CRAN.R-project.org/package=ashr).


## DEPENDENCIES
* R (>= 4.5.0)
* biomaRt (>=2.65.0)
* SummarizedExperiment (>= 1.39.1)
* Bioc.gff (>= 0.99.17)
* DESeq2 (>=1.49.4)
* GenomicRanges (>=1.61.5)
* IRanges (>=2.43.5)
* S4Vectors (0.47.4)
* ashr (>=2.2-63)
* DHARMa (>=0.4.7)
* emmeans (>=1.11.2-8)
* glmmTMB (>=1.1.12)
* eulerr (>=7.0.4)
* pbapply (>=1.7-4)

## INSTALLATION
Please ensure the dependencies listed above are installed using the 
following steps before installing `DOTSeq`:
```r
# Create a directory for R packages if not already
package_dir <- file.path(Sys.getenv("HOME"), "R/4.5")
dir.create(package_dir, showWarnings = TRUE, recursive = TRUE)
.libPaths(c(package_dir, .libPaths()))

# Install BiocManager if not already available
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = package_dir)

# Initialise usage of Bioconductor 3.22
options(repos = BiocManager::repositories(version = "3.22"))

# Install DOTSeq and required packages with automatic update confirmation
BiocManager::install("compgenom/DOTSeq", lib = package_dir, ask = FALSE)
```

## DOCUMENTATION
TL;DR: To understand how to use `DOTSeq` without going through the 
preprocessing steps, please refer to the 
[vignettes](https://github.com/compgenom/DOTSeq/tree/main/vignettes).

### Preprocessing Steps Required Before Running `DOTSeq`:

#### Step 1. Align Ribo-seq and RNA-seq reads
We use a publicly available HeLa cell cycle dataset from 
[Ly 2024](https://pubmed.ncbi.nlm.nih.gov/39443796/). 

```shell
# Clone DOTSeq repository
git clone https://github.com/compgenom/DOTSeq.git

# Assume FASTQ files are downloaded via SRA Toolkit and stored in:
# DOTSeq/inst/extdata/ly_2024
# STAR index will be generated in:
# DOTSeq/inst/extdata/hg38_star_index

# Download GENCODE annotation and transcript FASTA files
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz \
  -O DOTSeq/inst/extdata/gencode.v47.annotation.gtf.gz

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.pc_transcripts.fa.gz \
  -O DOTSeq/inst/extdata/gencode.v47.pc_transcripts.fa.gz

# Assume reference genome FASTA is downloaded
# Generate STAR index
STAR --runMode genomeGenerate \
  --runThreadN 32 \
  --genomeFastaFiles DOTSeq/inst/extdata/hg38.fa \
  --sjdbGTFfile DOTSeq/inst/extdata/gencode.v47.annotation.gtf.gz \
  --genomeDir DOTSeq/inst/extdata/hg38_star_index

# Trim and align reads
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
```

#### Step 2. Prepare ORF-level GTF and BED files
Step-by-step on how to prepare ORF-level annotation is available in the 
[vignettes](https://github.com/compgenom/DOTSeq/tree/main/vignettes). 
GTF files from GENCODE, Ensembl, or Araport should be used as input.

Alternatively, DOTSeq accept flattened annotation files generated using 
the [`RIBOSS`](https://github.com/lcscs12345/riboss) engine. 

```shell
# Generate ORF-level GTF using DOTSeq's Python script
python DOTSeq/inst/python_scripts/orf_to_gtf.py \
  --gtf DOTSeq/inst/extdata/gencode.v47.annotation.gtf.gz \
  --transcripts DOTSeq/inst/extdata/gencode.v47.pc_transcripts.fa.gz \
  --output DOTSeq/inst/extdata/dotseq
```

#### Step 3: Count reads
```shell
featureCounts -f -O -s 0 -T 16 \
-F GTF -a DOTSeq/inst/extdata/dotseq.gtf \
-o DOTSeq/inst/extdata/featureCounts.dotseq.out DOTSeq/inst/extdata/ly_2024/*/*/*Aligned.sortedByCoord.out.bam
```

#### Step 4: Run `DOTSeq`
Follow the analysis workflow as demonstrated in the 
[vignettes](https://github.com/compgenom/DOTSeq/tree/main/vignettes).

## CONTRIBUTING
We welcome contributions from the community! Whether it's fixing bugs, 
improving documentation, or suggesting new features, your input is valuable.

By participating in this project, you agree to abide by the terms outlined in 
this [Contributor Code of Conduct](https://www.contributor-covenant.org/version/1/0/0/code-of-conduct/).

To get started:
- Fork the repository
- Create a new branch for your feature or fix
- Submit a pull request with a clear description of your changes

If you have questions or ideas, feel free to open an issue or start 
a discussion.

## CONTACTS AND BUG REPORTS
- Chun Shen Lim: 
chunshen [dot] lim [at] otago [dot] ac [dot] nz
- Gabrielle Chieng: 
gabrielle [dot] chieng [at] postgrad [dot] otago [dot] ac [dot] nz
