#' Check if a BAM file is paired-end
#'
#' Quickly checks whether a BAM file contains paired-end reads by 
#' inspecting the first 1000 flags. Assumes that all reads in the 
#' file are consistently either paired-end or single-end.
#'
#' @param bam_file Character. Path to a BAM file.
#'
#' @return Logical. `TRUE` if any of the first 1000 reads have the 
#' 'paired' flag set (0x1), `FALSE` otherwise.
#' 
#' @importFrom Rsamtools ScanBamParam scanBam BamFile
#' 
#' @keywords internal
#' 
is_paired_end <- function(bam_file) {
    param <- ScanBamParam(what = "flag")
    # bam <- scanBam(BamFile(bam_file), param = param, yieldSize = 1000)[[1]]
    bf <- BamFile(bam_file, yieldSize = 1000)
    open(bf)
    bam <- scanBam(bf, param = param)[[1]]
    close(bf)
    
    # Check if any reads are paired
    libtype <- any(bitwAnd(bam$flag, 0x1) != 0)
    return(libtype)
}


#' Group BAM files by read type (paired-end or single-end)
#'
#' Splits a list of BAM file paths into paired-end and single-end 
#' groups based on the read flags.
#'
#' @param bam_files Character vector. Paths to BAM files.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{paired}{Character vector of paired-end BAM file paths}
#'   \item{single}{Character vector of single-end BAM file paths}
#' }
#'
#' @keywords internal
#' 
group_bam_files <- function(bam_files) {
    paired <- c()
    single <- c()
    for (bam in bam_files) {
        if (is_paired_end(bam)) {
            paired <- c(paired, bam)
        } else {
            single <- c(single, bam)
        }
    }
    return(list(paired = paired, single = single))
}


#' Filter BAM files to retain only reads overlapping exonic regions
#'
#' @description
#' This function filters BAM files to retain reads overlapping exonic 
#' regions #' defined in a TxDb object stored in the metadata of a 
#' GRanges object. Optionally, it restricts to coding genes only. 
#' Filtered BAMs are sorted and saved with `.exonic.sorted.bam` suffixes.
#'
#' @param gr A \code{GRanges} object with a \code{TxDb} object stored 
#' in its metadata slot under \code{metadata(gr)$txdb}.
#' @param bam_files A character vector of paths to BAM files to be 
#' filtered.
#' @param coding_genes_only Logical; if \code{TRUE}, restrict filtering 
#' to coding genes only (i.e., genes with CDS).
#' @param verbose Logical; if \code{TRUE}, print progress and runtime 
#' messages.
#'
#' @return This function does not return a value. It creates filtered 
#' and sorted BAM files for each input BAM.
#' 
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom AnnotationDbi loadDb
#' @importFrom GenomicFeatures cdsBy exonsBy transcripts
#' @importFrom BiocGenerics unlist
#' @importFrom Rsamtools scanBamHeader BamFile filterBam ScanBamParam
#' @importFrom Rsamtools indexBam sortBam
#' @importFrom GenomeInfoDb keepSeqlevels
#' 
#' @export
#' @examplesIf requireNamespace("TxDb.Dmelanogaster.UCSC.dm3.ensGene", quietly = TRUE) && requireNamespace("pasillaBamSubset", quietly = TRUE)
#' library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
#' library(pasillaBamSubset)
#' library(GenomeInfoDb)
#' 
#' txdb_chr4 <- keepSeqlevels(
#'     TxDb.Dmelanogaster.UCSC.dm3.ensGene, 
#'     "chr4", 
#'     pruning.mode = "coarse"
#' )
#' gr <- GRanges(seqnames = "chr4", ranges = IRanges(start = 233, end = 2300))
#' metadata(gr)$txdb <- txdb_chr4
#' 
#' getExonicReads(gr, bam_files = c(untreated1_chr4()))
#' 
#' # Get the directory of the BAM file and clean up output files
#' bam_dir <- dirname(untreated1_chr4())
#' output_files <- list.files(
#'     path = bam_dir,
#'     pattern = "*exonic.*",
#'     full.names = TRUE
#' )
#' file.remove(output_files)
#' 
getExonicReads <- function(gr, bam_files, coding_genes_only = TRUE, verbose = TRUE) {
    
    if (!file.exists(metadata(gr)$txdb)) {
        stop("The TxDb object has been removed. Please rerun getORFs().")
    }
    
    # Filter annotations to coding genes
    annotation <- loadDb(metadata(gr)$txdb)
    exons_by_genes <- exonsBy(annotation, by = "gene")
    
    if (coding_genes_only) {
        cds_by_genes <- cdsBy(annotation, by = "gene")
        common_genes <- intersect(names(exons_by_genes), names(cds_by_genes))
        exons_by_genes <- exons_by_genes[common_genes]
    }
    
    exons_by_genes <- unlist(exons_by_genes)
    
    # Get BAM header seqlevels
    for (bam in bam_files) {
        
        if (verbose) {
            start_filterbam <- Sys.time()
            message("starting filtering ", bam)
        }
        
        bamfile <- BamFile(bam)
        bam_header <- scanBamHeader(BamFile(bam))
        bam_seqlevels <- names(bam_header[[1]])
        
        # Get seqlevels from annotation
        annotation_seqlevels <- seqlevels(exons_by_genes)
        
        # Find common seqlevels
        common_seqlevels <- intersect(annotation_seqlevels, bam_seqlevels)
        
        # Keep only common seqlevels in annotation
        exons_for_bam <- keepSeqlevels(exons_by_genes, common_seqlevels, pruning.mode = "coarse")
        
        # Filter reads
        if (!file.exists(paste0(bam, ".bai"))) {
            indexBam(bam)
        }
        
        filtered_bam <- paste0(tools::file_path_sans_ext(bam), ".exonic.bam")
        filterBam(file = bam,
                  destination = filtered_bam,
                  indexDestination = FALSE,
                  param = ScanBamParam(which = exons_for_bam))
        sorted_bam <- paste0(tools::file_path_sans_ext(filtered_bam), ".sorted")
        sortBam(file = filtered_bam, destination = sorted_bam)
        invisible(file.remove(filtered_bam))
        
        if (verbose) {
            end_filterbam <- Sys.time()
            elapsed_filterbam <- runtime(end_filterbam, start_filterbam)
            
            if (!is.null(elapsed_filterbam$mins)) {
                message(sprintf("BAM filtering runtime: %d mins %.3f secs", elapsed_filterbam$mins, elapsed_filterbam$secs))
            } else {
                message(sprintf("BAM filtering runtime: %.3f secs", elapsed_filterbam$secs))
            }
        }
    }
}


#' Count reads from BAM files over genomic features
#'
#' Uses \code{\link[GenomicAlignments]{summarizeOverlaps}} to count reads 
#' overlapping genomic ranges, handling both single-end and paired-end 
#' BAM files.
#'
#' @param gr A \code{GRanges} object representing genomic features (e.g., 
#' ORFs).
#' @param bam_files Character vector. Paths to BAM files.
#' @param ignore.strand A named list with logical values for `single` 
#' and `paired` indicating whether to ignore strand information.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A matrix of read counts with features as rows and samples as 
#' columns.
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom Rsamtools BamFileList
#' 
#' @export
#' @examplesIf requireNamespace("pasillaBamSubset", quietly = TRUE) && requireNamespace("GenomicRanges", quietly = TRUE)
#' library(GenomicRanges)
#' library(pasillaBamSubset)
#' 
#' bam_list <- c(untreated1_chr4(), untreated3_chr4())
#' 
#' gr <- GRanges(seqnames = "chr4", 
#'     ranges = IRanges(start = 233, end = 2300))
#'     
#' countReads(gr = gr, bam_files = bam_list)
#' 
countReads <- function(
        gr, 
        bam_files, 
        ignore.strand = list(single = TRUE, paired = FALSE),
        verbose = TRUE
) {
    
    start_counts <- Sys.time()
    
    split_bams <- group_bam_files(bam_files)
    
    se_list <- list()
    single_length <- length(split_bams$single)
    paired_length <- length(split_bams$paired)
    
    if (single_length > 0) {
        if (verbose) {
            message("assigning reads from ", single_length, " single-end BAM files")
        }
        se_list$single <- summarizeOverlaps(
            features = gr,
            reads = BamFileList(split_bams$single),
            singleEnd = TRUE,
            fragments = FALSE,
            ignore.strand = ignore.strand$single
        )
    }
    
    if (paired_length > 0) {
        if (verbose) {
            message("assigning reads for ", paired_length, " paired-end BAM files")
        }
        
        se_list$paired <- summarizeOverlaps(
            features = gr,
            reads = BamFileList(split_bams$paired),
            singleEnd = FALSE,
            fragments = TRUE,
            ignore.strand = ignore.strand$paired
        )
    }
    
    # Combine results if both types are present
    if (length(se_list) == 2) {
        count_table <- cbind(assay(se_list$paired), assay(se_list$single))
    } else {
        count_table <- assay(se_list[[1]])
    }
    
    if (verbose) {
        end_counts <- Sys.time()
        elapsed_counts <- runtime(end_counts, start_counts)
        
        if (!is.null(elapsed_counts$mins)) {
            message(sprintf("read counting runtime: %d mins %.3f secs", elapsed_counts$mins, elapsed_counts$secs))
        } else {
            message(sprintf("read counting runtime: %.3f secs", elapsed_counts$secs))
        }
    }
    
    return(as.data.frame(count_table))
}


#' Assign strategy levels for RNA-seq and Ribo-seq
#'
#' This function identifies and assigns the reference and target levels 
#' for RNA-seq and Ribo-seq strategies from a given input data frame. 
#' It uses a flexible regular expression to match common representations 
#' of RNA-seq (e.g., "rna", "RNA", "RNA-seq", or "0") and assigns the 
#' remaining unique value(s) as Ribo-seq.
#'
#' @param input_df A data frame containing a column that encodes 
#' strategy labels.
#' @param strategy_col A string specifying the name of the column in 
#' `input_df` that contains strategy labels. Default is "strategy".
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{ribo_level}{The detected level corresponding to Ribo-seq.}
#'   \item{rna_level}{The detected level corresponding to RNA-seq.}
#' }
#'
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' df <- data.frame(strategy = c("RNA", "Ribo"))
#' assign_strategy_levels(df)
#' }
#' 
assign_strategy_levels <- function(input_df, strategy_col = "strategy") {
    strategy_vals <- as.character(input_df[[strategy_col]])
    
    # Define a flexible regex pattern: 
    # matches "rna", "RNA", "RNA-seq", and "0"
    rna_pattern <- "(?i)rna|^0$"
    
    # Find matching values to RNA-seq
    rna_level <- unique(strategy_vals[grepl(rna_pattern, strategy_vals)])
    
    # Handle multiple matches
    if (length(rna_level) > 1) {
        warning(
            "Multiple RNA-seq levels found. ", 
            "Using the first match as reference."
        )
        rna_level <- rna_level[1]
    } else if (length(rna_level) == 0) {
        stop(
            "No RNA-seq level found. Please use rna, RNA, RNA-seq, ", 
            "or 0 to represent RNA-seq."
        )
    }
    
    # Find the Ribo-seq level
    ribo_level <- setdiff(unique(strategy_vals), rna_level)
    
    if (length(ribo_level) > 1) {
        warning(
            "Multiple Ribo-seq levels found. Use ", 
            ribo_level[1], 
            " as target."
        )
        ribo_level <- ribo_level[1]
    } else if (length(ribo_level) == 0) {
        stop(
            "No Ribo-seq level found. ", 
            "Please check your data."
        )
    }
    
    return(list(ribo_level = ribo_level, rna_level = rna_level))
}


#' Annotate ORF Types from BED File
#'
#' @description
#' This function reads a BED file and converts it into a \code{GRanges} object.
#' It then finds overlaps with a reference set of ORFs (provided as a
#' \code{GRanges} object, typically derived from a GFF or DEXSeqDataSet), and
#' annotates each overlapping ORF with a label (e.g., score or type) from the
#' BED file.
#'
#' @param bed Character string. Path to a BED file (tab-delimited, no header)
#' with six columns:
#' \describe{
#'     \item{chr}{Chromosome name}
#'     \item{start}{Start coordinate (0-based)}
#'     \item{end}{End coordinate (exclusive)}
#'     \item{name}{Feature name}
#'     \item{score}{Annotation label to assign (e.g., ORF type)}
#'     \item{strand}{Strand information ("+" or "-")}
#' }
#'
#' @param gff_granges A \code{GRanges} object containing ORF coordinates to be
#' annotated.
#'
#' @return A \code{GRanges} object with an added metadata column 
#' \code{orf_type}, containing the BED score for overlapping ORFs. 
#' Non-overlapping ORFs will have \code{NA}.
#'
#' @importFrom utils read.table
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits mcols mcols<-
#' @importFrom IRanges IRanges
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' # Annotate a GRanges object using a BED file's score column
#' gff_granges <- parseBed("example.bed", granges)
#' head(orfs)
#' }
#' 
annotate_orf_type <- function(bed, gff_granges) {
    bedDf <- read.table(
        bed, 
        header = FALSE, 
        sep = "\t", 
        stringsAsFactors = FALSE
    )
    colnames(bedDf) <- c("chr", "start", "end", "name", "score", "strand")
    
    mORFs <- bedDf[bedDf$score == "mORF", ]
    
    # Identify duplicated names in both directions
    dups_name <- duplicated(mORFs$name) | 
        duplicated(mORFs$name, fromLast = TRUE)
    
    # Keep only rows with unique names
    bedDf <- bedDf[!dups_name, ]
    
    # Convert to GRanges
    bed_granges <- GRanges(
        seqnames = bedDf$chr,
        ranges = IRanges(start = bedDf$start + 1, end = bedDf$end),
        strand = bedDf$strand,
        name = bedDf$name,
        label = bedDf$score
    )
    
    # Initialize labels column
    mcols(gff_granges)$orf_type <- NA
    
    # Find overlaps
    hits <- findOverlaps(gff_granges, bed_granges)
    
    # Assign BED scores to matched ORFs
    mcols(gff_granges)[queryHits(hits), "orf_type"] <- mcols(bed_granges)[subjectHits(hits), "label"]
    
    # Return the updated GRanges object
    return(gff_granges)
}


#' Parse and validate sample metadata for DOTSeq analysis
#'
#' @description
#' This internal helper function reads and validates the sample metadata
#' provided either as a file path or a data frame. It ensures that the
#' required columns—\code{run}, \code{strategy}, \code{condition}, and
#' \code{replicate}—are present (case-insensitive match), and normalizes
#' column names for consistency. If a file path is provided, the function
#' reads the file and checks the header before loading the full table.
#' If a data frame is provided, it performs similar validation directly.
#'
#' @param condition_table A character string specifying the path to a
#' tab-delimited metadata file, or a data frame containing sample metadata.
#' The metadata must include the following columns (case-insensitive):
#' \code{run}, \code{strategy}, \code{condition}, and \code{replicate}.
#'
#' @return A data frame with normalized column names and row names set to
#' the \code{run} column. Used internally to construct the
#' \code{\link{DOTSeqDataSets-class}} object.
#'
#' @keywords internal
#' 
parse_condition_table <- function(condition_table) {
    
    # Define expected column names
    condCols <- c("run", "strategy", "condition", "replicate")
    
    # Helper to normalize column names
    normaliseNames <- function(x) tolower(trimws(x))
    
    if (is.character(condition_table) && file.exists(condition_table)) {
        # Read header line and normalize
        firstLine <- readLines(condition_table, n = 1)
        condHeader <- normaliseNames(strsplit(firstLine, "\t|,|\\s+")[[1]])
        
        missingCols <- setdiff(normaliseNames(condCols), condHeader)
        
        if (length(missingCols) == 0) {
            cond <- read.table(
                condition_table, 
                header = TRUE, 
                comment.char = "#", 
                stringsAsFactors = FALSE
            )
            names(cond) <- normaliseNames(names(cond))
        } else {
            stop(
                paste(
                    "Missing expected columns (case-insensitive match):", 
                    paste(missingCols, collapse = ", ")
                )
            )
        }
    } else if (is.data.frame(condition_table)) {
        condHeader <- normaliseNames(names(condition_table))
        missingCols <- setdiff(normaliseNames(condCols), condHeader)
        
        if (length(missingCols) == 0) {
            cond <- condition_table
            names(cond) <- normaliseNames(names(cond))
        } else {
            stop(
                paste(
                    "Data frame is missing expected columns", 
                    "(case-insensitive match):", 
                    paste(missingCols, collapse = ", ")
                )
            )
        }
    } else {
        stop("`condition_table` must be either a valid file path or a data frame.")
    }
    
    # Rename cond rownames
    rownames(cond) <- cond$run
    
    return(cond)
}


#' Match and align sample metadata with count table
#'
#' @description
#' This internal function matches sample identifiers between the count table
#' and the condition table, ensuring consistency in sample naming and 
#' ordering. It verifies that sufficient replicates exist, renames columns 
#' to include metadata, and prepares the count and condition tables for 
#' downstream differential translation analysis. It also assigns binary 
#' strategy labels (e.g., RNA = 0, Ribo = 1), sets factor levels, and 
#' optionally relevels the condition factor to a specified baseline.
#'
#' @param cnt A data frame containing count data. Columns must include sample
#' names matching the \code{run} column in \code{cond}. If \code{num_feat_cols}
#' is 6, the first six columns must be: \code{Geneid}, \code{Chr}, \code{Start},
#' \code{End}, \code{Strand}, \code{Length}.
#'
#' @param cond A data frame containing sample metadata. Must include columns:
#' \code{run}, \code{strategy}, \code{condition}, and \code{replicate}.
#'
#' @param num_feat_cols Integer specifying the number of feature columns in
#' \code{cnt} before sample columns. Typically 6 for featureCounts output,
#' or 0 for raw count matrices.
#'
#' @param baseline Optional character string specifying the reference level
#' for the \code{condition} factor.
#'
#' @return A list with two elements:
#' \describe{
#'     \item{\code{cnt}}{A reordered count table with renamed sample columns.}
#'     \item{\code{cond}}{A metadata data frame with binary strategy labels 
#'     and factor columns, ready for modeling.}
#' }
#'
#' @keywords internal
#' 
#' @importFrom stats relevel
#' 
match_runs <- function(cnt, cond, num_feat_cols = 0, baseline = NULL, verbose = TRUE) {
    # Find common identifiers
    common <- intersect(rownames(cond), names(cnt))
    
    if (length(common) < 4) {
        stop(
            "Run ID mismatch between count_table and count_table. ", 
            "Please check the column and row names. Only ", 
            length(common), " runs are matched."
        )
    }
    
    cond <- cond[common, , drop = FALSE]
    cond <- cond[order(cond$strategy, cond$replicate), ]
    
    cond_run <- cond$run
    cnt_run <- colnames(cnt[c((num_feat_cols + 1):ncol(cnt))])
    missing_samples <- setdiff(cnt_run, cond_run)
    
    if (length(cnt_run) < length(cond_run)) {
        warning(
            paste(missing_samples, collapse = ", "), 
            " are missing in count table."
        )
    } else if (length(cnt_run) > length(cond_run)) {
        warning(
            paste(missing_samples, collapse = ", "),
            " are missing in condition table."
        )
    }
    
    if (num_feat_cols == 6) {
        # Combine with metadata columns
        cntCols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
        cnt <- cnt[, c(cntCols, rownames(cond))]
    } else if (num_feat_cols == 0) {
        cnt <- cnt[, rownames(cond)]
    } else {
        stop(
            "count_table format currently not supported. ", 
            "Please use DOTSeqDataSetsFromSE or DOTSeqDataSetsFromFeaturCounts."
        )
    }
    
    # Rename header
    rownames(cond) <- apply(cond[, c("run", "condition", "replicate", "strategy")], 1, function(x) {
        paste0(
            trimws(x[1]), 
            ".", 
            trimws(x[2]), 
            ".", 
            trimws(x[3]), 
            ".", 
            trimws(x[4])
        )
    })
    # Find the indices of the columns that match cond$run
    libIndices <- match(cond$run, names(cnt))
    # Replace those column names with rownames(cond)
    if (ncol(cnt) == (length(rownames(cond)) + num_feat_cols)) {
        names(cnt)[libIndices] <- rownames(cond)
    } else {
        stop(
            "Number of samples in count_table and condition_table ", 
            "doesn't match."
        )
    }
    
    strategy_levels <- assign_strategy_levels(
        input_df = cond, 
        strategy_col = "strategy"
    )
    
    riboCond <- cond[cond$strategy == strategy_levels$ribo_level, ]
    riboCond$strategy <- 1
    rnaCond <- cond[cond$strategy == strategy_levels$rna_level, ]
    rnaCond$strategy <- 0
    
    # Combine and duplicate columns
    combined_cond <- rbind(rnaCond, riboCond)
    combined_cond <- combined_cond[, !(names(combined_cond) %in% "run")]
    
    # Convert to factors
    for (col in colnames(combined_cond)) {
        combined_cond[[col]] <- factor(combined_cond[[col]])
    }
    
    # Set baseline
    if (!is.null(baseline)) {
        if (isTRUE(verbose)) {
            message("setting condition", baseline, " as baseline")
        }
        combined_cond$condition <- relevel(
            combined_cond$condition, 
            ref = baseline
        )
    }
    
    if (num_feat_cols == 6) {
        # Combine with metadata columns
        cntCols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
        cnt <- cnt[, c(cntCols, rownames(combined_cond))]
    } else if (num_feat_cols == 0) {
        cnt <- cnt[, rownames(combined_cond)]
    }
    
    if (is.null(combined_cond) || nrow(combined_cond) < 4) {
        stop(
            "condition_table is malformed. ", 
            "It must have at least two replicates per conditions that ", 
            "match with run IDs in count_table."
        )
    }
    
    if (nrow(cnt)<2 && ncol(cnt)<4) {
        stop(
            "count_table. is malformed. ", 
            "It must have at least two replicates per conditions that ", 
            "match with run IDs in condition_table"
        )
    }
    
    return(list(cnt = cnt, cond = combined_cond))
}


#' Create DOU and DESeq2 datasets for differential translation analysis
#'
#' @description
#' This internal function constructs a \code{DOUData} object and a 
#' \code{DESeqDataSet} object from raw count data, sample metadata, and ORF 
#' annotations. It applies filtering based on count thresholds and singlet 
#' status, prepares metadata for modeling, and stores formulas and contrast 
#' specifications for downstream analysis.
#'
#' @param count_table A matrix or data frame of raw counts. Columns must 
#' match sample names in \code{condition_table}.
#'
#' @param condition_table A data frame containing sample metadata. Must 
#' include columns: \code{run}, \code{strategy}, \code{condition}, and 
#' \code{replicate}.
#'
#' @param annotation A \code{GRanges} object containing ORF annotations, 
#' typically parsed from a flattened GTF or BED file.
#'
#' @param reduced_formula A formula object specifying the reduced model 
#' used for estimating marginal means (e.g., \code{~ condition + strategy}).
#'
#' @param emm_specs A list of specifications for estimated marginal means 
#' contrasts, typically generated using \code{emmeans::contrast()}.
#'
#' @param deseq_formula A formula object specifying the design for DESeq2 
#' modeling (e.g., \code{~ condition * strategy}).
#'
#' @param min_count Integer specifying the minimum count threshold for 
#' filtering ORFs. Default is \code{1}.
#'
#' @param stringent Logical or \code{NULL}; determines the filtering strategy:
#' \describe{
#'     \item{\code{TRUE}}{
#'         Keep ORFs where all replicates in at least one condition pass 
#'         \code{min_count}.
#'     }
#'     \item{\code{FALSE}}{
#'         Keep ORFs where all replicates in at least one condition-strategy 
#'         group pass \code{min_count}.
#'     }
#'     \item{\code{NULL}}{
#'         Keep ORFs where total counts across all samples pass 
#'         \code{min_count}.
#'     }
#' }
#'
#' @param verbose Logical; if \code{TRUE}, prints progress and runtime 
#' messages. Default is \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'     \item{\code{sumExp}}{
#'         A \code{DOUData} object containing raw counts, metadata, and 
#'         filtering status.
#'     }
#'     \item{\code{dds}}{
#'         A \code{DESeqDataSet} object prepared for differential 
#'         translation efficiency analysis.
#'     }
#' }
#'
#' @keywords internal
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<-
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom S4Vectors metadata metadata<-
#' 
create_datasets <- function(
        count_table, 
        condition_table,
        annotation, 
        reduced_formula, 
        emm_specs, 
        deseq_formula,
        min_count = 1,
        stringent = TRUE,
        verbose = TRUE
) {
    
    start_parsing <- Sys.time()
    
    # Check required arguments
    required_args <- list(
        count_table = count_table,
        condition_table = condition_table,
        annotation = annotation,
        reduced_formula = reduced_formula,
        emm_specs = emm_specs,
        deseq_formula = deseq_formula
    )
    
    missing_args <- names(required_args)[vapply(required_args, is.null, logical(1))]
    
    if (length(missing_args) > 0) {
        stop(
            "Missing required arguments:",
            paste(missing_args, collapse = ", ")
        )
    }
    
    if (!is.matrix(count_table) && !is.data.frame(count_table)) {
        stop("count_table must be a matrix or data.frame.")
    }
    
    if (!is.data.frame(condition_table)) {
        stop("condition_table should be a data.frame.")
    }
    
    sumExp <- SummarizedExperiment(
        assays = list(counts = count_table),
        colData = condition_table,
        rowRanges = annotation
    )
    
    sumExp <- new(
        "DOUData",
        sumExp,
        formula = reduced_formula,
        specs = emm_specs
    )
    
    # Filtering logic
    if (stringent == TRUE) {
        # TRUE: Keep ORFs where all replicates in at least one condition pass min_count
        keep_list <- lapply(
            levels(condition_table$condition), 
            function(cond_name) {
                # Select both ribo and rna samples for this condition
                samples <- rownames(condition_table[condition_table$condition == cond_name, ])
                rowSums(count_table[, samples, drop = FALSE] >= min_count) == length(samples)
            }
        )
        # if (verbose) message("filtering: All replicates in at least one condition >= ", min_count)
        keep <- Reduce("|", keep_list)
    } else if (stringent == FALSE) {
        # FALSE: Keep ORFs where all replicates in at least one condition-strategy group pass min_count
        condition_table$group <- paste(
            condition_table$condition, 
            condition_table$strategy, 
            sep = "_"
        )
        keep_list <- lapply(levels(condition_table$group), function(group_name) {
            samples <- rownames(condition_table[condition_table$group == group_name, ])
            rowSums(count_table[, samples, drop = FALSE] >= min_count) == length(samples)
        })
        # if (verbose) message("filtering: All replicates in at least one condition-strategy group >= ", min_count)
        keep <- Reduce("|", keep_list)
    } else { # is.null(stringent)
        # NULL: Keep ORFs where total counts across all samples pass min_count
        # if (verbose) message("filtering: Total counts across all samples >= ", min_count)
        keep <- rowSums(count_table) >= min_count
    }
    
    # Add initial filter status to the raw object
    rowData(sumExp)$count_filter <- keep
    
    # Find singlets based on the *initially filtered* data (count_table[keep, ])
    filtered_genes <- rowData(sumExp)$gene_id[keep]
    singlets <- names(which(table(filtered_genes) == 1))
    
    # Status for singlet filter
    rowData(sumExp)$singlet_filter <- !rowData(sumExp)$gene_id %in% singlets
    
    # Create the final filter mask based on count_filter AND non-singlets
    rowData(sumExp)$is_kept <- rowData(sumExp)$count_filter & 
        rowData(sumExp)$singlet_filter
    
    # Creating DESeqDataSet
    dds <- DESeqDataSetFromMatrix(
        countData = assay(sumExp),
        colData = colData(sumExp),
        design = deseq_formula
    )
    
    # Store formulas
    # metadata(dds)$formula <- deseq_fmla$reduced_formula
    metadata(dds)$specs <- emm_specs
    
    if (is.null(sumExp) || !inherits(sumExp, "DOUData") || nrow(sumExp) == 0 || ncol(sumExp) == 0) {
        stop("Failed to create valid a DOUData object with non-empty data.")
    }
    
    if (is.null(dds) || !inherits(dds, "DESeqDataSet") || nrow(dds) == 0 || ncol(dds) == 0) {
        stop("Failed to create valid a DESeqDataSet objects with non-empty data.")
    }
    
    if (verbose) {
        message("SummarizedExperiment objects created successfully")
        
        end_parsing <- Sys.time()
        elapsed_parsing <- runtime(end_parsing, start_parsing)
        
        if (!is.null(elapsed_parsing$mins)) {
            message(
                sprintf(
                    "data parsing runtime: %d mins %.3f secs", 
                    elapsed_parsing$mins, 
                    elapsed_parsing$secs
                )
            )
        } else {
            message(
                sprintf(
                    "data parsing runtime: %.3f secs", 
                    elapsed_parsing$secs
                )
            )
        }
        start_dou <- Sys.time()
    }
    return(list(sumExp = sumExp, dds = dds))
}


#' Construct DOTSeqDatasets from featureCounts output for 
#' Differential ORF Translation Analysis
#'
#' @description
#' This function initialize and construct the 
#' \code{\link{DOTSeqDataSets-class}} object. This includes loading count 
#' and metadata tables, parsing ORF annotations, filtering ORFs based on 
#' count thresholds, and preparing objects for downstream differential 
#' translation analysis using beta-binomial and negative binomial GLM.
#'
#' @param count_table Path to a count table file or a data frame. Must 
#' contain columns: \code{Geneid}, \code{Chr}, \code{Start}, 
#' \code{End}, \code{Strand}, \code{Length}, plus one column per 
#' sample.
#'
#' @param condition_table Path to a sample metadata file or a data frame.
#' Must include columns: \code{run}, \code{strategy}, \code{condition},
#' \code{replicate}.
#'
#' @param flattened_gtf Path to a flattened GFF/GTF file containing ORF 
#' annotations.
#'
#' @param flattened_bed Path to a flattened BED file with ORF annotations.
#'
#' @param formula A formula object specifying the design.
#' Default is \code{~ condition * strategy}.
#'
#' @param target Character string specifying the non-reference condition 
#' level to extract the corresponding interaction term. Contrasted 
#' against the baseline condition. Default is \code{NULL}.
#'
#' @param baseline Character string specifying the desired reference 
#' level. Default is \code{NULL}.
#'
#' @param min_count Minimum count threshold for filtering ORFs.
#' Default is \code{1}.
#'
#' @param stringent Logical or \code{NULL}; determines the filtering 
#' strategy:
#' \describe{
#'     \item{\code{TRUE}}{
#'         Keep ORFs where all replicates in at least one condition 
#'         pass \code{min_count}.
#'     }
#'     \item{\code{FALSE}}{
#'         Keep ORFs where all replicates in at least one 
#'         condition-strategy group pass \code{min_count}.
#'     }
#'     \item{\code{NULL}}{
#'         Keep ORFs where total counts across replicates pass
#'         \code{min_count}.
#'     }
#' }
#'
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#' Default is \code{TRUE}.
#'
#' @return A \code{\link{DOTSeqDataSets-class}} object containing:
#' \describe{
#'     \item{DOU}{
#'         A \code{\link{DOUData-class}} object containing pre-filtered 
#'         raw counts (\code{assay} slot), sample metadata 
#'         (\code{colData} slot), and ORF-level annotation (\code{rowRanges}) 
#'         used for modeling Differential ORF Usage (DOU).
#'     }
#'     \item{DTE}{
#'         A \code{\link{DTEData-class}} object used for modeling 
#'         Differential Translation Efficiency (DTE). Stores all data above
#'         except for \code{rowRanges}.
#'     }
#' }
#' 
#' @rdname DOTSeqDataSets
#' 
#' @importFrom methods new
#' @importFrom SummarizedExperiment mcols mcols<-
#' @importFrom rtracklayer import
#' @importFrom utils read.table
#'
#' @export
#'
#' @examples
#' # Read in count matrix, condition table, and annotation files
#' dir <- system.file("extdata", package = "DOTSeq")
#'
#' cnt <- read.table(
#'     file.path(dir, "featureCounts.cell_cycle_subset.txt.gz"),
#'     header = TRUE,
#'     comment.char = "#"
#' )
#' names(cnt) <- gsub(".*(SRR[0-9]+).*", "\\1", names(cnt))
#'
#' gtf <- file.path(dir, "gencode.v47.orf_flattened_subset.gtf.gz")
#' bed <- file.path(dir, "gencode.v47.orf_flattened_subset.bed.gz")
#'
#' meta <- read.table(file.path(dir, "metadata.txt.gz"))
#' names(meta) <- c("run", "strategy", "replicate", "treatment", "condition")
#' cond <- meta[meta$treatment == "chx", ]
#' cond$treatment <- NULL
#' 
#' # Create a DOTSeqDataSets object
#' d <- DOTSeqDataSetsFromFeatureCounts(
#'     count_table = cnt,
#'     condition_table = cond,
#'     flattened_gtf = gtf,
#'     flattened_bed = bed
#' )
#' 
#' show(d)
#'
DOTSeqDataSetsFromFeatureCounts <- function(
        count_table,
        condition_table,
        flattened_gtf,
        flattened_bed,
        formula = ~ condition * strategy,
        target = NULL,
        baseline = NULL,
        min_count = 1,
        stringent = TRUE,
        verbose = TRUE
) {
    
    # Check required arguments
    required_args <- list(
        count_table = count_table,
        condition_table = condition_table,
        flattened_gtf = flattened_gtf,
        flattened_bed = flattened_bed
    )
    
    missing_args <- names(required_args)[vapply(required_args, is.null, logical(1))]
    
    if (length(missing_args) > 0) {
        stop(
            paste(
                "Missing required arguments:",
                paste(missing_args, collapse = ", ")
            )
        )
    }

    # Validate and reduce formula if needed
    fmla <- reduce_formula(formula, condition_table)
    reduced_formula <- fmla$reduced_formula
    emm_specs <- fmla$emm_specs

    deseq_fmla <- remove_random_effects(formula)
    deseq_fmla <- reduce_formula(deseq_fmla, condition_table)

    cntCols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

    if (is.character(count_table) && file.exists(count_table)) {
        # If count_table is a file path (character) and the file exists
        # Read the first line to get column names
        firstLine <- readLines(count_table, n = 1)
        cntHeader <- strsplit(firstLine, "\t|,|\\s+")[[1]]

        # Check if all expected columns are present
        if (all(cntCols %in% cntHeader)) {
            cond <- read.table(
                count_table, 
                header = TRUE, 
                comment.char = "#", 
                stringsAsFactors = FALSE
            )
        } else {
            stop(
                "Header must contain the following columns: ", 
                "Geneid, Chr, Start, End, Strand, Length"
            )
        }
    } else if (is.data.frame(count_table) && all(cntCols %in% names(count_table))) {
        # If count_table is already a data frame
        cnt <- count_table
    } else {
        stop("count_table must be either a valid file path or a data frame.")
    }

    cond <- parse_condition_table(condition_table)
    
    matched <- match_runs(
        cnt, 
        cond, 
        num_feat_cols = 6, 
        baseline = baseline, 
        verbose = verbose
    )
    cnt <- matched$cnt
    cond <- matched$cond

    dcounts <- cnt[c(1, 7:ncol(cnt))]
    colnames(dcounts) <- c("GeneID", rownames(cond))
    id <- as.character(dcounts[, 1])
    n <- id
    split(n, id) <- lapply(split(n, id), seq_along)
    rownames(dcounts) <- sprintf("%s%s%03.f", id, ":O", as.numeric(n))
    dcounts <- dcounts[, 2:ncol(dcounts)]
    dcounts <- dcounts[, rownames(cond)]

    # Parse flattened GTF and return a GRanges object
    gff_granges <- import(flattened_gtf)
    gff_granges <- gff_granges[gff_granges$type == "exon", ]

    # Ensure all required attributes are present in mcols
    if (!("gene_id" %in% names(mcols(gff_granges)))) {
        stop(
            "The flattened GTF file is missing the 'gene_id' ", 
            "attribute in the ninth column."
        )
    }

    # Prepare the rowData for the SummarizedExperiment
    # Set the ORF names on the GRanges object
    names(gff_granges) <- paste0(
        mcols(gff_granges)$gene_id, 
        ":O", 
        mcols(gff_granges)$exon_number
    )

    # dcounts and gff_granges must have matching names/order before filtering
    dcounts <- dcounts[names(gff_granges), , drop = FALSE]

    # Set final dcounts structure
    dcounts <- as.matrix(dcounts)
    storage.mode(dcounts) <- "integer"

    stopifnot(identical(rownames(dcounts), names(gff_granges)))
    stopifnot(identical(colnames(dcounts), rownames(cond)))

    if (verbose) {
        message("GTF and data parsing successful")
    }

    gr <- annotate_orf_type(flattened_bed, gff_granges)
    gr$orf_number <- gr$exon_number
    mcols(gr) <- mcols(gr)[, c("gene_id", "orf_number", "orf_type")] # "transcripts", 

    se_datasets <- create_datasets(
        count_table = dcounts, 
        condition_table = cond,
        annotation = gr, 
        reduced_formula = reduced_formula, 
        emm_specs = emm_specs, 
        deseq_formula = deseq_fmla$reduced_formula,
        min_count = min_count,
        stringent = stringent,
        verbose = verbose
    )
    
    return(new("DOTSeqDataSets", DOU = se_datasets$sumExp, DTE = se_datasets$dds))
}


#' Generate DOTSeq count matrix from BAM files and ORF annotations
#'
#' This function identifies ORFs from genome and annotation, then counts 
#' reads from BAM files over those ORFs using `summarizeOverlaps()`, 
#' handling both single-end and paired-end libraries.
#'
#' @param count_table A matrix of read counts with features as rows and 
#' samples as columns.
#' 
#' @param condition_table Path to a sample metadata file or a data frame.
#' Must include columns: \code{run}, \code{strategy}, \code{condition},
#' \code{replicate}.
#' 
#' @param annotation A GRanges object with ORF level annotation, 
#' typically obtained from \code{\link{getORFs}}.
#' 
#' @param formula A formula object specifying the design.
#' Default is \code{~ condition * strategy}.
#'
#' @param target Character string specifying the non-reference condition 
#' level to extract the corresponding interaction term. Contrasted 
#' against the baseline condition. Default is \code{NULL}.
#'
#' @param baseline Character string specifying the desired reference 
#' level. Default is \code{NULL}.
#'
#' @param min_count Minimum count threshold for filtering ORFs.
#' Default is \code{1}.
#'
#' @param stringent Logical or \code{NULL}; determines the filtering 
#' strategy:
#' \describe{
#'     \item{\code{TRUE}}{
#'         Keep ORFs where all replicates in at least one condition 
#'         pass \code{min_count}.
#'     }
#'     \item{\code{FALSE}}{
#'         Keep ORFs where all replicates in at least one 
#'         condition-strategy group pass \code{min_count}.
#'     }
#'     \item{\code{NULL}}{
#'         Keep ORFs where total counts across replicates pass
#'         \code{min_count}.
#'     }
#' }
#'
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#' Default is \code{TRUE}.
#'
#' @return A \code{\link{DOTSeqDataSets-class}} object containing:
#' \describe{
#'     \item{DOU}{
#'         A \code{\link{DOUData-class}} object containing pre-filtered 
#'         raw counts (\code{assay} slot), sample metadata 
#'         (\code{colData} slot), and ORF-level annotation (\code{rowRanges}) 
#'         used for modeling Differential ORF Usage (DOU).
#'     }
#'     \item{DTE}{
#'         A \code{\link{DTEData-class}} object used for modeling 
#'         Differential Translation Efficiency (DTE). Stores all data above
#'         except for \code{rowRanges}.
#'     }
#' }
#' 
#' @rdname DOTSeqDataSets
#' 
#' @importFrom methods new
#' 
#' @export
#' 
DOTSeqDataSetsFromSE <- function(
        count_table,
        condition_table,
        annotation,
        formula = ~ condition * strategy,
        target = NULL,
        baseline = NULL,
        min_count = 1,
        stringent = TRUE,
        verbose = TRUE
) {
    # Check required arguments
    required_args <- list(
        count_table = count_table,
        condition_table = condition_table,
        annotation = annotation
    )
    
    missing_args <- names(required_args)[vapply(required_args, is.null, logical(1))]
    
    if (length(missing_args) > 0) {
        stop(
            paste(
                "Missing required arguments:",
                paste(missing_args, collapse = ", ")
            )
        )
    }
    
    # Validate and reduce formula if needed
    fmla <- reduce_formula(formula, condition_table)
    reduced_formula <- fmla$reduced_formula
    emm_specs <- fmla$emm_specs
    
    deseq_fmla <- remove_random_effects(formula)
    deseq_fmla <- reduce_formula(deseq_fmla, condition_table)
    
    if (!is.data.frame(count_table)) {
        stop("count_table must be a data frame.")
    } 
    
    cnt <- count_table
    
    cond <- parse_condition_table(condition_table)
    
    matched <- match_runs(
        cnt, 
        cond, 
        num_feat_cols = 0, 
        baseline = baseline, 
        verbose = verbose
    )
    cnt <- as.matrix(matched$cnt)
    cond <- matched$cond
    
    se_datasets <- create_datasets(
        count_table = cnt, 
        condition_table = cond,
        annotation = annotation, 
        reduced_formula = reduced_formula, 
        emm_specs = emm_specs, 
        deseq_formula = deseq_fmla$reduced_formula,
        min_count = min_count,
        stringent = stringent,
        verbose = verbose
    )
    
    return(new("DOTSeqDataSets", DOU = se_datasets$sumExp, DTE = se_datasets$dds))
}
