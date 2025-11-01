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


#' Construct DOTSeqDatasets for Differential ORF Translation Analysis
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
#' @return A \code{DOTSeqDataSets} object containing:
#' \describe{
#'     \item{DOU}{
#'         A \code{\link{DOUData-class}} object containing pre-filtered 
#'         raw counts and sample metadata, used for modeling 
#'         Differential ORF Usage (DOU).
#'     }
#'     \item{DTE}{
#'         A \code{\link{DTEData-class}} object used for modeling 
#'         Differential Translation Efficiency (DTE).
#'     }
#' }
#' 
#' @rdname DOTSeqDataSets
#' 
#' @importFrom S4Vectors DataFrame SimpleList metadata metadata<-
#' @importFrom methods new
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom SummarizedExperiment rowData rowData<- mcols mcols<-
#' @importFrom rtracklayer import
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom utils read.table tail
#' @importFrom stats relevel
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
#' # Create a raw input list
#' raw <- list(
#'     count_table = cnt,
#'     condition_table = cond,
#'     flattened_gtf = gtf,
#'     flattened_bed = bed
#' )
#' 
#' # Create a DOTSeqDataSets object
#' d <- DOTSeqDataSets(
#'     count_table = cnt,
#'     condition_table = cond,
#'     flattened_gtf = gtf,
#'     flattened_bed = bed
#' )
#' 
#' show(d)
#'
DOTSeqDataSets <- function(
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
    start_parsing <- Sys.time()

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
    
    cond_run <- cond$run

    # Rename cond rownames
    rownames(cond) <- cond$run
    # Find common identifiers
    common <- intersect(rownames(cond), names(cnt))

    cond <- cond[common, , drop = FALSE]
    cond <- cond[order(cond$strategy, cond$replicate), ]
    
    cnt_run <- colnames(cnt[c(7:ncol(cnt))])
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

    # Combine with metadata columns
    cnt <- cnt[, c(cntCols, rownames(cond))]

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
    if (ncol(cnt) == (length(rownames(cond)) + 6)) {
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

    dcounts <- cnt[c(1, 7:ncol(cnt))]
    colnames(dcounts) <- c("GeneID", rownames(cond))
    id <- as.character(dcounts[, 1])
    n <- id
    split(n, id) <- lapply(split(n, id), seq_along)
    rownames(dcounts) <- sprintf("%s%s%03.f", id, ":O", as.numeric(n))
    dcounts <- dcounts[, 2:ncol(dcounts)]
    dcounts <- dcounts[, rownames(combined_cond)]

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
    stopifnot(identical(colnames(dcounts), rownames(combined_cond)))

    if (verbose) {
        message("GTF and data parsing successful")
    }

    gr <- annotate_orf_type(flattened_bed, gff_granges)
    gr$orf_number <- gr$exon_number
    mcols(gr) <- mcols(gr)[, c("gene_id", "orf_number", "orf_type")] # "transcripts", 

    sumExp <- SummarizedExperiment(
        assays = list(counts = dcounts),
        colData = combined_cond,
        rowRanges = gr
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
            levels(combined_cond$condition), 
            function(cond_name) {
                # Select both ribo and rna samples for this condition
                samples <- rownames(cond[combined_cond$condition == cond_name, ])
                rowSums(dcounts[, samples, drop = FALSE] >= min_count) == length(samples)
            }
        )
        # if (verbose) message("filtering: All replicates in at least one condition >= ", min_count)
        keep <- Reduce("|", keep_list)
    } else if (stringent == FALSE) {
        # FALSE: Keep ORFs where all replicates in at least one condition-strategy group pass min_count
        combined_cond$group <- paste(
            combined_cond$condition, 
            combined_cond$strategy, 
            sep = "_"
        )
        keep_list <- lapply(levels(combined_cond$group), function(group_name) {
            samples <- rownames(cond[combined_cond$group == group_name, ])
            rowSums(dcounts[, samples, drop = FALSE] >= min_count) == length(samples)
        })
        # if (verbose) message("filtering: All replicates in at least one condition-strategy group >= ", min_count)
        keep <- Reduce("|", keep_list)
    } else { # is.null(stringent)
        # NULL: Keep ORFs where total counts across all samples pass min_count
        # if (verbose) message("filtering: Total counts across all samples >= ", min_count)
        keep <- rowSums(dcounts) >= min_count
    }

    # Add initial filter status to the raw object
    rowData(sumExp)$count_filter <- keep

    # Find singlets based on the *initially filtered* data (dcounts[keep, ])
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
        design = deseq_fmla$reduced_formula
    )

    # Store formulas
    # metadata(dds)$formula <- deseq_fmla$reduced_formula
    metadata(dds)$specs <- emm_specs

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
    
    return(new("DOTSeqDataSets", DOU = sumExp, DTE = dds))
}

