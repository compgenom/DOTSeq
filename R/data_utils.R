#' Construct DOTSeq Dataset for Differential ORF Translation Analysis
#'
#' Internal function to initialize and construct the DOTSeq dataset object.
#' This includes loading count and metadata tables, parsing ORF annotations,
#' filtering ORFs based on count thresholds, and preparing objects for downstream
#' differential translation analysis using beta-binomial and negative binomial models.
#'
#' @param count_table Path to a count table file or a data frame. Must contain columns:
#'   \code{Geneid}, \code{Chr}, \code{Start}, \code{End}, \code{Strand}, \code{Length},
#'   plus one column per sample.
#' @param condition_table Path to a sample metadata file or a data frame. Must include columns:
#'   \code{run}, \code{strategy}, \code{condition}, \code{replicate}.
#' @param flattened_gtf Optional path to a flattened GFF/GTF file containing exon definitions.
#' @param bed Path to a BED file with ORF annotations.
#' @param formula A formula object specifying the design, e.g., \code{~ condition * strategy}.
#' @param target Character string specifying the non-reference condition level to extract the corresponding interaction term.
#'   Contrasted against the baseline condition (default: \code{NULL}).
#' @param baseline Character string specifying the desired reference level (default: \code{NULL}).
#' @param min_count Minimum count threshold for filtering ORFs (default: \code{1}).
#' @param stringent Logical or \code{NULL}; determines the filtering strategy:
#'   \describe{
#'     \item{\code{TRUE}}{Keep ORFs where all replicates in at least one condition pass \code{min_count}.}
#'     \item{\code{FALSE}}{Keep ORFs where all replicates in at least one condition-strategy group pass \code{min_count}.}
#'     \item{\code{NULL}}{Keep ORFs where total counts across replicates pass \code{min_count}.}
#'   }
#' @param verbose Logical; if \code{TRUE}, prints progress messages (default: \code{TRUE}).
#'
#' @return A named \code{list} containing:
#' \describe{
#'   \item{sumExp}{A \code{SummarizedExperiment} object containing pre-filtered raw counts and sample metadata, 
#'   used for modeling Differential ORF Usage (DOU) within the DOTSeq framework.}
#'   \item{dds}{A \code{DESeqDataSet} object used for modeling Differential Translation Efficiency (DTE) within the DOTSeq framework via DESeq2.}
#' }
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment assay colData colData<- rowData rowData<- mcols mcols<-
#' @importFrom DESeq2 DESeqDataSetFromMatrix 
#' @importFrom utils read.table tail
#' @importFrom stats relevel
#' 
#' @export
#' 
DOTSeqDataSet <- function(count_table, 
                          condition_table, 
                          flattened_gtf, 
                          bed, 
                          formula = ~ condition * strategy,
                          target = NULL,
                          baseline = NULL,
                          min_count = 1, 
                          stringent = TRUE, 
                          verbose = TRUE) {
  
  start_parsing <- Sys.time()
  
  # Validate and reduce formula if needed
  fmla <- reduce_formula(formula, cond)
  reduced_formula <- fmla$reduced_formula
  emm_specs <- fmla$emm_specs
  
  deseq_fmla <- remove_random_effects(formula)
  deseq_fmla <- reduce_formula(deseq_fmla, cond)
  
  cntCols <- c("Geneid", "Chr", "Start",  "End", "Strand", "Length")
  
  if (is.character(count_table) && file.exists(count_table)) {
    # If count_table is a file path (character) and the file exists
    # Read the first line to get column names
    firstLine <- readLines(count_table, n = 1)
    cntHeader <- strsplit(firstLine, "\t|,|\\s+")[[1]]
    
    # Check if all expected columns are present
    if (all(cntCols %in% cntHeader)) {
      cond <- read.table(count_table, header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
    } else {
      stop("Header must contain the following columns: Geneid, Chr, Start, End, Strand, Length")
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
      cond <- read.table(condition_table, header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
      names(cond) <- normaliseNames(names(cond))
    } else {
      stop(paste("Missing expected columns (case-insensitive match):", paste(missingCols, collapse = ", ")))
    }
    
  } else if (is.data.frame(condition_table)) {
    condHeader <- normaliseNames(names(condition_table))
    missingCols <- setdiff(normaliseNames(condCols), condHeader)
    
    if (length(missingCols) == 0) {
      cond <- condition_table
      names(cond) <- normaliseNames(names(cond))
    } else {
      stop(paste("Data frame is missing expected columns (case-insensitive match):", paste(missingCols, collapse = ", ")))
    }
    
  } else {
    stop("`condition_table` must be either a valid file path or a data frame.")
  }
  
  # Rename cond rownames
  rownames(cond) <- cond$run
  
  # Find common identifiers (sample names that are count columns)
  common_runs <- intersect(cond$run, names(cnt))
  
  if (length(common_runs) < 1) {
    stop("No common 'run' IDs found between `condition_table` and `count_table` sample columns.")
  }
  
  # Use the 'cond' object as the base metadata for sample duplication
  cond_base <- cond[match(common_runs, cond$run), ]
  cnt_base_samples <- cnt[, common_runs, drop = FALSE]
  
  # Find the common count columns
  common_runs <- intersect(rownames(cond), names(cnt))
  # Subset cnt to GeneID + common runs
  cnt <- cnt[, c(cntCols, common_runs)]
  
  if (all(colnames(cnt) %in% rownames(cond))) {
    cond <- cond[colnames(cnt), , drop = FALSE] # Enforce sample order
  }
  
  # Convert to factors
  for (col in colnames(cond)) {
    if (is.character(cond[[col]])) {
      cond[[col]] <- factor(cond[[col]])
    }
  }
  
  # Set baseline
  strategy_vals <- as.character(cond$strategy)
  
  # Define a flexible regex pattern: matches "rna", "RNA", "RNA-seq", and "0"
  rna_pattern <- "(?i)rna|^0$"
  
  # Find matching values
  rna_like <- unique(strategy_vals[grepl(rna_pattern, strategy_vals)])
  
  # Handle multiple matches to RNA-seq
  if (length(rna_like) > 1) {
    warning("Multiple RNA-seq levels found. Use ", rna_like[1], " as reference.")
    rna_like <- rna_like[1]
  } else if (length(rna_like) == 0) {
    stop("No RNA-seq level found. Please use rna, RNA, RNA-seq, or 0 to represent RNA-seq.")
  }
  
  # Relevel RNA-seq as baseline
  cond$strategy <- relevel(cond$strategy, ref = rna_like)
  
  # Relevel the baseline condition as baseline
  if (is.null(baseline)) {
    condition_vals <- as.character(cond$condition)
    numeric_vals <- suppressWarnings(as.numeric(condition_vals))
    
    baseline <- NULL  # initialize
    
    if (!any(is.na(numeric_vals))) {
      unique_numeric <- sort(unique(numeric_vals))
      
      # If binary, use the smaller value
      if (length(unique_numeric) == 2) {
        baseline <- as.character(unique_numeric[1])
      }
    } else {
      baseline <- as.character(tail(sort(unique(condition_vals)), 1))
    }
  }
  
  baseline <- as.character(baseline)
  
  # Relevel the factor
  if (isTRUE(verbose)) {
    message("setting condition", baseline, " as baseline")
    cond$condition <- relevel(factor(cond$condition), ref = baseline)
  }
  
  # Remove columns with fewer than 2 levels
  validCols <- sapply(cond, function(x) !(is.factor(x) && length(unique(x)) < 2))
  cond <- cond[, validCols]
  
  cnt <- cnt[order(cnt$Geneid, cnt$Start, cnt$End), ]
  dcounts <- cnt[, c("Geneid", rownames(cond))] # Subset by Geneid and combined sample names
  colnames(dcounts)[1] <- "GeneID"
  
  # Create ORF IDs
  id <- as.character(dcounts$GeneID)
  n <- id
  split(n, id) <- lapply(split(n, id), seq_along)
  orf_ids <- sprintf("%s:O%03.f", id, as.numeric(n))
  
  # Filter out rows starting with "_" and align counts
  dcounts <- dcounts[substr(id, 1, 1) != "_", ]
  orf_ids <- orf_ids[substr(id, 1, 1) != "_"]
  
  # Set final dcounts structure
  dcounts <- round(as.matrix(dcounts[, -1])) # Remove GeneID column and convert to matrix
  storage.mode(dcounts) <- "integer"
  
  rownames(dcounts) <- orf_ids
  # dcounts <- dcounts[, rownames(cond)] # Final alignment of samples
  
  
  if (!is.null(flattened_gtf)) {
    # This automatically returns a GRanges object
    gff_granges <- Bioc.gff::import(flattened_gtf)
    gff_granges <- gff_granges[gff_granges$type=="exon", ]
    
    # Ensure all required attributes are present in mcols (e.g., gene_id, exon_number)
    if (!("gene_id" %in% names(mcols(gff_granges)))) {
      stop("The flattened GTF file is missing the 'gene_id' attribute in the ninth column.")
    }
    
    # Prepare the rowData for the SummarizedExperiment
    # Set the ORF names on the GRanges object
    names(gff_granges) <- paste0(mcols(gff_granges)$gene_id,":O" , mcols(gff_granges)$exon_number)
    
    # dcounts and gff_granges must have matching names/order before filtering
    dcounts <- dcounts[names(gff_granges), , drop = FALSE] 
    stopifnot(identical(rownames(dcounts), names(gff_granges)))
    stopifnot(identical(colnames(dcounts), rownames(cond)))
    
    if (verbose) {
      message("GTF and data parsing successful")
    }
  }
  
  gr <- annotate_orf_type(bed, gff_granges)
  gr$orf_number <- gr$exon_number
  mcols(gr) <- mcols(gr)[, c("gene_id", "transcripts", "orf_number")]
  
  # Create a temporary SummarizedExperiment object for filtering
  sumExp <- SummarizedExperiment(
    assays = list(counts = dcounts),
    colData = cond,
    rowRanges = gr
  )
  
  # Filtering logic
  if (stringent == TRUE) {
    # TRUE: Keep ORFs where all replicates in at least one condition pass min_count
    keep_list <- lapply(levels(cond$condition), function(cond_name) {
      # Select both ribo and rna samples for this condition
      samples <- rownames(cond[cond$condition == cond_name, ]) 
      rowSums(dcounts[, samples, drop = FALSE] >= min_count) == length(samples)
    })
    # if (verbose) message("filtering: All replicates in at least one condition >= ", min_count)
    keep <- Reduce("|", keep_list)
  } else if (stringent == FALSE) {
    # FALSE: Keep ORFs where all replicates in at least one condition-strategy group pass min_count
    cond$group <- paste(cond$condition, cond$strategy, sep = "_")
    keep_list <- lapply(levels(cond$group), function(group_name) {
      samples <- rownames(cond[cond$group == group_name, ])
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
  rowData(sumExp)$is_kept <- rowData(sumExp)$count_filter & rowData(sumExp)$singlet_filter

  
  # Creating DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = assay(sumExp),
                                        colData = colData(sumExp),
                                        design = deseq_fmla$reduced_formula)
  
  # Store formulas
  metadata(sumExp)$formula <- reduced_formula
  metadata(sumExp)$emm_specs <- emm_specs
  metadata(dds)$formula <- deseq_fmla$reduced_formula
  metadata(dds)$emm_specs <- emm_specs
  
  
  if (verbose) {
    message("SummarizedExperiment objects created successfully")
  }

  return(list(sumExp = sumExp,
              dds = dds))
  
}


#' Annotate ORF Types from BED File
#'
#' This function reads a BED file and converts it into a `GRanges` object.
#' It then finds overlaps with a reference set of ORFs (provided as a `GRanges` object,
#' typically derived from a GFF or DEXSeqDataSet), and annotates each overlapping ORF
#' with a label (e.g., score or type) from the BED file.
#'
#' @param bed Character string. Path to a BED file (tab-delimited, no header) with six columns:
#'   \describe{
#'     \item{chr}{Chromosome name}
#'     \item{start}{Start coordinate (0-based)}
#'     \item{end}{End coordinate (exclusive)}
#'     \item{name}{Feature name}
#'     \item{score}{Annotation label to assign (e.g., ORF type)}
#'     \item{strand}{Strand information ("+" or "-")}
#'   }
#' @param gff_granges A `GRanges` object containing ORF coordinates to be annotated.
#'
#' @return A `GRanges` object with an added metadata column `orf_type`, containing
#'   the BED score for overlapping ORFs. Non-overlapping ORFs will have `NA`.
#'
#' @importFrom utils read.table
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits mcols mcols<-
#' @importFrom IRanges IRanges
#' 
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' # Annotate a GRange using a BED file's score 
#' gff_granges <- parseBed("example.bed", granges)
#' head(orfs)
#' }
annotate_orf_type <- function(bed, gff_granges) {
  bedDf <- read.table(bed, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(bedDf) <- c("chr", "start", "end", "name", "score", "strand")
  
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