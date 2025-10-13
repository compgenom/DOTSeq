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
#' @importFrom Bioc.gff import
#' @importFrom DESeq2 DESeqDataSetFromMatrix 
#' @importFrom utils read.table tail
#' @importFrom stats relevel
#' 
#' @export
#' 
DOTSeqDataSet <- function(
    count_table, 
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
  fmla <- reduce_formula(formula, condition_table)
  reduced_formula <- fmla$reduced_formula
  emm_specs <- fmla$emm_specs
  
  deseq_fmla <- remove_random_effects(formula)
  deseq_fmla <- reduce_formula(deseq_fmla, condition_table)
  
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
  # Find common identifiers
  common <- intersect(rownames(cond), names(cnt))

  cond <- cond[common, , drop = FALSE]
  cond <- cond[order(cond$strategy, cond$replicate), ]

  # Combine with metadata columns
  cnt <- cnt[, c(cntCols, rownames(cond))]
  
  # Rename header
  rownames(cond) <- apply(cond[, c("run","condition", "replicate", "strategy")], 1, function(x) {
    paste0(trimws(x[1]), ".", trimws(x[2]), ".", trimws(x[3]), ".", trimws(x[4]))
  })
  # Find the indices of the columns that match cond$run
  libIndices <- match(cond$run, names(cnt))
  # Replace those column names with rownames(cond)
  if (ncol(cnt) == (length(rownames(cond))+6)) {
    names(cnt)[libIndices] <- rownames(cond)
  } else {
    stop("Number of samples in count_table and condition_table doesn't match.")
  }
  
  riboCond <- cond[cond$strategy=="ribo",]
  riboCond$strategy <- 1
  rnaCond <- cond[cond$strategy=="rna",]
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
    combined_cond$condition <- relevel(combined_cond$condition, ref = baseline)
  }
  
  # Remove columns with fewer than 2 levels
  validCols <- sapply(combined_cond, function(x) !(is.factor(x) && length(unique(x)) < 2))
  combined_cond <- combined_cond[, validCols]
  
  
  cnt <- cnt[order(cnt$Geneid, cnt$Start, cnt$End), ]
  dcounts <- cnt[c(1, 7:ncol(cnt))]
  colnames(dcounts) <- c("GeneID", rownames(cond))
  id <- as.character(dcounts[,1])
  n <- id
  split(n,id) <- lapply(split(n ,id), seq_along )
  rownames(dcounts) <- sprintf("%s%s%03.f",id,":O",as.numeric(n))
  dcounts <- dcounts[,2:ncol(dcounts)]
  dcounts <- dcounts[, rownames(combined_cond)]
  
  
  # Parse flattened GTF and return a GRanges object
  gff_granges <- import(flattened_gtf)
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
  
  # Set final dcounts structure
  dcounts <- as.matrix(dcounts)
  storage.mode(dcounts) <- "integer"
  
  stopifnot(identical(rownames(dcounts), names(gff_granges)))
  stopifnot(identical(colnames(dcounts), rownames(combined_cond)))
  
  if (verbose) {
    message("GTF and data parsing successful")
  }
  
  gr <- annotate_orf_type(bed, gff_granges)
  gr$orf_number <- gr$exon_number
  mcols(gr) <- mcols(gr)[, c("gene_id", "transcripts", "orf_number", "orf_type")]

  # combined_cond$sample <- rownames(combined_cond)
  sumExp <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = dcounts), 
                                                       colData = combined_cond, 
                                                       rowRanges = gr)
  
  # Filtering logic
  if (stringent == TRUE) {
    # TRUE: Keep ORFs where all replicates in at least one condition pass min_count
    keep_list <- lapply(levels(combined_cond$condition), function(cond_name) {
      # Select both ribo and rna samples for this condition
      samples <- rownames(cond[combined_cond$condition == cond_name, ]) 
      rowSums(dcounts[, samples, drop = FALSE] >= min_count) == length(samples)
    })
    # if (verbose) message("filtering: All replicates in at least one condition >= ", min_count)
    keep <- Reduce("|", keep_list)
  } else if (stringent == FALSE) {
    # FALSE: Keep ORFs where all replicates in at least one condition-strategy group pass min_count
    combined_cond$group <- paste(combined_cond$condition, combined_cond$strategy, sep = "_")
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
  rowData(sumExp)$is_kept <- rowData(sumExp)$count_filter & rowData(sumExp)$singlet_filter
  
  # Creating DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = assay(sumExp),
    colData = colData(sumExp),
    design = deseq_fmla$reduced_formula
  )
  
  # Store formulas
  metadata(sumExp)$formula <- reduced_formula
  metadata(sumExp)$emm_specs <- emm_specs
  metadata(dds)$formula <- deseq_fmla$reduced_formula
  metadata(dds)$emm_specs <- emm_specs
  
  if (verbose) {
    message("SummarizedExperiment objects created successfully")

    end_parsing <- Sys.time()
    elapsed_parsing <- runtime(end_parsing, start_parsing)
    
    if (!is.null(elapsed_parsing$mins)) {
      message(sprintf("data parsing runtime: %d mins %.3f secs", elapsed_parsing$mins, elapsed_parsing$secs))
    } else {
      message(sprintf("data parsing runtime: %.3f secs", elapsed_parsing$secs))
    }
    start_dou <- Sys.time()
  }
  
  return(list(
    sumExp = sumExp,
    dds = dds
  ))
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
  
  mORFs <- bedDf[bedDf$score == "mORF", ]
  
  # Identify duplicated names in both directions
  dups_name <- duplicated(mORFs$name) | duplicated(mORFs$name, fromLast = TRUE)
  
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