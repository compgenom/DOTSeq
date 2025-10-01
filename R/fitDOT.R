#' Fit DOTSeq Differential Translation Models
#'
#' This function performs the complete DOTSeq analysis pipeline:
#' loading count data, aligning with sample metadata, filtering ORFs,
#' normalizing counts, calculating translational efficiency (TE), and
#' fitting beta-binomial and negative binomial generalized linear models for 
#' differential ORF translation using \code{DOTSeq::fitDOU} and \code{DESeq2::DESeq}.
#'
#' @param count_table Either a path to a count table file or a data frame.
#'   Must contain columns: Geneid, Chr, Start, End, Strand, Length,
#'   plus one column per sample.
#' @param condition_table Either a path to a condition metadata file or a data frame.
#'   Must include columns: run, strategy, condition, replicate.
#' @param flattened_gtf Optional path to a flattened GFF/GTF file containing exon definitions.
#' @param bed Path to a BED file with ORF annotations.
#' @param target Character string specifying the non-reference condition level to extract the corresponding interaction term from the model. 
#' This is contrasted against the baseline condition (default: \code{NULL}).
#' @param baseline Character string specifying the desired reference level.
#' @param formula A formula object specifying the design, e.g., \code{~ condition * strategy}.
#' @param pseudocount Numeric pseudo-count to avoid division by zero when computing TE (default: \code{1e-6}).
#' @param min_count Minimum count threshold for filtering ORFs (default: \code{1}).
#' @param stringent Logical or \code{NULL}; determines the filtering strategy:
#'   \describe{
#'     \item{\code{TRUE}}{Keep ORFs where all replicates in at least one condition pass \code{min_count}.}
#'     \item{\code{FALSE}}{Keep ORFs where all replicates in at least one condition-strategy group pass \code{min_count}.}
#'     \item{\code{NULL}}{Keep ORFs where total counts across replicates pass \code{min_count}.}
#'   }
#' @param parallel A list passed to \code{glmmTMBControl} to configure parallel optimization.
#'   If \code{NULL}, parallelism is disabled (default: \code{list(n = 4L, autopar = TRUE)}).
#' @param optimizers Logical; if \code{TRUE}, enables brute-force optimization using multiple optimizers
#'   in \code{glmmTMB}: \code{nlminb}, \code{bobyqa}, and \code{optim} (default: \code{FALSE}).
#' @param dispersion_modeling Optional string specifying the dispersion modeling strategy.
#'   Used to select between default, constant, or custom dispersion models.
#' @param dispformula Optional formula object for custom dispersion modeling.
#' @param lrt Logical; if \code{TRUE}, performs a likelihood ratio test to compare the full model (with interaction) against a reduced model 
#' (without interaction) to assess translation-specific effects (default: \code{FALSE}).
#' @param diagnostic Logical; if \code{TRUE}, enables model diagnostics including tests for overdispersion,
#'   zero inflation, and residual properties (default: \code{FALSE}).
#' @param verbose Logical; if \code{TRUE}, prints progress messages (default: \code{TRUE}).
#'
#' @return A named \code{list} with the following elements:
#' \describe{
#'   \item{raw_counts}{Raw counts matrix for all samples.}
#'   \item{norm_counts}{Normalized counts matrix for all samples.}
#'   \item{orfs}{Data frame of ORFs derived from the BED file matched to the DOTSeq object.}
#'   \item{dds}{DESeq2 object used for modeling differential translation efficiency.}
#'   \item{sumExp}{SummarizedExperiment object containing normalized counts and sample metadata.}
#'   \item{dxd}{DOTSeq object used for modeling exon/ORF-level counts.}
#'   \item{formula}{Conditional formula used in \code{DOTSeq::fitDOU}.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- fitDOT(
#'   count_table = "counts.txt",
#'   condition_table = "samples.txt",
#'   flattened_gtf = "flattened.gff",
#'   bed = "orfs.bed",
#'   pseudocount = 1e-6,
#'   min_count = 1,
#'   stringent = TRUE,
#'   verbose = TRUE
#' )
#' head(result$sumExp)
#' }
#'
#' @export
fitDOT <- function(count_table, condition_table, 
                   flattened_gtf, bed, 
                   # rna_suffix = ".rna", 
                   # ribo_suffix = ".ribo", 
                   target = NULL,
                   baseline = NULL,
                   formula = ~ condition * strategy,
                   dispersion_modeling = "strategy",
                   dispformula = NULL,
                   lrt = FALSE,
                   diagnostic = FALSE,
                   # sampleDelim = NULL,
                   pseudocount = 1e-6, 
                   min_count = 1, 
                   stringent = TRUE, 
                   parallel = list(n=4L, autopar=TRUE),
                   optimizers = FALSE,
                   verbose = TRUE) {
  
  if (verbose) {
    start_parsing <- Sys.time()
    message(" - Use ", parallel$n, " threads")
  }
  
  cnt_cols <- c("Geneid", "Chr", "Start",  "End", "Strand", "Length")
  
  if (is.character(count_table) && file.exists(count_table)) {
    # If count_table is a file path (character) and the file exists
    # Read the first line to get column names
    first_line <- readLines(count_table, n = 1)
    cnt_header <- strsplit(first_line, "\t|,|\\s+")[[1]]
    
    # Check if all expected columns are present
    if (all(cnt_cols %in% cnt_header)) {
      cond <- read.table(count_table, header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
    } else {
      stop("Header must contain the following columns: Geneid, Chr, Start, End, Strand, Length")
    }
    
  } else if (is.data.frame(count_table) && all(cnt_cols %in% names(count_table))) {
    # If count_table is already a data frame
    cnt <- count_table
  } else {
    stop("count_table must be either a valid file path or a data frame.")
  }
  
  # Define expected column names
  cond_cols <- c("run", "strategy", "condition", "replicate")
  
  # Helper to normalize column names
  normalise_names <- function(x) tolower(trimws(x))
  
  if (is.character(condition_table) && file.exists(condition_table)) {
    # Read header line and normalize
    first_line <- readLines(condition_table, n = 1)
    cond_header <- normalise_names(strsplit(first_line, "\t|,|\\s+")[[1]])
    
    missing_cols <- setdiff(normalise_names(cond_cols), cond_header)
    
    if (length(missing_cols) == 0) {
      cond <- read.table(condition_table, header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
      names(cond) <- normalise_names(names(cond))
    } else {
      stop(paste("Missing expected columns (case-insensitive match):", paste(missing_cols, collapse = ", ")))
    }
    
  } else if (is.data.frame(condition_table)) {
    cond_header <- normalise_names(names(condition_table))
    missing_cols <- setdiff(normalise_names(cond_cols), cond_header)
    
    if (length(missing_cols) == 0) {
      cond <- condition_table
      names(cond) <- normalise_names(names(cond))
    } else {
      stop(paste("Data frame is missing expected columns (case-insensitive match):", paste(missing_cols, collapse = ", ")))
    }
    
  } else {
    stop("`condition_table` must be either a valid file path or a data frame.")
  }
  
  # Rename cond rownames
  rownames(cond) <- cond$run
  # Find common identifiers
  common <- intersect(rownames(cond), names(cnt))
  # # Subset and reorder cnt columns to match cond rownames
  # cnt_aligned <- cnt[, common, drop = FALSE]
  # Subset cond to match cnt columns
  cond <- cond[common, , drop = FALSE]
  
  cond <- cond[order(cond$strategy, cond$replicate), ]
  
  # Combine with metadata columns
  # cnt <- cbind(cnt[names(cnt) %in% cnt_cols], cnt_aligned)
  cnt <- cnt[, c(cnt_cols, rownames(cond))]
  
  # Rename header
  rownames(cond) <- apply(cond[, c("run","condition", "replicate", "strategy")], 1, function(x) {
    paste0(trimws(x[1]), ".", trimws(x[2]), ".", trimws(x[3]), ".", trimws(x[4]))
  })
  # Find the indices of the columns that match cond$run
  lib_indices <- match(cond$run, names(cnt))
  # Replace those column names with rownames(cond)
  if (ncol(cnt) == (length(rownames(cond))+6)) {
    names(cnt)[lib_indices] <- rownames(cond)
  } else {
    stop("Number of samples in count_table and condition_table doesn't match.")
  }
  
  cond$strategy <- ifelse(grepl("ribo", cond$strategy, ignore.case = TRUE), 1,
                                 ifelse(grepl("rna", cond$strategy, ignore.case = TRUE), 0, NA))
  
  combined_cond <- cond[, !(names(cond) %in% "run")]
  
  # Convert to factors
  for (col in colnames(combined_cond)) {
    if (is.character(combined_cond[[col]])) {
      combined_cond[[col]] <- factor(combined_cond[[col]])
    }
  }
  
  # Set baseline
  if (!is.null(baseline)) {
    combined_cond$condition <- relevel(combined_cond$condition, ref = baseline)
  }
  
  # Remove columns with fewer than 2 levels
  valid_cols <- sapply(combined_cond, function(x) !(is.factor(x) && length(unique(x)) < 2))
  combined_cond <- combined_cond[, valid_cols]
  
  fmla <- DOTSeq:::reduce_formula(formula, combined_cond)
  
  if (verbose) {
    message(" - Start differential ORF usage analysis")
    message(" - Conditional formula: ", deparse(fmla))
  }
  cnt <- cnt[order(cnt$Geneid, cnt$Start, cnt$End), ]
  dcounts <- cnt[c(1, 7:ncol(cnt))]
  colnames(dcounts) <- c("GeneID", rownames(cond))
  id <- as.character(dcounts[,1])
  n <- id
  split(n,id) <- lapply(split(n ,id), seq_along )
  rownames(dcounts) <- sprintf("%s%s%03.f",id,":O",as.numeric(n))
  dcounts <- dcounts[,2:ncol(dcounts)]
  
  dcounts <- dcounts[substr(rownames(dcounts), 1, 1) != "_", ] #remove _ from beginning of gene name 
  dcounts <- dcounts[, rownames(combined_cond)]
  
  ## get genes and exon names out
  splitted <- strsplit(rownames(dcounts), ":")
  exons <- sapply(splitted, "[[", 2)
  genesrle <- sapply(splitted, "[[", 1)
  
  
  if (!is.null(flattened_gtf)) {
    if (verbose) {
      message(" - Parse the flattened GFF file")
    }
    aggregates <- read.delim(flattened_gtf, stringsAsFactors = FALSE, 
                             header = FALSE)
    colnames(aggregates) <- c("chr", "source", "class", "start", 
                              "end", "ex", "strand", "ex2", "attr")
    aggregates$strand <- gsub("\\.", "*", aggregates$strand)
    aggregates <- aggregates[which(aggregates$class == "exon"), # "exonic_part"
    ]
    aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
    aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", 
                              aggregates$attr)
    # trim the gene_ids to 255 chars in order to match with featurecounts
    longIDs <- sum(nchar(unique(aggregates$gene_id)) > 255)
    warning(paste0(longIDs, 
                   " aggregate geneIDs were found truncated in featureCounts output"), 
            call. = FALSE)
    aggregates$gene_id <- substr(aggregates$gene_id,1,255)
    
    # Sort as dcounts
    aggregates <- aggregates[order(aggregates$gene_id, aggregates$start, aggregates$end), ]
    
    
    transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", 
                        aggregates$attr)
    transcripts <- strsplit(transcripts, "\\+")
    exonids <- gsub(".*exon_number\\s(\\S+).*", "\\1", # exonic_part_number
                    aggregates$attr)
    
    exoninfo <- GRanges(
      seqnames = as.character(aggregates$chr),
      ranges = IRanges(start = aggregates$start, end = aggregates$end),
      strand = aggregates$strand
    )
    
    names(exoninfo) <- paste(aggregates$gene_id, exonids, sep = ":O")
    
    names(transcripts) <- names(exoninfo) 
    if (!all(rownames(dcounts) %in% names(exoninfo))) {
      stop("Count files do not correspond to the flattened annotation file")
    }
    matching <- match(rownames(dcounts), names(exoninfo))
    stopifnot(all(names(exoninfo[matching]) == rownames(dcounts)))
    stopifnot(all(names(transcripts[matching]) == rownames(dcounts)))
    
    # Check for NA
    if (any(is.na(matching))) {
      stop("Some dcounts rows do not match exoninfo names. Possible truncation or naming mismatch.")
    }
    
    # Extract metadata from count table
    cnt_meta <- cnt[, c("Geneid", "Chr", "Start", "End", "Strand")]
    
    # Extract metadata from flattened GFF
    gff_meta <- data.frame(
      Geneid = aggregates$gene_id,
      Chr = aggregates$chr,
      Start = aggregates$start,
      End = aggregates$end,
      Strand = aggregates$strand
    )
    
    cnt_meta$key <- paste(cnt_meta$Geneid, cnt_meta$Start, cnt_meta$End, cnt_meta$Strand, sep = "_")
    gff_meta$key <- paste(gff_meta$Geneid, gff_meta$Start, gff_meta$End, gff_meta$Strand, sep = "_")
    
    if (any(duplicated(cnt_meta$key)) || any(duplicated(gff_meta$key))) {
      stop("Duplicate keys found in count table or flattened GFF. Cannot guarantee reliable alignment.")
    }
    
    if (all(cnt_meta$key == gff_meta$key)) {
      if (verbose) {
        message(" - Create a DEXSeq object")
      }
      
      # Convert integer-like and character columns to factors if they are meant to be categorical
      for (col in colnames(combined_cond)) {
        if (is.numeric(combined_cond[[col]]) && all(combined_cond[[col]] == as.integer(combined_cond[[col]]))) {
          combined_cond[[col]] <- factor(combined_cond[[col]])
        } else if (is.character(combined_cond[[col]])) {
          combined_cond[[col]] <- factor(combined_cond[[col]])
        }
      }
      
      design = ~ sample + exon + condition:exon
      dxd <- DEXSeqDataSet(dcounts, combined_cond, design, exons, 
                           genesrle, exoninfo[matching], transcripts[matching])
      
      # Get ORF annotation
      ORF_df <- parseBed(bed, dxd)
      names(ORF_df)[names(ORF_df) == "exonBaseMean"] <- "orfBaseMean"
      names(ORF_df)[names(ORF_df) == "exonBaseVar"] <- "orfBaseVar"
      
      counts <- counts(dxd)
      
      if (stringent == TRUE) {
        keep_list <- lapply(unique(cond$condition), function(cond_name) {
          samples <- rownames(cond[cond$condition == cond_name, ])
          rowSums(dcounts[, samples, drop = FALSE] >= min_count) == length(samples)
        })
        
        if (verbose) {
          message(" - Keep ORFs where all replicates in at least one condition have counts >= ", min_count)
        }
        
        keep <- Reduce("|", keep_list)
        dxd <- dxd[keep, ]
      } else if (stringent == FALSE) {
        cond$group <- paste(cond$condition, cond$strategy, sep = "_")
        keep_list <- lapply(unique(cond$group), function(group_name) {
          samples <- rownames(cond[cond$group == group_name, ])
          rowSums(dcounts[, samples, drop = FALSE] >= min_count) == length(samples)
        })
        
        if (verbose) {
          message(" - Keep ORFs where all replicates in at least one condition-strategy group have counts >= ", min_count)
        }
        
        keep <- Reduce("|", keep_list)
        dxd <- dxd[keep, ]
        
      } else if (is.null(stringent)) {
        if (verbose) {
          message(" - Keep ORFs where the total counts across replicates are >= ", min_count)
        }
        dxd <- dxd[rowSums(featureCounts(dxd)) >= min_count,]
      }
      
      # Find the gene IDs that have only one ORF "(exon)"
      singlets <- names(which(table(rowData(dxd)$groupID) == 1))
      if (verbose) {
        message(" - Filter out ", length(singlets), " single ORF genes")
      }
      # `! ... %in% ...` means "not in this list"
      keep_exons <- !rowData(dxd)$groupID %in% singlets
      dxd <- dxd[keep_exons, ]
      if (verbose) {
        message(" - Number of ORFs passing filter: ", nrow(featureCounts(dxd)))
      }
      
      # Get normalised counts
      dxd <- estimateSizeFactors(dxd)
      norm_counts <- counts(dxd, normalized = TRUE)
      colnames(norm_counts) <- colData(dxd)$sample
      norm_counts <- norm_counts[, seq_len(ncol(norm_counts)/2)]
      
      # te <- calculateTE(norm_counts, sampleDelim = sampleDelim, rna_suffix = rna_suffix, ribo_suffix = ribo_suffix, pseudocount = pseudocount)
      
      # Run DOTSeq::fitDOU
      exonInfo <- rowData(dxd)
      colnames(exonInfo)[1:2] <- c("isoform_id", "gene_id")
      exonInfo$isoform_id <- rownames(exonInfo)
      
      # Get sampleAnnotation
      anno <- sampleAnnotation(dxd)
      anno$sizeFactor <- NULL
      
      sumExp <- SummarizedExperiment::SummarizedExperiment(assays = list(counts=featureCounts(dxd)), 
                                                           colData = anno, 
                                                           rowData = exonInfo)
      
      
      if (verbose) {
        end_parsing <- Sys.time()
        
        elapsed_parsing <- DOTSeq:::runtime(end_parsing, start_parsing)
        
        if (!is.null(elapsed_parsing$mins)) {
          message(sprintf(" - Data parsing runtime: %d mins %.3f secs", elapsed_parsing$mins, elapsed_parsing$secs))
        } else {
          message(sprintf(" - Data parsing runtime: %.3f secs", elapsed_parsing$secs))
        }
        
        message(" - Fit beta-binomial generalized linear models")
        start_dou <- Sys.time()
      }
      
      sumExp <- DOTSeq::fitDOU(object = sumExp,
                               formula = fmla,
                               target = target,
                               dispersion_modeling = dispersion_modeling,
                               dispformula = dispformula,
                               lrt = lrt,
                               diagnostic = diagnostic,
                               parallel = parallel,
                               optimizers = optimizers,
                               verbose = verbose)
      
      # L <- contrastMatrix(sumExp, fmla, baseline = baseline)
      
      
      if (verbose) {
        end_dou <- Sys.time()
        
        elapsed_dou <- DOTSeq:::runtime(end_dou, start_dou)
        
        if (!is.null(elapsed_parsing$mins)) {
          message(sprintf(" - DOU runtime: %d mins %.3f secs", elapsed_dou$mins, elapsed_dou$secs))
        } else {
          message(sprintf(" - DOU runtime: %.3f secs", elapsed_dou$secs))
        }
        
        message(" - Start differential translation efficiency analysis")
        start_dte <- Sys.time()
      }
      
      deseq_formula <- remove_random_effects(formula)
      deseq_formula <- DOTSeq:::reduce_formula(deseq_formula, combined_cond)
      
      if (verbose) {
        message(" - Design formula: ", deparse(deseq_formula))
      }
      
      # Create DESeqDataSet
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = dcounts,
                                       colData = anno,
                                       design = deseq_formula)
      
      # Run the DESeq2 analysis
      dds <- DESeq2::DESeq(dds)
      
      # # The interaction term directly tests for differential TE
      # results_name <- resultsNames(dds)[grep("condition.*strategy", resultsNames(dds))]
      # 
      # if (length(results_name) > 0) {
      #   dte_results <- DESeq2::results(dds, name = results_name)
      # } else {
      #   stop("Could not find interaction term in DESeq2 results. Check your design formula.")
      # }
      
      if (verbose) {
        end_dte <- Sys.time()
        
        elapsed_dte <- DOTSeq:::runtime(end_dte, start_dte)
        
        if (!is.null(elapsed_dte$mins)) {
          message(sprintf("\n - DTE runtime: %d mins %.3f secs", elapsed_dte$mins, elapsed_dte$secs))
        } else {
          message(sprintf(" - DTE runtime: %.3f secs", elapsed_dte$secs))
        }
        
        # DOU summary
        dou_res <- extract_results(sumExp)
        msg <- capture.output(print(table(dou_res$model_type)))
        msg <- msg[2:length(msg)]
        message(paste(" - DOU model fitting summary:\n  ", paste(msg, collapse = "\n")))
        
        # DTE summary
        dte_res <- results(dds)
        dte_success <- sum(!is.na(dte_res$padj))
        dte_total <- nrow(dte_res)
        dte_failed <- dte_total - dte_success
        
        msg_dte <- c(
          "nonTestable" = dte_failed,
          "DESeq2" = dte_success
        )
        
        msg_dte_str <- capture.output(print(as.table(msg_dte)))
        message(paste(" - DTE model fitting summary:\n", paste(msg_dte_str, collapse = "\n")))

      }
      
    } else {
      mismatches <- which(cnt_meta$key != gff_meta$key)
      print(colnames(data.frame(cnt_key = cnt_meta$key[mismatches], gff_key = gff_meta$key[mismatches])))
      stop("The orders of the count table and flatten GFF don't match!")
    }
    
  }
  
  return(list(
    raw_counts = dcounts,
    norm_counts = norm_counts,
    orfs = ORF_df,
    # te = te$absoluteTE,
    # occupancyShift = te$occupancyShift,
    dds = dds,
    sumExp = sumExp,
    dxd = dxd,
    formula = fmla
    # contrastMat = L
  ))
}



#' @title Remove random effects from a formula
#' @description This function takes an R formula object that may contain random effect terms and returns a new formula with all such terms removed.
#' @param formula An R formula object. Random effect terms are identified by the pattern `(1 | group)`.
#' @return A new R formula object containing only the fixed effect terms.
#' @examples
#' # Example formula with both fixed and random effects
#' my_formula <- ~ condition + strategy + (1 | group)
#' 
#' # Remove the random effects
#' new_formula <- remove_random_effects(my_formula)
#' 
#' # The result is a formula with only the fixed effects
#' print(new_formula)
#' # ~condition + strategy
#'
remove_random_effects <- function(formula) {
  # Convert the formula to a character string and remove the '~'
  formula_str <- as.character(formula)[2]
  # Split the string by the '+' sign to get individual terms
  terms <- unlist(strsplit(formula_str, "\\s*\\+\\s*"))
  # Identify and remove any terms that contain a random effect pattern
  clean_terms <- terms[!grepl("\\(1\\s*\\|.*?\\)", terms)]
  # Collapse the remaining terms back into a single formula string
  clean_formula_str <- paste(clean_terms, collapse = " + ")
  # Convert the string back to a formula object
  as.formula(paste("~", clean_formula_str))
}


#' Reduce a formula by removing terms with insufficient variation in the data
#'
#' This function takes a formula and a data frame, and returns a simplified formula
#' that includes only terms corresponding to variables with at least two levels.
#' Interaction terms using `*` or `:` are retained only if all involved variables are valid.
#' If the formula is reduced, a message is printed to inform the user.
#'
#' @param formula_input A formula object (e.g., `~ batch + condition * strategy`) or a character string representing a formula.
#' @param data A data frame containing the variables referenced in the formula.
#'
#' @return A reduced formula object that excludes terms with only one level in the data.
#' If no valid terms remain, the function will throw an error.
#'
#' @examples
#' df <- data.frame(
#'   condition = factor(c(0, 0, 1, 1)),
#'   strategy = factor(c("ribo", "ribo", "rna", "rna")),
#'   batch = factor(rep(1, 4)) # single-level factor
#' )
#' formula_input <- ~ batch + condition * strategy
#' reduce_formula(formula_input, df)
#' # Returns: ~ condition * strategy, with a message about reduction
reduce_formula <- function(formula_input, data) {
  # Convert to character if it's a formula
  formula_str <- if (inherits(formula_input, "formula")) {
    deparse(formula_input)
  } else if (is.character(formula_input)) {
    formula_input
  } else {
    stop("Input must be a formula or a character string.")
  }
  
  # Extract RHS of formula
  rhs <- strsplit(formula_str, "~")[[1]][2]
  terms <- strsplit(rhs, "\\+")[[1]]
  terms <- trimws(terms)
  
  # Identify valid terms
  valid_terms <- c()
  for (term in terms) {
    if (grepl("\\*|:", term)) {
      components <- trimws(unlist(strsplit(term, "\\*|:")))
      missing_vars <- setdiff(components, colnames(data))
      if (length(missing_vars) > 0) {
        stop(paste0(
          "Invalid formula: interaction term '", term, "' includes missing variable(s): ",
          paste(missing_vars, collapse = ", "), "."
        ))
      }
      valid_terms <- c(valid_terms, term)
    } else if (term %in% colnames(data)) {
      valid_terms <- c(valid_terms, term)
    }
  }
  
  if (length(valid_terms) == 0) {
    stop("No valid terms left in formula after reduction.")
  }
  
  new_formula_str <- paste("~", paste(valid_terms, collapse = " + "))
  reduced_formula <- as.formula(new_formula_str)
  
  if (new_formula_str != formula_str) {
    message(paste0(
      " - Formula has been reduced due to missing variables or terms with only one level.\n",
      "   Original formula: ", formula_str, "\n",
      "   Reduced formula:  ", new_formula_str
    ))
  }
  
  return(reduced_formula)
}


#' Format elapsed time between two timestamps
#'
#' @description
#' Computes the elapsed time between a start and end timestamp, returning a list
#' with minutes and seconds. If the duration is less than one minute, only seconds
#' are returned. Useful for logging runtimes in a human-readable format.
#'
#' @param end_time POSIXct. The end time of the process.
#' @param start_time POSIXct. The start time of the process.
#' @param units Character. Units for `difftime()` calculation. Default is `"mins"`.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{mins} (optional) — integer number of minutes
#'   \item \code{secs} — integer number of seconds
#' }
#'
#' @examples
#' start <- Sys.time()
#' Sys.sleep(3.5)
#' end <- Sys.time()
#' runtime(end, start)
runtime <- function(end_time, start_time, units = "mins") {
  elapsed <- as.numeric(difftime(end_time, start_time, units = units))
  
  total_seconds <- elapsed * 60
  mins <- total_seconds %/% 60
  secs <- total_seconds %% 60
  
  if (mins > 0) {
    return(list(mins = mins, secs = secs))
  } else {
    return(list(secs = secs))
  }
}


#' Extract Model Results into a Unified Data Frame
#'
#' @description
#' Extracts scalar results from a named list of model objects (typically from `glmmTMB` fits)
#' into a tidy data frame. Handles nested lists, missing or `NULL` models, and ensures
#' consistent columns across all entries.
#'
#' @param m A named list or SummarizedExperiment-like object containing fitted model objects,
#'   typically stored in `rowData(m$sumExp)[['fitDOUModels']]`.
#' @param verbose Logical. If \code{TRUE}, messages will be printed for skipped or failed models.
#'
#' @return A data frame where each row corresponds to an ORF and its associated model results.
#'   Columns include scalar parameters extracted from the model, \code{ORF_ID}, and \code{model_type}.
#'
#' @details
#' The function uses a recursive helper (`flatten_scalars`) to extract scalar values
#' from nested lists. It supports models of type \code{"glmmTMB"} and \code{"glmmTMB_joint"},
#' and gracefully handles \code{NULL} or unsupported model types by returning minimal rows.
#'
#' Missing columns across models are filled with \code{NA} to ensure a consistent structure.
#'
#' @import SummarizedExperiment
#' @export
extract_results <- function(sumExp, verbose = TRUE) {
  
  models_list <- rowData(sumExp)[["fitDOUModels"]]
  # Helper to recursively flatten and extract scalar values
  flatten_scalars <- function(x, prefix = NULL) {
    if (is.atomic(x) && length(x) == 1) {
      name <- if (is.null(prefix)) "value" else prefix
      return(setNames(list(x), name))
    } else if (is.list(x)) {
      result <- list()
      # Use `names(x)` and `seq_along` to handle unnamed elements gracefully
      list_names <- names(x)
      if (is.null(list_names)) {
        list_names <- paste0("elem", seq_along(x))
      }
      for (i in seq_along(x)) {
        n <- list_names[i]
        sub_prefix <- if (is.null(prefix)) n else paste0(prefix, ".", n)
        result <- c(result, flatten_scalars(x[[i]], sub_prefix))
      }
      return(result)
    } else {
      return(list()) # Return an empty list for non-atomic, non-list objects
    }
  }
  
  # Process a single model
  process_model <- function(model_obj, name) {
    # Check for NULL objects
    if (is.null(model_obj)) {
      if (verbose) message("Skipping ORF_ID ", name, " due to NULL model object.")
      return(data.frame(
        ORF_ID = name,
        model_type = "NULL", # Or "Failed", etc.
        stringsAsFactors = FALSE
      ))
    }
    
    model_type <- model_obj@type
    
    if ((model_type == "glmmTMB_joint") | (model_type == "glmmTMB")) {
      params <- model_obj@results
      param_values <- flatten_scalars(params)
      
      param_values$ORF_ID <- name
      param_values$model_type <- model_type
      
      df <- as.data.frame(param_values, stringsAsFactors = FALSE)
      
    } else {
      # Return minimal row with NA for failed or single-ORF models
      df <- data.frame(
        ORF_ID = name,
        model_type = model_type,
        stringsAsFactors = FALSE
      )
    }
    
    return(df)
  }
  
  # Apply to all models
  results_df_list <- lapply(names(models_list), function(name) {
    process_model(models_list[[name]], name)
  })
  
  # Get all unique column names
  all_cols <- unique(unlist(lapply(results_df_list, names)))
  
  # Fill missing columns with NA
  results_df_list_filled <- lapply(results_df_list, function(df) {
    missing <- setdiff(all_cols, names(df))
    for (col in missing) df[[col]] <- NA
    df[all_cols]
  })
  
  # Combine all rows
  results_df <- do.call(rbind, results_df_list_filled)
  rownames(results_df) <- NULL
  
  return(results_df)
}
