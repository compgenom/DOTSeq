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


#' Fit DOTSeq Differential Translation Models
#'
#' This function performs the complete DOTSeq analysis pipeline:
#' loading count data, aligning with sample metadata, filtering ORFs,
#' normalizing counts, calculating translational efficiency (TE), and
#' fitting beta-binomial and negative binomial generalized linear models for 
#' differential ORF translation using \code{DOTSeq::fitDOU} and \code{DESeq2::DESeq}.
#'
#' @param countTable Either a path to a count table file or a data frame.
#'   Must contain columns: Geneid, Chr, Start, End, Strand, Length,
#'   plus one column per sample.
#' @param conditionTable Either a path to a condition metadata file or a data frame.
#'   Must include columns: run, strategy, condition, replicate.
#' @param flattenedFile Optional path to a flattened GFF/GTF file containing exon definitions.
#' @param bed Path to a BED file with ORF annotations.
#' @param rnaSuffix Character suffix to identify RNA-seq columns (default: \code{".rna"}).
#' @param riboSuffix Character suffix to identify Ribo-seq columns (default: \code{".ribo"}).
#' @param target Character string specifying the non-reference condition level to extract the corresponding interaction term from the model. 
#' This is contrasted against the baseline condition (default: \code{NULL}).
#' @param baseline Character string specifying the desired reference level.
#' @param formula A formula object specifying the design, e.g., \code{~ condition * strategy}.
#' @param batchCol String; name of the column in \code{conditionTable} that specifies batch assignments, 
#'   e.g., \code{"batch"}. If \code{NULL}, batch effects are not modeled (default: \code{NULL}).
#' @param pseudoCnt Numeric pseudo-count to avoid division by zero when computing TE (default: \code{1e-6}).
#' @param minCount Minimum count threshold for filtering ORFs (default: \code{1}).
#' @param stringent Logical or \code{NULL}; determines the filtering strategy:
#'   \describe{
#'     \item{\code{TRUE}}{Keep ORFs where all replicates in at least one condition pass \code{minCount}.}
#'     \item{\code{FALSE}}{Keep ORFs where all replicates in at least one condition-strategy group pass \code{minCount}.}
#'     \item{\code{NULL}}{Keep ORFs where total counts across replicates pass \code{minCount}.}
#'   }
#' @param parallel A list passed to \code{glmmTMBControl} to configure parallel optimization.
#'   If \code{NULL}, parallelism is disabled (default: \code{list(n = 4L, autopar = TRUE)}).
#' @param optimizers Logical; if \code{TRUE}, enables brute-force optimization using multiple optimizers
#'   in \code{glmmTMB}: \code{nlminb}, \code{bobyqa}, and \code{optim} (default: \code{FALSE}).
#' @param dispersionStrategy Optional string specifying the dispersion modeling strategy.
#'   Used to select between default, constant, or custom dispersion models.
#' @param dispformula Optional formula object for custom dispersion modeling.
#' @param diagnostic Logical; if \code{TRUE}, enables model diagnostics including tests for overdispersion,
#'   zero inflation, and residual properties (default: \code{FALSE}).
#' @param verbose Logical; if \code{TRUE}, prints progress messages (default: \code{FALSE}).
#'
#' @return A named \code{list} with the following elements:
#' \describe{
#'   \item{rawCnts}{Raw counts matrix for all samples.}
#'   \item{normCnts}{Normalized counts matrix for all samples.}
#'   \item{orfs}{Data frame of ORFs derived from the BED file matched to the DOTSeq object.}
#'   \item{dds}{DESeq2 object used for modeling differential gene translation.}
#'   \item{sumExp}{SummarizedExperiment object containing normalized counts and sample metadata.}
#'   \item{dxd}{DOTSeq object used for modeling exon/ORF-level counts.}
#'   \item{formula}{Design formula used in \code{DOTSeq::fitDOU}.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- fitDOT(
#'   countTable = "counts.txt",
#'   conditionTable = "samples.txt",
#'   flattenedFile = "flattened.gff",
#'   bed = "orfs.bed",
#'   rnaSuffix = ".rna",
#'   riboSuffix = ".ribo",
#'   pseudoCnt = 1e-6,
#'   minCount = 1,
#'   stringent = TRUE,
#'   verbose = TRUE
#' )
#' head(result$sumExp)
#' }
#'
#' @export
fitDOT <- function(countTable, conditionTable, 
                   flattenedFile, bed, 
                   rnaSuffix = ".rna", 
                   riboSuffix = ".ribo", 
                   target = NULL,
                   baseline = NULL,
                   formula = ~ condition * strategy,
                   dispersionStrategy = "strategy",
                   dispformula = NULL,
                   diagnostic = FALSE,
                   # sampleDelim = NULL,
                   batchCol = NULL,
                   pseudoCnt = 1e-6, 
                   minCount = 1, 
                   stringent = TRUE, 
                   parallel = list(n=4L, autopar=TRUE),
                   optimizers = FALSE,
                   verbose = FALSE) {
  
  if (verbose) {
    start_parsing <- Sys.time()
    message(" - Use ", parallel$n, " threads")
  }
  
  cntCols <- c("Geneid", "Chr", "Start",  "End", "Strand", "Length")
  
  if (is.character(countTable) && file.exists(countTable)) {
    # If countTable is a file path (character) and the file exists
    # Read the first line to get column names
    firstLine <- readLines(countTable, n = 1)
    cntHeader <- strsplit(firstLine, "\t|,|\\s+")[[1]]
    
    # Check if all expected columns are present
    if (all(cntCols %in% cntHeader)) {
      cond <- read.table(countTable, header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
    } else {
      stop("Header must contain the following columns: Geneid, Chr, Start, End, Strand, Length")
    }
    
  } else if (is.data.frame(countTable) && all(cntCols %in% names(countTable))) {
    # If countTable is already a data frame
    cnt <- countTable
  } else {
    stop("countTable must be either a valid file path or a data frame.")
  }
  
  # Define expected column names
  condCols <- c("run", "strategy", "condition", "replicate")
  
  # Helper to normalize column names
  normaliseNames <- function(x) tolower(trimws(x))
  
  if (is.character(conditionTable) && file.exists(conditionTable)) {
    # Read header line and normalize
    firstLine <- readLines(conditionTable, n = 1)
    condHeader <- normaliseNames(strsplit(firstLine, "\t|,|\\s+")[[1]])
    
    missingCols <- setdiff(normaliseNames(condCols), condHeader)
    
    if (length(missingCols) == 0) {
      cond <- read.table(conditionTable, header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
      names(cond) <- normaliseNames(names(cond))
    } else {
      stop(paste("Missing expected columns (case-insensitive match):", paste(missingCols, collapse = ", ")))
    }
    
  } else if (is.data.frame(conditionTable)) {
    condHeader <- normaliseNames(names(conditionTable))
    missingCols <- setdiff(normaliseNames(condCols), condHeader)
    
    if (length(missingCols) == 0) {
      cond <- conditionTable
      names(cond) <- normaliseNames(names(cond))
    } else {
      stop(paste("Data frame is missing expected columns (case-insensitive match):", paste(missingCols, collapse = ", ")))
    }
    
  } else {
    stop("`conditionTable` must be either a valid file path or a data frame.")
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
  # cnt <- cbind(cnt[names(cnt) %in% cntCols], cnt_aligned)
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
    stop("Number of samples in countTable and conditionTable doesn't match.")
  }
  
  riboCond <- cond[cond$strategy=="ribo",]
  riboCond$strategy <- 1
  rnaCond <- cond[cond$strategy=="rna",]
  rnaCond$strategy <- 0
  
  # Combine and duplicate columns
  combinedCond <- rbind(rnaCond, riboCond)
  
  combinedCond <- combinedCond[, !(names(combinedCond) %in% "run")]
  
  # Robust batch column detection if batchCol is NULL
  if (is.null(batchCol)) {
    possibleBatchCols <- grep("^batch$", names(combinedCond), ignore.case = TRUE, value = TRUE)
    if (length(possibleBatchCols) > 0) {
      batchCol <- possibleBatchCols[1]
      if (verbose) message(" - Auto-detected batch column: ", batchCol)
    }
  } else if (batchCol==FALSE) {
    possibleBatchCols <- grep("^batch$", names(combinedCond), ignore.case = TRUE, value = TRUE)
    if (length(possibleBatchCols) > 0) {
      batchCol <- possibleBatchCols[1]
      combinedCond[[batchCol]] <- NULL
      if (verbose) message(" - Auto-removed batch column: ", batchCol)
    }
  }
  
  # Convert to factors
  for (col in colnames(combinedCond)) {
    if (is.character(combinedCond[[col]])) {
      combinedCond[[col]] <- factor(combinedCond[[col]])
    }
  }
  
  # Set baseline
  if (!is.null(baseline)) {
    combinedCond$condition <- relevel(combinedCond$condition, ref = baseline)
  }
  
  # Remove columns with fewer than 2 levels
  validCols <- sapply(combinedCond, function(x) !(is.factor(x) && length(unique(x)) < 2))
  combinedCond <- combinedCond[, validCols]
  
  fmla <- as.formula(formula)
  
  if (verbose) {
    message(" - Start differential ORF translation analysis")
    message(" - Design formula: ", deparse(fmla))
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
  dcounts <- dcounts[, rownames(combinedCond)]
  
  ## get genes and exon names out
  splitted <- strsplit(rownames(dcounts), ":")
  exons <- sapply(splitted, "[[", 2)
  genesrle <- sapply(splitted, "[[", 1)
  
  
  if (!is.null(flattenedFile)) {
    if (verbose) {
      message(" - Parse the flattened GFF file")
    }
    aggregates <- read.delim(flattenedFile, stringsAsFactors = FALSE, 
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
      for (col in colnames(combinedCond)) {
        if (is.numeric(combinedCond[[col]]) && all(combinedCond[[col]] == as.integer(combinedCond[[col]]))) {
          combinedCond[[col]] <- factor(combinedCond[[col]])
        } else if (is.character(combinedCond[[col]])) {
          combinedCond[[col]] <- factor(combinedCond[[col]])
        }
      }
      
      design = ~ sample + exon + condition:exon
      dxd <- DEXSeqDataSet(dcounts, combinedCond, design, exons, 
                           genesrle, exoninfo[matching], transcripts[matching])
      
      # Get ORF annotation
      orfDf <- parseBed(bed, dxd)
      names(orfDf)[names(orfDf) == "exonBaseMean"] <- "orfBaseMean"
      names(orfDf)[names(orfDf) == "exonBaseVar"] <- "orfBaseVar"
      
      counts <- counts(dxd)
      
      if (stringent == TRUE) {
        keep_list <- lapply(unique(cond$condition), function(cond_name) {
          samples <- rownames(cond[cond$condition == cond_name, ])
          rowSums(dcounts[, samples, drop = FALSE] >= minCount) == length(samples)
        })
        
        if (verbose) {
          message(" - Keep ORFs where all replicates in at least one condition have counts >= ", minCount)
        }
        
        keep <- Reduce("|", keep_list)
        dxd <- dxd[keep, ]
      } else if (stringent == FALSE) {
        cond$group <- paste(cond$condition, cond$strategy, sep = "_")
        keep_list <- lapply(unique(cond$group), function(group_name) {
          samples <- rownames(cond[cond$group == group_name, ])
          rowSums(dcounts[, samples, drop = FALSE] >= minCount) == length(samples)
        })
        
        if (verbose) {
          message(" - Keep ORFs where all replicates in at least one condition-strategy group have counts >= ", minCount)
        }
        
        keep <- Reduce("|", keep_list)
        dxd <- dxd[keep, ]
        
      } else if (is.null(stringent)) {
        if (verbose) {
          message(" - Keep ORFs where the total counts across replicates are >= ", minCount)
        }
        dxd <- dxd[rowSums(featureCounts(dxd)) >= minCount,]
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
      normCnts <- counts(dxd, normalized = TRUE)
      colnames(normCnts) <- colData(dxd)$sample
      normCnts <- normCnts[, seq_len(ncol(normCnts)/2)]
      
      # te <- calculateTE(normCnts, sampleDelim = sampleDelim, rnaSuffix = rnaSuffix, riboSuffix = riboSuffix, pseudoCnt = pseudoCnt)
      
      # Run DOTSeq::fitDOU
      exonInfo <- rowData(dxd)
      colnames(exonInfo)[1:2] <- c("isoform_id", "gene_id")
      exonInfo$isoform_id <- rownames(exonInfo)
      
      # Get sampleAnnotation and set effect1 and batch to none for RNA-seq samples
      anno <- sampleAnnotation(dxd)
      anno$sizeFactor <- NULL
      
      sumExp <- SummarizedExperiment::SummarizedExperiment(assays = list(counts=featureCounts(dxd)), 
                                                           colData = anno, 
                                                           rowData = exonInfo)
      
      
      if (verbose) {
        end_parsing <- Sys.time()
        elapsed_parsing <- as.numeric(difftime(end_parsing, start_parsing, units = "mins"))
        message(" - Data parsing runtime: ", elapsed_parsing)
        
        message(" - Fit beta-binomial generalised linear models")
        start_dou <- Sys.time()
      }
      sumExp <- DOTSeq::fitDOU(object = sumExp,
                               formula = fmla,
                               target = target,
                               dispersionStrategy= dispersionStrategy,
                               dispformula = dispformula,
                               diagnostic = diagnostic,
                               parallel = parallel,
                               optimizers = optimizers,
                               verbose = verbose)
      
      # L <- contrastMatrix(sumExp, fmla, baseline = baseline)
      
      if (verbose) {
        end_dou <- Sys.time()
        elapsed_dou <- as.numeric(difftime(end_dou, start_dou, units = "mins"))
        message(" - DOU runtime: ", elapsed_dou)
        
        message(" - Start differential gene translation analysis")
        start_dot <- Sys.time()
      }
      
      deseq_formula <- remove_random_effects(formula)
      
      if (verbose) {
        message(" - Design formula: ", deparse(deseq_formula))
      }
      
      # Try creating DESeqDataSet with full formula
      dds <- tryCatch({
        DESeq2::DESeqDataSetFromMatrix(countData = dcounts,
                                       colData = anno,
                                       design = deseq_formula)
      }, error = function(e) {
        message(" - Failed with full formula: ", e$message)
        
        # Remove 'batch' from formula and retry
        reduced_formula_str <- gsub("(?<=\\+|\\~|^)\\s*batch\\s*\\+?|\\+?\\s*batch(?=\\+|$)", "", deseq_formula_str[2], perl = TRUE)
        reduced_formula_str <- gsub("^\\s*\\+\\s*|\\s*\\+\\s*$", "", reduced_formula_str)  # clean up leading/trailing '+'
        reduced_formula <- as.formula(paste(deseq_formula_str[1], reduced_formula_str))
        
        if (verbose) {
          message(" - Retrying with reduced formula: ", deparse(reduced_formula))
        }
        
        DESeq2::DESeqDataSetFromMatrix(countData = dcounts,
                                       colData = anno,
                                       design = reduced_formula)
      })
      
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
        end_dot <- Sys.time()
        elapsed_dot <- as.numeric(difftime(end_dot, start_dot, units = "mins"))
        message(" - DOT runtime: ", elapsed_dot)
      }
      
    } else {
      mismatches <- which(cnt_meta$key != gff_meta$key)
      print(colnames(data.frame(cnt_key = cnt_meta$key[mismatches], gff_key = gff_meta$key[mismatches])))
      stop("The orders of the count table and flatten GFF don't match!")
    }
    
  }
  
  return(list(
    rawCnts = dcounts,
    normCnts = normCnts,
    orfs = orfDf,
    # te = te$absoluteTE,
    # occupancyShift = te$occupancyShift,
    dds = dds,
    sumExp = sumExp,
    dxd = dxd,
    formula = fmla
    # contrastMat = L
  ))
}




