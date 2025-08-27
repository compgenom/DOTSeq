#' Fit DOTSeq Differential Translation Models
#'
#' This function performs the complete DOTSeq analysis pipeline:
#' loading count data, aligning with sample metadata, filtering ORFs,
#' normalising counts, calculating translational efficiency (TE), and
#' fitting quasi-binomial and negative binomial generalised linear models for 
#' differential ORF translation using \code{satuRn::fitDTU} and \code{DESeq2::DESeq}.
#'
#' @param countTable Either a path to a count table file or a data frame.
#'   The file/data frame must have columns: Geneid, Chr, Start, End, Strand, Length,
#'   plus one column per sample.
#' @param conditionTable Either a path to a condition metadata file or a data frame.
#'   Must include columns: run, strategy, condition, replicate.
#' @param flattenedFile Optional path to a flattened GFF/GTF file containing exon definitions.
#' @param bed Path to a BED file with ORF annotations.
#' @param rnaSuffix Character suffix to identify RNA-seq columns (default: \code{".rna"}).
#' @param riboSuffix Character suffix to identify Ribo-seq columns (default: \code{".ribo"}).
#' @param formula A formula object specifying the design, e.g., ~ 0 + replicate + condition * strategy.
#' @param baseline Optional character specifying the baseline condition for contrasts.
#'   If NULL, the function attempts to infer it automatically (default: NULL).
#' @param sampleDelim Delimiter used in sample name. Use "." if sample IDs are unique, 
#'    e.g. SRR, DRR, or ERR run accessions (default: NULL).
#' @param batchCol String; name of the column in `conditionTable` that specifies batch assignments, 
#'    e.g., "batch". If NULL, batch effects are not modeled (default: NULL).
#' @param pseudoCnt Numeric pseudo-count to avoid division by zero when computing TE (default: 1e-6).
#' @param minCount Minimum count threshold for filtering ORFs (default: 1).
#' @param stringent Logical or NULL; determines the filtering strategy:
#'   - TRUE: keep ORFs where all replicates in at least one condition pass \code{minCount}.
#'   - FALSE: keep ORFs where all replicates in at least one condition-strategy group pass \code{minCount}.
#'   - NULL: keep ORFs where total counts across replicates pass \code{minCount}.
#' @param seed An optional integer. Use to set the seed for the random number
#'   generator to ensure reproducible results (default: NULL).
#' @param parallel A logical value indicating whether to use parallel
#'   computation. If `TRUE`, the function will distribute tasks across
#'   multiple cores as configured by the `BiocParallel` package. If `FALSE`,
#'   the function will run sequentially on a single core (default: FALSE).
#' @param verbose Logical; if TRUE, prints progress messages (default: FALSE).
#'
#' @return A named \code{list} with the following elements:
#' \describe{
#'   \item{rawCnts}{Raw counts matrix for all samples.}
#'   \item{normCnts}{Normalized counts matrix for all samples.}
#'   \item{orfs}{Data frame of ORFs derived from the BED file matched to the DOTSeq object.}
#'   \item{absoluteTE}{Matrix of translational efficiency values per ORF and sample.}
#'   \item{occupancyShift}{Matrix of log2-transformed ribo/rna proportions within genes.}
#'   \item{dds}{DESeq2 object usd for modelling differential gene translation analysis.}
#'   \item{sumExp}{SummarizedExperiment object containing normalized counts and sample metadata.}
#'   \item{dxd}{DOTSeq object used for modelling exon/ORF-level counts.}
#'   \item{formula}{Design formula used in \code{satuRn::fitDTU}.}
#'   \item{contrastMat}{Contrast matrix used for differential ORF translation testing.}
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
#'   seed = 42,
#'   verbose = TRUE
#' )
#' head(result$absoluteTE)
#' head(result$occupancyShift)
#' }
#'
#' @export
fitDOT <- function(countTable, conditionTable, 
                   flattenedFile, bed, 
                   rnaSuffix = ".rna", 
                   riboSuffix = ".ribo", 
                   baseline = NULL,
                   formula = NULL,
                   sampleDelim = NULL,
                   batchCol = NULL,
                   pseudoCnt = 1e-6, 
                   minCount = 1, 
                   stringent = NULL, 
                   parallel = FALSE,
                   seed = NULL,
                   verbose = FALSE) {
  
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
      if (verbose) cat(" - Auto-detected batch column:", batchCol, "\n")
    }
  } else if (batchCol==FALSE) {
    possibleBatchCols <- grep("^batch$", names(combinedCond), ignore.case = TRUE, value = TRUE)
    if (length(possibleBatchCols) > 0) {
      batchCol <- possibleBatchCols[1]
      combinedCond[[batchCol]] <- NULL
      if (verbose) cat(" - Auto-removed batch column:", batchCol, "\n")
    }
  }
  
  # Convert to factors
  for (col in colnames(combinedCond)) {
    if (is.character(combinedCond[[col]])) {
      combinedCond[[col]] <- factor(combinedCond[[col]])
    }
  }
  
  # Remove columns with fewer than 2 levels
  validCols <- sapply(combinedCond, function(x) !(is.factor(x) && length(unique(x)) < 2))
  combinedCond <- combinedCond[, validCols]
  
  fmla <- as.formula(formula)
  
  if (verbose) {
    cat(" - Start differential ORF translation analysis\n")
    cat(" - Design formula:", deparse(fmla), "\n")
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
      cat(" - Parse the flattened GFF file\n")
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
    exoninfo <- GRanges(as.character(aggregates$chr), IRanges(start = aggregates$start, 
                                                              end = aggregates$end), strand = aggregates$strand)
    names(exoninfo) <- paste(aggregates$gene_id, exonids, 
                             sep = ":O")
    
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
        cat(" - Create a DEXSeq object\n")
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
      
      counts <- counts(dxd)
      
      if (stringent == TRUE) {
        keep_list <- lapply(unique(cond$condition), function(cond_name) {
          samples <- rownames(cond[cond$condition == cond_name, ])
          rowSums(dcounts[, samples, drop = FALSE] >= minCount) == length(samples)
        })
        
        if (verbose) {
          cat(" - Keep ORFs where all replicates in at least one condition have counts >=", minCount, "\n")
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
          cat(" - Keep ORFs where all replicates in at least one condition-strategy group have counts >=", minCount, "\n")
        }
        
        keep <- Reduce("|", keep_list)
        dxd <- dxd[keep, ]
        
      } else if (is.null(stringent)) {
        if (verbose) {
          cat(" - Keep ORFs where the total counts across replicates are >=", minCount, "\n")
        }
        dxd <- dxd[rowSums(featureCounts(dxd)) >= minCount,]
      }
      
      # Find the gene IDs that have only one ORF "(exon)"
      singlets <- names(which(table(rowData(dxd)$groupID) == 1))
      if (verbose) {
        cat(" - Filter out", length(singlets), "single ORF genes\n")
      }
      # `! ... %in% ...` means "not in this list"
      keep_exons <- !rowData(dxd)$groupID %in% singlets
      dxd <- dxd[keep_exons, ]
      if (verbose) {
        cat(" - Number of ORFs passing filter:", nrow(featureCounts(dxd)), "\n")
      }
      
      # Get normalised counts
      dxd <- estimateSizeFactors(dxd)
      normCnts <- counts(dxd, normalized = TRUE)
      colnames(normCnts) <- colData(dxd)$sample
      normCnts <- normCnts[, seq_len(ncol(normCnts)/2)]
      
      orfDf <- parseBed(bed, dxd)
      names(orfDf)[names(orfDf) == "exonBaseMean"] <- "orfBaseMean"
      names(orfDf)[names(orfDf) == "exonBaseVar"] <- "orfBaseVar"
      
      te <- calculateTE(normCnts, sampleDelim = sampleDelim, rnaSuffix = rnaSuffix, riboSuffix = riboSuffix, pseudoCnt = pseudoCnt)
      
      # Run satuRn::fitDTU
      exonInfo <- rowData(dxd)
      colnames(exonInfo)[1:2] <- c("isoform_id", "gene_id")
      exonInfo$isoform_id <- rownames(exonInfo)
      
      # Get sampleAnnotation and set effect1 and batch to none for RNA-seq samples
      anno <- sampleAnnotation(dxd)
      
      sumExp <- SummarizedExperiment::SummarizedExperiment(assays = list(counts=featureCounts(dxd)), 
                                                           colData = anno, 
                                                           rowData = exonInfo)
      
      if (verbose) {
        cat(" - Fit quasi-binomial generalised linear models\n")
      }
      sumExp <- satuRn::fitDTU(object = sumExp,
                               formula = fmla,
                               parallel = parallel,
                               BPPARAM = BiocParallel::bpparam(RNGseed = seed),
                               verbose = TRUE)
      
      L <- contrastMatrix(sumExp, fmla, baseline = baseline)
      
      if (verbose) {
        cat(" - Start differential gene translation analysis\n")
      }
      
      # Convert formula to string and remove '0 +' or '-1 +' if present
      deseq_formula_str <- as.character(formula)
      if (grepl("(-1|0)\\s*\\+", deseq_formula_str[2])) {
        deseq_formula_str[2] <- gsub("(-1|0)\\s*\\+\\s*", "", deseq_formula_str[2])
      }
      deseq_formula <- as.formula(paste(deseq_formula_str[1], deseq_formula_str[2]))
      
      if (verbose) {
        cat(" - Design formula:", deparse(deseq_formula), "\n")
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
          cat(" - Retrying with reduced formula:", deparse(reduced_formula), "\n")
        }
        
        DESeq2::DESeqDataSetFromMatrix(countData = dcounts,
                                       colData = anno,
                                       design = reduced_formula)
      })
      
      # Run the DESeq2 analysis
      dds <- DESeq2::DESeq(dds, parallel = parallel)
      
      # # The interaction term directly tests for differential TE
      # results_name <- resultsNames(dds)[grep("condition.*strategy", resultsNames(dds))]
      # 
      # if (length(results_name) > 0) {
      #   te_results <- DESeq2::results(dds, name = results_name)
      # } else {
      #   stop("Could not find interaction term in DESeq2 results. Check your design formula.")
      # }
      
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
    te = te$absoluteTE,
    occupancyShift = te$occupancyShift,
    dds = dds,
    sumExp = sumExp,
    dxd = dxd,
    formula = fmla, 
    contrastMat = L
  ))
}



#' Plot and Extract Dispersion Estimates from DOU Models
#'
#' This function extracts and compares original and posterior dispersion estimates
#' from DOU (Differential ORF Usage) models fitted using `DOTSeq::fitDOT`. It also computes
#' mean proportions of ORF-level counts relative to their parent gene, and visualizes
#' the relationship between mean proportion and dispersion. Additionally, it plots
#' dispersion estimates from the DESeq2 model used in the DTE (Differential Translation Efficiency) analysis.
#'
#' @param m A list-like object returned by `DOTSeq::fitDOT`, containing:
#'   \itemize{
#'     \item \code{sumExp}: A `SummarizedExperiment` object with satuRn DTU model fits in the `rowData` under `"fitDTUModels"`, and a `"counts"` assay.
#'     \item \code{dds}: A `DESeqDataSet` object from DESeq2 used for DTE analysis.
#'   }
#'
#' @return A list containing:
#' \describe{
#'   \item{meanProportions}{A numeric vector of mean proportions of ORF counts relative to their parent gene.}
#'   \item{originalDispersions}{A numeric vector of dispersion estimates from the original GLM fit.}
#'   \item{posteriorDispersions}{A numeric vector of posterior (empirical Bayes regularized) dispersion estimates.}
#' }
#'
#' @details
#' The function assumes that row names of the count matrix follow the format `"geneID:orfID"`.
#' It calculates the proportion of each ORF's counts relative to the total counts of its parent gene,
#' and plots both original and posterior dispersion estimates against these mean proportions.
#'
#' The posterior dispersions are generally more stable and recommended for downstream hypothesis testing.
#' The function also calls `plotDispEsts()` to visualize dispersion estimates from the DESeq2 model.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom stats lowess
#' @importFrom graphics plot points lines legend
#' @importFrom DESeq2 plotDispEsts
#'
meanDispersion <- function(m) {
  # Access the list-column of StatModel objects
  models <- rowData(m$sumExp)[["fitDTUModels"]]
  # Extract the dispersion estimates from the original GLM fit (before empirical Bayes squeezing)
  original_dispersions <- vapply(models, function(x) x@params$dispersion, FUN.VALUE = numeric(1))
  # Extract the posterior dispersion estimates (after empirical Bayes squeezing)
  posterior_dispersions <- vapply(models, function(x) x@varPosterior, FUN.VALUE = numeric(1))
  # The posterior dispersion is the regularized, or "squeezed," dispersion value,
  # which is generally more stable and recommended for downstream analysis and hypothesis testing.
  # The original dispersion is the raw value estimated from the GLM for each ORF.
  
  counts <- SummarizedExperiment::assay(m$sumExp, "counts")
  counts <- as.data.frame(counts)
  counts$groupID <- sapply(strsplit(rownames(counts), ":"), `[`, 1)
  groupTotal <- rowsum(counts[, -ncol(counts)], group = counts$groupID)
  groupTotal$groupID <- rownames(groupTotal)
  
  # Match groupID and divide row-wise
  # Get indices of matching groupIDs
  idx <- match(counts$groupID, groupTotal$groupID)
  # Calculate proportions for numeric columns only
  numeric_cols <- setdiff(names(counts), "groupID")
  props <- counts[numeric_cols] / groupTotal[idx, numeric_cols]
  # # Combine with groupID if needed
  # result <- cbind(groupID = counts$groupID, proportions)
  mean_props <- rowMeans(props, na.rm = TRUE)
  
  # pdf("mean_dispersion-dou.pdf", 4, 4.5)
  plot(x = mean_props, y = original_dispersions, col = "black",  
       pch = 16, cex = 0.3, log = "xy", 
       xlab = "Mean proportion", ylab = "Dispersion", 
       main = "Mean-Dispersion Plot for DOU")
  points(x = mean_props, y = posterior_dispersions, col = "#6495ED",  
         pch = 16, cex = 0.3) # , xlab = "Mean proportion", ylab = "Posterior dispersion")
  lines(lowess(mean_props, original_dispersions), col = "red", lwd = 2)
  legend(
    "bottomleft",
    legend = c("Original dispersions", "Posterior dispersions", "LOWESS fit (original)"),
    col = c("black", "#6495ED", "red"),
    pch = c(16, 16, NA),
    lty = c(NA, NA, 1),
    lwd = c(NA, NA, 2),
    pt.cex = 0.5,
    bty = "n"
  )
  # dev.off()
  
  # Plot mean-dispersion plot using the DESeq2 object
  DESeq2::plotDispEsts(m$dds)
  
  return(list(meanProportions = mean_props, 
              originalDispersions = original_dispersions, 
              posteriorDispersions = posterior_dispersions))
}
