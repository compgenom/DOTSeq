#' Fit DOTSeq Differential Translation Model
#'
#' This function performs the complete DOTSeq analysis pipeline:
#' loading count data, aligning with sample metadata, filtering ORFs,
#' normalising counts, calculating translational efficiency (TE), and
#' fitting quasi-binomial generalised linear models for differential ORF
#' translation using \code{satuRn::fitDTU}.
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
#' @param sampleDelim Delimiter used in sample name (default: NULL).
#' @param batchCol String; name of the column in `conditionTable` that specifies batch assignments, e.g., "batch". 
#'    If NULL, batch effects are not modeled (default: NULL).
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
#' @return A named list with the following elements:
#' \describe{
#'   \item{rawCnts}{Raw counts matrix for all samples.}
#'   \item{normCnts}{Normalized counts matrix for all samples.}
#'   \item{orfs}{Data frame of ORFs derived from the BED file matched to the DOTSeq object.}
#'   \item{absoluteTE}{Matrix of translational efficiency values per ORF and sample.}
#'   \item{occupancyShift}{Matrix of log2-transformed ribo/rna proportions within genes.}
#'   \item{sumExp}{SummarizedExperiment object containing normalized counts and sample metadata.}
#'   \item{dxd}{DOTSeq object used for modeling exon/ORF-level counts.}
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
                   sampleDelim = NULL,
                   batchCol = NULL,
                   pseudoCnt = 1e-6, 
                   minCount = 1, 
                   stringent = NULL, 
                   parallel = FALSE,
                   verbose = FALSE) {
  
  # Define expected column names
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
  normaliseNamesnames <- function(x) tolower(trimws(x))
  
  if (is.character(conditionTable) && file.exists(conditionTable)) {
    # Read header line and normalize
    firstLine <- readLines(conditionTable, n = 1)
    condHeader <- normaliseNamesnames(strsplit(firstLine, "\t|,|\\s+")[[1]])
    
    missingCols <- setdiff(normaliseNamesnames(condCols), condHeader)
    
    if (length(missingCols) == 0) {
      cond <- read.table(conditionTable, header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
      names(cond) <- normaliseNamesnames(names(cond))
    } else {
      stop(paste("Missing expected columns (case-insensitive match):", paste(missingCols, collapse = ", ")))
    }
    
  } else if (is.data.frame(conditionTable)) {
    condHeader <- normaliseNamesnames(names(conditionTable))
    missingCols <- setdiff(normaliseNamesnames(condCols), condHeader)
    
    if (length(missingCols) == 0) {
      cond <- conditionTable
      names(cond) <- normaliseNamesnames(names(cond))
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
  rnaCond <- cond[cond$strategy=="rna",]
  
  # Combine and duplicate columns
  combinedCond <- rbind(rnaCond, riboCond)
  combinedCond <- combinedCond[, !(names(combinedCond) %in% c("run","strategy"))]
  
  # Specify the columns to move
  colsToFront <- c("condition", "replicate")
  # Get the remaining columns in their original order
  remainingCols <- setdiff(names(combinedCond), colsToFront)
  # Rearrange the data frame
  combinedCond <- combinedCond[, c(colsToFront, remainingCols)]
  
  numRNASmps <- nrow(rnaCond)
  numRiboSmps <- nrow(riboCond)
  
  # Extract batch column if present
  hasBatch <- !is.null(batchCol) && batchCol %in% colnames(combinedCond)
  if (hasBatch) {
    batch <- combinedCond[[batchCol]]
    combinedCond[[batchCol]] <- NULL
  }
  
  numCond <- ncol(combinedCond)
  
  # Duplicate condition columns
  combinedCond <- combinedCond[, rep(1:numCond, 2)]
  
  # Add modality column
  modality <- c(rep("mRNA", numRNASmps), rep("RIBO", numRiboSmps))
  combinedCond <- cbind(combinedCond[1:numCond], modality, combinedCond[(numCond+1):ncol(combinedCond)])
  
  # Standardise duplicated column names
  effectCols <- (numCond + 2):ncol(combinedCond)
  colnames(combinedCond)[effectCols] <- paste0("effect", seq_along(effectCols))
  
  # Overwrite RNA rows with first row values for EXTRA columns
  for (i in effectCols) {
    combinedCond[1:numRNASmps, i] <- combinedCond[1, i]
  }
  
  # Add batch back
  if (hasBatch) {
    combinedCond$batch <- batch
  }
  # Create formula dynamically
  extendedConds <- colnames(combinedCond)
  fmla <- as.formula(paste("~ 0 +", paste(extendedConds, collapse= "+")))
  if (verbose) {
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
      anno$modality <- as.factor(anno$modality)
      anno$effect1 <- ifelse(anno$modality == "RIBO", as.character(anno$condition), "none")
      anno$effect1 <- factor(anno$effect1)
      
      anno$batch <- ifelse(anno$modality == "RIBO", as.character(anno$condition), "none")
      anno$batch <- factor(anno$batch)
      
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
      
      L <- contrastMatrix(sumExp, fmla)
      
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
    absoluteTE = te$absoluteTE,
    occupancyShift = te$occupancyShift,
    sumExp = sumExp,
    dxd = dxd,
    formula = fmla,
    contrastMat = L
  ))
}
