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
#' Filters one or more BAM files to retain reads overlapping exonic regions
#' defined in a TxDb object stored in the metadata of a GRanges object.
#' Optionally restricts filtering to coding genes only (genes with CDS).
#' Filtered BAMs are sorted and saved with the suffix `.exonic.sorted.bam`
#' in a user-specified or temporary directory.
#'
#' @param gr A \code{GRanges} object with a \code{TxDb} SQLite file path stored
#' in its metadata under \code{metadata(gr)$txdb}.
#' @param seqlevels_style Character; the naming style for chromosome 
#' identifiers (e.g., \code{"UCSC"}, \code{"NCBI"}). This is applied to both 
#' the GRanges object and the TxDb annotation. Default is \code{"UCSC"}.
#' @param bam_files A character vector of paths to BAM files to be filtered.
#' @param bam_output_dir A writable directory where filtered BAM files will be
#' saved. Defaults to \code{tempdir()}.
#' @param coding_genes_only Logical; if \code{TRUE}, restrict filtering to
#' coding genes only (i.e., genes with CDS). Default is \code{TRUE}.
#' @param verbose Logical; if \code{TRUE}, print progress and runtime messages.
#' Default is \code{TRUE}.
#'
#' @return This function is called for its side effect of creating filtered and
#' sorted BAM files in \code{bam_output_dir}. It returns \code{NULL} invisibly.
#'
#' @details
#' The function uses \code{Rsamtools::filterBam()} to extract reads overlapping
#' exonic regions and \code{Rsamtools::sortBam()} to sort the filtered BAM.
#' The output files are named based on the input BAM file with the suffix
#' \code{.exonic.sorted.bam}.
#'
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom AnnotationDbi loadDb
#' @importFrom GenomicFeatures cdsBy exonsBy
#' @importFrom BiocGenerics unlist
#' @importFrom Rsamtools scanBamHeader BamFile filterBam ScanBamParam
#' @importFrom Rsamtools indexBam sortBam
#' @importFrom GenomeInfoDb keepSeqlevels
#'
#' @export
#'
#' @examplesIf requireNamespace("TxDb.Dmelanogaster.UCSC.dm3.ensGene", quietly = TRUE) && requireNamespace("pasillaBamSubset", quietly = TRUE) && requireNamespace("withr", quietly = TRUE)
#' library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
#' library(pasillaBamSubset)
#' library(AnnotationDbi)
#' library(GenomeInfoDb)
#' library(withr)
#'
#' # Save a subset of TxDb as an SQLite file
#' txdb_chr4 <- keepSeqlevels(
#'     TxDb.Dmelanogaster.UCSC.dm3.ensGene,
#'     "chr4",
#'     pruning.mode = "coarse"
#' )
#' txdb_path <- file.path(tempdir(), "dm3_chr4.sqlite")
#' saveDb(txdb_chr4, file = txdb_path)
#'
#' # Create a GRanges object with a link to the TxDb SQLite file
#' gr <- GRanges(seqnames = "chr4", ranges = IRanges(start = 233, end = 2300))
#' metadata(gr)$txdb <- txdb_path
#'
#' # Filter BAM file and save output in a temporary directory
#' temp_dir <- tempdir()
#' getExonicReads(gr,
#'     bam_files = untreated1_chr4(),
#'     bam_output_dir = temp_dir
#' )
#'
#' # Clean up
#' withr::defer(unlink(txdb_path))
#' withr::defer(unlink(list.files(temp_dir, pattern = "exonic", full.names = TRUE)))
#' 
getExonicReads <- function(gr, seqlevels_style = "UCSC", bam_files, bam_output_dir = tempdir(), coding_genes_only = TRUE, verbose = TRUE) {
    
    if (!file.exists(metadata(gr)$txdb)) {
        stop("The TxDb object has been removed. Please rerun getORFs().")
    }
    
    # Filter annotations to coding genes
    annotation <- loadDb(metadata(gr)$txdb)
    
    tryCatch(
        {
            seqlevelsStyle(gr) <- seqlevels_style
            seqlevelsStyle(annotation) <- seqlevels_style
        }, 
        error = function(e) {
            warning("Failed to set seqlevels style: ", conditionMessage(e))
        }
    )
    
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
        
        bam_header <- scanBamHeader(BamFile(bam))
        bam_seqlevels <- names(bam_header[[1]])
        
        # Get seqlevels from annotation
        annotation_seqlevels <- seqlevels(exons_by_genes)
        
        # Find common seqlevels
        common_seqlevels <- intersect(annotation_seqlevels, bam_seqlevels)
        
        # Keep only common seqlevels in annotation
        exons_for_bam <- keepSeqlevels(exons_by_genes, common_seqlevels, pruning.mode = "coarse")
        
        if (length(exons_for_bam) == 0) {
            warning(
                "No overlapping chromosomes found for BAM: ", bam,
                "\nCheck seqlevels_style: ", seqlevels_style
            )
            next
        }
        
        # Filter reads
        if (!file.exists(paste0(bam, ".bai"))) {
            indexBam(bam)
        }
        
        if (!dir.exists(bam_output_dir)) dir.create(bam_output_dir, recursive = TRUE)
        
        filtered_bam <- file.path(bam_output_dir, paste0(tools::file_path_sans_ext(basename(bam)), ".exonic.bam"))
        
        sorted_bam <- paste0(tools::file_path_sans_ext(filtered_bam), ".sorted")
        
        filterBam(
            file = bam,
            destination = filtered_bam,
            indexDestination = FALSE,
            param = ScanBamParam(which = exons_for_bam)
        )
        
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


#' Pad a sparse count matrix with missing cell barcode columns
#'
#' Align a sparse matrix by a desired set/order of column names (cell barcodes),
#' adding zero-filled columns for any barcodes not currently present. This is
#' primarily used when summing sparse matrices generated from different chunks
#' that may contain different subsets of cells.
#'
#' @param M A sparse matrix (typically a \code{dgCMatrix}) with non-NULL column
#'   names representing cell barcodes. Row names are preserved.
#' @param cols Character vector of target column names (barcodes). The returned
#'   matrix will have columns exactly in this order. Missing barcodes are added
#'   as all-zero columns.
#'
#' @return A sparse matrix with the same number of rows as \code{M} and columns
#'   exactly equal to \code{cols} (same order). If \code{cols} contains names not
#'   present in \code{M}, they are appended as zero columns before reordering.
#'
#' @details
#' This function assumes that \code{colnames(M)} uniquely identify columns and
#' can be used as a key for alignment across matrices. It will error if column
#' names are \code{NULL}.
#'
#' When repeatedly called in a loop, this can incur overhead due to repeated
#' sparse matrix reallocation and subsetting; consider accumulating triplets
#' (\code{i}, \code{j}, \code{x}) and building a single sparse matrix when
#' performance is critical.
#' 
#' @importFrom Matrix Matrix cbind2
#' @importFrom IRanges IRanges
#' 
#' @keywords internal
#' 
pad_cols <- function(M, cols) {
    # Pad sparse matrix M to include 'cols' (barcodes), adding zero columns where missing.
    if (is.null(colnames(M))) stop("Matrix has NULL column names; cannot align by barcodes.")
    missing <- setdiff(cols, colnames(M))
    if (length(missing) > 0L) {
        zero <- Matrix(0, nrow = nrow(M), ncol = length(missing),
                       dimnames = list(rownames(M), missing), sparse = TRUE)
        M <- cbind2(M, zero)
    }
    m <- M[, cols, drop = FALSE]
    return(m)
}


#' Convert paired-end alignments to fragment ranges
#'
#' Given a \code{GAlignmentPairs} object, construct a \code{GRanges} representing
#' the full inferred fragment span for each pair by taking the minimum start and
#' maximum end across the first and last mate. This is useful for counting
#' overlaps at the fragment level rather than at the read level.
#'
#' @param gap A \code{GAlignmentPairs} object (typically produced by
#'   \code{GenomicAlignments::readGAlignmentPairs()}).
#' @param ignore.strand Logical. If \code{TRUE}, returned ranges have strand
#'   set to \code{"*"} regardless of input. If \code{FALSE}, the strand from the
#'   first mate is used.
#'
#' @return A \code{GRanges} object with one range per input pair, with:
#' \itemize{
#'   \item \code{seqnames} from the first mate,
#'   \item \code{start = min(start(first), start(last))},
#'   \item \code{end   = max(end(first),   end(last))},
#'   \item \code{strand} depending on \code{ignore.strand}.
#' }
#'
#' @details
#' This function assumes mates are on the same reference sequence. For unusual
#' cases where mates map to different chromosomes, the resulting behavior may be
#' undefined and should be handled upstream (e.g., by filtering such pairs).
#'
#' @import GenomicRanges
#' @importFrom GenomicAlignments first last
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors Rle
#' @keywords internal
fragmentize_pairs <- function(gap, ignore.strand = TRUE) {
    fst <- first(gap)
    lst <- last(gap)
    GRanges(
        seqnames = seqnames(fst),
        ranges   = IRanges(start = pmin(start(fst), start(lst)),
                           end   = pmax(end(fst),   end(lst))),
        strand   = if (ignore.strand) Rle("*") else strand(fst)
    )
}


#' Tally per-chunk single-cell feature counts using per-seqname NCList indexes
#'
#' Vectorized overlap counting for a chunk of alignments with cell barcode (CB)
#' and optional UMI (UB) annotations. Reads are split by \code{seqnames} and
#' overlapped against pre-built \code{GNCList} indexes for the matching seqname.
#' Counts are accumulated into a sparse matrix of dimension
#' \code{n_features x n_cells}.
#'
#' @param read_gr \code{GRanges} of read (or fragment) intervals for a single
#'   chunk. Must be aligned to the same reference naming as the feature indexes.
#' @param cb Character vector of cell barcodes (one per read in \code{read_gr}).
#'   \code{NA} values are dropped.
#' @param ub Optional character vector of UMIs (one per read). If \code{NULL},
#'   UMI-based deduplication is disabled regardless of \code{dedup}.
#' @param mapq_vec Optional integer/numeric vector of MAPQ values (one per read).
#'   If provided, reads with \code{mapq < mapq_min} are dropped.
#'
#' @param gr_common \code{GRanges} of all features used to define the global row
#'   space. (Not directly used for overlap here, but defines global indexing.)
#' @param feature_ids Character vector of length \code{n_features} giving global
#'   feature identifiers used as row names of the returned matrix.
#' @param n_features Integer number of global features (rows).
#'
#' @param nclist_by_seq Named list of \code{GNCList} objects, one per seqname,
#'   built from the per-seqname split of \code{gr_common}. Names must correspond
#'   to seqnames present in \code{read_gr}.
#' @param subj_global_idx_by_seq Named list mapping each seqname to an integer
#'   vector of global feature indices (rows of \code{gr_common}) corresponding
#'   to the features used to build \code{nclist_by_seq[[seq]]}. This is used to
#'   translate \code{subjectHits()} (local per-seqname feature indices) into
#'   global row indices in \code{1:n_features}.
#'
#' @param ignore.strand Logical; passed to \code{findOverlaps()}.
#' @param mapq_min Minimum MAPQ threshold. Only used if \code{mapq_vec} is not
#'   \code{NULL}.
#' @param cells Optional character vector of allowed cell barcodes. If provided,
#'   only these cells are retained and output columns are returned in this order.
#' @param dedup Logical. If \code{TRUE} and \code{ub} is provided, deduplicate
#'   UMIs per \code{(cell, feature)} within the chunk (i.e., each unique UMI
#'   contributes at most 1 count to a \code{(cell, feature)} pair).
#' @param mode Character overlap mode placeholder. Currently unused in the
#'   implementation shown; reserved for future extensions (e.g., intersection
#'   logic).
#'
#' @return A sparse matrix (\code{dgCMatrix}) with \code{n_features} rows and one
#'   column per cell observed/retained in the chunk. Row names are \code{feature_ids}.
#'   Column names are either:
#' \itemize{
#'   \item \code{cells} (if provided; fixed order), or
#'   \item sorted unique barcodes observed in the chunk (if \code{cells} is \code{NULL}).
#' }
#' If no reads pass filtering or no overlaps are found, returns \code{NULL}.
#'
#' @details
#' \strong{Deduplication semantics.} When \code{dedup=TRUE} and \code{ub} is not
#' \code{NULL}, counts are deduplicated by unique triplets \code{(cell, feature, UMI)}
#' within the current chunk and seqname. If the same \code{(cell, feature, UMI)}
#' appears multiple times in the chunk, it is counted once.
#'
#' \strong{Column ordering.} If \code{cells} is supplied, columns are aligned to
#' this whitelist/order (via \code{match()}). Otherwise, columns are based on
#' the set of barcodes present within the chunk (sorted).
#'
#' \strong{Performance notes.} This function currently constructs and merges
#' per-seqname sparse matrices using repeated column union and padding. For very
#' large numbers of cells, it can be faster to accumulate triplets and build a
#' single sparse matrix per chunk.
#'
#' @seealso \code{\link[GenomicRanges]{findOverlaps}}, \code{\link[GenomicRanges]{GNCList}},
#'   \code{\link[Matrix]{sparseMatrix}}
#'
#' @import GenomicRanges
#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors queryHits subjectHits
#' 
#' @keywords internal
#' 
tally_chunk_vec <- function(read_gr, cb, ub, mapq_vec,
                            gr_common, feature_ids, n_features,
                            nclist_by_seq, subj_global_idx_by_seq,
                            ignore.strand, mapq_min, cells, dedup, mode) {
    
    # Filter reads by CB, MAPQ
    keep_read <- !is.na(cb)
    if (!is.null(mapq_vec)) keep_read <- keep_read & (mapq_vec >= mapq_min)
    if (!any(keep_read)) return(NULL)
    read_gr <- read_gr[keep_read]
    cb <- cb[keep_read]
    ub <- if (!is.null(ub)) ub[keep_read] else NULL
    mapq_vec <- if (!is.null(mapq_vec)) mapq_vec[keep_read] else NULL
    
    # Split reads by seqname to use the matching NCList
    reads_by_seq <- split(read_gr, seqnames(read_gr))
    seq_levels <- names(reads_by_seq)
    if (length(seq_levels) == 0L) return(NULL)
    
    # Accumulate i (feature idx), j (cell idx) vectors across seqnames
    i_all <- integer(0); j_all <- integer(0)
    
    for (seqn in seq_levels) {
        # Skip if features have no index for this seqname
        if (!(seqn %in% names(nclist_by_seq))) next
        
        reads_seq <- reads_by_seq[[seqn]]
        if (length(reads_seq) == 0L) next
        
        # Overlaps against NCList
        hits <- findOverlaps(reads_seq, nclist_by_seq[[seqn]], ignore.strand = ignore.strand)
        if (length(hits) == 0L) next
        q <- queryHits(hits)    # indices into reads_seq
        s <- subjectHits(hits)  # indices into features subset on this seqname
        
        # Map to global read indices: indices of reads_seq inside full read_gr after filtering
        reads_seq_global_idx <- which(as.character(seqnames(read_gr)) == seqn)
        q_global <- reads_seq_global_idx[q]
        # Map to global feature row indices (rows of gr_common)
        s_global <- subj_global_idx_by_seq[[seqn]][s]
        
        # Extract per-hit cell barcodes and UMIs
        cell_hits <- cb[q_global]
        umi_hits  <- if (!is.null(ub)) ub[q_global] else NULL
        
        # Optional restrict to whitelist of cells
        if (!is.null(cells)) {
            keep_c <- cell_hits %in% cells
            if (!any(keep_c)) next
            cell_hits <- cell_hits[keep_c]
            s_global  <- s_global[keep_c]
            umi_hits  <- if (!is.null(umi_hits)) umi_hits[keep_c] else NULL
        }
        
        if (length(cell_hits) == 0L) next
        
        # Vectorizsd UMI dedup per (cell, feature)
        if (dedup && !is.null(umi_hits)) {
            # Drop NA UMIs
            ok <- !is.na(umi_hits)
            if (!any(ok)) next
            cell_hits <- cell_hits[ok]; s_global <- s_global[ok]; umi_hits <- umi_hits[ok]
            keys <- paste(cell_hits, s_global, umi_hits, sep = "\t")
            dedup_idx <- !duplicated(keys)
            cell_hits <- cell_hits[dedup_idx]
            s_global  <- s_global[dedup_idx]
        }
        # One record per hit to count as 1
        # Factorize cells for columns
        cell_levels <- if (!is.null(cells)) cells else sort(unique(cell_hits))
        j <- match(cell_hits, cell_levels)
        i <- s_global
        
        # Append to accumulator vectors
        i_all <- c(i_all, i)
        # Build a sparse matrix once, using full union of cell_levels across seqnames.
        # Capture j values relative to 'cell_levels', use a local sparse then pad+sum.
        # Return a sparse matrix per chunk with column names = unique cells in chunk.
        
        # Build chunk-level sparse immediately for this seqname
        M_seq <- sparseMatrix(i = i, j = j, x = 1L,
                              dims = c(n_features, length(cell_levels)),
                              dimnames = list(feature_ids, cell_levels))
        # Combine per-seqname matrices into one chunk matrix
        if (!exists("M_chunk", inherits = FALSE)) {
            M_chunk <- M_seq
        } else {
            cols_union <- union(colnames(M_chunk), colnames(M_seq))
            M_chunk <- pad_cols(M_chunk, cols_union) + pad_cols(M_seq, cols_union)
        }
    }
    
    if (!exists("M_chunk", inherits = FALSE)) return(NULL)
    return(M_chunk)
}


#' Count single-cell reads/UMIs from BAM over genomic features (optimized)
#'
#' Vectorized, memory-efficient counting from BAMs with CB/UB tags.
#' Overlaps are computed against per-seqname NCList indexes; UMIs are deduplicated
#' per (cell, feature) if provided. Results are sparse matrices (features × cells).
#'
#' @param gr GRanges of genomic features (e.g., ORFs). Must be on the same reference
#'   naming style as BAMs or convertible via seqlevelsStyle.
#' @param bam Character vector of BAM paths (STARsolo/CellRanger or similar).
#' @param seqlevels_style Character; the naming style for chromosome 
#' identifiers (e.g., \code{"UCSC"}, \code{"NCBI"}). This is applied to both 
#' the GRanges object and the TxDb annotation. Default is \code{"UCSC"}.
#' @param tags Named list with `cell` and `umi` tag names (default CB/UB).
#' @param ignore.strand Logical: whether to ignore strand in overlaps.
#' @param yieldSize Integer chunk size for streaming from BAM.
#' @param mapq Minimum MAPQ to keep (default 10).
#' @param dedup Logical: if TRUE and UMI tag present, deduplicate UMIs per (cell, feature).
#' @param mode Overlap mode; currently supports "Union" (others placeholder for extension).
#' @param cells Optional character vector of barcodes to keep; if provided, columns are
#'   restricted to these and ordered accordingly.
#' @param verbose Logical for progress messages.
#' @param BPPARAM A BiocParallelParam. Default:
#'   \code{BiocParallel::MulticoreParam(workers = 12, stop.on.error = FALSE, progressbar = TRUE)}.
#'   
#' @return A sparse dgCMatrix (features × cells), with rownames taken from `names(gr)` if present.
#' 
#' @import GenomicRanges 
#' @importFrom GenomicAlignments first last readGAlignmentPairs readGAlignments
#' @importFrom Rsamtools ScanBamParam BamFile scanBamFlag
#' @importFrom S4Vectors queryHits subjectHits Rle mcols mcols<-
#' @importFrom IRanges IRanges
#' @importFrom Matrix Matrix sparseMatrix
#' @importFrom GenomeInfoDb seqinfo seqlevels seqlevelsStyle
#' @importFrom BiocParallel MulticoreParam bpnworkers bplapply
#' 
#' @export
#' 
countReadsSingleCell <- function(
        gr,
        bam,
        seqlevels_style = "UCSC",
        tags = list(cell = "CB", umi = "UB"),
        ignore.strand = TRUE,
        yieldSize = 1e6,
        mapq = 10,
        dedup = TRUE,
        mode = c("Union", "IntersectionStrict", "IntersectionNotEmpty"),
        cells = NULL,
        verbose = TRUE,
        BPPARAM = MulticoreParam(workers = 12, stop.on.error = FALSE, progressbar=TRUE)
) {
    mode <- match.arg(mode)
    stopifnot(is(gr, "GRanges"))
    stopifnot(length(bam) == 1L) # Enforce single BAM file
    stopifnot(file.exists(bam))
    paired <- length(group_bam_files(bam)$paired) > 0
    start_counts <- Sys.time()
    
    # Global Setup
    # Setup BAM/GR info
    bf <- BamFile(bam, yieldSize = yieldSize)
    # Harmonization and filtering of seqlevels
    try({
        seqlevelsStyle(gr) <- if (!is.null(seqlevelsStyle(bf))) {
            seqlevelsStyle(bf)
        } else {
            seqlevels_style
        }
    }, silent = TRUE)
    bam_si   <- seqinfo(bf)
    common_seqs <- intersect(seqlevels(gr), seqlevels(bam_si))
    if (length(common_seqs) == 0L) {
        if (verbose) message("No overlapping seqlevels between features and BAM; returning empty matrix.")
        feature_ids <- if (!is.null(names(gr))) names(gr) else paste0("feat_", seq_along(gr))
        return(Matrix(0, nrow = length(feature_ids), ncol = 0,
                      dimnames = list(feature_ids, character()), sparse = TRUE))
    }
    gr_common <- keepSeqlevels(gr, common_seqs, pruning.mode = "coarse")
    feature_ids <- if (!is.null(names(gr_common))) names(gr_common) else paste0("feat_", seq_along(gr_common))
    n_features  <- length(gr_common)
    # Decompose features into regions for parallelization
    gr_regions <- split(gr_common, seqnames(gr_common))
    region_names <- names(gr_regions)
    if (verbose) message(sprintf("Processing BAM %s across %d genomic regions in parallel.", basename(bam), length(region_names)))
    # Prebuild NCList indexes and global index map (Required for tally_chunk_vec)
    gr_by_seq <- split(gr_common, seqnames(gr_common))
    nclist_by_seq <- lapply(gr_by_seq, GNCList)
    subj_global_idx_by_seq <- lapply(names(gr_by_seq), function(seqn) {
        which(as.character(seqnames(gr_common)) == seqn)
    })
    names(subj_global_idx_by_seq) <- names(gr_by_seq)
    # Core Parallel Function
    process_single_region <- function(region_name) {
        # Get subset of features for this region
        gr_subset <- gr_regions[[region_name]]
        # Local setup of ScanBamParam for this specific region (TARGETED I/O)
        which_ranges <- gr_subset # Restricts I/O to only this region's coordinates
        param <- ScanBamParam(
            tag  = c(tags$cell, tags$umi),
            what = c("mapq"),
            which = which_ranges, # Critical for efficiency
            flag = scanBamFlag(isSecondaryAlignment = FALSE,
                               isSupplementaryAlignment = FALSE)
        )
        # Local BamFile object for the worker to stream the BAM
        bf_local <- BamFile(bam, yieldSize = yieldSize)
        open(bf_local)
        on.exit(close(bf_local), add = TRUE)
        # Initialize accumulator for this region's counts
        final_mat_region <- Matrix(0, nrow = n_features, ncol = 0,
                                   dimnames = list(feature_ids, character()), sparse = TRUE)
        k <- 0L
        repeat {
            # Read Chunk (Targeted I/O)
            if (paired) {
                ga <- readGAlignmentPairs(bf_local, use.names = FALSE, param = param)
                if (length(ga) == 0L) break
                rd_gr <- fragmentize_pairs(ga)
                tags_list <- mcols(first(ga))
            } else {
                ga <- readGAlignments(bf_local, use.names = FALSE, param = param)
                if (length(ga) == 0L) break
                rd_gr <- granges(ga)
                tags_list <- mcols(ga)
            }
            # Read seqlevel harmonization
            rd_gr <- keepSeqlevels(rd_gr, common_seqs, pruning.mode = "coarse")
            seqinfo(rd_gr) <- seqinfo(gr_common)[common_seqs]
            cb <- tags_list[[tags$cell]]
            ub <- if (!is.null(tags$umi) && (tags$umi %in% names(tags_list))) tags_list[[tags$umi]] else NULL
            mapq_vec <- mcols(ga)[["mapq"]]
            # Vectorized tally for the chunk (using global NCList/index info)
            chunk_mat <- tally_chunk_vec(
                read_gr = rd_gr, cb = cb, ub = ub, mapq_vec = mapq_vec,
                gr_common = gr_common, feature_ids = feature_ids, n_features = n_features,
                nclist_by_seq = nclist_by_seq, subj_global_idx_by_seq = subj_global_idx_by_seq,
                ignore.strand = ignore.strand, mapq_min = mapq,
                cells = cells, dedup = dedup, mode = mode
            )
            # Iterative Accumulation
            if (!is.null(chunk_mat) && ncol(chunk_mat) > 0L) {
                k <- k + 1L
                # The accumulation step is only merging small chunks into the region's final_mat_region
                if (ncol(final_mat_region) == 0L) {
                    final_mat_region <- chunk_mat
                } else {
                    cols_union <- union(colnames(final_mat_region), colnames(chunk_mat))
                    final_mat_region <- pad_cols(final_mat_region, cols_union) + pad_cols(chunk_mat, cols_union)
                }
            }
        }
        # Worker returns the final accumulated sub-matrix for this region
        return(final_mat_region)
    }
    # Parallel Execution and Final Assembly
    if (verbose) message(sprintf("Starting BiocParallel with %d workers...", bpnworkers(BPPARAM)))
    sub_matrices <- bplapply(
        region_names,
        FUN = process_single_region,
        BPPARAM = BPPARAM,
        BPOPTIONS = list( # BPOPTIONS is a list of options
            export = c( # 'export' is a named field within BPOPTIONS
                "gr_regions", "gr_common", "n_features", "feature_ids",
                "common_seqs", "tags", "yieldSize", "mapq", "dedup",
                "mode", "cells", "paired", "bam",
                "nclist_by_seq", "subj_global_idx_by_seq",
                "tally_chunk_vec", "pad_cols", "fragmentize_pairs"
            )
        )
    )
    # Filter out empty results
    sub_matrices <- Filter(function(X) ncol(X) > 0L, sub_matrices)
    # 3. Final Assembly (Serial, Fast)
    if (length(sub_matrices) == 0L) {
        mat <- Matrix(0, nrow = n_features, ncol = 0,
                      dimnames = list(feature_ids, character()), sparse = TRUE)
    } else {
        # Rbind the resulting matrices to reduce the total number of cell barcodes to align
        mat <- Reduce(function(A, B) {
            # Only need to align the cell barcodes (columns) before summing
            cols <- union(colnames(A), colnames(B))
            pad_cols(A, cols) + pad_cols(B, cols)
        }, sub_matrices)
    }
    if (verbose) {
        end_counts <- Sys.time()
        elapsed <- runtime(end_counts, start_counts)
        if (!is.null(elapsed$mins)) {
            message(sprintf("constructed matrix: %d features x %d cells", nrow(mat), ncol(mat)))
            message(sprintf("read counting runtime: %d mins %.3f secs", elapsed$mins, elapsed$secs))
        } else {
            message(sprintf("constructed matrix: %d features x %d cells", nrow(mat), ncol(mat)))
            message(sprintf("read counting runtime: %.3f secs", elapsed$secs))
        }
    }
    return(mat)
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
#' @import GenomicRanges
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
#' required columns (\code{run}, \code{strategy}, \code{condition}, and
#' \code{replicate}) are present (case-insensitive match), and normalizes
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
            "Please use DOTSeqDataSetsFromSummarizeOverlaps or DOTSeqDataSetsFromFeaturCounts."
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
#' @import SummarizedExperiment
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
#' @seealso \code{\link{DOTSeqDataSetsFromFeatureCounts}}
#' 
#' @rdname DOTSeqDataSets
#' 
#' @importFrom methods new
#' @import SummarizedExperiment
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


#' Construct DOTSeqDatasets from summarizeOverlaps read counts 
#' for Differential ORF Translation Analysis
#' 
#' @description
#' This function initialize and construct the 
#' \code{\link{DOTSeqDataSets-class}} object. This includes loading count 
#' and metadata tables, read a \code{GRanges} object with ORF-level 
#' annotation, filtering ORFs based on count thresholds, and preparing 
#' objects for downstream differential translation analysis using 
#' beta-binomial and negative binomial GLM.
#'
#' @param count_table A dataframe of read counts with features as rows and 
#' samples as columns. Generated from \code{\link{countReads}} based on 
#' \code{\link[GenomicAlignments]{summarizeOverlaps}}.
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
#' @seealso \code{\link{DOTSeqDataSetsFromFeatureCounts}}
#' 
#' @rdname DOTSeqDataSets
#' 
#' @importFrom methods new
#' 
#' @export
#' 
DOTSeqDataSetsFromSummarizeOverlaps <- function(
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
