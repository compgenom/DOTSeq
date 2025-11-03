#' Filter GTF file for transcript features with required metadata
#'
#' This function imports a GTF file and filters transcript features 
#' based on the presence of required metadata columns and an optional 
#' source filter.
#'
#' @param gtf_file Character. Path to the GTF file.
#' @param require_ids Character vector. Metadata column names that must 
#' be present and non-NA. Default: \code{c("hgnc_id", "protein_id", 
#' "ccdsid")}.
#' @param source_filter Character or NULL. If provided, filters 
#' transcripts by the 'source' column. Default: \code{NULL}.
#' @param verbose Logical. Whether to print summary messages. 
#' Default: \code{TRUE}.
#'
#' @return A filtered \code{GRanges} object containing transcript 
#' features.
#'
#' @importFrom rtracklayer import
#' @importFrom S4Vectors mcols mcols<-
#' @keywords internal
#' 
filter_gtf <- function(
        gtf_file,
        require_ids = c("hgnc_id", "protein_id", "ccdsid"),
        source_filter = NULL,
        verbose = TRUE
) {
    
    # Import only transcript features
    gtf <- import(gtf_file, feature.type = "transcript")
    input_length <- length(gtf)
    
    for (col in require_ids) {
        if (col %in% names(mcols(gtf))) {
            gtf <- gtf[!is.na(mcols(gtf)[[col]]), ]
        } else {
            warning("Column", col, "not found in gtf")
        }
    }
    
    if (!is.null(source_filter)) {
        if ("source" %in% names(gtf)) {
            gtf <- gtf[gtf$source == source_filter, ]
        } else {
            warning("Column 'source' not found in gtf")
        }
    }
    
    if (verbose) {
        message("selected ", length(gtf), " from ", input_length, " transcript IDs in the GTF file")
    }
    
    return(gtf)
}


#' Extract Genomic ORFs from Transcript Sequences
#'
#' Identifies open reading frames (ORFs) from transcript sequences and maps 
#' them to genomic coordinates using a GTF/GFF annotation file or a TxDb 
#' object. Supports input sequences as a FASTA file, a \code{DNAStringSet}, 
#' or a \code{BSgenome} object. Classifies small ORFs (sORFs) as upstream 
#' (uORF), downstream (dORF), or overlapping (oORF) relative to the main 
#' ORFs (mORFs).
#'
#' @param sequences Transcript sequences. Can be a character string (path 
#' to a FASTA file),
#' a \code{DNAStringSet} object, or a \code{BSgenome} object.
#' @param annotation Transcript annotation. Can be a character string (path 
#' to a GTF or GFF file)
#' or a \code{TxDb} object.
#' @param organism Character string specifying the organism name (used  
#' only when building a TxDb from a GTF/GFF file). Default is 
#' \code{"Homo sapiens"}.
#' @param require_ids Character vector. Metadata column names that must 
#' be present and non-NA. Default: \code{c("hgnc_id", "protein_id", 
#' "ccdsid")}.
#' @param source_filter Character or NULL. If provided, filters 
#' transcripts by the 'source' column. Default: \code{NULL}.
#' @param circ_seqs Character vector of circular sequences to exclude 
#' (e.g., \code{"chrM"}).
#' Default is \code{"chrM"}.
#' @param start_codons Character vector of start codons to search for 
#' (e.g., \code{"ATG"}).
#' Default is \code{"ATG"}.
#' @param stop_codons Character string of stop codons separated by 
#' \code{"|"} (e.g., \code{"TAA|TAG|TGA"}).
#' Default is \code{"TAA|TAG|TGA"}.
#' @param min_len Integer specifying the minimum ORF length in bases. 
#' Default is \code{0}.
#' @param longest_orf Logical. If \code{TRUE}, only the longest ORF per 
#' transcript is returned. Default is \code{TRUE}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages and 
#' timing information. Default is \code{TRUE}.
#'
#' @return A \code{GRanges} object containing genomic coordinates of 
#' ORFs, with metadata columns \code{gene_id} and \code{orf_type}. 
#' Main ORFs are labeled as \code{"mORF"}, and small ORFs are 
#' classified as \code{"uORF"}, \code{"dORF"}, or \code{"oORF"}.
#'
#' @details
#' - ORFs are identified in transcript space using \code{findORFsFasta()}.
#' - Coordinates are mapped to the genome using \code{mapFromTranscripts()} 
#' and exon annotations.
#' - Main ORFs are defined by overlap with annotated CDS regions.
#' - Small ORFs are classified relative to mORFs based on strand-aware 
#' genomic position.
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom BiocGenerics unlist
#' @importFrom BSgenome getSeq
#' @importFrom GenomicRanges reduce strand GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomicFeatures cdsBy exonsBy mapFromTranscripts transcripts genes
#' @importFrom txdbmaker makeTxDbFromGFF
#' @importFrom AnnotationDbi select keys
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom stats ave
#' @export
#'
#' @examplesIf requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) && requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE) && requireNamespace("GenomicFeatures", quietly = TRUE)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(GenomicFeatures)
#'
#' # Load genome and TxDb
#' genome <- BSgenome.Hsapiens.UCSC.hg38
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#'
#' # Get exons grouped by transcript
#' exons_by_tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
#'
#' # Select a single transcript for demonstration
#' tx1 <- head(exons_by_tx, 100)
#'
#' # Extract transcript sequence
#' tx_seqs <- extractTranscriptSeqs(genome, tx1)
#'
#' # Run getORFs on the transcript sequence
#' gr <- getORFs(
#'     sequences = tx_seqs,
#'     annotation = txdb
#' )
#' print(gr)
#' 
#' @references
#' Lawrence, M., Huber, W., Pagès, H., Aboyoun, P., Carlson, M., 
#' Gentleman, R., Morgan, M., Carey, V. (2013). Software for Computing 
#' and Annotating Genomic Ranges. PLoS Computational Biology, 9. 
#' DOI: 10.1371/journal.pcbi.1003118
#' 
#' Tjeldnes, H., Labun, K., Torres Cleuren, Y. et al. 
#' ORFik: a comprehensive R toolkit for the analysis of translation. 
#' BMC Bioinformatics 22, 336 (2021). DOI: 10.1186/s12859-021-04254-w
#' 
getORFs <- function(
        sequences,
        annotation,
        organism = "Homo sapiens", 
        require_ids = c("hgnc_id", "protein_id", "ccdsid"),
        source_filter = NULL,
        circ_seqs = NULL,
        start_codons = "ATG",
        stop_codons = "TAA|TAG|TGA",
        min_len = 0,
        longest_orf = TRUE,
        verbose = TRUE
) {
    
    gtf <- NULL
    
    # Load annotation
    if (is.character(annotation)) {
        annotation <- path.expand(annotation)
        if (!file.exists(annotation)) stop("File does not exist: ", annotation)
        
        # Check if it's a single string and a valid file path
        if (length(annotation) != 1 || !nzchar(annotation) || !file.exists(annotation)) {
            stop("Invalid file path: ", annotation)
        }
        
        if (verbose) {
            start_annotate <- Sys.time()
            message("invoking makeTxDbFromGFF:")
        }
        
        # Apply filter to GTF
        gtf <- filter_gtf(gtf_file = annotation, require_ids = require_ids, source_filter = source_filter, verbose = verbose)
        
        format <- if (grepl("\\.gtf", annotation, ignore.case = TRUE)) "gtf" else "gff"
        
        annotation <- makeTxDbFromGFF(
            file = annotation, 
            format = format, 
            organism = organism, 
            circ_seqs = circ_seqs
        )
        
        if (verbose) {
            end_annotate <- Sys.time()
            elapsed_annotate <- runtime(end_annotate, start_annotate)
            
            format_upper <- toupper(format)
            
            if (!is.null(elapsed_annotate$mins)) {
                message(sprintf("%s parsing runtime: %d mins %.3f secs", format_upper, elapsed_annotate$mins, elapsed_annotate$secs))
            } else {
                message(sprintf("%s parsing runtime: %.3f secs", format_upper, elapsed_annotate$secs))
            }
            
            start_seq <- Sys.time()
        }
        
        # Extract exons grouped by transcript
        exons_by_tx <- exonsBy(annotation, by = "tx", use.names = TRUE)
        
    } else if (inherits(annotation, "TxDb")) {
        start_seq <- Sys.time()
        exons_by_tx <- exonsBy(annotation, by = "tx", use.names = TRUE)
    } else if (inherits(annotation, "GRangesList")) {
        exons_by_tx <- annotation
    } else {
        stop("Unsupported 'annotation' input type: must be file path or TxDb object.")
    }
    
    
    # Load sequences
    if (is.character(sequences)) {
        sequences <- path.expand(sequences)
        if (!file.exists(sequences)) stop("File does not exist: ", sequences)
        
        seqs <- readDNAStringSet(sequences)
        names(seqs) <- vapply(strsplit(names(seqs), "[\\|| ]"), function(x) x[1], character(1))
        
    } else if (inherits(sequences, "DNAStringSet")) {
        if (is.null(names(sequences))) stop("DNAStringSet must have named sequences.")
        seqs <- sequences
        
    } else if (inherits(sequences, "BSgenome")) {
        seqs <- getSeq(sequences)
        
    } else {
        stop("Unsupported 'sequences' input type: must be file path, DNAStringSet, or BSgenome.")
    }
    
    if (length(seqs) == 0) {
        stop("No sequences found. Please check your input sequences.")
    }
    
    if (verbose) {
        end_seq <- Sys.time()
        elapsed_seq <- runtime(end_seq, start_seq)
        message("read ", length(seqs), " fasta sequences successfully")
        
        if (!is.null(elapsed_seq$mins)) {
            message(sprintf("sequence reading runtime: %d mins %.3f secs", elapsed_seq$mins, elapsed_seq$secs))
        } else {
            message(sprintf("sequence reading runtime: %.3f secs", elapsed_seq$secs))
        }
        
        start_orfs <- Sys.time()
        
        message(
            "finding ORFs with start codons: ", 
            paste(start_codons, collapse = ", "), 
            " on transcript sequences..."
        )
    }
    
    if (!is.null(gtf)) {
        seqs <- seqs[names(seqs) %in% gtf$transcript_id]
    }
    
    # Find ORFs
    gr <- findORFsFasta(
        sequences = seqs, 
        start_codons = start_codons, 
        stop_codons = stop_codons, 
        min_len = min_len, 
        longest_orf = longest_orf,
        plus_strand_only = TRUE
    )
    
    if (verbose) {
        end_orfs <- Sys.time()
        elapsed_orfs <- runtime(end_orfs, start_orfs)
        
        if (!is.null(elapsed_orfs$mins)) {
            message(sprintf("ORF finding runtime: %d mins %.3f secs", elapsed_orfs$mins, elapsed_orfs$secs))
        } else {
            message(sprintf("ORF finding runtime: %.3f secs", elapsed_orfs$secs))
        }
        
        start_map <- Sys.time()
        
        message("starting transcript-to-genomic coordinate conversion")
    }
    
    # Flatten ORFs
    gr_orfs <- reduce(gr, with.revmap = TRUE)
    
    gr_tx <- transcripts(annotation, use.names = TRUE)
    # Split transcripts by strand
    tx_plus  <- gr_tx[strand(gr_tx) == "+"]
    tx_minus <- gr_tx[strand(gr_tx) == "-"]
    # Split ORFs by transcript strand (based on seqnames)
    orfs_plus  <- gr_orfs[as.character(seqnames(gr_orfs)) %in% names(tx_plus)]
    orfs_minus <- gr_orfs[as.character(seqnames(gr_orfs)) %in% names(tx_minus)]
    strand(orfs_minus) <- "-"# A workaround to do minus strand coordinate mapping
    
    exons_by_tx_plus <- exons_by_tx[as.character(names(exons_by_tx)) %in% names(tx_plus)]
    exons_by_tx_minus <- exons_by_tx[as.character(names(exons_by_tx)) %in% names(tx_minus)]
    
    # Map separately
    orfs_plus  <- mapFromTranscripts(orfs_plus, exons_by_tx_plus)
    orfs_minus <- mapFromTranscripts(orfs_minus, exons_by_tx_minus)
    flattened_orfs <- reduce(c(orfs_plus, orfs_minus)) # All flattened ORFs
    
    # Flatten mORFs
    # Extract cds grouped by transcript
    cds_by_tx <- cdsBy(annotation, by = "tx", use.names = TRUE)
    gr_cds <- unlist(cds_by_tx, use.names = TRUE)
    
    # Find overlaps
    hits <- findOverlaps(flattened_orfs, gr_cds, ignore.strand = FALSE)
    
    # overlapping_orfs <- flattened_orfs[unique(queryHits(hits))]
    overlapping_cds  <- gr_cds[unique(subjectHits(hits))]
    
    # ID mapping
    tx2gene <- select(
        annotation, 
        keys = keys(annotation, "TXID"), 
        columns = c("TXNAME", "GENEID"), 
        keytype = "TXID"
    )
    
    # Create a named vector for fast lookup
    txname_to_gene <- setNames(tx2gene$GENEID, tx2gene$TXNAME)
    
    # Match transcript names in overlapping_cds to gene IDs
    cds_txnames <- names(overlapping_cds)
    
    # Attach gene IDs to overlapping_cds
    mcols(overlapping_cds)$gene_id <- txname_to_gene[cds_txnames]
    
    # Split morfs by gene_id
    overlapping_cds_by_gene <- split(overlapping_cds, mcols(overlapping_cds)$gene_id)
    
    # Reduce each group to a single range
    flattened_overlapping_cds <- range(overlapping_cds_by_gene)
    
    # Add gene_id back as metadata
    mcols(flattened_overlapping_cds)$gene_id <- names(flattened_overlapping_cds)
    flattened_overlapping_cds <- unlist(flattened_overlapping_cds)
    gr_morfs <- unique(flattened_overlapping_cds)
    
    if (length(gr_morfs) > 0) {
        mcols(gr_morfs)$gene_id <- names(gr_morfs)
    } else {
        stop("No mORFs found after filtering and overlap. Check input data or filtering criteria.")
    }
    
    mcols(gr_morfs)$gene_id <- names(gr_morfs)
    names(gr_morfs) <- NULL
    
    # Classify sORFs
    gr_genes <- genes(annotation)
    
    hits <- findOverlaps(flattened_orfs, gr_genes, ignore.strand = FALSE)
    
    # Extract gene IDs from gr_genes
    gene_ids <- mcols(gr_genes)$gene_id
    
    # Attach gene IDs to sORFs
    flattened_orfs <- flattened_orfs[queryHits(hits)]
    if (length(flattened_orfs) > 0) {
        mcols(flattened_orfs)$gene_id <- gene_ids[subjectHits(hits)]
    } else {
        stop("No sORFs found overlapping with genes.")
    }
    
    # Drop duplicates
    flattened_orfs <- unique(flattened_orfs)
    
    # Convert to data.frame
    flattened_orfs <- as.data.frame(flattened_orfs)
    morf_df <- as.data.frame(gr_morfs)
    
    orf_df <- merge(
        flattened_orfs,
        morf_df,
        by = c("gene_id", "seqnames", "strand"),
        suffixes = c("_flattened_orf", "_morf")
    )
    
    orf_df$orf_type <- with(orf_df, ifelse(
        strand == "+",
        ifelse(
            start_flattened_orf == start_morf, "mORF",
            ifelse(
                end_flattened_orf < start_morf, "uORF",
                ifelse(
                    start_flattened_orf > end_morf, "dORF",
                    "oORF"
                )
            )
        ),
        ifelse(
            end_flattened_orf == end_morf, "mORF",
            ifelse(
                start_flattened_orf > end_morf, "uORF",
                ifelse(
                    end_flattened_orf < start_morf, "dORF",
                    "oORF"
                )
            )
        )
    ))
    
    # TODO: Omit overlapping ORFs (oORFs) for now 
    # - It is a challenge to assign reads to the correct ORFs
    # - These oORFs are as large as mORFs
    # # Overwrite orf_type for oORFs only
    # orf_df$orf_type <- with(orf_df, ifelse(
    #     orf_type == "oORF" & strand == "+" & start_flattened_orf < start_morf, "uoORF",
    #     ifelse(
    #         orf_type == "oORF" & strand == "+" & start_flattened_orf >= start_morf, "doORF",
    #         ifelse(
    #             orf_type == "oORF" & strand == "-" & end_flattened_orf > end_morf, "uoORF",
    #             ifelse(
    #                 orf_type == "oORF" & strand == "-" & end_flattened_orf <= end_morf, "doORF",
    #                 orf_type  # keep original classification for non-oORFs
    #             )
    #         )
    #     )
    # ))
    
    # Create GRanges
    flattened_gr <- GRanges(
        seqnames = orf_df$seqnames,
        ranges = IRanges(start = orf_df$start_flattened_orf, end = orf_df$end_flattened_orf),
        strand = orf_df$strand
    )
    
    # Add gene_id and orf_type to metadata
    mcols(flattened_gr)$gene_id <- orf_df$gene_id
    mcols(flattened_gr)$orf_type <- orf_df$orf_type
    
    flattened_gr <- flattened_gr[flattened_gr$orf_type != "oORF", ]
    flattened_gr <- flattened_gr[flattened_gr$orf_type != "mORF", ]
    
    # Merge with mORF GRanges
    mcols(gr_morfs)$orf_type <- "mORF"
    flattened_gr <- c(gr_morfs, flattened_gr)
    flattened_gr <- sort(flattened_gr)
    
    # Add ORF IDs
    # Get the gene_id vector
    gene_ids <- mcols(flattened_gr)$gene_id
    
    # Create a factor for grouping (more efficient for ave)
    gene_id_factor <- as.factor(gene_ids)
    
    # Calculate the sequential number for each row within its gene_id group
    orf_seq_num <- as.integer(ave(seq_along(flattened_gr), gene_id_factor, FUN = seq_along))
    
    # Format the number string and create the full ORF ID
    orf_num_str <- sprintf("%03d", orf_seq_num)
    full_orf_id <- paste0(gene_ids, ":O", orf_num_str)
    
    # Assign the full ID to the names
    names(flattened_gr) <- full_orf_id
    
    # Assign just the number part to the metadata column
    mcols(flattened_gr)$orf_number <- orf_num_str
    
    if (verbose) {
        end_map <- Sys.time()
        elapsed_map <- runtime(end_map, start_map)
        
        if (!is.null(elapsed_map$mins)) {
            message(sprintf("coordinate mapping runtime: %d mins %.3f secs", elapsed_map$mins, elapsed_map$secs))
        } else {
            message(sprintf("coordinate mapping runtime: %.3f secs", elapsed_map$secs))
        }
        
        message("flattened ", length(gr), " ORFs to ", length(flattened_gr))
    }
    return(flattened_gr)
}


utils::globalVariables(c(".", ".I", ".N", ".SD", ":="))

#' Get longest ORF per stop site
#'
#' Rule: if seqname, strand and stop site is equal, take longest one.
#' Else keep.
#' If IRangesList or IRanges, seqnames are groups, if GRanges or GRangesList
#' seqnames are the seqlevels (e.g. chromosomes/transcripts)
#'
#' @author Haakon Tjeldnes et al.
#' 
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}/IRangesList, GRanges/IRanges of ORFs
#' @return a \code{\link[GenomicRanges]{GRangesList}}/IRangesList, GRanges/IRanges
#' (same as input)
#' 
#' @importFrom data.table data.table
#' @importFrom GenomicRanges seqnames width width<-
#' 
#' @family ORFHelpers
#' 
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' ORF1 = GRanges("1", IRanges(10,21), "+")
#' ORF2 = GRanges("1", IRanges(1,21), "+") # <- longest
#' grl <- GRangesList(ORF1 = ORF1, ORF2 = ORF2)
#' longestORFs(grl) # get only longest
#' }
#' 
#' @references
#' Tjeldnes, H., Labun, K., Torres Cleuren, Y. et al. 
#' ORFik: a comprehensive R toolkit for the analysis of translation. 
#' BMC Bioinformatics 22, 336 (2021). DOI: 10.1186/s12859-021-04254-w
#' 
longestORFs <- function(grl) {
    if(length(grl) == 0) return(grl) # if empty
    
    if (is(grl, "GRangesList")) { # only for GRangesList
        stops <- stopSites(grl, is.sorted = TRUE)
        widths <- widthPerGroup(grl, FALSE)
        seqnames <- seqnamesPerGroup(grl, FALSE)
        strands <- strandPerGroup(grl, FALSE)
    } else { # GRanges, IRanges or IRangesList
        stops <- unlist(end(grl), use.names = FALSE)
        widths <- unlist(width(grl), use.names = FALSE)
        
        if (is(grl, "GRanges")) { # GRanges
            seqnames <- as.character(seqnames(grl))
            strands <- as.character(strand(grl))
            stops[strands == "-"] <- start(grl)[strands == "-"]
        } else { # IRanges or IRangesList
            strands <- rep("+", length(widths))
            if (is(grl, "IRanges")) {
                seqnames <- rep.int(1, length(widths))
            } else if (is(grl, "IRangesList")) {
                seqnames <- rep.int(seq.int(length(grl)), lengths(grl))
            }
        }
    }
    dt <- data.table(seqnames, strands, stops, widths)
    longestORFs <- dt[, .I[which.max(widths)], by = .(seqnames, strands, stops)]$V1
    if (is(grl, "IRangesList")) {
        ir <- unlist(grl, use.names = FALSE)
        ir <- ir[longestORFs]
        irl <- split(ir, seqnames[longestORFs])
        names(irl) <- names(grl)
        return(irl)
    }
    return(grl[longestORFs])
}


#' Get logical list of strands
#'
#' Helper function to get a logical list of True/False,
#'  if GRangesList group have + strand = T, if - strand = F
#' Also checks for * strands, so a good check for bugs
#' 
#' @author Haakon Tjeldnes et al.
#' 
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} or GRanges object
#' @return a logical vector
#' 
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' gr <- GRanges(Rle(c("chr2", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#'               IRanges(1:10, width = 10:1),
#'               Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)))
#' strandBool(gr)
#' }
#' 
#' @references
#' Tjeldnes, H., Labun, K., Torres Cleuren, Y. et al. 
#' ORFik: a comprehensive R toolkit for the analysis of translation. 
#' BMC Bioinformatics 22, 336 (2021). DOI: 10.1186/s12859-021-04254-w
#' 
strandBool <- function(grl) {
    if (is(grl, "GRanges") | is(grl, "GAlignments") | is(grl, "GAlignmentPairs")) {
        posIndices <- as.character(strand(grl)) == "+"
    } else {
        posIndices <- strandPerGroup(grl, FALSE) == "+"
    }
    
    sums <- sum(posIndices) + sum(!posIndices)
    if (is.na(sums)) {
        stop(
            "could not get strands from grl object",
            " most likely NULL object was passed."
        )
    }
    if (sums != length(grl)) {
        stop("grl contains * strands, set them to either + or -")
    }
    return(posIndices)
}


#' Get list of strands per granges group
#' 
#' @author Haakon Tjeldnes et al.
#' 
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
#' @param keep.names a boolean, keep names or not, default: (TRUE)
#' @return a vector named/unnamed of characters
#' @importFrom IRanges heads PartitioningByEnd
#' @importFrom GenomicRanges start strand
#' 
#' @keywords internal
#' 
#' @references
#' Tjeldnes, H., Labun, K., Torres Cleuren, Y. et al. 
#' ORFik: a comprehensive R toolkit for the analysis of translation. 
#' BMC Bioinformatics 22, 336 (2021). DOI: 10.1186/s12859-021-04254-w
#' 
strandPerGroup <- function(grl, keep.names = TRUE) {
    # validGRL(class(grl))
    starts <- start(PartitioningByEnd(grl))
    res <- strand(grl@unlistData)[starts]
    if (!keep.names) {
        return(as.character(res))
    }
    return(res)
}


#' Get the stop sites from a GRangesList of orfs grouped by orfs
#'
#' In ATGTTTTGC, get the position of the C.
#' 
#' @author Haakon Tjeldnes et al.
#' 
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} object
#' @param asGR a boolean, return as GRanges object
#' @param keep.names a logical (FALSE), keep names of input.
#' @param is.sorted a speedup, if you know the ranges are sorted
#' @return if asGR is False, a vector, if True a GRanges object
#' @family ORFHelpers
#' 
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' stopSites(grl, is.sorted = FALSE)
#' }
#' 
#' @references
#' Tjeldnes, H., Labun, K., Torres Cleuren, Y. et al. 
#' ORFik: a comprehensive R toolkit for the analysis of translation. 
#' BMC Bioinformatics 22, 336 (2021). DOI: 10.1186/s12859-021-04254-w
#' 
stopSites <- function(
        grl, 
        asGR = FALSE, 
        keep.names = FALSE,
        is.sorted = FALSE
) {
    
    if (!is.sorted) {
        grl <- sortPerGroup(grl)
    }
    posIds <- strandBool(grl)
    
    stopSites <- rep(NA, length(grl))
    stopSites[posIds] <- lastExonEndPerGroup(grl[posIds], FALSE)
    stopSites[!posIds] <- lastExonStartPerGroup(grl[!posIds], FALSE)
    
    if (asGR) {
        stopSites <- GRanges(
            seqnames = seqnamesPerGroup(grl, FALSE),
            ranges = IRanges(stopSites, width = 1),
            strand = strandPerGroup(grl, FALSE),
            seqinfo = seqinfo(grl)
        )
    }
    if (keep.names) {
        names(stopSites) <- names(grl)
    }
    return(stopSites)
}


#' Sort a GRangesList
#'
#' A faster, more versatile reimplementation of
#' \code{\link[GenomicRanges]{sort.GenomicRanges}} for GRangesList,
#' needed since the original works poorly for more than 10k groups.
#' This function sorts each group, where "+" strands are
#' increasing by starts and "-" strands are decreasing by ends.
#'
#' Note: will not work if groups have equal names.
#' 
#' @author Haakon Tjeldnes et al.
#' 
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
#' @param ignore.strand a boolean, (default FALSE): should minus strands be
#' sorted from highest to lowest ends. If TRUE: from lowest to highest ends.
#' @param quick.rev default: FALSE, if TRUE, given that you know all ranges are
#' sorted from min to max for both strands, it will only reverse coordinates for
#' minus strand groups, and only if they are in increasing order. Much quicker
#' @return an equally named GRangesList, where each group is
#'  sorted within group.
#'  
#' @keywords internal
#' 
#' @references
#' Tjeldnes, H., Labun, K., Torres Cleuren, Y. et al. 
#' ORFik: a comprehensive R toolkit for the analysis of translation. 
#' BMC Bioinformatics 22, 336 (2021). DOI: 10.1186/s12859-021-04254-w
#' 
sortPerGroup <- function(grl, ignore.strand = FALSE, quick.rev = FALSE){
    if (quick.rev) {
        return(reverseMinusStrandPerGroup(grl))
    }
    if (!ignore.strand) {
        indicesPos <- strandBool(grl)
        
        grl[indicesPos] <- gSort(grl[indicesPos])
        grl[!indicesPos] <- gSort(
            grl[!indicesPos], 
            decreasing = TRUE,
            byStarts = FALSE)
        return(grl)
    }
    return(gSort(grl))
}


#' Reverse minus strand
#'
#' Reverse minus strand per group in a GRangesList
#' Only reverse if minus strand is in increasing order
#' 
#' @author Haakon Tjeldnes et al.
#' 
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
#' @param onlyIfIncreasing logical, default (TRUE), only reverse if decreasing
#' @return a \code{\link[GenomicRanges]{GRangesList}}
#' @keywords internal
#' 
#' @references
#' Tjeldnes, H., Labun, K., Torres Cleuren, Y. et al. 
#' ORFik: a comprehensive R toolkit for the analysis of translation. 
#' BMC Bioinformatics 22, 336 (2021). DOI: 10.1186/s12859-021-04254-w
#' 
reverseMinusStrandPerGroup <- function(grl, onlyIfIncreasing = TRUE) {
    minus <- !strandBool(grl)
    if (onlyIfIncreasing) {
        minGrl <- grl[minus & numExonsPerGroup(grl, FALSE) > 1]
        if (length(minGrl) == 0) return(grl)
        decreasing <- start(minGrl[[1]])[1] > start(minGrl[[1]])[2]
        if (decreasing) return(grl)
    }
    oldGrl <- rev(grl)
    oldGrl[rev(minus)]@unlistData@ranges <- rev(grl[minus]@unlistData@ranges)
    return(rev(oldGrl))
}


#' Get list of the number of exons per group
#'
#' Can also be used generaly to get number of GRanges object
#'  per GRangesList group
#'  
#' @author Haakon Tjeldnes et al.
#' 
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
#' @param keep.names a logical, keep names or not, default: (TRUE)
#' @return an integer vector of counts
#' 
#' @references
#' Tjeldnes, H., Labun, K., Torres Cleuren, Y. et al. 
#' ORFik: a comprehensive R toolkit for the analysis of translation. 
#' BMC Bioinformatics 22, 336 (2021). DOI: 10.1186/s12859-021-04254-w
#' 
numExonsPerGroup <- function(grl, keep.names = TRUE) {
    # validGRL(class(grl))
    return(lengths(grl, keep.names))
}


utils::globalVariables(c("group", "grnames"))

#' Sort a GRangesList, helper.
#'
#' A helper for sortPerGroup.
#' A faster, more versatile reimplementation of GenomicRanges::sort()
#' Normally not used directly.
#' Groups first each group, then either decreasing or increasing
#' (on starts if byStarts == T, on ends if byStarts == F)
#' 
#' @author Haakon Tjeldnes et al.
#' 
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
#' @param decreasing should the first in each group have max(start(group))
#'   ->T or min-> default(F) ?
#' @param byStarts a logical T, should it order by starts or ends F.
#' @importFrom data.table as.data.table :=
#' @importFrom BiocGenerics unlist
#' @return an equally named GRangesList, where each group is sorted within
#' group.
#' @keywords internal
#' 
#' @references
#' Tjeldnes, H., Labun, K., Torres Cleuren, Y. et al. 
#' ORFik: a comprehensive R toolkit for the analysis of translation. 
#' BMC Bioinformatics 22, 336 (2021). DOI: 10.1186/s12859-021-04254-w
#' 
gSort <- function(grl, decreasing = FALSE, byStarts = TRUE) {
    if (length(grl) == 0) return(GRangesList())
    
    DT <- as.data.table(grl)
    DT$group_name <- NULL
    group <- NULL # for not getting warning
    if (decreasing) {
        if (byStarts) {
            DT <- DT[order(group, -start)]
        } else {
            DT <- DT[order(group, -end)]
        }
    } else {
        if (byStarts) {
            DT <- DT[order(group, start)]
        } else {
            DT <- DT[order(group, end)]
        }
    }
    # TODO: test naming, this is still not perfect
    testName <- names(unlist(grl[1], use.names = FALSE)[1])
    if (!is.null(testName)) {
        DT[, grnames := names(unlist(grl, use.names = FALSE))]
    }
    
    asgrl <- makeGRangesListFromDataFrame(
        DT, split.field = "group",
        names.field = if(is.null(testName)) NULL else "grnames",
        keep.extra.columns = TRUE)
    
    names(asgrl) <- names(grl)
    
    return(asgrl)
}


#' Get last end per granges group
#' 
#' @author Haakon Tjeldnes et al.
#' 
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
#' @param keep.names a boolean, keep names or not, default: (TRUE)
#' @return a Rle(keep.names = T), or integer vector(F)
#' 
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' lastExonEndPerGroup(grl)
#' }
#' 
#' @references
#' Tjeldnes, H., Labun, K., Torres Cleuren, Y. et al. 
#' ORFik: a comprehensive R toolkit for the analysis of translation. 
#' BMC Bioinformatics 22, 336 (2021). DOI: 10.1186/s12859-021-04254-w
#' 
lastExonEndPerGroup <- function(grl,keep.names = TRUE) {
    # validGRL(class(grl))
    if (keep.names) {
        return(end(lastExonPerGroup(grl)))
    } else {
        return(as.integer(end(lastExonPerGroup(grl))))
    }
}


#' Get last start per granges group
#' 
#' @author Haakon Tjeldnes et al.
#' 
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
#' @param keep.names a boolean, keep names or not, default: (TRUE)
#' @return a Rle(keep.names = T), or integer vector(F)
#' 
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' lastExonStartPerGroup(grl)
#' }
#' 
#' @references
#' Tjeldnes, H., Labun, K., Torres Cleuren, Y. et al. 
#' ORFik: a comprehensive R toolkit for the analysis of translation. 
#' BMC Bioinformatics 22, 336 (2021). DOI: 10.1186/s12859-021-04254-w
#' 
lastExonStartPerGroup <- function(grl, keep.names = TRUE) {
    # validGRL(class(grl))
    if (keep.names) {
        return(start(lastExonPerGroup(grl)))
    } else {
        return(as.integer(start(lastExonPerGroup(grl))))
    }
}


#' Get last exon per GRangesList group
#'
#' grl must be sorted, call ORFik:::sortPerGroup if needed
#' 
#' @author Haakon Tjeldnes et al.
#' 
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
#' @return a GRangesList of the last exon per group
#' @importFrom IRanges tails
#' 
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' lastExonPerGroup(grl)
#' }
#' 
#' @references
#' Tjeldnes, H., Labun, K., Torres Cleuren, Y. et al. 
#' ORFik: a comprehensive R toolkit for the analysis of translation. 
#' BMC Bioinformatics 22, 336 (2021). DOI: 10.1186/s12859-021-04254-w
#' 
lastExonPerGroup <- function(grl) {
    # validGRL(class(grl))
    return(tails(grl, 1L))
}


#' Get list of widths per granges group
#' 
#' @author Haakon Tjeldnes et al.
#' 
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
#' @param keep.names a boolean, keep names or not, default: (TRUE)
#' 
#' @importFrom GenomicRanges width width<-
#' 
#' @return an integer vector (named/unnamed) of widths
#' 
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' widthPerGroup(grl)
#' }
#' 
#' @references
#' Tjeldnes, H., Labun, K., Torres Cleuren, Y. et al. 
#' ORFik: a comprehensive R toolkit for the analysis of translation. 
#' BMC Bioinformatics 22, 336 (2021). DOI: 10.1186/s12859-021-04254-w
#' 
widthPerGroup <- function(grl, keep.names = TRUE) {
    # validGRL(class(grl))
    if (keep.names) {
        return(sum(width(grl)))
    } else {
        return(as.integer(sum(width(grl))))
    }
}


#' Get list of seqnames per granges group
#' 
#' @author Haakon Tjeldnes et al.
#' 
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
#' @param keep.names a boolean, keep names or not, default: (TRUE)
#' @importFrom IRanges heads
#' @importFrom GenomicRanges seqnames
#' 
#' @return a character vector or Rle of seqnames(if seqnames == T)
#' 
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' seqnamesPerGroup(grl)
#' }
#' 
#' @references
#' Tjeldnes, H., Labun, K., Torres Cleuren, Y. et al. 
#' ORFik: a comprehensive R toolkit for the analysis of translation. 
#' BMC Bioinformatics 22, 336 (2021). DOI: 10.1186/s12859-021-04254-w
#' 
seqnamesPerGroup <- function(grl, keep.names = TRUE) {
    # validGRL(class(grl))
    if (keep.names) {
        return(heads(seqnames(grl), 1L))
    } else {
        return(as.character(heads(seqnames(grl), 1L)))
    }
}


#' Find ORFs in FASTA sequences using ORFik's C++ engine
#'
#' This function identifies ORFs in DNA sequences from FASTA files, 
#' `DNAStringSet`, or `BSgenome` objects. It supports both linear and 
#' circular genomes, and can detect ORFs on the sense (+) or both strands.
#'
#' This method is optimized for prokaryotic genomes or transcript 
#' sequences. Direct use for eukaryotic whole genomes is not suitable 
#' due to the presence of splicing. See also \pkg{ORFik}'s \code{findMapORFs}.
#'
#' Each FASTA header is treated independently, and the name (up to the 
#' first space) is used as the `seqnames` in the returned `GRanges` 
#' object. Circular genome support is included, and ORFs that span the 
#' start/end boundary are handled.
#'
#' Note: Ensure your FASTA file is valid and headers are formatted as 
#' `>name info`, with the name first and no hidden characters, as this 
#' affects coordinate parsing.
#'
#' @author Haakon Tjeldnes et al. (original), 
#' Chun Shen Lim (modification).
#' 
#' @param sequences Character path to a FASTA file, or a `DNAStringSet` 
#' or `BSgenome` object.
#' @param start_codons Character string of start codons 
#' (e.g., "ATG|GTG")
#' @param stop_codons Character string of stop codons 
#' (e.g., "TAA|TAG")
#' @param min_len Integer. Minimum ORF length in bases. Default is \code{0}.
#' @param longest_orf Logical. If TRUE, return only the longest ORF per 
#' sequence.
#' @param is_circular Logical. Whether the genome is circular 
#' (e.g., bacterial genomes).
#' @param plus_strand_only Logical. If TRUE, scan only the forward strand; 
#' if FALSE, scan both strands.
#'
#' @return A `GRanges` object containing the ORFs found.
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom GenomeInfoDb isCircular
#' @importFrom BSgenome getSeq
#' @importFrom Rcpp sourceCpp
#' @useDynLib DOTSeq, .registration = TRUE
#' 
#' @references
#' Tjeldnes, H., Labun, K., Torres Cleuren, Y. et al. 
#' ORFik: a comprehensive R toolkit for the analysis of translation. 
#' BMC Bioinformatics 22, 336 (2021). DOI: 10.1186/s12859-021-04254-w
#' 
findORFsFasta <- function(
        sequences,
        start_codons = "ATG",
        stop_codons = "TAA",
        min_len = 0,
        longest_orf = TRUE,
        is_circular = FALSE,
        plus_strand_only = TRUE
) {
    
    # Load sequences
    seqs <- if (is.character(sequences)) {
        sequences <- path.expand(sequences)
        if (!file.exists(sequences)) stop("File does not exist: ", sequences)
        readDNAStringSet(sequences)
    } else if (inherits(sequences, "DNAStringSet")) {
        if (is.null(names(sequences))) stop("DNAStringSet must have named sequences.")
        sequences
    } else if (inherits(sequences, "BSgenome")) {
        BSgenome::getSeq(sequences)
    } else {
        stop("Unsupported input type: must be file path, DNAStringSet, or BSgenome.")
    }
    
    gr <- findORFsFastaCpp(
        as.character(seqs, use.names = TRUE), 
        start_codons, 
        stop_codons,
        min_len, 
        is_circular,
        plus_strand_only
    )
    
    if (longest_orf) {
        gr <- longestORFs(gr)
    }
    if (is_circular) {
        isCircular(gr) <- rep(TRUE, length(seqlevels(gr)))
    }
    return(gr)
}


#' Extract Significant Genes Based on LFSR Threshold
#'
#' @description
#' Identifies genes with significant differential ORF usage (DOU)
#' based on a local false sign rate (LFSR) threshold. Extracts gene
#' IDs from ORF-level results and filters those with LFSR below the
#' specified threshold.
#'
#' @seealso \code{\link{plotDOT}}
#'
#' @param results A data frame containing ORF-level DOU results.
#' Must include columns \code{orf_id} and the specified \code{padj_col}.
#'
#' @param padj_col Character string specifying the column name for
#' LFSR values. Default is \code{"lfsr"}.
#'
#' @param padj_threshold Numeric threshold for filtering significant
#' ORFs. Default is \code{0.05}.
#'
#' @return A character vector of Ensembl gene IDs corresponding to
#' significant ORFs.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' sig_genes <- get_significant_genes(
#'   results_df,
#'   padj_col = "lfsr",
#'   padj_threshold = 0.05
#' )
#' }
#' 
get_significant_genes <- function(
    results,
    padj_col = "lfsr",
    padj_threshold = 0.05
) {
    results$gene_id <- sub("[.:].*", "", results$orf_id)
    
    gene_ids <- results[
        results[[padj_col]] < padj_threshold,
    ]$gene_id
    
    return(na.omit(gene_ids))
}



#' Internal wrapper for BioMart functions
#' @keywords internal
#' 
bm_use_ensembl <- function(biomart, dataset, host = NULL) {
    
    if (!requireNamespace("biomaRt", quietly=TRUE)) {
        stop(
            "ID mapping require the 'biomaRt' package. ", 
            "Please install it by running: BiocManager::install('biomaRt')"
        )
    }
    
    if (is.null(host)) {
        biomaRt::useEnsembl(
            biomart = biomart, 
            dataset = dataset
        )
    } else {
        biomaRt::useEnsembl(
            biomart = biomart, 
            dataset = dataset, 
            host = host
        )
    }
}

bm_use_ensembl_genomes <- function(biomart, dataset, host = NULL) {
    if (is.null(host)) {
        biomaRt::useEnsemblGenomes(
            biomart = biomart, 
            dataset = dataset
        )
    } else {
        biomaRt::useEnsemblGenomes(
            biomart = biomart, 
            dataset = dataset, 
            host = host
        )
    }
}

bm_get <- function(attributes, filters, values, mart) {
    biomaRt::getBM(
        attributes = attributes,
        filters = filters,
        values = values,
        mart = mart
    )
}


#' Retrieve Gene Symbols from Ensembl or Ensembl Genomes
#'
#' @description
#' Queries Ensembl or Ensembl Genomes BioMart databases to retrieve
#' gene symbols, descriptions, and optionally Gene Ontology (GO)
#' terms. Supports multiple organism groups including vertebrates,
#' plants, fungi, protists, metazoa, and bacteria.
#'
#' If the specified \code{symbol_col} returns only \code{NA} values,
#' the function automatically falls back to using the
#' \code{description} field instead.
#'
#' @param ensembl_ids A character vector of Ensembl gene IDs to query.
#'
#' @param dataset A string specifying the dataset name (e.g.,
#' \code{"hsapiens_gene_ensembl"}, \code{"athaliana_eg_gene"}).
#'
#' @param symbol_col A string specifying the attribute to use as the
#' gene symbol. Common options include \code{hgnc_symbol}, 
#' \code{"external_gene_name"}, \code{description}. The default is 
#' \code{"external_gene_name"}, which is widely used in vertebrate 
#' datasets such as human and mouse.
#'
#' @param include_go Logical; if \code{TRUE}, includes GO annotations
#' (\code{go_id}, \code{name_1006}, \code{namespace_1003}) in the output.
#'
#' @param mart_source A string indicating the BioMart source. One of
#' \code{"ensembl"}, \code{"plants"}, \code{"fungi"}, \code{"protists"},
#' \code{"metazoa"}, or \code{"bacteria"}.
#'
#' @param host Optional. A custom host URL (e.g., for archived Ensembl
#' versions).
#'
#' @return A data frame containing gene symbols for the input Ensembl
#' IDs. If \code{symbol_col} is unavailable, the \code{description}
#' field is used instead and renamed to match \code{symbol_col}.
#'
#' @export
#' @examples
#' # Human gene example
#' mapIDs(
#'     c("ENSG00000139618"),
#'     dataset = "hsapiens_gene_ensembl",
#'     mart_source = "ensembl"
#' )
#' 
#' # Arabidopsis gene example
#' # mapIDs(
#' #    c("AT1G01010"),
#' #    dataset = "athaliana_eg_gene",
#' #    symbol_col = "tair_symbol",
#' #    mart_source = "plants"
#' # )
#'
#' # Plasmodium falciparum gene example with fallback
#' # mapIDs(
#' #    c("PF3D7_0100100"),
#' #    dataset = "pfalciparum_eg_gene",
#' #    mart_source = "protists"
#' # )
#'
#' @references
#' Durinck S, Spellman P, Birney E, Huber W (2009). Mapping identifiers 
#' for the integration of genomic datasets with the R/Bioconductor 
#' package biomaRt. Nature Protocols, 4, 1184–1191.
#' DOI: 10.1038/nprot.2009.97
#'
#' Durinck S, Moreau Y, Kasprzyk A, Davis S, De Moor B, Brazma A, 
#' Huber W (2005). BioMart and Bioconductor: a powerful link between 
#' biological databases and microarray data analysis. Bioinformatics, 
#' 21, 3439–3440. DOI: 10.1093/bioinformatics/bti525
#'
mapIDs <- function(
        ensembl_ids,
        dataset,
        symbol_col = "external_gene_name",
        include_go = FALSE,
        mart_source = "ensembl",
        host = NULL
) {
    
    if (missing(ensembl_ids) || !is.character(ensembl_ids) || length(ensembl_ids) == 0) {
        stop("'ensembl_ids' must be a non-empty character vector.")
    }
    
    if (missing(dataset) || !is.character(dataset) || length(dataset) != 1) {
        stop("'dataset' must be a string.")
    }
    
    if (!is.logical(include_go) || length(include_go) != 1) {
        stop("'include_go' must be TRUE or FALSE.")
    }
    
    valid_source <- c(
        "ensembl",
        "plants",
        "fungi",
        "protists",
        "metazoa",
        "bacteria"
    )
    
    mart_source <- match.arg(mart_source, choices = valid_source)
    
    biomart_map <- list(
        ensembl = "genes",
        plants = "plants_mart",
        fungi = "fungi_mart",
        protists = "protists_mart",
        metazoa = "metazoa_mart",
        bacteria = "bacteria_mart"
    )
    
    # Connect to BioMart
    mart <- if (mart_source == "ensembl") {
        bm_use_ensembl(
            biomart = biomart_map[[mart_source]], 
            dataset = dataset,
            host = host
        )
    } else {
        bm_use_ensembl_genomes(
            biomart = biomart_map[[mart_source]], 
            dataset = dataset, 
            host = host
        )
    }
    
    # Define attributes
    base_attributes <- c("ensembl_gene_id", symbol_col)
    if (include_go) {
        base_attributes <- c(
            base_attributes, 
            "go_id", 
            "name_1006", 
            "namespace_1003"
        )
    }
    
    # Try primary symbol_col
    genes <- bm_get(
        attributes = base_attributes,
        filters = "ensembl_gene_id",
        values = ensembl_ids,
        mart = mart
    )
    
    # Fallback to description if all symbol_col are NA
    if (all(is.na(genes[[symbol_col]]))) {
        message(
            sprintf(
                "'%s' returned all NA; falling back to 'description'", 
                symbol_col
            )
        )
        fallback_attributes <- gsub(
            symbol_col,
            "description", 
            base_attributes
        )
        genes <- bm_get(
            attributes = fallback_attributes,
            filters = "ensembl_gene_id",
            values = ensembl_ids,
            mart = mart
        )
        colnames(genes)[colnames(genes) == "description"] <- symbol_col
    }
    
    genes
}

