#' Extract Genomic ORFs from Transcript Sequences
#'
#' This function identifies open reading frames (ORFs) from transcript 
#' sequences and maps them to genomic coordinates using either a GTF/GFF 
#' annotation file or a TxDb object. It supports input sequences as a 
#' FASTA file, a DNAStringSet object, or a BSgenome object.
#'
#' @param sequences Input transcript sequences. Can be a character string 
#' (path to a FASTA file), a \code{DNAStringSet} object, or a 
#' \code{BSgenome} object.
#' @param annotation Transcript annotation. Can be a character string 
#' (path to a GTF or GFF file) or a \code{TxDb} object.
#' @param organism Character string specifying the organism name (used only 
#' when building a TxDb from a GTF/GFF file). Default is 
#' \code{"Homo sapiens"}.
#' @param circ_seqs Character vector of circular sequences to exclude 
#' (e.g., \code{"chrM"}). Default is \code{"chrM"}.
#' @param start_codon Character string specifying the start codon(s) 
#' (e.g., \code{"ATG"}). Default is \code{"ATG"}.
#' @param stop_codon Character string specifying stop codon(s), separated 
#' by \code{"|"} (e.g., \code{"TAA|TAG|TGA"}). Default is 
#' \code{"TAA|TAG|TGA"}.
#' @param min_len Integer specifying the minimum ORF length in bases. 
#' Default is \code{0}.
#' @param longest_orf Logical. If \code{TRUE}, only the longest ORF per 
#' region is returned. Default is \code{TRUE}.
#'
#' @return A \code{GRanges} object containing genomic coordinates of ORFs.
#'
#' @details
#' The function uses \code{findORFsFasta()} to identify ORFs from transcript 
#' sequences, filters for ORFs on the positive strand, and optionally merges 
#' overlapping ORFs using \code{reduce()}. It then maps transcript-relative 
#' coordinates to genomic coordinates using exon annotations extracted from 
#' a TxDb object.
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom BiocGenerics unlist
#' @importFrom BSgenome getSeq
#' @importFrom GenomicRanges reduce strand
#' @importFrom GenomicFeatures cdsBy exonsBy mapToTranscripts
#' @importFrom txdbmaker makeTxDbFromGFF
#'
#' @export
#'
#' @examples
#' # Use a BSgenome and TxDb
#' if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
#'     requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE) &&
#'     requireNamespace("GenomicFeatures", quietly = TRUE)) {
#'
#'     library(BSgenome.Hsapiens.UCSC.hg38)
#'     library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#'     library(GenomicFeatures)
#'
#'     # Load genome and TxDb
#'     genome <- BSgenome.Hsapiens.UCSC.hg38
#'     txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#'
#'     # Get exons grouped by transcript
#'     exons_by_tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
#'
#'     # Select a single transcript for demonstration
#'     tx1 <- head(exons_by_tx, 1)
#'
#'     # Extract transcript sequence
#'     tx_seqs <- extractTranscriptSeqs(genome, tx1)
#'
#'     # Run getORFs on the transcript sequence
#'     orfs <- getORFs(
#'         sequences = tx_seqs,
#'         annotation = txdb
#'     )
#'     print(orfs)
#' }
#' 
getORFs <- function(
        sequences,
        annotation,
        organism = "Homo sapiens", 
        circ_seqs = "chrM",
        start_codon = "ATG",
        stop_codon = "TAA|TAG|TGA",
        min_len = 0,
        longest_orf = TRUE,
        verbose = TRUE
) {
    
    if (verbose) {
        start_seq <- Sys.time()
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
    
    if (verbose) {
        end_seq <- Sys.time()
        
        elapsed_seq <- runtime(end_seq, start_seq)
        
        message("parsed ", length(seqs), " sequences successfully")
        
        if (!is.null(elapsed_seq$mins)) {
            message(sprintf("parsing runtime: %d mins %.3f secs", elapsed_seq$mins, elapsed_seq$secs))
        } else {
            message(sprintf("parsing runtime: %.3f secs", elapsed_seq$secs))
        }
        
        start_orfs <- Sys.time()
        
        message(
            "finding ORFs with start codons: ", 
            paste0(start_codon, collapse = ", "), 
            " on transcript sequences"
        )
    }
    
    # Find ORFs
    gr <- findORFsFasta(
        sequences = seqs, 
        start_codon = start_codon, 
        stop_codon = stop_codon, 
        min_len = min_len, 
        longest_orf = longest_orf,
        plus_strand_only = TRUE
    )
    # gr <- gr[strand(gr) == "+"]
    
    if (verbose) {
        end_orfs <- Sys.time()
        
        elapsed_orfs <- runtime(end_orfs, start_orfs)
        
        message("found ", length(gr), " ORFs")
        
        if (!is.null(elapsed_orfs$mins)) {
            message(sprintf("ORF finding runtime: %d mins %.3f secs", elapsed_orfs$mins, elapsed_orfs$secs))
        } else {
            message(sprintf("ORF finding runtime: %.3f secs", elapsed_orfs$secs))
        }
        
        # start_map <- Sys.time()
        # 
        # message("starting transcript-to-genomic coordinate conversion")
    }
    
    # Flatten ORFs
    gr_orfs <- reduce(gr, with.revmap = TRUE)
    
    # Load annotation
    if (is.character(annotation)) {
        annotation <- path.expand(annotation)
        if (!file.exists(annotation)) stop("File does not exist: ", annotation)
        
        format <- if (grepl("\\.gtf", annotation, ignore.case = TRUE)) "gtf" else "gff"
        
        if (verbose) {
            message("invoking makeTxDbFromGFF:")
        }
        annotation <- makeTxDbFromGFF(
            file = annotation, 
            format = format, 
            organism = organism, 
            circ_seqs = circ_seqs
        )
        # Extract exons grouped by transcript
        exons_by_tx <- exonsBy(annotation, by = "tx", use.names = TRUE)
        
    } else if (inherits(annotation, "TxDb")) {
        exons_by_tx <- exonsBy(annotation, by = "tx", use.names = TRUE)
        
    # } else if (inherits(annotation, "GRangesList")) {
    #     exons_by_tx <- annotation
        
    } else {
        stop("Unsupported 'annotation' input type: must be file path or TxDb object.")
    }
    
    # Map ORFs to genome
    # flattened_orfs <- mapFromTranscripts(gr_orfs, exons_by_tx, ignore.strand = TRUE)
    
    gr_exons <- unlist(exons_by_tx, use.names = TRUE)
    
    # Extract cds grouped by transcript
    cds_by_tx <- cdsBy(annotation, by = "tx", use.names = TRUE)
    gr_cds <- unlist(cds_by_tx, use.names = TRUE)
    
    # if (verbose) {
    #     end_map <- Sys.time()
    #     
    #     elapsed_map <- runtime(end_map, start_map)
    #     
    #     message("flattened ", length(gr), " ORFs to ", length(flattened_orfs))
    #     
    #     if (!is.null(elapsed_map$mins)) {
    #         message(sprintf("coordinate conversion runtime: %d mins %.3f secs", elapsed_map$mins, elapsed_map$secs))
    #     } else {
    #         message(sprintf("coordinate conversion runtime: %.3f secs", elapsed_map$secs))
    #     }
    # }
    
    return(list(gr_cds = gr_cds, gr_exons = gr_exons, gr_orfs = gr_orfs))
}


#' Get longest ORF per stop site
#'
#' Rule: if seqname, strand and stop site is equal, take longest one.
#' Else keep.
#' If IRangesList or IRanges, seqnames are groups, if GRanges or GRangesList
#' seqnames are the seqlevels (e.g. chromosomes/transcripts)
#'
#' @author Haakon Tjeldnes et al.
#' 
#' @param grl a \code{\link{GRangesList}}/IRangesList, GRanges/IRanges of ORFs
#' @return a \code{\link{GRangesList}}/IRangesList, GRanges/IRanges
#' (same as input)
#' 
#' @importFrom data.table data.table .I
#' @importFrom GenomicRanges seqnames width width<-
#' 
#' @family ORFHelpers
#' @examples
#' \dontrun{
#' ORF1 = GRanges("1", IRanges(10,21), "+")
#' ORF2 = GRanges("1", IRanges(1,21), "+") # <- longest
#' grl <- GRangesList(ORF1 = ORF1, ORF2 = ORF2)
#' longestORFs(grl) # get only longest
#' }
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
    longestORFs <- dt[, .I[which.max(widths)],
                      by = .(seqnames, strands, stops)]$V1
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
#' @param grl a \code{\link{GRangesList}} or GRanges object
#' @return a logical vector
#' @examples
#' \dontrun{
#' gr <- GRanges(Rle(c("chr2", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#'               IRanges(1:10, width = 10:1),
#'               Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)))
#' strandBool(gr)
#' }
#'
strandBool <- function(grl) {
    if (is(grl, "GRanges") | is(grl, "GAlignments") | is(grl, "GAlignmentPairs")) {
        posIndices <- as.character(strand(grl)) == "+"
    } else {
        posIndices <- strandPerGroup(grl, FALSE) == "+"
    }
    
    sums <- sum(posIndices) + sum(!posIndices)
    if (is.na(sums)) {
        stop("could not get strands from grl object",
             " most likely NULL object was passed.")
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
#' @param grl a \code{\link{GRangesList}}
#' @param keep.names a boolean, keep names or not, default: (TRUE)
#' @return a vector named/unnamed of characters
#' @importFrom IRanges heads
#' @examples
#' \dontrun{
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' strandPerGroup(grl)
#' }
#' 
strandPerGroup <- function(grl, keep.names = TRUE) {
    # validGRL(class(grl))
    starts <- start(IRanges::PartitioningByEnd(grl))
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
#' @param grl a \code{\link{GRangesList}} object
#' @param asGR a boolean, return as GRanges object
#' @param keep.names a logical (FALSE), keep names of input.
#' @param is.sorted a speedup, if you know the ranges are sorted
#' @return if asGR is False, a vector, if True a GRanges object
#' @family ORFHelpers
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
stopSites <- function(grl, asGR = FALSE, keep.names = FALSE,
                      is.sorted = FALSE) {
    if (!is.sorted) {
        grl <- sortPerGroup(grl)
    }
    posIds <- strandBool(grl)
    
    stopSites <- rep(NA, length(grl))
    stopSites[posIds] <- lastExonEndPerGroup(grl[posIds], FALSE)
    stopSites[!posIds] <- lastExonStartPerGroup(grl[!posIds], FALSE)
    
    if (asGR) {
        stopSites <- GRanges(seqnames = seqnamesPerGroup(grl, FALSE),
                             ranges = IRanges(stopSites, width = 1),
                             strand = strandPerGroup(grl, FALSE),
                             seqinfo = seqinfo(grl))
    }
    if (keep.names) {
        names(stopSites) <- names(grl)
    }
    return(stopSites)
}


#' Get last end per granges group
#' 
#' @author Haakon Tjeldnes et al.
#' 
#' @param grl a \code{\link{GRangesList}}
#' @param keep.names a boolean, keep names or not, default: (TRUE)
#' @return a Rle(keep.names = T), or integer vector(F)
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
#' @param grl a \code{\link{GRangesList}}
#' @param keep.names a boolean, keep names or not, default: (TRUE)
#' @return a Rle(keep.names = T), or integer vector(F)
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
#' @param grl a \code{\link{GRangesList}}
#' @return a GRangesList of the last exon per group
#' @importFrom IRanges tails
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
lastExonPerGroup <- function(grl) {
    # validGRL(class(grl))
    return(tails(grl, 1L))
}


#' Get list of widths per granges group
#' 
#' @author Haakon Tjeldnes et al.
#' 
#' @param grl a \code{\link{GRangesList}}
#' @param keep.names a boolean, keep names or not, default: (TRUE)
#' 
#' @importFrom GenomicRanges width width<-
#' 
#' @return an integer vector (named/unnamed) of widths
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
#' @param grl a \code{\link{GRangesList}}
#' @param keep.names a boolean, keep names or not, default: (TRUE)
#' @importFrom IRanges heads
#' @importFrom GenomicRanges seqnames
#' 
#' @return a character vector or Rle of seqnames(if seqnames == T)
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
#' due to the presence of splicing; for spliced ORFs.
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
#' @seealso \code{link[ORFik]{findMapORFs}}
#' 
#' @param sequences Character path to a FASTA file, or a `DNAStringSet` 
#' or `BSgenome` object.
#' @param start_codon Character string of start codons 
#' (e.g., "ATG|GTG")
#' @param stop_codon Character string of stop codons 
#' (e.g., "TAA|TAG")
#' @param min_len Integer. Minimum ORF length in bases. Default is \code{0}.
#' @param is_circular Logical. Whether the genome is circular 
#' (e.g., bacterial genomes).
#' @param strand Character. One of `"+"` or `"both"`. `"+"` scans 
#' only the forward strand; any other value scans both strands.
#'
#' @return A `GRanges` object containing the ORFs found.
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom GenomeInfoDb isCircular
#' @importFrom BSgenome getSeq
#' @importFrom Rcpp sourceCpp
#' @useDynLib DOTSeq, .registration = TRUE
#' 
findORFsFasta <- function(
        sequences,
        start_codon = "ATG",
        stop_codon = "TAA",
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
        start_codon, 
        stop_codon,
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
#' @import biomaRt
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

