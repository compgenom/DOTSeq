#' Parse BED File and Match to ORFs
#'
#' This function reads a BED file and converts it into a `GRanges` object.
#' It then finds overlaps with ORFs in the DOTSeq object and returns
#' a data frame of exact matches, including BED score annotations.
#'
#' @param bed Path to a BED file (tab-delimited, no header required) with columns:
#'   chr, start, end, name, score, strand.
#' @param dxd A `DEXSeqDataSet` object containing ORF ranges.
#'
#' @return A data frame of ORFs with exact BED matches. Includes BED score in
#'   the `labels` column.
#'
#' @examples
#' # Parse a BED file and match to ORFs
#' orfDf <- parseBed("example.bed", dxd)
#' head(orfDf)
#'
#' @export
parseBed <- function(bed, dxd) {
  bedDf <- read.table(bed, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(bedDf) <- c("chr", "start", "end", "name", "score", "strand")
  
  # Convert to GRanges
  bedGR <- GRanges(
    seqnames = bedDf$chr,
    ranges = IRanges(start = bedDf$start + 1, end = bedDf$end),
    strand = bedDf$strand,
    name = bedDf$name,
    label = bedDf$score
  )
  
  # Find overlaps
  orfRanges <- rowRanges(dxd)
  hits <- findOverlaps(orfRanges, bedGR)
  
  # Extract matched pairs
  orfs <- orfRanges[queryHits(hits)]
  beds  <- bedGR[subjectHits(hits)]
  
  # Filter for exact matches
  # exactHits <- which(
  #   seqnames(orfs) == seqnames(beds) &
  #     start(orfs) == start(beds) &
  #     end(orfs) == end(beds) &
  #     strand(orfs) == strand(beds)
  # )
  # Coerce to character/integer vectors
  orfs_chr   <- as.character(seqnames(orfs))
  beds_chr   <- as.character(seqnames(beds))
  orfs_start <- start(orfs)
  beds_start <- start(beds)
  orfs_end   <- end(orfs)
  beds_end   <- end(beds)
  orfs_str   <- as.character(strand(orfs))
  beds_str   <- as.character(strand(beds))
  
  # Compare element-wise
  exactHits <- which(
    orfs_chr == beds_chr &
      orfs_start == beds_start &
      orfs_end == beds_end &
      orfs_str == beds_str
  )
  
  
  # Put ORF label back
  exactBeds <- beds[exactHits]
  exactORFs <- orfs[exactHits]
  
  # Get corresponding BED scores and add BED score to overlapping exons
  mcols(exactORFs)$labels <- mcols(exactBeds)$label
  
  orfDf <- as.data.frame(exactORFs)
  rownames(orfDf) <- paste0(orfDf$groupID, ":", orfDf$featureID)
  
  return(orfDf = orfDf)
}
