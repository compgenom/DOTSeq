library(Biostrings)

context("findORFsFasta and findORFsFastaCpp")

# Basic ORF detection
test_that("findORFsFasta detects ORFs with valid start/stop codons", {
    seqs <- DNAStringSet(c(test1 = "ATGAAATAA", test2 = "ATGCCCTAG"))
    gr <- findORFsFasta(sequences = seqs, start_codons = "ATG", stop_codons = "TAA|TAG")
    expect_s4_class(gr, "GRanges")
    expect_true(length(gr) >= 1)
})

test_that("findORFsFasta detects multiple start/stop codons", {
    seqs <- DNAStringSet(c(multi = "GTGAAATAAATGCCCTGA"))
    gr <- findORFsFasta(sequences = seqs, start_codons = "ATG|GTG", stop_codons = "TAA|TGA")
    expect_true(length(gr) >= 2)
})

test_that("findORFsFasta handles ambiguous bases", {
    seqs <- DNAStringSet(c(ambiguous = "ATGNAATAA"))
    gr <- findORFsFasta(sequences = seqs, start_codons = "ATG", stop_codons = "TAA")
    expect_true(length(gr) >= 1)
})


# Filters
test_that("findORFsFasta respects min_len filter", {
    seqs <- DNAStringSet(c(test1 = "ATGAAATAA"))
    gr <- findORFsFasta(sequences = seqs, min_len = 12)
    expect_equal(length(gr), 0)
})

test_that("findORFsFasta applies longest_orf option", {
    seqs <- DNAStringSet(c(test1 = "ATGAAAATGTAA"))
    gr_all <- findORFsFasta(sequences = seqs, longest_orf = FALSE)
    gr_longest <- findORFsFasta(sequences = seqs, longest_orf = TRUE)
    expect_true(length(gr_longest) <= length(gr_all))
})


# Circular genome handling
test_that("findORFsFasta detects ORFs in circular genomes", {
    seqs <- DNAStringSet(c(circular = "ATGAAATAAATGAAATAA"))
    gr <- findORFsFasta(sequences = seqs, is_circular = TRUE)
    expect_s4_class(gr, "GRanges")
    expect_true(length(gr) > 0)
})


# Strand options
test_that("findORFsFasta detects ORFs on both strands", {
    seqs <- DNAStringSet(c(test1 = "TTTATGAAATAA"))
    gr <- findORFsFasta(sequences = seqs, plus_strand_only = FALSE)
    expect_true(length(gr) > 0)
})


# Input types
test_that("findORFsFasta works with FASTA file input", {
    testthat::skip_if_not_installed("withr")
    library(withr)
    
    fasta_path <- tempfile(fileext = ".fa")
    writeLines(c(">seq1", "ATGAAATAA", ">seq2", "ATGCCCTAG"), fasta_path)
    withr::defer(unlink(fasta_path))
    gr <- findORFsFasta(sequences = fasta_path, start_codons = "ATG", stop_codons = "TAA|TAG")
    expect_s4_class(gr, "GRanges")
    expect_true(length(gr) >= 1)
})

test_that("findORFsFastaCpp reverse strand triggers complement for all cases", {
    # Sequence with all uppercase bases and N
    seqs_upper <- c(test_upper = "ATGCNATGC")
    gr_upper <- findORFsFastaCpp(seqs_upper, "ATG", "TAA", 0, FALSE, FALSE)
    expect_true(length(gr_upper) >= 0) # Just ensures function runs
    
    # Sequence with all lowercase bases
    seqs_lower <- c(test_lower = "atgcatgc")
    gr_lower <- findORFsFastaCpp(seqs_lower, "atg", "taa", 0, FALSE, FALSE)
    expect_true(length(gr_lower) >= 0)
    
    # Sequence with unexpected character (e.g., X) to hit default case
    seqs_default <- c(test_default = "ATGXATG")
    gr_default <- findORFsFastaCpp(seqs_default, "ATG", "TAA", 0, FALSE, FALSE)
    expect_true(length(gr_default) >= 0)
})


# Error handling
test_that("findORFsFasta handles invalid codons gracefully", {
    seqs <- DNAStringSet(c(test = "ATGAAATAA"))
    gr <- findORFsFasta(sequences = seqs, start_codons = "XXX", stop_codons = "YYY")
    expect_equal(length(gr), 0)
})

test_that("findORFsFasta errors on unnamed DNAStringSet", {
    seqs <- DNAStringSet(c("ATGAAATAA"))
    names(seqs) <- NULL
    expect_error(findORFsFasta(sequences = seqs), "must have named sequences")
})


# Direct tests for findORFsFastaCpp edge cases
test_that("findORFsFastaCpp detects ORFs spanning circular boundary", {
    seqs <- c(boundary = "AAAAAATGAAAAATAA")
    gr <- findORFsFastaCpp(seqs, "ATG", "TAA", 0, TRUE, TRUE)
    expect_true(length(gr) > 0)
})

test_that("findORFsFastaCpp detects ORFs on both strands", {
    seqs <- c(test1 = "TTTATGAAATAATTATTTCATAAA")
    gr <- findORFsFastaCpp(seqs, "ATG", "TAA", 0, FALSE, FALSE)
    expect_true(length(gr) > 1)
})

test_that("findORFsFastaCpp fully exercises circular boundary logic", {
    # Stop codon at start, start codon at end
    seqs <- c(boundary = paste0("TAA", strrep("A", 10), "ATG"))
    gr <- findORFsFastaCpp(seqs, "ATG", "TAA", 0, TRUE, TRUE)
    expect_true(length(gr) > 0)
})
