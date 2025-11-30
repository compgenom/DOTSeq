library(IRanges)

# Create mock GRanges and GRangesList for testing
gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
                   ranges = IRanges(c(7, 14), width = 3),
                   strand = c("+", "+"))
gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
                    ranges = IRanges(c(4, 1), c(9, 3)),
                    strand = c("-", "-"))
gr_mixed <- GRanges(seqnames = c("chr3", "chr3"),
                    ranges = IRanges(c(1, 10), width = 5),
                    strand = c("+", "-"))

grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)

test_that("strandBool returns correct logical vector", {
    expect_equal(strandBool(gr_plus), c(TRUE, TRUE))
    expect_equal(strandBool(gr_minus), c(FALSE, FALSE))
    expect_equal(strandBool(grl), c(TRUE, FALSE))
})

test_that("strandPerGroup returns strand per group", {
    strands <- strandPerGroup(grl)
    expect_equal(as.character(strands), c("+", "-"))
    expect_s4_class(strands, "Rle")
})

test_that("stopSites returns correct positions", {
    stops <- stopSites(grl, is.sorted = FALSE)
    expect_length(stops, 2)
    expect_true(all(stops > 0))
    
    # Test as GRanges output
    stops_gr <- stopSites(grl, asGR = TRUE)
    expect_s4_class(stops_gr, "GRanges")
    expect_equal(length(stops_gr), 2)
})

test_that("sortPerGroup sorts correctly by strand", {
    sorted_grl <- sortPerGroup(grl)
    expect_s4_class(sorted_grl, "GRangesList")
    expect_equal(names(sorted_grl), names(grl))
    
    # Check that plus strand is sorted ascending
    expect_true(all(start(sorted_grl[[1]]) <= end(sorted_grl[[1]])))
    
    # For minus strand, ends should be in decreasing order
    expect_true(all(diff(end(sorted_grl[[2]])) <= 0))
})

# Mock GRanges and GRangesList for testing
gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
                   ranges = IRanges(c(7, 14), width = 3),
                   strand = c("+", "+"))
gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
                    ranges = IRanges(c(1, 4), width = 3),
                    strand = c("-", "-"))
grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)

test_that("numExonsPerGroup returns correct counts", {
    counts <- numExonsPerGroup(grl)
    expect_equal(counts, c(tx1 = 2L, tx2 = 2L))
})

test_that("lastExonPerGroup returns last exon per group", {
    last_exons <- lastExonPerGroup(grl)
    expect_s4_class(last_exons, "GRangesList")
    expect_equal(length(last_exons), 2)
})

test_that("lastExonEndPerGroup returns correct ends", {
    ends <- lastExonEndPerGroup(grl, keep.names = FALSE)
    expect_type(ends, "integer")
    expect_true(all(ends > 0))
})

test_that("lastExonStartPerGroup returns correct starts", {
    starts <- lastExonStartPerGroup(grl, keep.names = FALSE)
    expect_type(starts, "integer")
    expect_true(all(starts > 0))
})

test_that("widthPerGroup returns total widths", {
    widths <- widthPerGroup(grl, keep.names = FALSE)
    expect_type(widths, "integer")
    expect_true(all(widths > 0))
})

test_that("seqnamesPerGroup returns correct seqnames", {
    seqnames <- seqnamesPerGroup(grl, keep.names = FALSE)
    expect_type(seqnames, "character")
    expect_equal(seqnames, c("chr1", "chr2"))
})

test_that("reverseMinusStrandPerGroup reverses minus strand groups", {
    # Add a minus strand group with increasing order
    gr_minus_inc <- GRanges(seqnames = c("chr2", "chr2"),
                            ranges = IRanges(c(1, 4), width = 3),
                            strand = c("-", "-"))
    grl2 <- GRangesList(tx1 = gr_plus, tx2 = gr_minus_inc)
    
    reversed <- reverseMinusStrandPerGroup(grl2)
    expect_s4_class(reversed, "GRangesList")
    expect_equal(names(reversed), names(grl2))
})
