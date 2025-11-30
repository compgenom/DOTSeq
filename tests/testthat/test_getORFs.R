test_that("getORFs works with TxDb and BSgenome integration", {
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
    skip_if_not_installed("GenomicFeatures")
    
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(GenomicFeatures)
    
    genome <- BSgenome.Hsapiens.UCSC.hg38
    
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    
    # Get exons grouped by transcript
    exons_by_tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
    
    # Select a small subset for speed
    tx_subset <- head(exons_by_tx, 50) #[25:40]
    
    # Extract transcript sequences
    tx_seqs <- extractTranscriptSeqs(genome, tx_subset)
    
    txdb_output_dir <- tempdir()
    
    # Run getORFs
    gr <- getORFs(
        sequences = tx_seqs,
        annotation = txdb,
        txdb_output_dir = txdb_output_dir,
        start_codons = "ATG",
        stop_codons = "TGA",
        min_len = 0,
        longest_orf = TRUE
    )
    
    # Expectations
    expect_s4_class(gr, "GRanges")
    expect_true(length(gr) > 0)
    expect_true(all(grepl("^ATG", as.character(Biostrings::getSeq(genome, gr)))))
    
    # Check file exists
    sqlite_files <- list.files(txdb_output_dir, pattern = "\\.sqlite$", full.names = TRUE)
    expect_true(length(sqlite_files) > 0)
    
    # Clean up
    unlink(sqlite_files)
    
    # invalid sequences input
    expect_error(
        getORFs(
            sequences = NULL, 
            annotation = txdb,
            start_codons = "ATG",
            stop_codons = "TGA",
            min_len = 0,
            longest_orf = TRUE
        ),
        "Unsupported 'sequences' input type: must be file path, DNAStringSet, or BSgenome"
    )
    
    expect_error(
        getORFs(
            sequences = "not_a_DNAStringSet", 
            annotation = txdb,
            start_codons = "ATG",
            stop_codons = "TGA",
            min_len = 0,
            longest_orf = TRUE
        ),
        "File does not exist"
    )
    
    # invalid annotation input
    expect_error(
        getORFs(
            sequences = tx_seqs, 
            annotation = NULL,
            start_codons = "ATG",
            stop_codons = "TGA",
            min_len = 0,
            longest_orf = TRUE
        ),
        "Unsupported 'annotation' input type"
    )
    
    expect_error(
        getORFs(
            sequences = tx_seqs, 
            annotation = "not_a_TxDb",
            start_codons = "ATG",
            stop_codons = "TGA",
            min_len = 0,
            longest_orf = TRUE
        ),
        "File does not exist"
    )
    
    # Run getORFs
    dir <- system.file("extdata", package = "DOTSeq")
    gtf <- file.path(dir, "gencode.v47.orf_flattened_subset.gtf.gz")
    expect_error(
        suppressWarnings(
            getORFs(
                sequences = tx_seqs,
                annotation = gtf,
                start_codons = "ATG",
                stop_codons = "TGA",
                min_len = 0,
                longest_orf = TRUE
            )
        ),
        "some genes have no \"ID\" attribute"
    )
    
})
