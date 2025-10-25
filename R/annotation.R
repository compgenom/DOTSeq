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
    results$gene_id <- sub("\\..*", "", results$orf_id)
    
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
#' \doi{10.1038/nprot.2009.97}
#'
#' Durinck S, Moreau Y, Kasprzyk A, Davis S, De Moor B, Brazma A, 
#' Huber W (2005). BioMart and Bioconductor: a powerful link between 
#' biological databases and microarray data analysis. Bioinformatics, 
#' 21, 3439–3440. \doi{10.1093/bioinformatics/bti525}
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

