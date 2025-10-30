# Declare global variables to use NSE in ggplot2 functions
utils::globalVariables(c("condition", "usage", "orf_id"))

#' Retrieve and format adjusted p-value for a specific ORF and contrast
#'
#' This helper function extracts the adjusted p-value (e.g., LFSR or padj) 
#' for a given ORF and contrast from a results data frame, and formats it 
#' to three decimal places. If the value is missing or NA, it returns "N/A".
#'
#' @param rowdata A data frame containing differential results, including 
#' columns for `orf_id`, `contrast`, and the specified p-value column.
#' 
#' @param orf A character string specifying the ORF ID to look up.
#' 
#' @param contrast_name A character string specifying the contrast name 
#' (e.g., "Treatment - Control").
#' 
#' @param dou_signif_col A character string specifying the column name to 
#' extract the adjusted p-value from. Defaults to `"lfsr"`.
#'
#' @return A character string representing the formatted p-value 
#' (e.g., `"0.023"`).
#' 
#' @keywords internal
#' 
get_lfsr_annotation <- function(
        rowdata, 
        orf, 
        contrast_name, 
        dou_signif_col = "lfsr"
) {
    
    if (!dou_signif_col %in% colnames(rowdata)) {
        stop(sprintf("Column '%s' not found in rowdata.", dou_signif_col))
    }
    
    lfsr_value <- rowdata[rowdata$orf_id == orf & rowdata$contrast == contrast_name, dou_signif_col]
    
    return(as.numeric(lfsr_value))
}


#' Plot ORF usage across conditions with significance annotations
#'
#' Generates a faceted bar plot of ORF usage across experimental 
#' conditions, with optional significance annotations based on adjusted 
#' p-values (e.g., LFSR). Requires 
#' \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2} and 
#' \href{https://CRAN.R-project.org/package=ggsignif}{ggsignif}.
#' 
#' @param data A \code{\link{DOUData-class}} object with an 
#' `interaction` slot storing interaction-specific post hoc 
#' contrast results.
#' 
#' @param gene_id A character string specifying the gene ID of interest
#' 
#' @param id_mapping Optional data frame with gene symbols (e.g., 
#' from biomaRt). Used to label heatmap rows with gene symbols.
#' 
#' @param levels Optional character vector specifying the order of 
#' conditions on the x-axis.
#' 
#' @param dou_signif_thresh Numeric threshold for significance 
#' annotation (default is `0.05`).
#'
#' @param colors Color palette for 
#' \code{\link[ggplot2]{scale_fill_brewer}}. Default is 
#' \code{Set2}.
#' 
#' @return A \code{ggplot} object visualizing ORF usage across 
#' conditions, with significance annotations for contrasts passing 
#' `dou_signif_thresh`. ORFs with non-significant contrasts are 
#' not annotated.
#' 
#' @keywords internal
#' 
#' @importFrom SummarizedExperiment rowRanges strand 
#' @importFrom SummarizedExperiment mcols mcols<- 
#' 
plot_orf_usage <- function(
        data, 
        gene_id = NULL,
        id_mapping = NULL,
        levels = NULL, 
        dou_signif_thresh = 0.05,
        colors = "Set2"
) {
    
    if (!requireNamespace("ggplot2", quietly = TRUE) || 
        !requireNamespace("ggsignif", quietly = TRUE)) {
        stop(
            "Model diagnostics require the 'ggplot2' and 'ggsignif' packages. ", 
            "Please install them by running: install.packages(c('ggplot2', 'ggsignif'))"
        )
    }
    
    annot_df <- getContrasts(data, type = "interaction")
    if (is.null(annot_df)) {
        stop("No interaction slot found in ", data)
    }
    
    rowranges <- rowRanges(data)
    if (is.null(rowranges)) {
        stop("No rowRanges slot found in ", data)
    }
    
    orfs <- data.frame(
        orf_id = rownames(mcols(rowranges)),
        orf_number = mcols(rowranges)$orf_number,
        orf_type = mcols(rowranges)$orf_type,
        strand = as.character(strand(rowranges))
    )
    
    # Remove input gene_id version if present
    gene_id_clean <- gsub("\\..*$", "", gene_id)
    fig_title <- gene_id_clean
    if (!is.null(id_mapping)) {
        # Extract the symbol column name
        symbol_col <- colnames(id_mapping)[2]
        
        # Match by Ensembl ID with and without version
        match_ensembl_ver <- id_mapping[id_mapping$ensembl_gene_id == gene_id, ]
        match_ensembl <- id_mapping[id_mapping$ensembl_gene_id == gene_id_clean, ]
        
        # Match by gene symbol
        match_symbol <- id_mapping[grepl(gene_id, id_mapping[[symbol_col]], ignore.case = TRUE), ]
        
        # Combine matches
        id_mapping <- rbind(match_ensembl_ver, match_ensembl, match_symbol)
        
        id_mapping <- unique(id_mapping)
        
        if (nrow(id_mapping) == 0) {
            stop("No match found for gene ID or symbol:", gene_id)
        } else {
            fig_title <- unique(id_mapping[[symbol_col]])
            gene_id_clean <- unique(id_mapping$ensembl_gene_id)
        }
    }
    
    
    # Calculate ORF usage from fitted models
    usage_df <- calculateUsage(
        data = data,
        gene_id = gene_id_clean
    )

    usage_df$usage <- as.numeric(usage_df$usage)
    if (!is.null(levels)) {
        usage_df$condition <- factor(usage_df$condition, levels = levels)
    }
    
    usage_df <- merge(orfs, usage_df, by = "orf_id")
    if (unique(usage_df$strand) == "+") {
        usage_df <- usage_df[order(-rank(usage_df$orf_type), usage_df$orf_number), ]
        usage_df$orf_number <- paste0(usage_df$orf_type, " (O", usage_df$orf_number, ")")
        usage_df$orf_number <- factor(usage_df$orf_number, levels = unique(usage_df$orf_number))
    } else {
        usage_df <- usage_df[order(-rank(usage_df$orf_type), -rank(usage_df$orf_number)), ]
        usage_df$orf_number <- paste0(usage_df$orf_type, " (O", usage_df$orf_number, ")")
        usage_df$orf_number <- factor(usage_df$orf_number, levels = unique(usage_df$orf_number))
    }
    
    annot_df$orf_id <- as.character(annot_df$orf_id)
    annot_df <- merge(orfs, annot_df, by = "orf_id")
    annot_df$orf_number <- paste0(usage_df$orf_type, " O", annot_df$orf_number, ")")
    annot_df$contrast <- as.character(annot_df$contrast)
    contrasts_to_annotate <- unique(annot_df$contrast)
    
    p <- ggplot2::ggplot(usage_df, ggplot2::aes(x = condition, y = usage, fill = condition)) +
        ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
        ggplot2::facet_wrap(~ orf_number) +
        ggplot2::theme_bw() + 
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), 
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "none"
        ) +
        ggplot2::labs(
            title = fig_title, 
            y = "ORF usage", 
            x = "Condition", 
            fill = "Condition"
        ) + 
        ggplot2::scale_fill_brewer(palette = colors)
    
    for (orf in unique(usage_df$orf_id)) {
        max_usage <- max(usage_df$usage[usage_df$orf_id == orf])
        y_base <- max_usage + 0.01
        y_increment <- 0.08
        y_pos <- y_base
        
        for (contrast in contrasts_to_annotate) {
            lfsr <- get_lfsr_annotation(annot_df, orf, contrast)
            if (as.numeric(lfsr) > dou_signif_thresh) next
            lfsr_formatted <- format.pval(lfsr, digits = 3, eps = .Machine$double.eps)
            # message(
            #     "ORF ID: ", orf, 
            #     "\nContrasts: ", contrast, 
            #     "\nLFSR: ", lfsr_formatted
            # )
            
            conds <- unlist(strsplit(contrast, " - "))
            
            p <- p + ggsignif::geom_signif(
                annotation = lfsr_formatted,
                comparisons = list(conds),
                y_position = y_pos,
                tip_length = 0.01,
                data = subset(usage_df, orf_id == orf)
            )
            
            y_pos <- y_pos + y_increment
        }
    }
    
    return(ggplot2::ggplotGrob(p))
}


#' Plot Venn Diagram of DTE and DOU Significance Overlap
#'
#' @description
#' Generates a Venn diagram (Euler diagram) showing the overlap of
#' significantly differentially translated ORFs (DTE) and differentially
#' used ORFs (DOU). ORFs are classified as significant in DTE only,
#' DOU only, or both. Requires 
#' \href{https://CRAN.R-project.org/package=eulerr}{eulerr}.
#'
#' @param results A data frame containing DTE and DOU results. Must include
#' row names corresponding to ORF identifiers, and columns for adjusted
#' p-values from both tests.
#'
#' @param dou_signif_col Character string specifying the column name for DOU
#' significance values. Should correspond to local false sign rate 
#' (lfsr). Default is \code{"lfsr"}.
#'
#' @param dte_signif_col Character string specifying the column name for DTE
#' adjusted p-values. Default is \code{"padj"}.
#'
#' @param dou_signif_thresh Numeric threshold for DOU LFSR significance.
#' Default is \code{0.05}.
#'
#' @param dte_signif_thresh Numeric threshold for DTE adjusted p-value
#' significance. Default is \code{0.05}.
#'
#' @return A Venn diagram (Euler plot) showing the number of ORFs 
#' significant in DTE only, DOU only, or both.
#'
#' @details
#' Significance is determined using a threshold of 0.05 on the adjusted
#' p-values. The plot uses a color-blind friendly palette and includes
#' counts for each region.
#'
#' @importFrom graphics plot
#'
#' @keywords internal
#'
#' @references
#' Larsson J, Gustafsson P (2018). “A Case Study in Fitting Area-Proportional
#' Euler Diagrams with Ellipses Using eulerr.” In Proceedings of 
#' International Workshop on Set Visualization and Reasoning, volume 2116, 
#' 84–91. \url{https://ceur-ws.org/Vol-2116/paper7.pdf}
#'
plot_venn <- function(
        results,
        dou_signif_col = "lfsr",
        dte_signif_col = "padj",
        dou_signif_thresh = 0.05,
        dte_signif_thresh = 0.05
) {
    
    if (!requireNamespace("eulerr", quietly=TRUE)) {
        stop(
            "Plotting Venn diagrams require the 'eulerr' package. ", 
            "Please install it by running: install.packages('eulerr')"
        )
    }
    
    signif_sig <- !is.na(results[[dte_signif_col]]) &
        results[[dte_signif_col]] < dte_signif_thresh
    fdr_sig <- !is.na(results[[dou_signif_col]]) &
        results[[dou_signif_col]] < dou_signif_thresh
    
    signif_set <- rownames(results)[signif_sig]
    fdr_set <- rownames(results)[fdr_sig]
    
    # Create Euler fit
    fit <- eulerr::euler(c(
        "DTE" = length(setdiff(signif_set, fdr_set)),
        "DOU" = length(setdiff(fdr_set, signif_set)),
        "DTE&DOU" = length(intersect(signif_set, fdr_set))
    ))
    
    # Define color-blind friendly palette
    venn_colors <- c(
        "DTE" = "#0072B2", 
        "DOU" = "#E69F00", 
        "DTE&DOU" = "#CC79A7"
    )
    
    # pdf("dte_dot_venn_cycling_arrest.pdf", 2.5, 2)
    venn <- plot(
        fit,
        fills = list(fill = venn_colors, alpha = 0.6),
        labels = list(font = 2),
        quantities = TRUE,
        main = "Differentially translated ORFs"
    )
    
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, fig = c(0, 1, 0, 1))
    
    return(venn)
}


#' Plot Composite Scatter and Marginal Plots for DTE vs DOU
#'
#' @description
#' Visualizes the relationship between differential translation efficiency 
#' (DTE) and differential ORF usage (DOU) using a scatter plot with 
#' marginal histograms or density curves. Points can be colored by 
#' significance (DTE, DOU, both) or by ORF type (uORF, mORF, dORF). 
#' This plot helps assess overlap and divergence between DTE and DOU 
#' signals across ORFs.
#'
#' @param results A data frame containing DTE and DOU post hoc
#' contrast results. Must include columns for DTE and DOU
#' estimates and their significance values.
#'
#' @param rowdata Optional data frame containing ORF metadata,
#' with row names corresponding to \code{orf_id}. Must
#' include an \code{orf_type} column with values "uORF",
#' "mORF", or "dORF" if \code{color_by = "orf_type"}.
#'
#' @param color_by Character string specifying how to color
#' points:
#' \itemize{
#'     \item{\code{"significance"}: Colors by DTE-only,
#'     DOU-only, both significant}
#'     \item{\code{"orf_type"}: Colors by ORF type (requires
#'     \code{rowdata})}
#' }
#' Default is \code{"significance"}.
#'
#' @param marginal_plot_type Character string specifying the
#' type of marginal plot: \code{"histogram"} or \code{"density"}. 
#' Default is \code{"density"}.
#'
#' @param dou_estimates_col Column name for DOU posterior
#' estimates. Default is \code{"posterior"}.
#'
#' @param dou_signif_col Column name for DOU significance values
#' (e.g., local false sign rate). Default is \code{"lfsr"}.
#'
#' @param dte_estimates_col Column name for DTE estimates
#' (e.g., log2 fold-change). Default is \code{"log2FoldChange"}.
#'
#' @param dte_signif_col Column name for DTE adjusted p-values.
#' Default is \code{"padj"}.
#'
#' @param dou_signif_thresh Numeric threshold for DOU significance. 
#' Default is \code{0.05}.
#'
#' @param dte_signif_thresh Numeric threshold for DTE significance. 
#' Default is \code{0.05}.
#'
#' @param flip_sign Logical; if \code{TRUE}, flips the sign of
#' DOU estimates to align directionality with DTE. Default
#' is \code{FALSE}.
#'
#' @param lhist Integer; number of bins for marginal histograms.
#' Default is \code{20}.
#'
#' @param colors A named list of colors used for plotting. Must include:
#' \itemize{
#'     \item{\code{dte}, \code{dou}, \code{both},
#'     \code{none}: for significance-based coloring}
#'     \item{\code{uorf}, \code{morf}, \code{dorf}: for ORF
#'     type-based coloring}
#' }
#' Each value should be a valid color string or result of
#' \code{adjustcolor()}.
#' 
#' @param legend_position Position of the legend. Options include
#' \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"},
#' \code{"left"}, \code{"topleft"}, \code{"top"},
#' \code{"topright"}, \code{"right"}, \code{"center"}.
#' 
#' @return A composite plot consisting of:
#' \describe{
#'     \item{Scatter plot}{
#'         Displays DTE (log2 fold-change) vs DOU
#'         (log-odds change) estimates, colored by
#'         significance or ORF type.
#'     }
#'     \item{Marginal plots}{
#'         Show the distribution of DTE and DOU estimates
#'         along the top and right margins, using either
#'         histograms or density curves.
#'     }
#' }
#' Additionally, the function returns the Spearman correlation 
#' between DTE and DOU estimates.
#'
#' @details
#' The function uses base R graphics and a custom layout matrix to 
#' arrange the scatter and marginal plots. Coloring by significance 
#' requires DTE and DOU significance columns. Coloring by ORF type 
#' requires a metadata table with \code{orf_type} values.
#'
#' @importFrom graphics plot points legend par layout barplot
#' @importFrom graphics plot.new plot.window hist text lines
#' @importFrom grDevices adjustcolor
#' @importFrom stats cor.test density na.omit
#'
#' @keywords internal
#' 
plot_composite <- function(
        results,
        rowdata = NULL,
        color_by = c("significance", "orf_type"),
        marginal_plot_type = c("histogram", "density"),
        dou_estimates_col = "posterior",
        dou_signif_col = "lfsr",
        dte_estimates_col = "log2FoldChange",
        dte_signif_col = "padj",
        dou_signif_thresh = 0.05,
        dte_signif_thresh = 0.05,
        flip_sign = FALSE,
        lhist = 20,
        colors = list(
            dte = adjustcolor("#0072B2", alpha.f = 0.6),
            dou = adjustcolor("#E69F00", alpha.f = 0.6),
            both = adjustcolor("#CC79A7", alpha.f = 0.6),
            none = adjustcolor("grey80", alpha.f = 0.4),
            uorf = adjustcolor("#D73027", alpha.f = 0.6),
            morf = adjustcolor("#4575B4", alpha.f = 0.6),
            dorf = adjustcolor("#A6A6A6", alpha.f = 0.6)
        ),
        legend_position = "bottomright"
) {
    
    layout(1)
    par(mfrow = c(1, 1))
    
    color_by <- match.arg(color_by)
    marginal_plot_type <- match.arg(marginal_plot_type)
    
    old_par <- par(no.readonly = TRUE)
    on.exit({
        try(par(old_par), silent = TRUE)
        try(layout(1), silent = TRUE)
        try(par(mfrow = c(1, 1)), silent = TRUE)
    }, add = TRUE)
    
    if (nrow(results) == 0) stop("Please provide results as input")
    
    results <- as.data.frame(results)
    
    if (isTRUE(flip_sign)) {
        results[[dou_estimates_col]] <- results[[dou_estimates_col]] * -1
    }
    
    if (!is.null(rowdata)) {
        rowdata <- as.data.frame(rowdata)
        results <- merge(results, rowdata, by.x = "orf_id", by.y = "row.names")
    }
    
    results <- na.omit(results)
    
    # Significance flags
    signif_sig <- !is.na(results[[dte_signif_col]]) & 
        results[[dte_signif_col]] < dte_signif_thresh
    fdr_sig <- !is.na(results[[dou_signif_col]]) & 
        results[[dou_signif_col]] < dou_signif_thresh
    both_sig <- signif_sig & fdr_sig
    signif_only <- signif_sig & !fdr_sig
    fdr_only <- fdr_sig & !signif_sig
    
    # Colors
    col_dte <- colors$dte
    col_dou <- colors$dou
    col_both <- colors$both
    col_none <- colors$none
    
    col_uorfs <- colors$uorf
    col_morfs <- colors$morf
    col_dorfs <- colors$dorf
    
    # Layout
    layMat <- matrix(c(1, 4, 3, 2), ncol = 2)
    layout(
        layMat, 
        widths = c(4 / 5, 1 / 5), 
        heights = c(1 / 5, 4 / 5)
    )
    ospc <- 0.5
    pext <- 4
    bspc <- 1
    par(mar = c(pext, pext, bspc, bspc), oma = rep(ospc, 4))
    
    xlim <- range(results[[dte_estimates_col]], na.rm = TRUE)
    ylim <- range(results[[dou_estimates_col]], na.rm = TRUE)
    
    # Marginal plots
    if (marginal_plot_type == "histogram") {
        par(mar = c(0, pext, 0, 0))
        xhist <- hist(
            results[[dte_estimates_col]], 
            breaks = seq(xlim[1], xlim[2], length.out = lhist), 
            plot = FALSE
        )
        barplot(
            xhist$density, 
            axes = FALSE, 
            space = 0, 
            col = "gray", 
            border = "black"
        )
        
        par(mar = c(pext, 0, 0, 0))
        yhist <- hist(
            results[[dou_estimates_col]], 
            breaks = seq(ylim[1], ylim[2], length.out = lhist), 
            plot = FALSE)
        barplot(
            yhist$density, 
            axes = FALSE, 
            space = 0, 
            col = "gray", 
            border = "black", 
            horiz = TRUE)
    } else if (marginal_plot_type == "density") {
        # Top: DTE distribution
        par(mar = c(0, pext, 0, 0))
        plot.new()
        plot.window(xlim = xlim, ylim = c(0, 1))
        dens_list <- lapply(c("uORF", "mORF", "dORF"), function(orf_type) {
            subset <- results[results$orf_type == orf_type, ]
            density(subset[[dte_estimates_col]], from = xlim[1], to = xlim[2])
        })
        ymax <- max(vapply(dens_list, function(d) max(d$y), numeric(1)))
        
        for (i in seq_along(dens_list)) {
            lines(
                dens_list[[i]]$x,
                dens_list[[i]]$y / ymax,
                col = c(col_uorfs, col_morfs, col_dorfs)[i],
                lwd = 2
            )
        }
        
        # Right: DOU distribution
        par(mar = c(pext, 0, 0, 0))
        plot.new()
        plot.window(xlim = c(0, 1), ylim = ylim)
        dens_list <- lapply(c("uORF", "mORF", "dORF"), function(orf_type) {
            subset <- results[results$orf_type == orf_type, ]
            density(subset[[dou_estimates_col]], from = ylim[1], to = ylim[2])
        })
        ymax <- max(vapply(dens_list, function(d) max(d$y), numeric(1)))
        
        for (i in seq_along(dens_list)) {
            lines(
                dens_list[[i]]$y / ymax,
                dens_list[[i]]$x,
                col = c(col_uorfs, col_morfs, col_dorfs)[i],
                lwd = 2
            )
        }
    }
    
    # Placeholder
    par(mar = c(0, 0, 0, 0))
    plot.new()
    
    # Main scatter plot
    par(mar = c(pext, pext, 0, 0))
    plot(
        results[[dte_estimates_col]],
        results[[dou_estimates_col]],
        col = col_none,
        pch = 16,
        xlab = "log2 fold-change in ORF TE",
        ylab = "log-odds change in ORF usage",
        xlim = xlim,
        ylim = ylim
    )
    
    if (color_by == "significance") {
        points(
            results[[dte_estimates_col]][signif_only], 
            results[[dou_estimates_col]][signif_only], 
            col = col_dte
        )
        points(
            results[[dte_estimates_col]][fdr_only], 
            results[[dou_estimates_col]][fdr_only], 
            col = col_dou
        )
        points(
            results[[dte_estimates_col]][both_sig], 
            results[[dou_estimates_col]][both_sig], 
            col = col_both
        )
        legend(
            legend_position, 
            legend = c("DTE", "DOU", "Both"), 
            col = c(col_dte, col_dou, col_both), 
            pch = 1, bty = "n", 
            inset = 0.04
        )
    } else if (!is.null(rowdata) && "orf_type" %in% colnames(results)) {
        is_uorfs <- results$orf_type == "uORF"
        is_morfs <- results$orf_type == "mORF"
        is_dorfs <- results$orf_type == "dORF"
        points(
            results[[dte_estimates_col]][is_morfs], 
            results[[dou_estimates_col]][is_morfs], 
            col = col_morfs
        )
        points(
            results[[dte_estimates_col]][is_uorfs], 
            results[[dou_estimates_col]][is_uorfs], 
            col = col_uorfs
        )
        points(
            results[[dte_estimates_col]][is_dorfs], 
            results[[dou_estimates_col]][is_dorfs], 
            col = col_dorfs
        )
        legend(
            legend_position, 
            legend = c("uORF", "mORF", "dORF"), 
            col = c(col_uorfs, col_morfs, col_dorfs), 
            pch = 1, 
            bty = "n", 
            inset = 0.04
        )
    }
    
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, fig = c(0, 1, 0, 1))
    
    # Spearman correlation
    correlation_results <- cor.test(
        results[[dte_estimates_col]],
        results[[dou_estimates_col]],
        method = "spearman"
    )
    
    return(correlation_results)
}


#' Volcano Plot for Differential ORF Usage (DOU)
#'
#' @description
#' Generates a volcano plot to visualize differential ORF usage (DOU) 
#' results. The x-axis represents log-odds changes in ORF usage 
#' (effect sizes), and the y-axis shows the negative log10-transformed 
#' local false sign rate (LFSR). Points can be colored either by 
#' significance in DTE/DOU or by ORF type (uORF, mORF, dORF). Optional 
#' gene labeling is supported via an ID mapping table. The y-axis is 
#' capped at the nearest multiple of \code{dou_signif_ceil} for 
#' cleaner plotting.
#'
#' @param results A data frame containing DOU and optionally DTE
#' results. Must include columns for DOU estimates and LFSR.
#' If \code{color_by = "orf_type"}, must also include an
#' \code{orf_type} column with values "uORF", "mORF", or "dORF".
#'
#' @param rowdata Optional data frame containing ORF metadata,
#' with row names corresponding to \code{orf_id}. Used to
#' merge ORF type information.
#'
#' @param id_mapping Optional data frame with gene symbols. Used
#' to label points with gene symbols. Must include
#' \code{ensembl_gene_id} and symbol columns.
#'
#' @param color_by Character string specifying how to color points:
#' \itemize{
#'     \item{\code{"significance"}: Colors by DTE-only,
#'     DOU-only, both significant}
#'     \item{\code{"orf_type"}: Colors by ORF type (requires
#'     \code{orf_type} column)}
#' }
#' Default is \code{"significance"}.
#'
#' @param dou_estimates_col Column name for DOU effect size estimates. 
#' Default is \code{"posterior"}.
#'
#' @param dou_signif_col Column name for DOU significance values (LFSR). 
#' Default is \code{"lfsr"}.
#'
#' @param dte_estimates_col Column name for DTE effect size estimates. 
#' Default is \code{"log2FoldChange"}.
#'
#' @param dte_signif_col Column name for DTE adjusted p-values.
#' Default is \code{"padj"}.
#' 
#' @param dte_signif_thresh Numeric threshold for DTE adjusted p-value
#' significance. Default is \code{0.05}.
#' 
#' @param flip_sign Logical. If \code{TRUE}, flips the sign of DOU 
#' estimates for plotting. Default is \code{FALSE}.
#'
#' @param dou_estimates_thresh Numeric threshold for DOU
#' effect size significance. Default is \code{1}.
#'
#' @param dou_signif_thresh Numeric threshold for DOU LFSR significance.
#' Default is \code{0.05}.
#'
#' @param dou_signif_ceil Numeric value to define the rounding ceiling
#' for -log10(LFSR). The maximum y-axis value will be rounded up 
#' to the nearest multiple of this value. Default is \code{10}.
#'
#' @param dte_signif_thresh Numeric threshold for DTE adjusted
#' p-value significance. Default is \code{0.05}.
#'
#' @param extreme_thresh Optional numeric threshold for
#' labeling extreme points based on -log10(LFSR).
#'
#' @param top_hits Optional numeric. If provided, labels the
#' top N most significant points.
#'
#' @param legend_position Position of the legend. Options include
#' \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"},
#' \code{"left"}, \code{"topleft"}, \code{"top"},
#' \code{"topright"}, \code{"right"}, \code{"center"}.
#'
#' @param colors A named list of colors used for plotting. Must
#' include:
#' \itemize{
#'     \item{\code{dte}, \code{dou}, \code{both},
#'     \code{none}: for significance-based coloring}
#'     \item{\code{uorf}, \code{morf}, \code{dorf}: for ORF
#'     type-based coloring}
#' }
#' Each value should be a valid color string or result of
#' \code{adjustcolor()}.
#'
#' @param verbose Logical. If \code{TRUE}, prints messages about
#' plot scaling and thresholds. Default is \code{TRUE}.
#'
#' @return A volcano plot is displayed in a new graphics device.
#' Points are colored by significance or ORF type, depending
#' on the \code{color_by} argument. Optionally, top or
#' extreme points are labeled.
#'
#' @importFrom graphics plot points abline legend text
#' @importFrom grDevices adjustcolor
#' @importFrom utils head
#'
#' @keywords internal
#' 
plot_volcano <- function(
        results,
        rowdata = NULL,
        id_mapping = NULL,
        color_by = c("significance", "orf_type"),
        dou_estimates_col = "posterior",
        dou_signif_col = "lfsr",
        dte_estimates_col = "log2FoldChange",
        dte_signif_col = "padj",
        dte_signif_thresh = 0.05,
        flip_sign = FALSE,
        dou_estimates_thresh = 1,
        dou_signif_thresh = 0.05,
        dou_signif_ceil = 10,
        extreme_thresh = NULL,
        top_hits = NULL,
        legend_position = "topright",
        colors = list(
            dte = adjustcolor("#0072B2", alpha.f = 0.6),
            dou = adjustcolor("#E69F00", alpha.f = 0.6),
            both = adjustcolor("#CC79A7", alpha.f = 0.6),
            none = adjustcolor("grey80", alpha.f = 0.6),
            uorf = adjustcolor("#D73027", alpha.f = 0.6),
            morf = adjustcolor("#4575B4", alpha.f = 0.6),
            dorf = adjustcolor("#A6A6A6", alpha.f = 0.6)
        ),
        verbose = TRUE
) {
    color_by <- match.arg(color_by)
    
    results <- as.data.frame(results)
    
    if (!is.null(rowdata)) {
        results <- merge(results, rowdata, by.x = "orf_id", by.y = "row.names")
    }
    
    if (!is.null(id_mapping)) {
        genes_unique <- id_mapping[!duplicated(id_mapping$ensembl_gene_id), ]
        results$ensembl_gene_id <- sub("[.:].*", "", results$orf_id)
        results <- merge(
            genes_unique, 
            results, 
            by = "ensembl_gene_id", 
            all.x = TRUE
        )
    }
    
    results <- na.omit(results)
    
    if (isTRUE(flip_sign)) {
        results[[dou_estimates_col]] <- results[[dou_estimates_col]] * -1
    }
    
    results$loglfsr <- -log10(pmax(results[[dou_signif_col]], .Machine$double.eps))
    loglfsr_ceil <- ceiling(
        max(
            results$loglfsr[is.finite(results$loglfsr)], 
            na.rm = TRUE
        ) / dou_signif_ceil
    ) * dou_signif_ceil
    
    if (verbose) {
        message(
            "capping volcano plot y-axis at -log10(", dou_signif_col, ") = ", 
            loglfsr_ceil
        )
    }
    
    results$loglfsr[results$loglfsr > loglfsr_ceil] <- loglfsr_ceil
    
    ylim_max <- max(results$loglfsr, na.rm = TRUE) * 1.1
    ylim_min <- 0
    y_range_step <- (ylim_max - ylim_min) * 0.02
    
    xlim_max <- max(results[[dou_estimates_col]], na.rm = TRUE) * 1.1
    xlim_min <- min(results[[dou_estimates_col]], na.rm = TRUE) * 1.1
    
    par(mfrow = c(1, 1), mar = c(5, 4, 2, 2) + 0.1, fig = c(0, 1, 0, 1))
    
    if (color_by == "significance") {
        signif_sig <- !is.na(results[[dte_signif_col]]) & 
            results[[dte_signif_col]] < dte_signif_thresh
        fdr_sig <- !is.na(results[[dou_signif_col]]) & 
            results[[dou_signif_col]] < dou_signif_thresh
        both_sig <- signif_sig & fdr_sig
        signif_only <- signif_sig & !fdr_sig
        fdr_only <- fdr_sig & !signif_sig
        col_dte <- colors$dte
        col_dou <- colors$dou
        col_both <- colors$both
        col_none <- colors$none
        
        plot(
            results[[dou_estimates_col]],
            results$loglfsr,
            pch = 20,
            col = col_none,
            xlab = "log-odds change in ORF usage",
            ylab = expression(paste("-log"[10], "(LFSR)")),
            xlim = c(xlim_min, xlim_max),
            ylim = c(ylim_min, ylim_max)
        )
        
        points(
            results[[dou_estimates_col]][signif_only], 
            results$loglfsr[signif_only], 
            col = col_dte
        )
        points(
            results[[dou_estimates_col]][fdr_only], 
            results$loglfsr[fdr_only], 
            col = col_dou
        )
        points(
            results[[dou_estimates_col]][both_sig], 
            results$loglfsr[both_sig], 
            col = col_both
        )
        
        legend(
            legend_position, 
            legend = c("DTE", "DOU", "Both"),
            col = c(col_dte, col_dou, col_both), 
            pch = 1, 
            bty = "n", 
            inset = c(0.02, 0.05)
        )
    } else if (color_by == "orf_type" && "orf_type" %in% colnames(results)) {
        col_uorfs <- colors$uorf
        col_morfs <- colors$morf
        col_dorfs <- colors$dorf
        
        is_uorfs <- results$orf_type == "uORF"
        is_morfs <- results$orf_type == "mORF"
        is_dorfs <- results$orf_type == "dORF"
        
        plot(
            results[[dou_estimates_col]][is_dorfs],
            results$loglfsr[is_dorfs],
            pch = 1,
            col = col_dorfs,
            xlab = "log-odds change in ORF usage",
            ylab = expression(paste("-log"[10], "(LFSR)")),
            xlim = c(xlim_min, xlim_max),
            ylim = c(ylim_min, ylim_max)
        )
        
        points(
            results[[dou_estimates_col]][is_uorfs], 
            results$loglfsr[is_uorfs], 
            col = col_uorfs, 
            pch = 1
        )
        points(
            results[[dou_estimates_col]][is_morfs], 
            results$loglfsr[is_morfs], 
            col = col_morfs, 
            pch = 1
        )
        
        legend(
            legend_position, 
            legend = c("uORF", "mORF", "dORF"),
            col = c(col_uorfs, col_morfs, col_dorfs), 
            pch = 1, 
            bty = "n", 
            inset = c(0.02, 0.05)
        )
    }
    
    # Labeling
    if (!is.null(extreme_thresh)) {
        idx_to_label <- which(results$loglfsr > extreme_thresh)
    } else if (is.numeric(top_hits)) {
        idx_to_label <- head(
            order(results$loglfsr, decreasing = TRUE), 
            top_hits
        )
    } else {
        idx_to_label <- NULL
    }
    
    if (!is.null(idx_to_label)) {
        x_offset <- rep(c(-0.1, 0, 0.1), length.out = length(idx_to_label))
        y_offset <- rep(
            c(-y_range_step, 0, y_range_step), 
            length.out = length(idx_to_label)
        )
        
        labels <- if (!is.null(id_mapping)) {
            results[[2]][idx_to_label]
        } else {
            results$orf_id[idx_to_label]
        }
        
        text(
            x = results[[dou_estimates_col]][idx_to_label] + x_offset,
            y = results$loglfsr[idx_to_label] + y_offset,
            labels = labels,
            pos = 3,
            cex = 0.7,
            xpd = TRUE
        )
    }
    
    abline(h = -log10(dou_signif_thresh), col = "black", lty = 2)
    abline(
        v = c(-dou_estimates_thresh, dou_estimates_thresh), 
        col = "black", 
        lty = 2
    )
    
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, fig = c(0, 1, 0, 1))
}


#' Prepare Data and Metadata for DOU Heatmap
#'
#' @description
#' Prepares a matrix of differential ORF usage (DOU) estimates for heatmap
#' visualization, comparing mORFs with either uORFs or dORFs. Filters
#' significant ORFs based on local false sign rate (LFSR), clusters genes
#' based on DOU profiles, and retrieves gene symbols for labeling.
#'
#' @param results A data frame from \code{testDOT} containing DOU estimates
#' and significance values. Must include columns for ORF IDs, gene IDs,
#' effect size estimates, and LFSR values.
#'
#' @param rowdata A data frame containing ORF-level metadata, including ORF
#' type (e.g., mORF, uORF).
#'
#' @param id_mapping Optional data frame with gene symbols (e.g., from 
#' biomaRt). Used to label heatmap rows with gene symbols.
#'
#' @param estimates_col Character string specifying the column name for DOU
#' effect size estimates. Default is \code{"posterior"}.
#' 
#' @param estimates_thresh Numeric threshold for effect size.
#' Default is \code{1}.
#' 
#' @param signif_col Character string specifying the column name for 
#' significance values (LFSR). Default is \code{"lfsr"}.
#'
#' @param signif_thresh Numeric threshold for filtering significant ORFs
#' based on LFSR. Default is \code{0.05}.
#' 
#' @param rank_by Character string specifying how genes should be 
#' ranked for visualization. Options are:
#' - `"score"`: ranks genes by a composite score defined as 
#' `abs(effect size) × (1 - lfsr)`, which prioritizes genes with large 
#' and confident differential ORF usage.
#' - `"significance"`: ranks genes by the median significance values 
#' across ORFs, prioritizing statistical confidence regardless of 
#' effect size. Default is `"significance"`.
#' 
#' @param top_hits Optional numeric. If specified, limits the heatmap to the
#' top N genes ranked by effect size magnitude and significance.
#'
#' @param flip_sign Logical; if \code{TRUE}, flips the sign of DOU estimates
#' to align directionality with differential translation.
#' Default is \code{TRUE}.
#'
#' @param sorf_type Character string specifying the short ORF type to compare
#' with mORFs. Accepts \code{"uORF"} or \code{"dORF"}.
#' Default is \code{"uORF"}.
#' @param colors A named list of three colors for heatmap gradient:
#' \itemize{
#'     \item{\code{low}: color for low values (e.g., "blue")}
#'     \item{\code{middle}: color for mid values (e.g., "white")}
#'     \item{\code{high}: color for high values (e.g., "red")}
#' }
#' 
#' @return A list containing:
#' \describe{
#'     \item{ordered_matrix}{
#'         Matrix of DOU estimates for significant genes, ordered by
#'         hierarchical clustering.
#'     }
#'     \item{row_dend_clean}{
#'         Dendrogram object for clustered gene rows.
#'     }
#'     \item{highlight_df}{
#'         Data frame indicating which ORF-gene tiles to highlight in
#'         the heatmap.
#'     }
#'     \item{gene_labels}{
#'         Vector of gene symbols or Ensembl IDs for heatmap row labeling.
#'     }
#'     \item{color_palette}{
#'         Color palette used for heatmap visualization.
#'     }
#'     \item{color_breaks}{
#'         Breaks used for color scaling in the heatmap.
#'     }
#'     \item{abs_max}{
#'         Maximum absolute value used for symmetric color scaling.
#'     }
#' }
#'
#' @details
#' This function supports visualization of differential ribosome loading
#' across ORFs within genes. It compares mORFs with a specified short ORF
#' type (e.g., uORFs), filters significant ORFs using local false sign rate
#' (LFSR), and clusters genes based on DOU profiles. Gene symbols are
#' optionally retrieved from Ensembl and used for labeling. The output is
#' designed for downstream heatmap plotting.
#'
#' @importFrom stats dist hclust as.dendrogram order.dendrogram complete.cases
#' @importFrom stats setNames median
#' @importFrom utils head
#' @importFrom grDevices colorRampPalette
#'
#' @keywords internal
#'
cistronic_data <- function(
        results,
        rowdata,
        id_mapping = NULL,
        estimates_col = "posterior",
        estimates_thresh = 1,
        signif_col = "lfsr",
        signif_thresh = 0.05,
        rank_by = "significance",
        top_hits = NULL,
        flip_sign = FALSE,
        sorf_type = "uORF",
        colors = list(
            low = "blue", 
            middle = "white", 
            high = "red"
        )
) {

    results <- na.omit(results)
    results$gene_id <- sub("[.:].*", "", results$orf_id)
    
    sig_genes <- results[results[[signif_col]] < signif_thresh, ]$gene_id
    sig_res <- results[results$gene_id %in% sig_genes, ]
    
    rowdata$DOUResults <- NULL
    rowdata$orf_id <- rownames(rowdata)
    rowdata <- as.data.frame(rowdata)[, c("orf_id", "orf_type")]
    
    sig_res <- merge(sig_res, rowdata, by.x = "orf_id")[, c(
        "gene_id",
        "orf_type",
        estimates_col,
        signif_col
    )]
    sig_res <- sig_res[sig_res$orf_type %in% c("mORF", sorf_type), ]
    
    sig_res$capped_lfsr <- pmax(sig_res[[signif_col]], .Machine$double.eps)

    if (rank_by == "score") {
        gene_scores <- aggregate(
            abs(sig_res[[estimates_col]]) * (1 - sig_res$capped_lfsr) ~ 
                sig_res$gene_id,
            FUN = median
        )
    } else if (rank_by == "significance") {
        gene_scores <- aggregate(
            sig_res$capped_lfsr ~ sig_res$gene_id,
            FUN = median
        )
        names(gene_scores) <- c("gene_id", "score")
        gene_scores$score <- -log10(gene_scores$score)
    }
    names(gene_scores) <- c("gene_id", "score")
    sig_res$capped_lfsr <- NULL
    
    sig_res <- merge(gene_scores, sig_res, by = "gene_id")
    sorf_res <- sig_res[sig_res$orf_type == sorf_type, ][, c(
        "gene_id",
        estimates_col
    )]
    sorf_res <- aggregate(
        sorf_res[[estimates_col]] ~ sorf_res$gene_id,
        FUN = median
    )
    names(sorf_res) <- c("gene_id", estimates_col)
    
    morf_res <- sig_res[sig_res$orf_type == "mORF", ][, c(
        "gene_id",
        estimates_col,
        "score"
    )]
    sig_mat <- merge(sorf_res, morf_res, by = "gene_id")
    names(sig_mat) <- c("gene_id", sorf_type, "mORF", "score")
    sig_mat <- sig_mat[!is.na(sig_mat$gene_id), ]
    sig_mat <- sig_mat[order(sig_mat$score, decreasing = TRUE), ]
    rownames(sig_mat) <- sig_mat$gene_id
    
    if (is.numeric(top_hits)) {
        sig_filt <- sig_mat[
            (abs(sig_mat[[sorf_type]]) > estimates_thresh) & (abs(sig_mat$mORF) > estimates_thresh), 
        ]
        if (nrow(sig_filt) > top_hits) {
            sig_mat <- head(sig_filt, top_hits)
        } else {
            sig_mat <- head(sig_mat, top_hits)
        }
    }
    
    sig_mat <- sig_mat[, !(names(sig_mat) %in% c("gene_id", "score"))]
    sig_mat[] <- lapply(sig_mat, as.numeric)
    sig_mat <- as.matrix(sig_mat)
    sig_mat_clean <- sig_mat[complete.cases(sig_mat), ]
    
    row_cluster_clean <- hclust(dist(sig_mat_clean))
    row_dend_clean <- as.dendrogram(row_cluster_clean)
    
    min_val <- min(sig_mat_clean, na.rm = TRUE)
    max_val <- max(sig_mat_clean, na.rm = TRUE)
    abs_max <- max(abs(min_val), abs(max_val))
    
    color_breaks <- seq(-abs_max, abs_max, length.out = 101)
    color_palette <- colorRampPalette(
        c(colors$low, colors$middle, colors$high)
    )(100)
    
    row_order <- order.dendrogram(row_dend_clean)
    ordered_matrix <- sig_mat_clean[row_order, ]
    sig_res_filtered <- sig_res[
        sig_res$orf_type %in%
            colnames(ordered_matrix) &
            sig_res$gene_id %in% rownames(ordered_matrix),
    ]
    group_to_row <- setNames(
        seq_len(nrow(ordered_matrix)),
        rownames(ordered_matrix)
    )
    label_to_col <- setNames(
        seq_len(ncol(ordered_matrix)),
        colnames(ordered_matrix)
    )
    highlight_df <- data.frame(
        row = group_to_row[sig_res_filtered$gene_id],
        col = label_to_col[sig_res_filtered$orf_type],
        color = "black"
    )
    highlight_df <- highlight_df[complete.cases(highlight_df), ]
    
    ensembl_ids <- sub("[.:].*", "", rownames(ordered_matrix))
    
    if (!is.null(id_mapping)) {
        genes_unique <- id_mapping[!duplicated(id_mapping$ensembl_gene_id), ]
        genes_unique_sorted <- genes_unique[
            match(ensembl_ids, genes_unique$ensembl_gene_id),
        ]
        gene_labels <- genes_unique_sorted[[2]]
    } else {
        gene_labels <- ensembl_ids
    }
    
    if (isTRUE(flip_sign)) {
        ordered_matrix <- ordered_matrix * -1
    }
    
    return(list(
        ordered_matrix = ordered_matrix,
        row_dend_clean = row_dend_clean,
        highlight_df = highlight_df,
        gene_labels = gene_labels,
        color_palette = color_palette,
        color_breaks = color_breaks,
        abs_max = abs_max
    ))
}


#' Plot Heatmap for Differential ORF Usage (DOU)
#'
#' @description
#' Generates a heatmap of differential ORF usage (DOU) estimates with
#' hierarchical clustering and optional tile highlighting for significant ORFs.
#' The function uses base R graphics and expects all required inputs to be
#' passed via a named list produced by \code{cistronic_data()}.
#'
#' @param paired_data A named list containing heatmap data and metadata,
#' typically the output from \code{cistronic_data()}. Must include the
#' following elements:
#' \describe{
#'     \item{ordered_matrix}{
#'         Matrix of DOU estimates for significant genes, ordered by
#'         hierarchical clustering.
#'     }
#'     \item{row_dend_clean}{
#'         Dendrogram object for clustered gene rows.
#'     }
#'     \item{highlight_df}{
#'         Data frame indicating which ORF-gene tiles to highlight in
#'         the heatmap.
#'     }
#'     \item{gene_labels}{
#'         Vector of gene symbols or Ensembl IDs for heatmap row labeling.
#'     }
#'     \item{color_palette}{
#'         Color palette used for heatmap visualization.
#'     }
#'     \item{color_breaks}{
#'         Breaks used for color scaling in the heatmap.
#'     }
#'     \item{abs_max}{
#'         Maximum absolute value used for symmetric color scaling.
#'     }
#' }
#'
#' @return A heatmap with a dendrogram and color key.
#'
#' @details
#' The heatmap displays log-odds changes in ORF usage across genes, comparing
#' mORFs with short ORFs (e.g., uORFs or dORFs). Significant ORFs are
#' highlighted using rectangles. The color scale is symmetric and centered at
#' zero. The function uses base R plotting functions and is intended for
#' interactive or scripted visualization.
#'
#' @importFrom graphics axis image mtext rect par layout plot
#'
#' @keywords internal
#'
plot_heatmap <- function(paired_data) {
    ordered_matrix <- paired_data$ordered_matrix
    row_dend_clean <- paired_data$row_dend_clean
    highlight_df <- paired_data$highlight_df
    gene_labels <- paired_data$gene_labels
    color_palette <- paired_data$color_palette
    color_breaks <- paired_data$color_breaks
    abs_max <- paired_data$abs_max
    
    old_par <- par(no.readonly = TRUE)
    on.exit({
        try(par(old_par), silent = TRUE)
        try(layout(1), silent = TRUE)
        try(par(mfrow = c(1, 1)), silent = TRUE)
    }, add = TRUE)
    
    n_rows <- nrow(ordered_matrix)
    n_cols <- ncol(ordered_matrix)
    
    # Dynamic scaling
    cex_row <- ifelse(n_rows > 50, 0.5, 0.8)
    cex_col <- ifelse(n_cols > 50, 0.5, 0.8)
    right_margin <- max(5, min(5, n_rows / 5))
    layout_widths <- c(2, max(6, n_cols / 10))
    
    layout(matrix(c(1, 2), nrow = 1), widths = layout_widths)
    
    ## Panel 1: Dendrogram
    par(mar = c(5, 1.5, 4, 0))
    # par(mar = c(4, 1, 3, 1))
    plot(row_dend_clean, horiz = TRUE, axes = FALSE, yaxs = "i")
    
    ## Panel 2: Heatmap
    par(mar = c(5, 0.5, 4, right_margin))
    image(
        x = seq_len(ncol(ordered_matrix)),
        y = seq_len(nrow(ordered_matrix)),
        z = t(ordered_matrix),
        col = color_palette,
        breaks = color_breaks,
        axes = FALSE,
        xlab = ""
    )
    
    for (i in seq_len(nrow(highlight_df))) {
        rect(
            xleft = highlight_df$col[i] - 0.5,
            xright = highlight_df$col[i] + 0.5,
            ybottom = highlight_df$row[i] - 0.5,
            ytop = highlight_df$row[i] + 0.5,
            border = highlight_df$color[i],
            lwd = 1.5
        )
    }
    
    axis(
        1,
        at = seq_len(ncol(ordered_matrix)),
        labels = colnames(ordered_matrix),
        las = 2,
        cex.axis = cex_col
    )
    axis(
        4,
        at = seq_len(nrow(ordered_matrix)),
        labels = gene_labels,
        las = 2,
        cex.axis = cex_row
    )
    
    mtext(
        "Differential ORF usage",
        side = 3,
        line = 1.5,
        cex = 1.2,
        adj = 0.5,
        font = 2
    )
    
    ## Dynamic Color Key Placement (Top-Left)
    dend_width <- layout_widths[1] / sum(layout_widths)
    legend_width <- 0.15
    legend_height <- 0.08
    legend_left <- 0.02
    legend_right <- legend_left + legend_width
    legend_bottom <- 0.92
    legend_top <- legend_bottom + legend_height
    
    par(
        fig = c(legend_left, legend_right, legend_bottom, legend_top),
        new = TRUE,
        mar = c(0, 0, 2, 0)
    )
    image(
        z = t(matrix(seq(-abs_max, abs_max, length.out = 100), nrow = 1)),
        col = color_palette,
        axes = FALSE
    )
    
    tick_vals <- pretty(c(-abs_max, abs_max), n = 5)
    tick_pos <- (tick_vals - (-abs_max)) / (2 * abs_max)
    axis(1, at = tick_pos, labels = tick_vals, las = 1, cex.axis = 0.7)
    mtext("log-odds", side = 1, line = 2, cex = 0.8)
    
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, fig = c(0, 1, 0, 1))
}



#' Reset and Recover Graphics Device
#'
#' @description
#' A wrapper function that resets the graphics layout and optionally
#' opens a new graphics device before executing a plotting function.
#' It is designed to prevent layout persistence and recover from
#' common base R graphics errors such as "invalid graphics state" or
#' "figure margins too large", which can occur when plotting multiple
#' figures sequentially.
#'
#' @param plot_fn A function object that performs plotting. This
#' function will be executed within a clean graphics context.
#'
#' @param force_new_device Logical; if \code{TRUE}, attempts to close
#' and reopen the graphics device before plotting. This helps isolate
#' plots and prevent layout or overlay issues. Default is \code{TRUE}.
#'
#' @return Invisibly returns the result of \code{plot_fn()}, if
#' successful. If an error occurs and recovery is triggered, the
#' function attempts to replot in a fresh graphics device.
#'
#' @details
#' This function is useful when plotting multiple complex figures
#' (e.g., composite plots with \code{layout()}) in sequence. It resets
#' the layout using \code{layout(1)} and \code{par(mfrow = c(1, 1))}
#' before plotting, and catches common graphics errors to retry the
#' plot in a clean device.
#'
#' @importFrom graphics layout par 
#' @importFrom grDevices dev.cur dev.new dev.off
#' 
#' @keywords internal
reset_graphics <- function(plot_fn, force_new_device = TRUE) {
    tryCatch({
        if (force_new_device && dev.cur() > 1) {
            msg <- sprintf(
                "resetting graphics device (force_new_device = %s)",
                force_new_device
            )
            message(msg)
            try(dev.off(), silent = TRUE)
            try(dev.new(), silent = TRUE)
        }
        layout(1)
        par(mfrow = c(1, 1))
        invisible(plot_fn())
    }, error = function(e) {
        if (grepl(
            "invalid graphics state|figure margins too large",
            e$message,
            ignore.case = TRUE
        )) {
            if (dev.cur() > 1) {
                try(dev.off(), silent = TRUE)
            }
            try(dev.new(), silent = TRUE)
            layout(1)
            par(mfrow = c(1, 1))
            invisible(plot_fn())
        } else {
            message("Plotting failed: ", e$message)
        }
    })
}


#' Generate Differential ORF Translation (DOT) Visualization Suite
#'
#' @description
#' Generates a suite of visualizations to explore Differential ORF Usage 
#' (DOU) and Translation Efficiency (DTE) results. Supports volcano 
#' plots, Venn diagrams, composite scatter plots, heatmaps, and usage 
#' plots. Integrates Ensembl gene symbols and highlights significant 
#' ORFs based on empirical Bayes shrinkage (via the \pkg{ashr} package's 
#' \code{\link[ashr]{ash}} function).
#' 
#' @param plot_type Type of plot to generate. Options include 
#' \code{"volcano"}, \code{"venn"}, \code{"composite"}, \code{"heatmap"}, 
#' and \code{"usage"}. Default is \code{"volcano"}.
#'
#' @param results A data frame containing DOU and DTE estimates and 
#' significance values. Required for all plots except \code{"usage"}. 
#' Must include ORF-level identifiers and columns 
#' specified in \code{dou_params} and \code{dte_params}.
#'
#' @param data A \code{\link{DOUData-class}} object required for all
#' plots except for \code{"venn"}. Should contain \code{DOUResults} in 
#' \code{rowData()} and post hoc results in the \code{interaction} slot.
#'
#' @param id_mapping Optional input controlling gene symbol annotation.
#' Can be one of:
#' \itemize{
#'     \item \code{FALSE} (default): disables annotation.
#'     \item \code{TRUE}: triggers automatic annotation via BioMart 
#'     using \code{mapIDs()}.
#'     \item a \code{data.frame}: a precomputed mapping of Ensembl IDs 
#'     to gene symbols, which can be reused across plots.
#' }
#' Used in volcano and heatmap plots to annotate genes with symbols. 
#' If \code{TRUE}, annotation is attempted and fallback to Ensembl IDs 
#' occurs if symbol coverage is low. Default is \code{FALSE}
#'
#' @param include_go Logical; if \code{TRUE}, includes GO annotations 
#' in \code{id_mapping}. Default is \code{FALSE}.
#'
#' @param gene_id Character string specifying the gene ID for usage 
#' plots. Default is \code{NULL}.
#'
#' @param dou_params A named list of parameters for DOU filtering 
#' and display:
#' \itemize{
#'     \item{\code{est_col}: Column name for DOU effect size (e.g., 
#'     \code{"posterior"})}
#'     \item{\code{est_thresh}: Threshold for effect size significance}
#'     \item{\code{signif_col}: Column name for DOU significance 
#'     (e.g., \code{"lfsr"})}
#'     \item{\code{signif_thresh}: Threshold for significance}
#'     \item{\code{signif_ceil}: Ceiling for -log10 significance 
#'     axis}
#'     \item{\code{extreme_thresh}: Optional threshold for labeling 
#'     extreme points}
#' }
#'
#' @param dte_params A named list of parameters for DTE filtering:
#' \itemize{
#'     \item{\code{est_col}: Column name for DTE effect size (e.g., 
#'     \code{"log2FoldChange"})}
#'     \item{\code{signif_col}: Column name for DTE significance (e.g., 
#'     \code{"padj"})}
#'     \item{\code{signif_thresh}: Threshold for significance}
#' }
#'
#' @param plot_params A named list controlling labeling of top hits:
#' \itemize{
#'     \item{\code{top_hits}: Number of top hits to label or show in 
#'     volcano plot and heatmap}
#'     \item{\code{color_by}:  How to color points. Options:
#'         \itemize{
#'             \item{\code{"significance"}: Based on DTE-only, DOU-only, 
#'             or both}
#'             \item{\code{"orf_type"}: Based on ORF type (requires 
#'             \code{rowdata})}
#'         }
#'.    }
#'     \item{\code{rank_by}: Ranking method (\code{"significance"} or 
#'     \code{"score"})}
#'     \item{\code{legend_position}: Character string specifying the
#'     position of the legend in the volcano and composite scatter plots. 
#'     Options include: \code{"bottomright"}, \code{"bottom"}, 
#'     \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, 
#'     \code{"top"}, \code{"topright"}, \code{"right"}, \code{"center"}.
#'     Default is \code{"topright"}.}
#'     \item{\code{order_by}: Optional character vector 
#'     specifying the order of conditions on the \code{"usage"} 
#'     plot x-axis.}
#'     \item{\code{flip_sign}: Logical; if \code{TRUE}, flips the sign of DOU 
#' estimates to align directionality with DTE. Default is \code{TRUE}.}
#' }
#'
#' @param annotation_params A named list for BioMart annotation:
#' \itemize{
#'     \item{\code{sorf_type}: Short ORF type (\code{"uORF"} or 
#'     \code{"dORF"})}
#'     \item{\code{dataset}: Ensembl dataset name}
#'     \item{\code{symbol_col}: Column name for gene symbols}
#'     \item{\code{mart_source}: BioMart source (\code{"ensembl"}, 
#'     \code{"plants"}, etc.)}
#' }
#'
#' @param colors A named list of colors used across plots:
#' \itemize{
#'     \item{\code{dte}, \code{dou}, \code{both}, \code{none}: for 
#'     significance-based coloring}
#'     \item{\code{uorf}, \code{morf}, \code{dorf}: for ORF 
#'     type-based coloring}
#'     \item{\code{low}, \code{middle}, \code{high}: for heatmap 
#'     gradient}
#'     \item{\code{usage}: for condition-specific coloring}
#' }
#' 
#' @param force_new_device Logical; if \code{TRUE}, detects graphics
#' error and resets graphics state unconditionally. 
#' Default is \code{TRUE}.
#'
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#' Default is \code{TRUE}.
#'
#' @return A data frame containing gene symbols retrieved from Ensembl, 
#' used for labeling and heatmap visualization.
#'
#' @export
#' @details
#' This function orchestrates multiple visualization components to 
#' explore differential translation across ORFs. It uses 
#' \code{\link{testDOU}} output to identify significant ORFs, 
#' retrieves gene symbols via  \code{\link[biomaRt]{getBM}}, and generates 
#' plots to summarize DOU and DTE relationships. The composite scatter 
#' plot includes marginal distributions by ORF type, helping to 
#' visualize the overlap and divergence between DTE and DOU signals. 
#' The volcano plot highlights extreme and top-ranked ORFs, while 
#' the heatmap summarizes DOU across top genes.
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom graphics par layout
#' @importFrom grid grid.newpage grid.draw
#' @importFrom utils modifyList
#' 
#' @export
#' @examples
#' # Example ORF-level results
#' results_df <- data.frame(
#'   orf_id = c(
#'     "ENSG00000139618.19:O001",
#'     "ENSG00000139618.19:O002",
#'     "ENSG00000157764.15:O003"
#'   ),
#'   lfsr = c(0.01, 0.2, 0.03),
#'   padj = c(0.02, 0.01, 0.1)
#' )
#'
#' plotDOT(plot_type = "venn", results = results_df)
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
#' Larsson J, Gustafsson P (2018). “A Case Study in Fitting
#' Area-Proportional Euler Diagrams with Ellipses Using eulerr.”
#' In Proceedings of International Workshop on Set Visualization
#' and Reasoning, volume 2116, 84–91.
#' \url{https://ceur-ws.org/Vol-2116/paper7.pdf}
#'
plotDOT <- function(
        plot_type = "volcano",
        results = NULL,
        data = NULL,
        id_mapping = FALSE,
        include_go = FALSE,
        gene_id = NULL,
        dou_params = list(
            est_col = "posterior",
            est_thresh = 1,
            signif_col = "lfsr",
            signif_thresh = 0.05,
            signif_ceil = 10,
            extreme_thresh = NULL
        ),
        dte_params = list(
            est_col = "log2FoldChange",
            signif_col = "padj",
            signif_thresh = 0.05
        ),
        plot_params = list(
            top_hits = 20,
            color_by = "significance",
            rank_by = "significance",
            legend_position = "topright",
            order_by = NULL,
            flip_sign = FALSE
        ),
        annotation_params = list(
            sorf_type = "uORF",
            dataset = "hsapiens_gene_ensembl",
            symbol_col = "hgnc_symbol",
            mart_source = "ensembl"
        ),
        colors = list(
            dte = adjustcolor("#0072B2", alpha.f = 0.6),
            dou = adjustcolor("#E69F00", alpha.f = 0.6),
            both = adjustcolor("#CC79A7", alpha.f = 0.6),
            none = adjustcolor("grey80", alpha.f = 0.6),
            uorf = adjustcolor("#D73027", alpha.f = 0.6),
            morf = adjustcolor("#4575B4", alpha.f = 0.6),
            dorf = adjustcolor("#A6A6A6", alpha.f = 0.6),
            low = "blue", 
            middle = "white", 
            high = "red",
            usage = "Set2"
        ),
        force_new_device = TRUE,
        verbose = TRUE
) {
    
    # Validate plot_types
    valid_types <- c("venn", "composite", "volcano", "heatmap", "usage")
    plot_type <- match.arg(plot_type, choices = valid_types)
    
    
    # Validate dou_params
    default_dou_params <- list(
        est_col = "posterior",
        est_thresh = 1,
        signif_col = "lfsr",
        signif_thresh = 0.05,
        signif_ceil = 10,
        extreme_thresh = NULL
    )
    
    valid_dou_params <- list(
        est_col = c("betahat", "posterior"),
        signif_col = c("waldpval", "waldpadj", "qvalue", "lfdr", "lfsr")
    )
    
    # Merge defaults with user overrides
    dou_params <- modifyList(default_dou_params, dou_params)
    dou_params$est_col <- match.arg(
        dou_params$est_col, 
        choices = valid_dou_params$est_col
    )
    dou_params$signif_col <- match.arg(
        dou_params$signif_col, 
        choices = valid_dou_params$signif_col
    )
    
    
    # Validate dte_params
    default_dte_params <- list(
        est_col = "log2FoldChange",
        signif_col = "padj",
        signif_thresh = 0.05
    )
    
    valid_dte_params <- list(
        signif_col = c("pvalue", "padj")
    )
    
    dte_params$signif_col <- match.arg(
        dte_params$signif_col, 
        choices = valid_dte_params$signif_col
    )
    
    
    # Validate plot_params
    default_plot_params <- list(
        top_hits = 20,
        color_by = "significance",
        rank_by = "significance",
        legend_position = "topright",
        order_by = NULL,
        flip_sign = FALSE
    )
    
    valid_plot_params <- list(
        color_by = c("significance", "orf_type"),
        rank_by = c("significance", "score"),
        legend_position = c(
            "bottomright", "bottom", "bottomleft", "left", 
            "topleft", "top", "topright", "right", "center"
        )
    )
    
    plot_params <- modifyList(default_plot_params, plot_params)
    plot_params$color_by <- match.arg(
        plot_params$color_by, 
        choices = valid_plot_params$color_by
    )
    plot_params$rank_by <- match.arg(
        plot_params$rank_by, 
        choices = valid_plot_params$rank_by
    )
    plot_params$legend_position <- match.arg(
        plot_params$legend_position, 
        choices = valid_plot_params$legend_position
    )
    
    marginal_plot_type <- if (plot_params$color_by == "significance") "histogram" else "density"
    
    
    # Validate annotation_params
    default_annotation_params <- list(
        sorf_type = "uORF",
        dataset = "hsapiens_gene_ensembl",
        symbol_col = "hgnc_symbol",
        mart_source = "ensembl"
    )
    
    valid_source <- list(
        sorf_type = c("uORF", "dORF"),
        mart_source = c(
            "ensembl",
            "plants",
            "fungi",
            "protists",
            "metazoa",
            "bacteria"
        )
    )
    
    annotation_params <- modifyList(default_annotation_params, annotation_params)
    annotation_params$sorf_type <- match.arg(
        annotation_params$sorf_type, 
        choices = annotation_params$sorf_type
    )
    annotation_params$mart_source <- match.arg(
        annotation_params$mart_source, 
        choices = annotation_params$mart_source
    )
    
    
    if (!is.null(results)) {
        results <- as.data.frame(results)
    } else if (plot_type %in% c("venn", "composite", "volcano", "heatmap")) {
        stop("A results dataframe is required for plotting: ", plot_type)
    }
    
    if (!is.null(data) && (inherits(data, "DOUData"))) {
        rowdata <- SummarizedExperiment::rowData(data)
    } else {
        rowdata <- NULL
    }
    
    if (any(c("heatmap", "volcano") %in% plot_type) && (identical(id_mapping, TRUE))) {
        if (verbose) message("retrieving gene annotation from BioMart")
        id_mapping_attempt <- tryCatch({
            mapIDs(
                ensembl_ids = get_significant_genes(results),
                dataset = annotation_params$dataset,
                symbol_col = annotation_params$symbol_col,
                mart_source = annotation_params$mart_source,
                include_go = include_go
            )
        }, error = function(e) {
            if (verbose) message("Failed to retrieve gene annotation: ", e$message)
            NULL
        })
        if (!is.null(id_mapping_attempt)) {
            retrieved_symbols <- id_mapping_attempt[[2]][nzchar(id_mapping_attempt[[2]])]
            if (verbose) message("retrieved ", length(retrieved_symbols), " gene symbols")
            frac_retrieved <- length(retrieved_symbols) / nrow(id_mapping_attempt)
            prop_retrieved <- round(frac_retrieved * 100, 3)
            if (frac_retrieved < 0.5) {
                warning(
                    "Retrieved only ", 
                    prop_retrieved, 
                    "% gene symbols; falling back to gene IDs"
                )
                id_mapping <- NULL
                
            } else {
                id_mapping <- id_mapping_attempt
            }
        }
    } else if (any(c("heatmap", "volcano") %in% plot_type) && isFALSE(id_mapping)) {
        id_mapping <- NULL 
    } else if (!isFALSE(id_mapping) && !inherits(id_mapping, "data.frame")) {
        stop("'id_mapping' must be FALSE, TRUE, or a data.frame")
    }
    
    if (plot_type == "venn") {
        reset_graphics(function() {
            grid::grid.newpage()
            grid::grid.draw(plot_venn(
                results = results,
                dou_signif_col = dou_params$signif_col,
                dte_signif_col = dte_params$signif_col,
                dou_signif_thresh = dou_params$signif_thresh,
                dte_signif_thresh = dte_params$signif_thresh
            ))
            if (verbose) message("Venn diagram plotted")
        }, force_new_device = force_new_device)
    }
    
    if (plot_type == "composite") {
        if (plot_params$color_by == "orf_type" && is.null(rowdata)) {
            warning("DOUData is missing; falling back to color_by = 'significance'.")
            plot_params$color_by <- "significance"
            marginal_plot_type <- "histogram"
        }
        reset_graphics(function() {
            correlation_results <- plot_composite(
                results = results,
                rowdata = rowdata,
                color_by = plot_params$color_by,
                marginal_plot_type = marginal_plot_type,
                dou_estimates_col = dou_params$est_col,
                dou_signif_col = dou_params$signif_col,
                dte_estimates_col = dte_params$est_col,
                dte_signif_col = dte_params$signif_col,
                flip_sign = plot_params$flip_sign,
                lhist = 20,
                colors = colors,
                legend_position = plot_params$legend_position
            )
            message(
                "Spearman correlation between DOU and DTE estimates: ",
                round(correlation_results$estimate[["rho"]], 3),
                ", p-value: ",
                format.pval(correlation_results$p.value, digits = 3, eps = .Machine$double.eps)
            )
            if (verbose) message("composite plot colored by significance plotted")
        }, force_new_device = force_new_device)
    }
    
    if (plot_type == "volcano") {
        if (plot_params$color_by == "orf_type" && is.null(rowdata)) {
            warning("DOUData is missing; falling back to color_by = 'significance'.")
            plot_params$color_by <- "significance"
        }
        reset_graphics(function() {
            plot_volcano(
                results = results,
                rowdata = rowdata,
                id_mapping = id_mapping,
                color_by = plot_params$color_by,
                dou_estimates_col = dou_params$est_col,
                dou_signif_col = dou_params$signif_col,
                dte_estimates_col = dte_params$est_col,
                dte_signif_col = dte_params$signif_col,
                dou_estimates_thresh = dou_params$est_thresh,
                dou_signif_thresh = dou_params$signif_thresh,
                dou_signif_ceil = dou_params$signif_ceil,
                extreme_thresh = dou_params$extreme_thresh,
                top_hits = plot_params$top_hits,
                colors = colors,
                legend_position = plot_params$legend_position,
                verbose = verbose
            )
            if (verbose) message("volcano plot colored by significance plotted")
        }, force_new_device = force_new_device)
    }
    
    if (plot_type == "heatmap") {
        if (verbose) message("plotting heatmap for the top ", plot_params$top_hits, " DOU genes")
        paired_data <- tryCatch({
            cistronic_data(
                results = results,
                rowdata = rowdata,
                id_mapping = id_mapping,
                estimates_col = dou_params$est_col,
                estimates_thresh = dou_params$est_thresh,
                signif_col = dou_params$signif_col,
                signif_thresh = dou_params$signif_thresh,
                rank_by = plot_params$rank_by,
                flip_sign = plot_params$flip_sign,
                sorf_type = annotation_params$sorf_type,
                top_hits = plot_params$top_hits,
                colors = colors
            )
        }, error = function(e) {
            if (verbose) message("Failed to prepare heatmap data: ", e$message)
            NULL
        })
        
        if (!is.null(paired_data)) {
            reset_graphics(function() {
                plot_heatmap(paired_data)
                if (verbose) message("heatmap plotted")
            }, force_new_device = force_new_device)
        }
    }
    
    if (plot_type == "usage") {
        
        if (is.null(gene_id)) stop("Please provide gene_id.")
        
        if (!is.character(gene_id)) stop("'gene_id' must be a character vector.")
        
        reset_graphics(function() {
            grid::grid.newpage()
            if (isFALSE(id_mapping)) {
                p <- plot_orf_usage(
                    data = data, 
                    gene_id = gene_id, 
                    levels = plot_params$order_by, 
                    dou_signif_thresh = dou_params$signif_thresh,
                    colors = colors$usage
                )
            } else if (!inherits(id_mapping, "data.frame")) {
                p <- plot_orf_usage(
                    data = data, 
                    gene_id = gene_id, 
                    id_mapping = id_mapping,
                    levels = plot_params$order_by, 
                    dou_signif_thresh = dou_params$signif_thresh,
                    colors = colors$usage
                )
            }
            grid::grid.draw(p)
            if (verbose) message("ORF usage plot for ", gene_id, " plotted")
        }, force_new_device = force_new_device)
    }
    return(invisible(id_mapping))
}
