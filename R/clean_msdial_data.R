
#' Combine pos and neg files from MSDial and filter peaks according to user parameters
#'
#' @eval recurrent_params("filter_blk", "filter_blk_threshold", "filter_mz", "filter_rsd", "filter_rsd_threshold", "filter_rmd", "filter_rmd_range", "filter_blk_ghost_peaks", "threshold_rt", "threshold_mz", "user_neg_neutral_refs,user_pos_neutral_refs,user_pos_adducts_refs,user_neg_adducts_refs")
#' @param compute_pearson_correlation Compute Pearson correlation between peaks to detect clusters.
#' @param pearson_correlation_threshold Ignore links having a Pearson correlation < threshold (default: 0.8).
#' @param pearson_p_value Ignore links having a non-significative Pearson correlation (default: 0.05).
#'
#' @section Architecture needed by clean_msdial_data in the project directory:
#' \itemize{
#'   \item pos
#'   \describe{
#'     \item{Normalized-<...>.txt}{Positive peaks info exported from MSDial}
#'     \item{peaks}{Positive peaks files exported from MSDial}
#'   }
#'   \item neg
#'   \describe{
#'     \item{Normalized-<...>.txt}{Negative peaks info exported from MSDial}
#'     \item{peaks}{Negative peaks files exported from MSDial}
#'   }
#' }
#'
#' @section Output files of clean_msdial_data in the project directory:
#' \itemize{
#'   \item pos
#'   \describe{
#'     \item{filtered_peaks}{Positive peaks files remaining after cleaning, copied from the folder pos/peaks}
#'   }
#'   \item neg
#'   \describe{
#'     \item{filtered_peaks}{Negative peaks files remaining after cleaning, copied from the folder neg/peaks}
#'   }
#'   \item intermediary_data
#'   \describe{
#'     \item{samples.csv}{Information about samples}
#'     \item{adducts.graphml}{Graph of adducts relations between peaks}
#'     \item{clusters.csv}{MSDial data with clusters information}
#'     \item{clusters.graphml}{Graph containing clusters and MSDial data}
#'     \item{MS_peaks-final.csv}{Peaks remaining after all filters have been applied}
#'     \item{links-*.csv}{Files containing information about links between peaks.}
#'   }
#' }
#'
#' @export
clean_msdial_data <- function(filter_blk = TRUE,
                              filter_blk_threshold = 0.8,
                              filter_blk_ghost_peaks = TRUE,
                              filter_mz = TRUE,
                              filter_rsd = TRUE,
                              filter_rsd_threshold = 30,
                              filter_rmd = TRUE,
                              filter_rmd_range = c(50, 3000),
                              threshold_mz = 0.05,
                              threshold_rt = 0.1,
                              user_pos_adducts_refs = NA,
                              user_neg_adducts_refs = NA,
                              user_pos_neutral_refs = NA,
                              user_neg_neutral_refs = NA,
                              compute_pearson_correlation = FALSE,
                              pearson_correlation_threshold = 0.8,
                              pearson_p_value = 0.05) {

    check_input_parameters_msdial_data(filter_blk                  = filter_blk,
                                       filter_blk_threshold        = filter_blk_threshold,
                                       filter_mz                   = filter_mz,
                                       filter_rsd                  = filter_rsd,
                                       filter_rsd_threshold        = filter_rsd_threshold,
                                       filter_blk_ghost_peaks      = filter_blk_ghost_peaks,
                                       filter_rmd                  = filter_rmd,
                                       filter_rmd_range            = filter_rmd_range,
                                       threshold_mz                = threshold_mz,
                                       threshold_rt                = threshold_rt,
                                       compute_pearson_correlation = compute_pearson_correlation,
                                       pearson_threshold           = pearson_correlation_threshold,
                                       pearson_p_value             = pearson_p_value,
                                       references_adduct_pos       = user_pos_adducts_refs,
                                       references_adduct_neg       = user_neg_adducts_refs,
                                       references_neutral_pos      = user_pos_neutral_refs,
                                       references_neutral_neg      = user_neg_neutral_refs)

    check_architecture_for_clean_msdial_data(filter_blk)

    print_message("*** Treating ", get("analysis_directory", envir = mscleanrCache), " ***")

    params <- as.list(environment())
    for (ref in c("user_pos_adducts_refs", "user_neg_adducts_refs", "user_pos_neutral_refs", "user_neg_neutral_refs")) {
        params[[ref]] <- ifelse(is.na(params[[ref]]), "package", "personalized")
    }
    export_params(params)

    # MSDIAL
    msdial <- import_msdial_data(filter_blk,
                                 filter_blk_threshold,
                                 filter_mz,
                                 filter_rsd,
                                 filter_rsd_threshold,
                                 filter_blk_ghost_peaks,
                                 filter_rmd,
                                 filter_rmd_range,
                                 threshold_mz)
    msdial_peak_data <- msdial$peak_data[!grepl('NA', rownames(msdial$peak_data)),]  # Fix for occasional small bug
    print_peaks_status(msdial_peak_data, "MSDial peaks after filtering:")
    msdial_links <- msdial$links
    print_message("MSDial links: ", nrow(msdial_links))
    print_message(nrow(msdial$identified_peaks), " peaks identified by MSDial")
    export_data(msdial$identified_peaks, "msdial_identified_peaks")
    msdial_detected_adducts <- msdial$detected_adducts
    export_data(msdial_detected_adducts, "msdial_detected_adducts")
    samples <- msdial$samples
    export_data(samples, "samples")
    rm(msdial)


    # PEARSON (separated pos/neg)
    if (compute_pearson_correlation) {
        cor_links <- data.frame()
        for (source in get("analysis_modes", envir = mscleanrCache)) {
            tmp_peak_data <- msdial_peak_data[msdial_peak_data$source == source, samples$Column_name]
            if(nrow(samples) < 5) {
                print_warning("Less than 5 samples used, correlation links will not be as robust.")
                tmp_cols <- data.frame(replicate(5, tmp_peak_data[, 1]))
                names(tmp_cols) <- c("dummy.1", "dummy.2", "dummy.3", "dummy.4", "dummy.5")
                tmp_peak_data <- cbind(tmp_peak_data, tmp_cols)
            }

            # Computing correlations
            suppressWarnings(
                cor_res <- Hmisc::rcorr(t(data.frame(tmp_peak_data,
                                                     row.names = msdial_peak_data[msdial_peak_data$source == source,]$id)),
                                        type="pearson")
            )
            cor_res <- merge(matrix_to_df(cor_res$r, "r"), matrix_to_df(cor_res$P, "P"))  # Var1 Var2 r P
            names(cor_res) <- c("id.1", "id.2", "correlation", "P")

            # Correlation and p.value filters
            cor_res <- cor_res[cor_res$correlation >= pearson_correlation_threshold & cor_res$P <= pearson_p_value,]

            # RT filter
            cor_res <- merge(cor_res, msdial_peak_data[, c("id", "Average.Rt.min.", "Adduct.type")], by.x = "id.1", by.y = "id")
            names(cor_res)[names(cor_res) == "Average.Rt.min."] <- "RT1"
            names(cor_res)[names(cor_res) == "Adduct.type"] <- "Adduct.type.1"
            cor_res <- merge(cor_res, msdial_peak_data[, c("id", "Average.Rt.min.", "Adduct.type")], by.x = "id.2", by.y = "id")
            names(cor_res)[names(cor_res) == "Average.Rt.min."] <- "RT2"
            names(cor_res)[names(cor_res) == "Adduct.type"] <- "Adduct.type.2"
            cor_res$delta_RT <- abs(cor_res$RT1 - cor_res$RT2)
            cor_res <- cor_res[cor_res$delta_RT <= threshold_rt,]
            print_message("Correlation links found in ", source, ": ", nrow(cor_res))
            cor_links <- rbind(cor_links, cor_res)
        }

        # merge with msdial_links
        cor_links$simple.nature <- "Pearson correlation"
        cor_links <- cor_links[, c("id.1", "id.2", "simple.nature", "Adduct.type.1", "Adduct.type.2", "correlation")]
        msdial_links <- rbind(msdial_links, cor_links)
    }


    # CLUSTERS MSDIAL [+ PEARSON]
    # Keeping only links between remaining peaks after filtering
    msdial_links <- msdial_links[msdial_links$id.1 %in% msdial_peak_data$id & msdial_links$id.2 %in% msdial_peak_data$id,]
    if (nrow(msdial_links) == 0) stop_script("No links available, can't construct clusters.")
    g <- igraph::graph_from_data_frame(msdial_links,
                                       vertices = msdial_peak_data[,c("id", "Average.Rt.min.", "Average.Mz",
                                                                      "Adduct.type", "Alignment.ID", "source")],
                                       directed = FALSE)
    msdial_clusters <- igraph::cluster_louvain(g)
    igraph::V(g)$cluster.msdial <- msdial_clusters$membership
    final_data <- merge(msdial_peak_data,
                        data.frame(id = igraph::V(g)$name,
                                   cluster.msdial = igraph::V(g)$cluster.msdial),
                        by = "id",
                        all.x = TRUE)

    # Adding cluster size
    cluster_sizes <- dplyr::count(final_data, .data$cluster.msdial)
    names(cluster_sizes) <- c("cluster.msdial", "cluster.msdial.size")
    final_data <- merge(final_data, cluster_sizes)
    print_message("Clusters detected with MSDial data: ", length(levels(factor(msdial_clusters$membership))))

    export_data(final_data, "clusters_msdial")


    # ADDUCTS / NEUTRAL LOSSES
    msc_data <- get_adducts_nl_links(final_data,
                                     samples$Column_name,
                                     threshold_rt,
                                     threshold_mz,
                                     msdial_detected_adducts,
                                     user_pos_neutral_refs,
                                     user_neg_neutral_refs,
                                     user_pos_adducts_refs,
                                     user_neg_adducts_refs)

    export_data(msc_data$computed_massdiff, "computed_massdiff")
    export_data(msc_data$filtered_massdiff, "filtered_computed_massdiff")
    msc_global <- msc_data$links
    if(nrow(msc_global) > 0) msc_global$source <- "Adducts/Neutral losses"
    print_message("Adducts and neutral losses links: ", nrow(msc_global))
    rm(msc_data)


    # FUSION WITH MSDIAL LINKS
    final_links <- merge(msdial_links, msdial_peak_data[, c("id", "Average.Mz", "Average.Rt.min.")], by.x = "id.1", by.y = "id")
    final_links <- merge(final_links, msdial_peak_data[, c("id", "Average.Mz", "Average.Rt.min.")], by.x = "id.2", by.y = "id")
    names(final_links) <- c("CpdID.2", "CpdID.1", "simple.nature", "Adduct.1", "Adduct.2", "Correlation",
                            "Mass.1", "RT.1", "Mass.2", "RT.2")
    final_links$source <- "MSDial"
    final_links$Adduct.1 <- as.character(final_links$Adduct.1)
    final_links$Adduct.2 <- as.character(final_links$Adduct.2)

    cols <- c("CpdID.1", "CpdID.2", "simple.nature", "Adduct.1", "Adduct.2", "Correlation", "RT.1", "RT.2", "Mass.1", "Mass.2",
              "Mean.1", "Mean.2", "N.1", "N.2", "source")
    if (nrow(msc_global) > 0) final_links <- gtools::smartbind(final_links, msc_global)
    else {
        # adding missing columns
        for (col_name in cols) {
            if (!col_name %in% names(final_links)) {
                tmp <- data.frame(tmp = NA)
                names(tmp) <- c(col_name)
                final_links <- cbind(final_links, tmp)
            }
        }
    }

    final_links <- merge(final_links, final_data[, c("id", "cluster.msdial")], by.x = "CpdID.1", by.y = "id")
    final_links <- merge(final_links, final_data[, c("id", "cluster.msdial")], by.x = "CpdID.2", by.y = "id",
                         suffixes = c(".1", ".2"))
    final_links <- final_links[, c("cluster.msdial.1", "cluster.msdial.2", cols)]

    # Weighting links based on probabilities computed from MSDial detections
    weighted_links <- merge(final_links[final_links$simple.nature == "pol / adduct",],
                            msdial_detected_adducts[, c("adduct", "p")],
                            by.x = "Adduct.1", by.y = "adduct", all.x = TRUE)
    weighted_links <- merge(weighted_links,
                            msdial_detected_adducts[, c("adduct", "p")],
                            by.x = "Adduct.2", by.y = "adduct", all.x = TRUE,
                            suffixes = c(".1", ".2"))
    weighted_links$weight <- weighted_links$p.1 * weighted_links$p.2
    final_links <- gtools::smartbind(weighted_links, final_links[final_links$simple.nature != "pol / adduct",])
    export_data(final_links, "links_ms_pre_selection")


    # GRAPH OF ADDUCTS
    # (1) Creation of a graph with each node representing an adduct_peak (peak identified with a specific adduct)
    filtered_links <- final_links[final_links$simple.nature == "pol / adduct",]
    filtered_links$node.1 <- paste(filtered_links$CpdID.1, filtered_links$Adduct.1, sep = "::")
    filtered_links$node.2 <- paste(filtered_links$CpdID.2, filtered_links$Adduct.2, sep = "::")
    g <- igraph::graph_from_data_frame(filtered_links[, c("node.1", "node.2", "weight")], directed = FALSE)
    export_data(g, "adduct_links_graph")
    rm(filtered_links)

    final_peaks_adducts <- data.frame()
    while (length(igraph::V(g)) != nrow(final_peaks_adducts)) {
        # (2) Selection of yet untreated node(s) with highest strength
        current_nodes <- igraph::V(g)[!igraph::V(g)$name %in% final_peaks_adducts$name]
        print_message(length(current_nodes) - nrow(final_peaks_adducts), " nodes to treat")
        max_nodes <- names(current_nodes[igraph::strength(g, v = current_nodes) == max(igraph::strength(g, v = current_nodes))])

        # (3) Deleting nodes with same peak but a different adduct
        for (node_name in max_nodes) {
            print_message(">> ", length(igraph::V(g)), " nodes, ", length(igraph::E(g)), " edges")
            peaks_ids     <- lapply(strsplit(igraph::V(g)$name, "::"), `[[`, 1)
            peaks_adducts <- lapply(strsplit(igraph::V(g)$name, "::"), `[[`, 2)
            node <- igraph::V(g)[igraph::V(g)$name == node_name]
            if (length(node) > 0) {
                print_message("Treating ", peaks_ids[[node]], " (", peaks_adducts[[node]], ")")
                final_peaks_adducts <- rbind(final_peaks_adducts,
                                             data.frame(id     = peaks_ids[[node]],
                                                        adduct = peaks_adducts[[node]],
                                                        name   = igraph::V(g)[node]$name,
                                                        stringsAsFactors = FALSE))
                # deleting nodes with same peak but different adduct
                g <- g - igraph::V(g)[peaks_ids == peaks_ids[[node]] & peaks_adducts != peaks_adducts[[node]]]
            }
        }
        # Back to step 2 until length(igraph::V(g)) == nrow(final_peaks_adducts)
    }
    export_data(g, "final_links_graph")

    # Correction of adducts in final_data
    final_peaks_adducts$name <- NULL
    final_peaks_adducts$adduct <- as.character(final_peaks_adducts$adduct)
    final_data$Adduct.type <- as.character(final_data$Adduct.type)
    final_data <- merge(final_data, final_peaks_adducts, by = "id", all.x = TRUE)
    final_data[!is.na(final_data$adduct),]$Adduct.type <- final_data[!is.na(final_data$adduct),]$adduct
    final_data$adduct <- NULL

    # Correction in final_links
    final_links <- merge(final_links, final_peaks_adducts, by.x = "CpdID.1", by.y = "id", all.x = TRUE)
    final_links$Adduct.1 <- as.character(final_links$Adduct.1)
    final_links$adduct   <- as.character(final_links$adduct)
    final_links <- final_links[final_links$Adduct.1 == final_links$adduct | final_links$simple.nature != "pol / adduct",]
    final_links$adduct <- NULL
    final_links <- merge(final_links, final_peaks_adducts, by.x = "CpdID.2", by.y = "id", all.x = TRUE)
    final_links$Adduct.2 <- as.character(final_links$Adduct.2)
    final_links$adduct <- as.character(final_links$adduct)
    final_links <- final_links[final_links$Adduct.2 == final_links$adduct | final_links$simple.nature != "pol / adduct",]
    final_links$adduct <- NULL

    # Adding strength in final_data
    final_data <- merge(final_data,
                        data.frame(id = gsub("(.+)::.+", "\\1", igraph::V(g)$name),
                                   strength = igraph::strength(g),
                                   stringsAsFactors = FALSE),
                        by = "id", all.x = TRUE)
    final_data[is.na(final_data$strength),]$strength <- 0  # Fix: some strength are NA

    export_data(final_links, "links_ms_post_selection")
    export_data(final_peaks_adducts, "final_adducts_data")


    # DEALING WITH MERS
    # Deleting heteromers and neutral losses if not linked with an adduct
    tmp.adducts <- c(final_links[final_links$simple.nature == "pol / adduct",]$CpdID.1,
                     final_links[final_links$simple.nature == "pol / adduct",]$CpdID.2)
    tmp.nl <- final_links[final_links$simple.nature == "neutral loss",]$CpdID.2
    tmp.heteromers <- final_links[final_links$simple.nature == "heteromers",]$CpdID.1
    final_data <- final_data[final_data$id %in% tmp.adducts | !(final_data$id %in% tmp.nl | final_data$id %in% tmp.heteromers),]
    final_links <- final_links[final_links$CpdID.1 %in% final_data$id & final_links$CpdID.2 %in% final_data$id,]

    # Renaming dimers
    dimers <- final_links[final_links$simple.nature == "homomers",]$CpdID.1
    if(nrow(final_data[final_data$id %in% dimers & final_data$source == "pos",]) > 0) {
        final_data[final_data$id %in% dimers & final_data$source == "pos",]$Adduct.type <- "[2M+H]+"
    }
    if(nrow(final_data[final_data$id %in% dimers & final_data$source == "neg",]) > 0) {
        final_data[final_data$id %in% dimers & final_data$source == "neg",]$Adduct.type <- "[2M-H]-"
    }
    rm(dimers)


    # CLUSTERS MSDIAL + ADDUCTS/NL
    g <- igraph::graph_from_data_frame(final_links[, c("CpdID.1", "CpdID.2", "simple.nature", "Correlation", "source")],
                                       vertices = final_data[c("id", "Average.Rt.min.", "Average.Mz",
                                                               "Adduct.type", "Alignment.ID", "source")],
                                       directed = FALSE)
    rel_clusters <- igraph::cluster_louvain(g)
    igraph::V(g)$cluster <- rel_clusters$membership
    final_data <- merge(final_data,
                        data.frame(id = igraph::V(g)$name,
                                   cluster = igraph::V(g)$cluster),
                        by = "id",
                        all.x = TRUE)
    export_data(g, "clusters_graph")

    # Adding cluster size
    cluster_sizes <- dplyr::count(final_data, .data$cluster)
    names(cluster_sizes) <- c("cluster", "cluster.size")
    final_data <- merge(final_data, cluster_sizes)
    final_links <- merge(final_links, final_data[, c("id", "cluster")], by.x = "CpdID.1", by.y = "id")
    final_links <- merge(final_links, final_data[, c("id", "cluster")], by.x = "CpdID.2", by.y = "id", suffixes = c(".1", ".2"))
    print_message("Clusters detected with MSDial and Adducts/Neutral losses data: ", length(levels(factor(rel_clusters$membership))))


    # CLEANING
    final_links <- final_links[, c("source", "cluster.1", "cluster.2", "cluster.msdial.1", "cluster.msdial.2",
                                   "CpdID.1", "CpdID.2", "simple.nature", "Adduct.1", "Adduct.2", "RT.1", "RT.2",
                                   "Mass.1", "Mass.2", "Mean.1", "Mean.2", "N.1", "N.2", "Correlation")]
    export_data(final_links, "links_ms_final")

    ord.columns <- c("cluster", "cluster.size", "cluster.msdial", "cluster.msdial.size", "strength", "id", "source",
                     "Alignment.ID", "Average.Rt.min.", "Average.Mz", "Adduct.type", "Post.curation.result",
                     "Fill..", "S.N.average", "Spectrum.reference.file.name", "MS1.isotopic.spectrum", "MS.MS.spectrum")
    final_data <- final_data[, c(ord.columns, names(final_data)[!names(final_data) %in% ord.columns])]
    export_data(final_data, "clusters_final")

    print_message("Cleaning done.")
}
