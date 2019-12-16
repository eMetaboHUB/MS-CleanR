
#' Annotates peaks based on files extracted from MSFinder.
#'
#' @eval recurrent_params("compound_levels", "biosoc_levels")
#' @param levels_scores A list of levels names and their corresponding multiplier to adapt final annotation scores.
#'
#' @section Architecture needed by launch_msfinder_annotation in the project_directory:
#' \itemize{
#'   \item pos
#'   \itemize{
#'     \item msf_X
#'     \describe{
#'       \item{Formula<...>.txt}{Formulas for possible identifications in the X biosoc level, exported from MSFinder}
#'       \item{Structure<...>.txt}{Structures for possible identifications in the X biosoc level, exported from MSFinder}
#'     }
#'   }
#'   \item neg
#'   \itemize{
#'     \item msf_X
#'     \describe{
#'       \item{Formula<...>.txt}{Formulas for possible identifications in the X biosoc level, exported from MSFinder}
#'       \item{Structure<...>.txt}{Structures for possible identifications in the X biosoc level, exported from MSFinder}
#'     }
#'   }
#'   \item intermediary_data
#'   \describe{
#'     \item{samples.csv}{Information on samples}
#'     \item{final_peaks_data.csv}{Peaks remaining after all filters have been applied}
#'     \item{links_ms-final.csv}{Links between peaks indicated by MSDial (global) and adducts/neutral losses detection (by cluster)}
#'   }
#' }
#'
#' @section Output files of launch_msfinder_annotation in the project directory:
#' \describe{
#'   \item{annotated_data-cleaned.csv}{Final auto-annotated data with only relevant peaks}
#' }
#' \itemize{
#'   \item intermediary_data
#'   \describe{
#'     \item{identifying_data.csv}{Union of MSDial and MSFinder data, with clusters info}
#'     \item{annotated_data.csv}{Final auto-annotated data, complete with all possibilities}
#'     \item{annotated_data-manual_check.csv}{Final auto-annotated data for clusters needing manual checking (only created in this case)}
#'   }
#' }
#'
#' @export
launch_msfinder_annotation <- function(compound_levels = NULL,  # c() = NULL
                                       biosoc_levels = c("generic"),
                                       levels_scores = NULL) {

    check_input_parameters_launch_msfinder(biosoc_levels, compound_levels, levels_scores)
    check_architecture_for_launch_msfinder_annotation(biosoc_levels)

    print_message("*** Treating ", get("project_directory", envir = mscleanrCache), " ***")
    samples     <- import_data("samples")
    final_data  <- import_data("final_peaks_data")
    final_links <- import_data("links_ms_final")

    # Samples check
    for(s in samples$Column_name) {
        if(!(s %in% names(final_data))) stop_script("Sample '", s, "' not found in peaks data file.")
    }

    # Merging all MSFinder files
    msfinder_data <- data.frame()
    for(source in c("pos", "neg")) {
        for(level in c(biosoc_levels, "generic")) {
            suppressWarnings(
                msf_tmp <- import_msfinder_data(source, level)
            )
            # Dealing with potentially missing columns
            if (!is.null(msf_tmp) & length(names(msfinder_data)) > 0) {
                main_names <- names(msfinder_data)
                tmp_names  <- names(msf_tmp)
                for (name in main_names) if (!name %in% tmp_names)  msf_tmp[name]       <- NA
                for (name in tmp_names)  if (!name %in% main_names) msfinder_data[name] <- NA
            }
            if (!is.null(msf_tmp)) {
                suppressWarnings(  # Column class mismatch for 'Classyfire_subclass'. Converting column to class 'character'.
                    msfinder_data <- gtools::smartbind(msfinder_data, msf_tmp)
                )
            }
        }
    }

    if(nrow(msfinder_data) == 0) stop_script("No usable MSFinder data found in ",
                                             get("project_directory", envir = mscleanrCache),
                                             ".")

    # Deleting duplicates present in several levels
    msfinder_data <- msfinder_data %>% dplyr::distinct(.data$id,
                                                       .data$Formula,
                                                       .data$Structure,
                                                       .data$SMILES,
                                                       .data$InChIKey,
                                                       .keep_all = TRUE)


    # Merging with MSDial data
    msd_main_cols <- c("cluster", "cluster.size", "source", "id", "Alignment.ID", "Average.Rt.min.",
                       "Average.Mz", "Adduct.type", "MS.MS.spectrum")
    msf_main_cols <- c("level", "rank.formula", "Formula", "rank.structure", "Structure", "Total.score")
    msf_other_cols <- names(msfinder_data)[!names(msfinder_data) %in% c(msf_main_cols,
                                                                        "source", "id", "Alignment.ID", "Databases.formula")]

    identifying_data <- merge(final_data, msfinder_data, by = c("id", "source", "Alignment.ID"), all.x = TRUE)
    identifying_data <- identifying_data[, c(msd_main_cols, msf_main_cols, msf_other_cols, samples$Column_name)]
    # sorting
    biosoc_levels <- biosoc_levels[biosoc_levels != "generic"]
    identifying_data$level <- ordered(identifying_data$level, levels = c(biosoc_levels, "generic"))
    identifying_data <- identifying_data[with(identifying_data, order(cluster, Alignment.ID, level, -xtfrm(Total.score))),]
    export_data(identifying_data, "identifying_data")


    # AUTOMATIC ANNOTATION
    identifying_data$Total.score        <- as.numeric(identifying_data$Total.score)
    identifying_data$rowid              <- 1:nrow(identifying_data)
    identifying_data$annotation         <- NA
    identifying_data$annotation_result  <- NA
    identifying_data$annotation_warning <- NA

    # Clusters of size 1 or 2 without identification and with links: Ignored
    # print_message("*** Ignoring small clusters without identification and with links ***")
    # for(ignored_size in c(1, 2)) {
    #     small_clusters <- dplyr::count(identifying_data[identifying_data$cluster.size == ignored_size
    #                                                     & is.na(identifying_data$level),],
    #                                    .data$cluster)
    #     for(cluster_id in small_clusters[small_clusters$n == ignored_size,]$cluster) {
    #         cluster_links <- final_links[   final_links$simple.nature != "found in higher mz's MsMs"
    #                                      &  final_links$simple.nature != "similar chromatogram in higher mz"
    #                                      & (final_links$cluster.1 == cluster_id | final_links$cluster.2 == cluster_id)
    #                                      &  final_links$cluster.2 != final_links$cluster.1,]
    #         if(nrow(cluster_links) > 0) {
    #             identifying_data[identifying_data$cluster == cluster_id,]$annotation_result <- "Ignored"
    #         }
    #     }
    # }

    # Clusters with a pair [M+H]+ / [M-H]-
    print_message("*** Annotating clusters with [M+H]+ / [M-H]- couples ***")
    m_pairs <- final_links[final_links$Adduct.1 %in% c("[M+H]+", "[M-H]-") & final_links$Adduct.2 %in% c("[M+H]+", "[M-H]-"),]
    m_pairs <- m_pairs[m_pairs$CpdID.1 %in% final_data$id & m_pairs$CpdID.2 %in% final_data$id,]
    m_pairs <- m_pairs[m_pairs$cluster.1 == m_pairs$cluster.2 & !is.na(m_pairs$cluster.1),]
    freq <- dplyr::count(m_pairs, .data$cluster.1)
    for(cluster_id in freq[freq$n == 1,]$cluster.1) {
        identifying_data <- annotate_cluster(identifying_data,
                                             cluster_id,
                                             couple_ids = c(m_pairs[m_pairs$cluster.1 == cluster_id,]$CpdID.1,
                                                            m_pairs[m_pairs$cluster.1 == cluster_id,]$CpdID.2),
                                             compound_levels = compound_levels,
                                             biosoc_levels = biosoc_levels)
    }

    # Clusters with several pairs [M+H]+ / [M-H]-, we only consider the pair having the highest mass
    for(cluster_id in freq[freq$n > 1,]$cluster.1) {
        ids_1 <- unlist(m_pairs[m_pairs$cluster.1 == cluster_id,]$CpdID.1)
        ids_2 <- unlist(m_pairs[m_pairs$cluster.1 == cluster_id,]$CpdID.2)
        max_mass <- max(final_data[final_data$id %in% c(ids_1, ids_2),]$Average.Mz)
        couple_max <- m_pairs[m_pairs$cluster.1 == cluster_id & (m_pairs$Mass.1 == max_mass | m_pairs$Mass.2 == max_mass),]
        if(length(couple_max$CpdID.1) == 1) {
            identifying_data <- annotate_cluster(identifying_data,
                                                 cluster_id,
                                                 couple_ids = c(couple_max$CpdID.1, couple_max$CpdID.2),
                                                 compound_levels = compound_levels,
                                                 biosoc_levels = biosoc_levels)
        } else {
            # If several pairs with highest mass, we consider all peaks in these pairs
            identifying_data <- annotate_cluster(identifying_data,
                                                 cluster_id,
                                                 couple_ids = c(ids_1, ids_2),
                                                 compound_levels = compound_levels,
                                                 biosoc_levels = biosoc_levels)
        }
    }

    # Reste
    print_message("*** Annotating remaining clusters ***")
    for(cluster_id in levels(factor(identifying_data[is.na(identifying_data$annotation_result),]$cluster))) {
        identifying_data <- annotate_cluster(identifying_data,
                                             cluster_id,
                                             compound_levels = compound_levels,
                                             biosoc_levels = biosoc_levels)
    }


    # Adapting scores
    identifying_data$Final.score <- as.numeric(identifying_data$Total.score)
    for (level in names(levels_scores)) identifying_data <- update_annotations_scores(identifying_data, level, levels_scores[[level]])

    msf_main_cols <- c(msf_main_cols, "Final.score")  # adding Final.score
    identifying_data <- identifying_data[, c(c("annotation_result", "annotation", "annotation_warning"),
                                             msd_main_cols[msd_main_cols != "MS.MS.spectrum"],
                                             msf_main_cols,
                                             msf_other_cols,
                                             samples$Column_name,
                                             "MS.MS.spectrum"),]
    identifying_data <- identifying_data[with(identifying_data, order(cluster, Alignment.ID, -xtfrm(Final.score))),]
    export_data(identifying_data, "annotated_data", empty_na = TRUE)


    # Exporting final annotations
    if(nrow(identifying_data[!is.na(identifying_data$annotation_warning),]) > 0) {
        export_data(identifying_data[!is.na(identifying_data$annotation_warning),],
                    "annotated_data-manual_check",
                    empty_na = TRUE)
    }

    final_annotations <- identifying_data[which(identifying_data$annotation),]
    export_data(final_annotations[!names(final_annotations) %in% c("annotation", "id", "rank.formula",
                                                                   "rank.structure", "cluster", "cluster.size")],
                "annotated_data-cleaned",
                empty_na = TRUE)


    # Normalization of remaining peaks
    final_annotations[samples$Column_name] <- scale(final_annotations[samples$Column_name],
                                                    center=FALSE,
                                                    scale=colSums(final_annotations[samples$Column_name]))
    final_annotations[samples$Column_name] <- round(final_annotations[samples$Column_name] * 1000, 2)

    export_data(final_annotations[!names(final_annotations) %in% c("annotation", "id", "rank.formula",
                                                                   "rank.structure", "cluster", "cluster.size")],
                "annotated_data-normalized",
                empty_na = TRUE)

}
