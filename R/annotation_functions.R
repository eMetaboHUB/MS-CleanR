#' Annotate a single cluster
#'
#' @eval recurrent_params("compound_levels", "biosoc_levels")
#' @param identifying_data A data.frame containing possible identifications.
#' @param cluster_id A string indicating the cluster id.
#' @param couple_ids A 2-tuple indicating specific compounds to consider for annotation
#' @return A data.frame with the selected identification annotated.
annotate_cluster <- function(identifying_data, cluster_id, couple_ids = NULL, compound_levels = NULL, biosoc_levels = NULL) {
    print_message("Annotating cluster ", cluster_id)
    in_cluster <- identifying_data$cluster == cluster_id
    possibilities <- identifying_data[in_cluster
                                      & !is.na(identifying_data$level)
                                      & !is.nan(identifying_data$Total.score),]
    if (is.null(couple_ids)) {
        mode_1 <- possibilities[possibilities$source == "pos",]
        mode_2 <- possibilities[possibilities$source == "neg",]
    } else {
        mode_1 <- possibilities[possibilities$id %in% couple_ids[1],]
        mode_2 <- possibilities[possibilities$id %in% couple_ids[2],]
    }

    # No ID
    if (nrow(mode_1) == 0 & nrow(mode_2) == 0) {
        identifying_data[in_cluster,]$annotation_result <- "Unknown compound"
        if (is.null(couple_ids)) {
            max.mz <- max(identifying_data[in_cluster,]$Average.Mz, na.rm=TRUE)
            identifying_data[in_cluster & identifying_data$Average.Mz == max.mz,]$annotation <- TRUE
        } else {
            max.mz <- max(identifying_data[identifying_data$id %in% couple_ids,]$Average.Mz, na.rm=TRUE)
            identifying_data[identifying_data$id %in% couple_ids & identifying_data$Average.Mz == max.mz,]$annotation <- TRUE
        }
    } else {
        if ("Compound_level" %in% names(mode_1)) {
            col_merge <- c("level", "Formula", "Structure", "SMILES", "Compound_level")
        } else {
            col_merge <- c("level", "Formula", "Structure", "SMILES")
        }
        common_ids <- merge(mode_1, mode_2, by = col_merge)
        if (nrow(common_ids) > 0) {
            # ID identique dans 2 modes
            id_type <- "Double ID"
            most_probable <- find_structure(common_ids, compound_levels = compound_levels, biosoc_levels = biosoc_levels)
        } else {
            # ID dans un seul mode
            id_type <- "Simple ID"
            most_probable <- find_structure(rbind(mode_1, mode_2), compound_levels = compound_levels, biosoc_levels = biosoc_levels)
        }
        identifying_data[in_cluster,]$annotation_result <- id_type

        if (length(most_probable) > 1) {
            identifying_data[in_cluster,]$annotation_warning <- paste0("Manual check recommanded : several final ",
                                                                       id_type,
                                                                       " possibilities, the first one was selected.")
            most_probable <- c(most_probable[1])
        }
        identifying_data[identifying_data$rowid %in% most_probable,]$annotation <- TRUE
    }

    return(identifying_data)
}



#' Find the most probable structure by iterating over all possible levels of annotations.
#'
#' @eval recurrent_params("possibilities", "compound_levels", "biosoc_levels")
#' @return A vector containing the row id(s) of the most probable structure(s).
find_structure <- function(possibilities, compound_levels = NULL, biosoc_levels = NULL) {
    # (1) Structures with a compound level
    for(cat in compound_levels) {
        for(lvl in biosoc_levels) {
            tmp <- possibilities[possibilities$level == lvl & possibilities$Compound_level == cat,]
            if(nrow(tmp) > 0) return(get_most_probable(tmp))
        }
    }

    # (2) Structures without a compound level
    for(lvl in biosoc_levels) {
        tmp <- possibilities[possibilities$level == lvl & !(possibilities$Compound_level %in% compound_levels),]
        if(nrow(tmp) > 0) return(get_most_probable(tmp))
    }

    # (3) If nothing found in personal databases, searching in last generic level
    return(get_most_probable(possibilities[possibilities$level == "generic",]))
}



#' Get the most probable structure for the whole input data.frame.
#'
#' @eval recurrent_params("possibilities")
#' @return A vector containing the row id(s) of the most probable structure(s).
get_most_probable <- function(possibilities) {
    # ID simple ou double
    simple_id <- length(grep("Total.score", names(possibilities))) == 1

    # Max score then min rank
    if(simple_id) {
        possibilities$sort_score <- possibilities$Total.score
        possibilities$sort_rank  <- as.numeric(possibilities$rank.structure)
    } else {
        possibilities$sort_score <- possibilities$Total.score.x + possibilities$Total.score.y
        possibilities$sort_rank  <- (as.numeric(possibilities$rank.structure.x) + as.numeric(possibilities$rank.structure.y)) / 2
    }
    possibilities <- possibilities[possibilities$sort_score == max(possibilities$sort_score, na.rm=TRUE),]
    possibilities <- possibilities[possibilities$sort_rank  == min(possibilities$sort_rank, na.rm=TRUE),]

    if (simple_id) final <- c(possibilities$rowid)
    else {
        # get rowid max score for each couple
        possibilities$max.rowid <- ifelse(possibilities$Total.score.x > possibilities$Total.score.y,
                                          possibilities$rowid.x,
                                          possibilities$rowid.y)
        final <- levels(factor(possibilities$max.rowid))
    }
    return(final)
}



#' Adapt the final annotations scores.
#'
#' @param data A data.frame containing annotations.
#' @param focus_level A string indicating the biosource or compound level to consider.
#' @param coef_mult A numeric value indicating the coefficient by which to multiply the original MSFinder score.
#' @return A data.frame containing annotations with multiplied scores.
update_annotations_scores <- function(data, focus_level, coef_mult) {
    filter_biosoc_lvl <- !is.na(data$level) & data$level == focus_level
    filter_cpd_lvl <- !is.na(data$Compound_level) & data$Compound_level == focus_level
    filter <- !is.na(data$Final.score) & (filter_biosoc_lvl | filter_cpd_lvl)
    data[filter,]$Final.score <- data[filter,]$Final.score * coef_mult
    return(data)
}
