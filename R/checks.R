
#' Check project architecture before running \code{\link{clean_msdial_data}}.
#'
#' @eval recurrent_params("filter_blk")
check_architecture_for_clean_msdial_data <- function(filter_blk) {
    check_main_architecture()

    clean_for_current_run("final_folder", create_dir = TRUE)
    clean_for_current_run("inter_folder", create_dir = TRUE)

    for (source in get("analysis_modes", envir = mscleanrCache)) {
        check_path("normalized_ms_data",
                   source = source,
                   na_msg = paste0("Can't find normalized MS data for ", source, " mode."))
        if (filter_blk) {
            check_path("height_ms_data",
                       source = source,
                       na_msg = paste0("Can't find height MS data for ", source, " mode (necessary for blank filtering)."))}
    }
}



#' Check project architecture before running \code{\link{keep_top_peaks}}.
#'
#' @param export_filtered_peaks Copy peaks files corresponding to the peaks remaining after filtering in a new folder.
check_architecture_for_keep_top_peaks <- function(export_filtered_peaks) {
    check_main_architecture()
    check_path("clusters_final")
    check_path("samples")
    check_path("final_adducts_data")
    if (export_filtered_peaks) {
        for (source in get("analysis_modes", envir = mscleanrCache)) {
            check_path("msdial_peaks", source = source)
            clean_for_current_run("msdial_filtered_peaks", source = source, create_dir = TRUE)
        }
    }
}



#' Check project architecture before running \code{\link{launch_msfinder_annotation}}.
#'
#' @param msfinder_biosoc_levels MSFinder annotation biosource levels to check (only used when \code{launch_msfinder_annotation = TRUE}).
check_architecture_for_launch_msfinder_annotation <- function(msfinder_biosoc_levels = NULL) {
    check_main_architecture()

    check_path("final_peaks_data")
    check_path("links_ms_final")
    check_path("samples")

    clean_for_current_run("final_folder", delete = FALSE, create_dir = TRUE)
    clean_for_current_run("annotated_data-manual_check")
    clean_for_current_run("annotated_data-cleaned")
    clean_for_current_run("annotated_data-normalized")
    for (source in c("pos", "neg")) clean_for_current_run(paste0("msp_", source))

    if (length(msfinder_biosoc_levels) == 0) print_warning("No levels provided for MSFinder annotation, only 'generic' will be used.")

    for (source in get("analysis_modes", envir = mscleanrCache)) {
        for (info in c("Formula", "Structure")) {
            for (lvl in c(msfinder_biosoc_levels, "generic")) {
                check_path("msfinder_data", source = source, msfinder_info = info, msfinder_lvl = lvl,
                           stop = FALSE,
                           na_msg = paste0("Can't find MSFinder ", info, " file for level ", lvl, " and ", source, " mode."))
            }
        }
    }
}



#' Check the project main architecture
#'
#' @param posneg_check Whether to check the existence of pos/neg directories.
check_main_architecture <- function(posneg_check = TRUE) {
    if (!exists("analysis_directory", envir = mscleanrCache)) {
        stop_script("No project directory chosen. Please do so with choose_project_directory().")
    }
    check_path(filetype = NULL, na_msg = "Project directory can't be NA.")
    if (posneg_check) {
        modes <- c()
        if (check_path("source_folder", source = "pos", stop = FALSE)) modes <- c(modes, "pos")
        if (check_path("source_folder", source = "neg", stop = FALSE)) modes <- c(modes, "neg")
        if (length(modes) == 0) stop_script("A pos or a neg directory must be present.")
        assign("analysis_modes", modes, envir = mscleanrCache)
    }
}



#' Check before running \code{\link{launch_msfinder_annotation}}.
#'
#' @param min_score Input parameter to check.
check_for_convert_csv_to_msp <- function(min_score) {
    check_main_architecture(posneg_check = FALSE)
    check_path("annotated_data-normalized")
    check_path("samples")
    check_positive_num(min_score)
    for (source in c("pos", "neg")) clean_for_current_run(paste0("msp_", source))
}



#' Check if input parameters are correct before running \code{clean_msdial_data}.
#'
#' @param filter_blk,filter_blk_threshold,filter_mz,filter_rsd,filter_rsd_threshold,filter_blk_ghost_peaks,filter_rmd,filter_rmd_range,threshold_mz,threshold_rt,compute_pearson_correlation,pearson_threshold,pearson_p_value,references_adduct_pos,references_adduct_neg,references_neutral_pos,references_neutral_neg Input parameters to check.
check_input_parameters_msdial_data <- function(filter_blk,
                                               filter_blk_threshold,
                                               filter_mz,
                                               filter_rsd,
                                               filter_rsd_threshold,
                                               filter_blk_ghost_peaks,
                                               filter_rmd,
                                               filter_rmd_range,
                                               threshold_mz,
                                               threshold_rt,
                                               compute_pearson_correlation,
                                               pearson_threshold,
                                               pearson_p_value,
                                               references_adduct_pos,
                                               references_adduct_neg,
                                               references_neutral_pos,
                                               references_neutral_neg) {
    check_boolean(filter_blk,                  "filter_blk")
    check_boolean(filter_mz,                   "filter_mz")
    check_boolean(filter_rsd,                  "filter_rsd")
    check_boolean(filter_rmd,                  "filter_rmd")
    check_boolean(filter_blk_ghost_peaks,      "filter_blk_ghost_peaks")
    check_boolean(compute_pearson_correlation, "compute_pearson_correlation")

    check_positive_num(threshold_mz,     "threshold_mz")
    check_positive_num(threshold_rt,     "threshold_rt")

    check_probability(filter_blk_threshold, "filter_blk_threshold")
    check_probability(pearson_threshold,    "pearson_threshold")
    check_probability(pearson_p_value,      "pearson_p_value")

    check_positive_int(filter_rsd_threshold, "filter_rsd_threshold")
    check_positive_int(filter_rmd_range[1],  "filter_rmd_range first value")
    check_positive_int(filter_rmd_range[2],  "filter_rmd_range last value")

    check_references(references_adduct_pos,  "positive adducts")
    check_references(references_adduct_neg,  "negative adducts")
    check_references(references_neutral_pos, "positive neutral losses")
    check_references(references_neutral_neg, "negative neutral losses")

    if (filter_rmd_range[1] >= filter_rmd_range[2]) stop_script("filter_rmd_range first value must be lower than filter_rmd_range last value.")
}



#' Check if input parameters are correct before running \code{keep_top_peaks}.
#'
#' @param selection_criterion,n,export_filtered_peaks Input parameters to check.
check_input_parameters_keep_top_peaks <- function(selection_criterion, n, export_filtered_peaks) {
    check_boolean(export_filtered_peaks, "export_filtered_peaks")
    check_positive_int(n, "n")
    if (!selection_criterion %in% c("intensity", "degree", "both")) {
        stop_script('selection_criterion must be "intensity", "degree" or "both"')
    }
}



#' Check if input parameters are correct before running \code{launch_msfinder_annotation}.
#'
#' @param msfinder_biosoc_levels,msfinder_compound_levels,levels_scores Input parameters to check.
check_input_parameters_launch_msfinder <- function(msfinder_biosoc_levels, msfinder_compound_levels, levels_scores) {
    suppressWarnings(  # is.na(list(...)) returns a warning: la condition a une longueur > 1 et seul le premier élément est utilisé
        if (length(msfinder_biosoc_levels) == 0 | is.na(msfinder_biosoc_levels)) {
            stop_script("No MSFinder levels provided.")
        }
    )
    for (lvl in msfinder_biosoc_levels) if (is.null(lvl) | is.na(lvl) | lvl == "") stop_script("An MSFinder level can't be empty.")
    for (lvl in names(levels_scores)) {
        if (!(lvl %in% msfinder_biosoc_levels | lvl %in% msfinder_compound_levels)) {
            print_warning("Level ", lvl, " present in scores but not in biosource or compound levels.")
        }
        check_positive_num(levels_scores[[lvl]], "levels_scores")
    }
}



#' Create a directory or delete a file for a new run.
#'
#' Rename the old directory or file if existing.
#'
#' @eval recurrent_params("filetype")
#' @param delete A boolean indicating whether the file or directory indicated by the path needs to be deleted.
#' @param create_dir A boolean indicating whether the directory indicated by the path needs to be created.
#' @param ... Additional parameters to pass to \code{\link{get_project_file_path}}.
clean_for_current_run <- function(filetype, delete = TRUE, create_dir = FALSE, ...) {
    path <- get_project_file_path(filetype, ...)
    if (file.exists(path) & delete) {
        if (!get("overwrite", envir = mscleanrCache) & !get("shiny_running", envir = mscleanrCache)) {
            assign("overwrite", get_confirmation_for_overwriting(path), envir = mscleanrCache)
        }
        if (get("overwrite", envir = mscleanrCache)) {
            print_warning("Deleting ", path)
            unlink(path, recursive = TRUE)
        }
        else stop_script(path, " already exists.", generic_msg = FALSE)
    }
    if (create_dir) dir.create(path, showWarnings = FALSE)
}



#' Check if a path exist.
#'
#' @eval recurrent_params("filetype")
#' @param stop If \code{TRUE}, stops the script if a path is not found.
#' @param na_msg A string with info to print if path doesn't exist.
#' @param ... Additional parameters passed on to \code{\link{get_project_file_path}}.
check_path <- function(filetype, stop = TRUE, na_msg = NA, ...) {
    path <- get_project_file_path(filetype, ...)
    if (stop) {
        if (is.na(path))        stop_script(na_msg)
        if (!file.exists(path)) stop_script("Cannot find ", path)
    } else {
             if (is.na(path))        print_warning(na_msg)
        else if (!file.exists(path)) {
            print_warning("Cannot find ", path)  # file.exists(NA) causes an exception
            return(FALSE)
        }
    }
    return(TRUE)
}



#' Check if value is boolean.
#'
#' @param x Value to check.
#' @param n Name of value.
check_boolean <- function(x, n) if (!is.logical(x)) stop_script(n, " must be a boolean value.")



#' Check if value is a positive integer.
#'
#' @param x Value to check.
#' @param n Name of value.
check_positive_int <- function(x, n) if (x %% 1 != 0 | x <= 0) stop_script(n, " must be a positive integer value.")



#' Check if value is a positive numerical value.
#'
#' @param x Value to check.
#' @param n Name of value.
#' @param no_NA if TRUE, \code{x} can't be NA.
check_positive_num <- function(x, n, no_NA = TRUE) {
    if ((no_NA & is.na(x)) | !(is.na(x) | (is.numeric(x) & x >= 0))) stop_script(n, " must be a numerical value >= 0.")
}



#' Check if value is a probability (numerical value between 0 and 1).
#'
#' @param x Value to check.
#' @param n Name of value.
check_probability <- function(x, n) {
    if (!is.na(x) & (x < 0 | x > 1)) stop_script(n, " must be a numerical value between 0 and 1.")
}



#' Check if value is a correct reference data.frame
#'
#' @param x Value to check.
#' @param n Name of value.
check_references <- function(x, n) {
    if (!is.na(x)) {
        if (!is.data.frame(x) | ncol(x) != 2) {
            stop_script("References for ", n, " masses needs to be a 2-columns data.frame.")
        }
        if (nrow(x) == 0) print_warning("The data.frame of references for ", n, " masses is empty.")
    }
}
