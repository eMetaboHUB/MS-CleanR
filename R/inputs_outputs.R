
#' Get path for a project file or directory.
#'
#' @eval recurrent_params("filetype", "source")
#' @param msfinder_info A string indicating which type of MSFinder file to get path from: "Structure" or "Formula".
#' @param msfinder_lvl A string indicating which level of MSFinder file to get path from.
#' @param peak_id A string indicating the wanted peak id.
#' @return The file or folder full path.
get_project_file_path <- function(filetype, source = NA, msfinder_info = NA, msfinder_lvl = NA, peak_id = NA) {
    # base
    if (is.null(filetype)) return(get("analysis_directory", envir = mscleanrCache))
    # useful to create a centralized interface between code and file system

    # pos / neg
    wd <- file.path(get("analysis_directory", envir = mscleanrCache), source)
    if (filetype == "source_folder")         return(wd)
    if (filetype == "normalized_ms_data")    return(get_pattern_file(wd, "^Normalized"))
    if (filetype == "height_ms_data")        return(get_pattern_file(wd, "^Height"))
    if (filetype == "msdial_peaks")          return(file.path(wd, "peaks"))
    if (filetype == "msdial_filtered_peaks") return(file.path(wd, "filtered_peaks"))
    if (filetype == "msfinder_data")         return(get_pattern_file(file.path(wd, paste0("msf_", msfinder_lvl)),
                                                                     paste0("^", msfinder_info)))
    if (endsWith(filetype, "_peak")) {
        peak_pattern <- paste0("^AlignmentID ", peak_id, "_")
        original_dir <- get_project_file_path("msdial_peaks",source = source)
        if (filetype == "msdial_peak") return(get_pattern_file(original_dir, peak_pattern))
        # else, msdial_filtered_peak
        filtered_dir <- get_project_file_path("msdial_filtered_peaks", source = source)
        tmp <- get_pattern_file(filtered_dir, peak_pattern)
        if (!is.na(tmp)) return(tmp)  # found in filtered peaks folder
        # else, search for file name in original peaks folder
        return(file.path(filtered_dir, get_pattern_file(original_dir, peak_pattern, only_name = TRUE)))
    }

    # intermediary_data
    wd <- file.path(get("analysis_directory", envir = mscleanrCache), "intermediary_data")
    if (filetype == "inter_folder")                return(wd)
    if (filetype == "samples")                     return(file.path(wd, "samples.csv"))
    if (filetype == "msdial_detected_adducts")     return(file.path(wd, "adducts-detected_by_MSDial.csv"))
    if (filetype == "adduct_links_graph")          return(file.path(wd, "adducts-initial.graphml"))
    if (filetype == "final_links_graph")           return(file.path(wd, "adducts-filtered.graphml"))
    if (filetype == "computed_massdiff")           return(file.path(wd, "adducts_massdiff-total.csv"))
    if (filetype == "filtered_computed_massdiff")  return(file.path(wd, "adducts_massdiff-filtered.csv"))
    if (filetype == "final_adducts_data")          return(file.path(wd, "adducts-final_selection.csv"))
    if (filetype == "clusters_msdial")             return(file.path(wd, "MS_peaks-clusters_msdial.csv"))
    if (filetype == "clusters_final")              return(file.path(wd, "MS_peaks-clusters_final.csv"))
    if (filetype == "clusters_graph")              return(file.path(wd, "MS_peaks-clusters.graphml"))
    if (filetype == "final_peaks_data")            return(file.path(wd, "MS_peaks-final_selection.csv"))
    if (filetype == "links_ms_pre_selection")      return(file.path(wd, "links-pre_selection.csv"))
    if (filetype == "links_ms_post_selection")     return(file.path(wd, "links-post_selection.csv"))
    if (filetype == "links_ms_final")              return(file.path(wd, "links-clusters_final.csv"))
    if (filetype == "msdial_identified_peaks")     return(file.path(wd, "annotated_MS_peaks-MSDial.csv"))
    if (filetype == "annotated_data")              return(file.path(wd, "annotated_MS_peaks.csv"))
    if (filetype == "identifying_data")            return(file.path(wd, "annotation_possibilities.csv"))
    if (filetype == "deleted_blk")                 return(file.path(wd, "deleted_blanks.csv"))
    if (filetype == "deleted_ghosts")              return(file.path(wd, "deleted_blank_ghosts.csv"))
    if (filetype == "deleted_mz")                  return(file.path(wd, "deleted_mz.csv"))
    if (filetype == "deleted_rsd")                 return(file.path(wd, "deleted_rsd.csv"))
    if (filetype == "deleted_rmd")                 return(file.path(wd, "deleted_rmd.csv"))

    # CAD/PDA
    if (startsWith(filetype, "CAD_") | startsWith(filetype, "PDA_")) {
        mode <- substr(filetype, 1, 3)
        subfiletype <- substring(filetype, 5)
        wd <- file.path(get("analysis_directory", envir = mscleanrCache), mode)
        if (subfiletype == "manual_data")          return(file.path(wd, "raw_data.csv"))
        if (subfiletype == "cleaned")              return(file.path(wd, "reformatted_peaks.csv"))
        if (subfiletype == "spectra_plot")         return(file.path(wd, "spectra_plot.pdf"))
        if (subfiletype == "roi_plot")             return(file.path(wd, "ROI_plot.pdf"))
        if (subfiletype == "filled")               return(file.path(wd, "grouped_filled_peaks.csv"))
        if (subfiletype == "features")             return(file.path(wd, "featured_peaks.csv"))
    }

    # final data
    wd <- file.path(get("analysis_directory", envir = mscleanrCache), "final_data")
    if (filetype == "final_folder")                return(wd)
    if (filetype == "params")                      return(file.path(wd, "parameters.csv"))
    # if (filetype == "params_clean")                return(file.path(wd, "parameters_clean_msdial_data.csv"))
    # if (filetype == "params_keep")                 return(file.path(wd, "parameters_keep_top_peaks.csv"))
    # if (filetype == "params_launch")               return(file.path(wd, "parameters_launch_msfinder_annotation.csv"))
    # if (filetype == "params_convert")              return(file.path(wd, "parameters_convert_csv_to_msp.csv"))
    if (filetype == "annotated_data-cleaned")      return(file.path(wd, "annotated_MS_peaks-cleaned.csv"))
    if (filetype == "annotated_data-normalized")   return(file.path(wd, "annotated_MS_peaks-normalized.csv"))
    if (filetype == "annotated_data-manual_check") return(file.path(wd, "annotated_MS_peaks-manual_check.csv"))
    if (filetype == "final_CAD")                   return(file.path(wd, "CAD-cleaned.csv"))
    if (filetype == "final_PDA")                   return(file.path(wd, "PDA-cleaned.csv"))
    if (filetype == "msp_pos")                     return(file.path(wd, "peaks_pos.msp"))
    if (filetype == "msp_neg")                     return(file.path(wd, "peaks_neg.msp"))
    if (filetype == "msp_metadata_pos")            return(file.path(wd, "peaks_metadata_pos.csv"))
    if (filetype == "msp_metadata_neg")            return(file.path(wd, "peaks_metadata_neg.csv"))

    stop_script("Unrecognized filetype: ", filetype)
}



#' Get path for a file found thanks to a pattern.
#'
#' @param dir_path A file.path indicating the directory in which to look for the file.
#' @param pattern A string indicating the pattern used to find the file.
#' @param only_name A boolean indicating whether to return the full path or only the name of the file.
#' @return A string with the found file path or NA if not found.
get_pattern_file <- function(dir_path, pattern, only_name = FALSE) {
    tmp <- dir(dir_path, pattern = pattern)
    if (length(tmp) == 0) return(NA)
    if (length(tmp) > 1)  print_warning("Several files matched in directory ", dir_path, ".\nLoading file ", tmp[1])
    if (only_name) return(tmp[1])
    else           return(file.path(dir_path, tmp[1]))
}



#' Export data for users.
#'
#' @eval recurrent_params("filetype", "source")
#' @param data_to_export A data.frame to export to the user file system.
#' @param empty_na A boolean indicating if NA must be replaced by an empty string in the output file.
export_data <- function(data_to_export, filetype, source = NA, empty_na = FALSE) {
    if (endsWith(filetype, "_graph")) {
        igraph::write_graph(data_to_export, get_project_file_path(filetype), format = "graphml")
    } else {
        path <- get_project_file_path(filetype, source = source)
        suppressMessages({
            if (empty_na) vroom::vroom_write(data_to_export, path, delim = ",", na = "", progress = FALSE)
            else          vroom::vroom_write(data_to_export, path, delim = ",", progress = FALSE)
        })
    }
}



#' Import data exported by the script.
#'
#' @eval recurrent_params("filetype")
#' @param ... Additional parameters passed on to \code{\link{get_project_file_path}}.
import_data <- function(filetype, ...) {
    suppressMessages(
        return(as.data.frame(vroom::vroom(get_project_file_path(filetype, ...),
                                          delim = ",",
                                          na = c("NA", "", "N/A"),
                                          progress = FALSE)))
    )
}



# FIXME
# Export parameters used for the current run.
#
# @param ... Parameters used for the current run of the script.
# export_params <- function(...) {
#     clean_list <- function(x) {
#         if (length(x) == 1) return(x)
#         if (length(x) == 0) return("")
#         if (is.null(names(x))) paste0("[", paste0(x, collapse = ","), "]")
#         else                   paste0("[",
#                                       paste0(lapply(names(x), FUN = function(y) { paste0(y, ":", x[[y]]) }), collapse = ","),
#                                       "]")
#     }
#
#     calling_func <- deparse(sys.calls()[[sys.nframe()-1]])
#     id_func <- unlist(strsplit(calling_func[1], "(", fixed = TRUE))[1]
#     params <- as.data.frame(lapply(list(...)[[1]], FUN = clean_list))
#     params[, ] <- lapply(params[, ], as.character) # converts all column to character, avoiding loss of logical values (sometimes converted to numeric otherwise)
#     params <- t(params)
#
#     current_params <- data.frame(FUNCTION  = id_func,
#                                  PARAMETER = c("Package version", rownames(params)),
#                                  VALUE     = c(as.character(utils::packageVersion("mscleanr")), params[,1]),
#                                  stringsAsFactors = FALSE)
#
#     # Better formatting for lists if only 1 element
#     for (p in c("compound_levels", "biosoc_levels", "levels_scores")) {
#         if (p %in% current_params$PARAMETER) {
#             val <- as.character(current_params[current_params$PARAMETER == p,]$VALUE)
#             if (!startsWith(val, "[")) {
#                 current_params[current_params$PARAMETER == p,]$VALUE <- paste0("[", val, "]")
#             }
#         }
#     }
#
#     if (file.exists(get_project_file_path("params"))) {
#         existing_params <- import_data("params")
#         existing_params <- existing_params[existing_params$FUNCTION != id_func,]
#         current_params  <- rbind(existing_params, current_params)
#     }
#     current_params$FUNCTION <- factor(current_params$FUNCTION,
#                                       levels = c("Package version",
#                                                  "clean_msdial_data",
#                                                  "keep_top_peaks",
#                                                  "launch_msfinder_annotation",
#                                                  "convert_csv_to_msp",
#                                                  "launch_cad_fusion"))
#     export_data(current_params[order(current_params$FUNCTION),], "params")
# }
