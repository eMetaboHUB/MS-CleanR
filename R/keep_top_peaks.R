
#' Filter peaks based on a given criterion.
#'
#' Keeps \code{n} peaks by cluster.
#' If there are ties, more than \code{n} peaks can be selected ( see \code{\link[dplyr]{top_n}} for more information).
#'
#' @eval recurrent_params("project_directory", "overwrite")
#' @param selection_criterion "intensity" or "degree" or "both".
#' @param n Keep only n.top peaks by clusters if intensity or degree filters.
#' @param export_filtered_peaks Copy peaks files corresponding to the peaks remaining after filtering in a new folder.
#'
#' @export
keep_top_peaks <- function(project_directory,
                           selection_criterion,
                           n = 1,
                           export_filtered_peaks = TRUE,
                           overwrite = FALSE) {

    check_input_parameters_keep_top_peaks(selection_criterion, n, export_filtered_peaks, overwrite)
    check_architecture_for_keep_top_peaks(project_directory, export_filtered_peaks, overwrite)

    samples             <- import_data(project_directory, "samples")
    final_data          <- import_data(project_directory, "clusters_final")
    final_peaks_adducts <- import_data(project_directory, "final_adducts_data")

    # Intensity filter
    final_data_cluster <- final_data %>% dplyr::group_by(.data$cluster)
    if (selection_criterion == "intensity" | selection_criterion == "both") {
        if("QC" %in% samples$Script_class) {
            final_data_intensity <- final_data_cluster %>% dplyr::top_n(n, .data$avg_QC)
        }
        else if("Standard" %in% samples$Script_class) {
            final_data_intensity <- final_data_cluster %>% dplyr::top_n(n, .data$avg_Standard)
        }
        else {
            final_data_intensity <- final_data_cluster %>% dplyr::top_n(n, .data$avg_Samples)
        }
        final_data_intensity <- data.frame(final_data_intensity)
    }

    # Strength filter
    if (selection_criterion == "degree" | selection_criterion == "both") {
        final_data_degree <- data.frame(final_data_cluster %>% dplyr::top_n(n, .data$strength))
    }

         if (selection_criterion == "intensity") final_data <- final_data_intensity
    else if (selection_criterion == "degree")    final_data <- final_data_degree
    else if (selection_criterion == "both")      final_data <- unique(rbind(final_data_intensity, final_data_degree))
    print_peaks_status(final_data, "MSDial peaks after peaks filtering:")

    final_data$cluster.msdial <- NULL
    final_data$cluster.msdial.size <- NULL
    export_data(final_data, project_directory, "final_peaks_data")


    # Exporting peaks files
    if (export_filtered_peaks) {
        for (mode in c("pos", "neg")) {
            for (peak_id in final_data[final_data$source == mode, "Alignment.ID"]) {
                peak_path <- get_project_file_path(project_directory, "msdial_peak", source = mode, peak_id = peak_id)
                if (is.na(peak_path)) print_warning("Can't find file corresponding to peak ", peak_id)
                else {
                    new_peak_path <- get_project_file_path(project_directory, "msdial_filtered_peak",
                                                           source = mode, peak_id = peak_id)
                    file.copy(from = peak_path, to = new_peak_path)

                    if (paste0(mode, "_", peak_id) %in% final_peaks_adducts$id) {
                        # Editing file for new adduct
                        print_message("Adduct modification in mat file for peak ", mode, " ", peak_id)
                        peak_data <-  readLines(new_peak_path, -1)
                        peak_data[5] <- paste0("PRECURSORTYPE: ",
                                               final_peaks_adducts[final_peaks_adducts$id == paste0(mode, "_", peak_id),
                                                                   "adduct"][1])
                        writeLines(peak_data, new_peak_path)
                    }
                }
            }
        }
    }

    print_message("Filtering top peaks done.")
}
