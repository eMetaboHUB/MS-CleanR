
#' Return documentation info for recurrent parameters.
#'
#' @param ... A list of parameters whose documentation is needed.
recurrent_params <- function(...) {
    known_doc <- get_doc()
    output_doc <- c()
    for (param_name in list(...)) output_doc <- c(output_doc, known_doc[[param_name]])
    return(output_doc)
}


#' Contains documentation data.
get_doc <- function() {
    return(list(
        "biosoc_levels"         = "@param biosoc_levels A vector containing the biosource levels to consider, in the given order.",
        "compound_levels"       = "@param compound_levels A vector containing the compound levels to consider, in the given order.",
        "filetype"              = "@param filetype A string indicating the needed directory or file.",
        "filter_blk"            = "@param filter_blk A boolean indicating whether or not to delete rows with too much noise.",
        "filter_blk_threshold"  = "@param filter_blk_threshold A numerical threshold for noise filtering: rows with ratio mean(blank columns)/mean(qc columns) >= \\code{filter_blk} are deleted (if there are no QC columns in the sample, the mean of standard columns is used, or the mean of all non-blank samples if needed).",
        "filter_mz"             = "@param filter_mz A boolean indicating whether or not to delete rows with masses ending in .8 or .9 (masses not found in natural products).",
        "filter_rsd"            = "@param filter_rsd A boolean indicating whether or not to delete rows with too much relative standard deviation in each class.",
        "filter_rsd_threshold"  = "@param filter_rsd_threshold A numerical threshold for relative standard deviation filtering: rows with relative standard deviation >= \\code{filter_rsd_threshold} in each class are deleted. Only used if \\code{filter_rsd} is \\code{TRUE}.",
        "level"                 = '@param level A string indicating the biosource level to consider.',
        "overwrite"             = "@param overwrite A boolean indicating whether to stop the script if there is an existing analysis or overwrite it.",
        "possibilities"         = "@param possibilities A data.frame of possible structures.",
        "project_directory"     = "@param project_directory The path of the project directory.",
        "source"                = '@param source A string indicating which mode to read from: "pos" or "neg".',
        "threshold_mz"          = "@param threshold_mz A numerical value indicating the mass tolerance in Dalton for the detection of adducts and neutral losses.",
        "threshold_rt"          = "@param threshold_rt A numerical value indicating the retention time tolerance.",
        "user_neg_neutral_refs,user_pos_neutral_refs,user_pos_adducts_refs,user_neg_adducts_refs" = "@param user_neg_neutral_refs,user_pos_neutral_refs,user_pos_adducts_refs,user_neg_adducts_refs An optional 2-column data.frame containing information about neutral losses, positive adducts or negative adducts (one column for the name and one column for the mass difference with the base compound). If no data.frame is provided, the package default list is used."
    ))
}



#' Contains output files and folders roles.
get_files_roles <- function() {
    return(list(
        "normalized_ms_data"          = "Normalized peaks info exported from MSDial.",
        "msdial_peaks"                = "Individual peaks files exported from MSDial.",
        "msdial_filtered_peaks"       = "Individual peaks files remaining after cleaning.",
        "msfinder_data"               = "",
        "cad_data"                    = "",
        "pda_data"                    = "",
        "inter_folder"                = "",
        "parameters"                  = "",
        "samples"                     = "Information about samples.",
        "msdial_detected_adducts"     = "",
        "adduct_links_graph"          = "Graph of adducts relations between peaks.",
        "final_links_graph"           = "",
        "computed_massdiff"           = "",
        "filtered_computed_massdiff"  = "",
        "clusters_msdial"             = "MSDial data with clusters information.",
        "clusters_final"              = "",
        "clusters_graph"              = "Graph containing clusters and MSDial data.",
        "links_ms_pre_selection"      = "Links between peaks before any processing.",
        "links_ms_post_selection"     = "Links between peaks after selecting the best adduct for each peak.",
        "links_ms_final"              = "Links between peaks after new phase of clustering and deletion of neutral losses and heteromers.",
        "final_peaks_data"            = "Data for peaks remaining after all filters have been applied.",
        "msdial_identified_peaks"     = "",
        "annotated_data"              = "Final auto-annotated data, complete with all possibilities.",
        "annotated_data-manual_check" = "Final auto-annotated data for clusters needing manual checking (only created in this case).",
        "identifying_data"            = "Union of MSDial and MSFinder data, with clusters info.",
        "final_folder"                = "",
        "annotated_data-cleaned"      = "Final auto-annotated data with only relevant peaks.",
        "annotated_data-normalized"   = "Final auto-annotated data with only relevant peaks and normalized intensities (on 1000)."
    ))
}
