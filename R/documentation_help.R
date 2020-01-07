
#' Return documentation info for recurrent parameters.
#'
#' @param ... A list of parameters whose documentation is needed.
recurrent_params <- function(...) {
    output_doc <- c()
    for (param_name in list(...)) {
        output_doc <- c(output_doc, paste0("@param ", param_name, " ", get(param_name, envir = mscleanrDocParams)))
    }
    return(output_doc)
}


# Contains documentation data
mscleanrDocParams <- new.env()
assign("biosoc_levels",
       "A vector containing the biosource levels to consider, in the given order.",
       envir = mscleanrDocParams)
assign("compound_levels",
       "A vector containing the compound levels to consider, in the given order.",
       envir = mscleanrDocParams)
assign("filetype",
       "A string indicating the needed directory or file.",
       envir = mscleanrDocParams)
assign("filter_blk",
       "A boolean indicating whether or not to delete rows with too much noise.",
       envir = mscleanrDocParams)
assign("filter_blk_threshold",
       "A numerical threshold for noise filtering: rows with ratio mean(blank columns)/mean(qc columns) >= \\code{filter_blk} are deleted (if there are no QC columns in the sample, the mean of standard columns is used, or the mean of all non-blank samples if needed).",
       envir = mscleanrDocParams)
assign("filter_mz",
       "A boolean indicating whether or not to delete rows with masses ending in .8 or .9 (masses not found in natural products).",
       envir = mscleanrDocParams)
assign("filter_rsd",
       "A boolean indicating whether or not to delete rows with too much relative standard deviation in each class.",
       envir = mscleanrDocParams)
assign("filter_rsd_threshold",
       "A numerical threshold for relative standard deviation filtering: rows with relative standard deviation >= \\code{filter_rsd_threshold} in each class are deleted. Only used if \\code{filter_rsd} is \\code{TRUE}.",
       envir = mscleanrDocParams)
assign("filter_rmd",
       "A boolean indicating whether or not to delete rows with a Relative Mass Defect outside of the range provided in \\code{filter_rmd_range}.",
       envir = mscleanrDocParams)
assign("filter_rmd_range",
       "A range of 2 integers indicating the acceptable Relative Mass Defects in ppm (only used if \\code{filter_rmd} is TRUE).",
       envir = mscleanrDocParams)
assign("filter_blk_ghost_peaks",
       "A boolean indicating whether or not to delete blank ghost peaks (only used if \\code{filter_blk} is TRUE, see publication for more information).",
       envir = mscleanrDocParams)
assign("level",
       "A string indicating the biosource level to consider.",
       envir = mscleanrDocParams)
assign("possibilities",
       "A data.frame of possible structures.",
       envir = mscleanrDocParams)
assign("source",
       'A string indicating which mode to read from: "pos" or "neg".',
       envir = mscleanrDocParams)
assign("threshold_mz",
       "A numerical value indicating the mass tolerance in Dalton for the detection of adducts and neutral losses.",
       envir = mscleanrDocParams)
assign("threshold_rt",
       "A numerical value indicating the retention time tolerance.",
       envir = mscleanrDocParams)
assign("user_neg_neutral_refs,user_pos_neutral_refs,user_pos_adducts_refs,user_neg_adducts_refs",
       "An optional 2-column data.frame containing information about neutral losses, positive adducts or negative adducts (one column for the name and one column for the mass difference with the base compound). If no data.frame is provided, the package default list is used.",
       envir = mscleanrDocParams)



# Contains output files and folders roles.
mscleanrDocRoles <- new.env()
assign("normalized_ms_data",
       "Normalized peaks info exported from MSDial.",
       envir = mscleanrDocRoles)
assign("msdial_peaks",
       "Individual peaks files exported from MSDial.",
       envir = mscleanrDocRoles)
assign("msdial_filtered_peaks",
       "Individual peaks files remaining after cleaning.",
       envir = mscleanrDocRoles)
assign("samples",
       "Information about samples.",
       envir = mscleanrDocRoles)
assign("adduct_links_graph",
       "Graph of adducts relations between peaks.",
       envir = mscleanrDocRoles)
assign("clusters_msdial",
       "MSDial data with clusters information.",
       envir = mscleanrDocRoles)
assign("clusters_graph",
       "Graph containing clusters and MSDial data.",
       envir = mscleanrDocRoles)
assign("links_ms_pre_selection",
       "Links between peaks before any processing.",
       envir = mscleanrDocRoles)
assign("links_ms_post_selection",
       "Links between peaks after selecting the best adduct for each peak.",
       envir = mscleanrDocRoles)
assign("links_ms_final",
       "Links between peaks after new phase of clustering and deletion of neutral losses and heteromers.",
       envir = mscleanrDocRoles)
assign("final_peaks_data",
       "Data for peaks remaining after all filters have been applied.",
       envir = mscleanrDocRoles)
assign("annotated_data",
       "Final auto-annotated data, complete with all possibilities.",
       envir = mscleanrDocRoles)
assign("annotated_data-manual_check",
       "Final auto-annotated data for clusters needing manual checking (only created in this case).",
       envir = mscleanrDocRoles)
assign("identifying_data",
       "Union of MSDial and MSFinder data, with clusters info.",
       envir = mscleanrDocRoles)
assign("annotated_data-cleaned",
       "Final auto-annotated data with only relevant peaks.",
       envir = mscleanrDocRoles)
assign("annotated_data-normalized",
       "Final auto-annotated data with only relevant peaks and normalized intensities (on 1000).",
       envir = mscleanrDocRoles)
