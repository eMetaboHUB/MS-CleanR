
#' Convert the final CSV file post annotations to MSP format.
#'
#' @eval recurrent_params("project_directory")
#' @param min_score Peaks must have a final score >= \code{min_score} to be exported to the msp file.
#' @export
convert_csv_to_msp <- function(project_directory, min_score = 20) {
    csv_data <- import_data(project_directory, "annotated_data-normalized")
    samples <- import_data(project_directory, "samples")

    mandatory <- c("Structure", "Formula", "source", "Average.Rt.min.", "PRECURSORMZ", "PRECURSORTYPE", "Final.score", "InChIKey", "SMILES", "MSMS.count", "MS.MS.spectrum")
    # optional <- c("Ontology", "Classyfire_subclass", "level", "Compound_level", "Internal_id")
    csv_data <- csv_data[csv_data$annotation_result != "Unknown compound"
                         & is.na(csv_data$annotation_warning)
                         & csv_data$Final.score >= min_score,
                         mandatory]
    print_message(nrow(csv_data), " peaks to convert in MSP.")

    if (nrow(csv_data) > 0) {
        # Name en 1er
        names(csv_data) <- toupper(names(csv_data))
        names(csv_data)[names(csv_data) == "STRUCTURE"] <- "NAME"
        names(csv_data)[names(csv_data) == "AVERAGE.RT.MIN."] <- "RETENTIONTIME"
        names(csv_data)[names(csv_data) == "SOURCE"] <- "IONMODE"
        csv_data[csv_data$IONMODE == "neg", "IONMODE"] <- "Negative"
        csv_data[csv_data$IONMODE == "pos", "IONMODE"] <- "Positive"

        output <- ""
        for (i in 1:nrow(csv_data)) {
            output <- paste0(output, dataframe_row_to_txt(csv_data[i,]), "\n\n")
        }

        output_file <- file(get_project_file_path(project_directory, "msp"))
        writeLines(output, output_file)
        close(output_file)
        print_message("Peaks converted, see file", get_project_file_path(project_directory, "msp"),".")
    }
}


dataframe_row_to_txt <- function(df_row) {
    output <- ""
    for (n in names(df_row)[!names(df_row) %in% c("MSMS.COUNT", "MS.MS.SPECTRUM")]) {
        output <- paste0(output, n, ": ", as.character(df_row[1, n]), "\n")
    }
    output <- paste0(output,
                     "MSTYPE: MS2\n",
                     "Num Peaks: ", as.character(df_row[1, "MSMS.COUNT"]), "\n")
    for (msms in strsplit(df_row$MS.MS.SPECTRUM, " ")) {
        for (details in strsplit(msms, ":")) {
            output <- paste0(output, details[1], "\t", details[2], "\n")
        }
    }
    return(output)
}
