
#' Convert the final CSV file post annotations to MSP format.
#'
#' @param all All peaks are exported to the MSP files.
#' @param min_score If \code{all} is \code{FALSE}, peaks must have a final score >= \code{min_score} to be exported to the MSP files.
#' @export
convert_csv_to_msp <- function(all = FALSE, min_score = 20) {
    check_for_convert_csv_to_msp(min_score)

    csv_data <- import_data("annotated_data-normalized")
    samples <- import_data("samples")

    mandatory <- c("Structure", "Formula", "source", "Average.Rt.min.", "Average.Mz", "Adduct.type", "Final.score", "InChIKey", "SMILES", "MSMS.count", "MS.MS.spectrum")
    # optional <- c("Ontology", "Classyfire_subclass", "level", "Compound_level", "Internal_id")
    if (all) {
        csv_data <- csv_data[mandatory]
    } else {
        csv_data <- csv_data[csv_data$annotation_result != "Unknown compound"
                             & is.na(csv_data$annotation_warning)
                             & csv_data$Final.score >= min_score,
                             mandatory]
    }
    print_message(nrow(csv_data), " peaks to convert in MSP.")

    if (nrow(csv_data) > 0) {
        # Name en 1er
        names(csv_data) <- toupper(names(csv_data))
        names(csv_data)[names(csv_data) == "STRUCTURE"]       <- "NAME"
        names(csv_data)[names(csv_data) == "AVERAGE.RT.MIN."] <- "RETENTIONTIME"
        names(csv_data)[names(csv_data) == "AVERAGE.MZ"]      <- "PRECURSORMZ"
        names(csv_data)[names(csv_data) == "ADDUCT.TYPE"]     <- "PRECURSORTYPE"
        names(csv_data)[names(csv_data) == "SOURCE"]          <- "IONMODE"
        csv_data[csv_data$IONMODE == "neg", "IONMODE"] <- "Negative"
        csv_data[csv_data$IONMODE == "pos", "IONMODE"] <- "Positive"

        for (source in c("pos", "neg")) {
            if (source == "pos") rows <- csv_data[csv_data$IONMODE == "Positive",]
            else                 rows <- csv_data[csv_data$IONMODE == "Negative",]

            if (nrow(rows) > 0) {
                output <- ""
                for (i in 1:nrow(rows)) {
                    output <- paste0(output, dataframe_row_to_txt(rows[i,]), "\n\n")
                }

                output_file <- file(get_project_file_path(paste0("msp_", source)))
                writeLines(output, output_file)
                close(output_file)
            }
        }

        print_message("Peaks converted, see MSP files in ", get_project_file_path("final_folder"),".")
    }
}


dataframe_row_to_txt <- function(df_row) {
    output <- ""
    for (n in names(df_row)[!names(df_row) %in% c("MSMS.COUNT", "MS.MS.SPECTRUM")]) {
        if (!is.na(df_row[1,n]) & !is.null(df_row[1, n]) & df_row[1, n] != "") {
            output <- paste0(output, n, ": ", as.character(df_row[1, n]), "\n")
        } else {
            output <- paste0(output, n, ": NA\n")
        }
    }
    msms_peaks <- strsplit(df_row$MS.MS.SPECTRUM, " ")
    output <- paste0(output,
                     "MSTYPE: MS2\n",
                     "Num Peaks: ", length(msms_peaks[[1]]), "\n")
    # "Num Peaks: ", as.character(df_row[1, "MSMS.COUNT"]), "\n")
    for (msms in msms_peaks) {
        for (details in strsplit(msms, ":")) {
            output <- paste0(output, details[1], "\t", details[2], "\n")
        }
    }
    return(output)
}
