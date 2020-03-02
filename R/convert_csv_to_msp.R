
#' Convert the final CSV file post annotations to MSP format.
#'
#' @param all All peaks are exported to the MSP files.
#' @param min_score If \code{all} is \code{FALSE}, peaks must have a final score >= \code{min_score} to be exported to the MSP files.
#' @export
convert_csv_to_msp <- function(all = FALSE, min_score = 20) {
    check_for_convert_csv_to_msp(min_score)
    export_params(as.list(environment()))

    samples  <- import_data("samples")

    csv_data <- import_data("annotated_data-cleaned")
    csv_data <- csv_data[c("annotation_result", "Structure", "Formula", "source", "Average.Rt.min.", "Average.Mz", "Adduct.type",
                           "Final.score", "InChIKey", "SMILES", "MSMS.count", "MS.MS.spectrum")]

    metadata <- import_data("annotated_data-normalized")
    metadata_mandatory <- c("cluster", "selected_feature", "annotation_result", "annotation_warning", "source", "Average.Rt.min.",
                            "Average.Mz", "Adduct.type", "Formula", "Structure", "Total.score", "Final.score", "Title", "PRECURSORMZ",
                            "PRECURSORTYPE", "Theoretical.mass", "Mass.error", "Formula.score", "InChIKey", "SMILES")
    metadata_optional <- names(metadata)[c(names(metadata) %in% c("level", "Compound_level", "Ontology",
                                                                  "Higher_biosoc", "Family", "Biosoc",
                                                                  "Classyfire_class", "Classyfire_subclass",
                                                                  "Nb_external_dbs", "Links"))]
    for (class in unique(samples[samples$Script_class != "Blank",]$Script_class)) {
        metadata[paste0("avg_", class)] <- rowMeans(metadata[samples[samples$Class == "FD", "Column_name"]])
        metadata_mandatory <- c(metadata_mandatory, paste0("avg_", class))
    }
    metadata <- metadata[c(metadata_mandatory, metadata_optional)]

    if (!all) {
        csv_data <- csv_data[csv_data$annotation_result != "Unknown compound" & csv_data$Final.score >= min_score,]
        metadata <- metadata[metadata$annotation_result != "Unknown compound" & metadata$Final.score >= min_score,]
    }
    csv_data$annotation_result <- NULL
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

                export_data(metadata[metadata$source == source,], paste0("msp_metadata_", source))
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
