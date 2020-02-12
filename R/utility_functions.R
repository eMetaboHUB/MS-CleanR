
#' Increment the final number in strings.
#'
#' @param strings_to_increment A list of strings to increment.
#' @param patterns A list of strings of patterns indicating which strings needs to be incremented.
#' @return A list of strings identical to \code{strings_to_increment} but with strings names corresponding to \code{patterns} having had their final number incremented.
increment_strings <- function(strings_to_increment, patterns) {
    new_names <- data.frame(old = strings_to_increment, new = NA, stringsAsFactors = FALSE)
    for(pattern in patterns) {
        start_pattern <- paste0("^", pattern)
        if (nrow(new_names[new_names$old == pattern,]) > 0) new_names[new_names$old == pattern, "old"] <- paste0(pattern, ".0")
        new_names[grepl(start_pattern, new_names$old), "new"] <- stringr::str_match(new_names[grepl(start_pattern, new_names$old), "old"],
                                                                                    ".+?(\\d+)$")[,2]
        new_names[grepl(start_pattern, new_names$old), "new"] <- paste(pattern,
                                                                       as.numeric(new_names[grepl(start_pattern, new_names$old), "new"]) + 1,
                                                                       sep = ".")
    }
    new_names[is.na(new_names$new), "new"] <- new_names[is.na(new_names$new), "old"]
    return(new_names$new)
}


#' Extract data concatenated in a single column.
#'
#' @param column_to_divide A 1-column data.frame with data to divide in several columns, with format "key1=value1,key2=value2,...".
#' @return A data.frame containing X columns.
extract_concatenated_data_from_column <- function(column_to_divide) {
    titles <- stringr::str_match_all(column_to_divide[1,1], "(^|,)(\\w+)=")[[1]][,3]
    if (length(titles) == 0) return(NULL)
    else {
        for(i in 2:length(titles)) {
            tmp <- data.frame(stringr::str_split_fixed(column_to_divide[[1]], paste0(",", titles[i], "="), n = 2))
            column_to_divide <- cbind(column_to_divide, tmp$X1)
            names(column_to_divide)[ncol(column_to_divide)] <- titles[i-1]
            if(i == 2) {
                column_to_divide[, titles[1]] <- gsub(paste0(titles[1], "=(.?)"),
                                                      "\\1",
                                                      as.character(column_to_divide[, titles[1]]))
            }
            column_to_divide[1] <- tmp$X2
        }
        column_to_divide <- cbind(column_to_divide, tmp$X2)
        names(column_to_divide)[ncol(column_to_divide)] <- titles[i]
        column_to_divide[1] <- NULL
        column_to_divide[column_to_divide == ""] <- NA
        return(column_to_divide)
    }
}



#' Transform a matrix to a data.frame.
#'
#' @param matrix_in A matrix to transform.
#' @param value_name A string indicating how to name the value column.
#' @return A 3-columns data.frame
matrix_to_df <- function(matrix_in, value_name = "Freq") {
    tmp <- matrix_in
    tmp[lower.tri(tmp, diag=TRUE)] <- NA # put NA
    tmp <- as.data.frame(as.table(tmp)) # as a dataframe
    tmp <- stats::na.omit(tmp) # remove NA
    names(tmp)[names(tmp) == "Freq"] <- value_name
    return(tmp)
}



#' Print number of peaks in data in a message for the user.
#'
#' @param peaks A data.frame contaning peaks data.
#' @param msg A string containing the message content.
print_peaks_status <- function(peaks, msg) {
    print_message(msg,
                  nrow(peaks[!is.na(peaks$source) & peaks$source == "pos",]), "positive,",
                  nrow(peaks[!is.na(peaks$source) & peaks$source == "neg",]), "negative,",
                  nrow(peaks[ is.na(peaks$source),]),                         "NA,",
                  nrow(peaks),                                                "total"
            )
}



#' Extract characters from a string from the right.
#'
#' @param x A string.
#' @param n An integer with the number of characters to get.
#' @return A string of length n.
substrRight <- function(x, n = 1) substr(x, nchar(x)-n+1, nchar(x))



#' Print a message in green.
#'
#' @param ... Strings to print.
print_message <- function(...) {
    if (get("shiny_running", envir = mscleanrCache)) cat(..., "\n")
    else                                             message(crayon::green(...))
}



#' Print a warning in yellow.
#'
#' @param ... Strings to print.
print_warning <- function(...) {
    if (get("shiny_running", envir = mscleanrCache)) cat("/!\\", ..., "\n")
    else                                             message(crayon::yellow(...))
}



#' Stop the script after printing an error.
#'
#' @param generic_msg Print the generic help message after the error.
#' @param ... Strings to print.
stop_script <- function(..., generic_msg = TRUE) {
    msg <- paste0(...)
    if (generic_msg) msg <- paste0(msg, "\nCheck your data and your parameters. If the problem persists, you can contact guillaume.marti@univ-tlse3.fr.")
    stop(crayon::red(msg), call. = FALSE)
}



#' Ask the user confirmation for overwriting if the path already exists.
#'
#' @param path A string indicating the path meant to be overwritten.
get_confirmation_for_overwriting <- function(path) {
    user_choice <- readline(print_warning(path, " already exists.\n",
                                          "Do you want to overwrite it? [y/n]\n",
                                          "Set 'overwrite = TRUE' to automatically overwrite the existing analysis."))
    return(tolower(user_choice) == "y")
}



#' Concatenate all MSP files present in the folder at \code{path}.
#'
#' @param path The path of the folder containing MSP files to concatenate.
#' @param output_file The name of the concatenated file.
#'
#' @export
concat_msp_files <- function(path, output_file = "_all_.msp") {
    if (!file.exists(path)) stop_script("Cannot find ", path)
    files <- list.files(path = path, pattern = "\\.msp$")
    if (length(files) > 0) {
        output_conn <- file(file.path(path, output_file), open = "a")
        for (input_f in files) {
            text <- readLines(file.path(path, input_f))
            writeLines(text, output_conn)
        }
        close(output_conn)
    }
    print_message(length(files), "MSP files concatenated in", file.path(path, output_file))
}

