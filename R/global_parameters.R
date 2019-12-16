
# Environment for the global variables
mscleanrCache <- new.env()
assign("overwrite",     FALSE, envir = mscleanrCache)
assign("shiny_running", FALSE, envir = mscleanrCache)


#' Declare the project directory for the analysis process.
#'
#' @param project_directory The path of the project directory.
#' @param overwrite A boolean indicating whether to stop the script if there is an existing analysis or overwrite it.
#' @export
choose_project_directory <- function(project_directory, overwrite = FALSE) {
    check_boolean(overwrite, "overwrite")
    assign("analysis_directory", project_directory, envir = mscleanrCache)
    assign("overwrite",          overwrite,         envir = mscleanrCache)
}

