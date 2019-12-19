
#' Launch the Graphical User Interface (GUI)
#'
#' @export
runGUI <- function() {

    # datasets loading
    utils::data("mass_neutral_loss_neg", package="mscleanr", envir=environment())
    mass_neutral_loss_neg <- get("mass_neutral_loss_neg", envir=environment())  # redundant, but makes the syntax checker happy
    utils::data("mass_neutral_loss_pos", package="mscleanr", envir=environment())
    mass_neutral_loss_pos <- get("mass_neutral_loss_pos", envir=environment())  # redundant, but makes the syntax checker happy
    utils::data("mass_adducts_neg", package="mscleanr", envir=environment())
    mass_adducts_neg <- get("mass_adducts_neg", envir=environment())  # redundant, but makes the syntax checker happy
    utils::data("mass_adducts_pos", package="mscleanr", envir=environment())
    mass_adducts_pos <- get("mass_adducts_pos", envir=environment())  # redundant, but makes the syntax checker happy
    utils::data("mass_isotopes", package="mscleanr", envir=environment())
    mass_isotopes <- get("mass_isotopes", envir=environment())  # redundant, but makes the syntax checker happy

    # personalized elements
    referenceInput <- function(id, text) {
        shiny::tabPanel(text,
            shiny::fluidRow(
                shiny::column(8,
                    shiny::fileInput(id,
                                     paste0("(OPTIONAL) Upload a CSV reference file for ", tolower(text), "."),
                                     accept = ".csv",
                                     placeholder = "  No file selected",
                                     width = "100%")
                ),
                shiny::column(4, shiny::tableOutput(id))
            )
        )
    }

    levelsInput <- function(id, label, placeholder, blank_if_none = TRUE) {
        shiny::fluidRow(
            shiny::column(8,
                shiny::textInput(id,
                                 paste0("Indicate the ",
                                        label,
                                        ", separated by commas",
                                        ifelse(blank_if_none, " (leave blank if none).", ".")),
                                 placeholder = placeholder,
                                 width = "90%")),
            shiny::column(4, shiny::tableOutput(id)))
    }

    dataset_panel <- function(title, id) shiny::tabPanel(title,
                                                         shiny::downloadButton(paste0("dl_", id)),
                                                         shiny::dataTableOutput(id))

    mainButton <- function(id, label) shiny::actionButton(id, label,  width = "100%", class = "btn btn-success")

    warning_perso <- function(condition, text) {
        shiny::conditionalPanel(condition = condition,
                                shiny::div(class = "panel panel-warning",
                                           shiny::div(class = "panel-heading",
                                                      shiny::h3(class = "panel-title", text))))
    }

    warning_wd <- function() warning_perso("document.getElementById('button_clean').disabled",
                                           "Please select the project directory.")

    fun_waiter <- function(h1)  waiter::Waiter$new(html = shiny::div(shiny::h1(h1, style = "color: white"),
                                                                     waiter::spin_double_bounce()),
                                                   color = "#4CAF50")


    ui <- shiny::tagList(
        shinyjs::useShinyjs(),
        waiter::use_waiter(),

        shiny::navbarPage(
            theme = shinythemes::shinytheme("paper"),
            "mscleanr",

            shiny::tabPanel(
                "Project directory",

                shiny::sidebarPanel(
                    shiny::h5("Project directory"),
                    "Select the directory containing the data you want to analyze.",
                    shiny::br(),
                    shiny::h5("Overwrite"),
                    "If checked, any existing file resulting from a previous analysis will be overwritten without asking."
                ),

                shiny::mainPanel(
                    warning_wd(),
                    #                   project_directory,
                    shinyFiles::shinyDirButton("project_directory",
                                               "Select the project directory",
                                               "Select the project directory",
                                               class = "btn btn-success",
                                               style = "width: 100%"),
                    shiny::verbatimTextOutput("project_directory_path"),
                    #                   overwrite = TRUE)
                    shiny::checkboxInput("overwrite",
                                         "Overwrite existing results?",
                                         value = TRUE)
                )
            ),


            shiny::tabPanel(
                "Clean MS-DIAL data",

                shiny::sidebarPanel(
                    shiny::h6("Combine positive and negative files from MSDial and filter peaks according to user parameters."),
                    shiny::br(),
                    shiny::h5("Filters"),
                    "Check which filters you want to use to clean your MS data.",
                    shiny::h5("Deltas"),
                    "Indicates the acceptable retention time and mass differences to consider that peaks are related.",
                    shiny::h5("Clusterisation options"),
                    "You can choose to use the Pearson correlation between peaks as a supplementary data used
                    during clusterisation",
                    shiny::h5("(Optional) Reference files"),
                    "Optionally, you can use your own files for adducts and neutral losses.
                    See the documentation for more information."
                ),

                shiny::mainPanel(
                    warning_wd(),
                    shiny::fluidRow(
                        #                   filter_blk = TRUE,
                        #                   filter_mz = TRUE,
                        #                   filter_rsd = TRUE,
                        shiny::column(4,
                            shiny::checkboxGroupInput("filters",
                                                      "What filters to use?",
                                                      choiceNames = c("Blank ratio", "Mass filter", "RSD filter"),
                                                      choiceValues = c("filter_blk", "filter_mz", "filter_rsd"),
                                                      selected = c("filter_blk", "filter_mz", "filter_rsd"))
                        ),
                        #                   filter_blk_threshold = 0.8,
                        shiny::column(4,
                            shiny::conditionalPanel(
                                condition = "input.filters.includes('filter_blk')",  # JS expression
                                shiny::numericInput("filter_blk_threshold",
                                                    "Minimum blank ratio",
                                                    value = 0.8,
                                                    min = 0,
                                                    max = 1,
                                                    step = 0.05)
                            )
                        ),
                        #                   filter_rsd_threshold = 30,
                        shiny::column(4,
                            shiny::conditionalPanel(
                                condition = "input.filters.includes('filter_rsd')",  # JS expression
                                shiny::numericInput("filter_rsd_threshold",
                                                    "Maximum RSD",
                                                    value = 30,
                                                    min = 0,
                                                    max = 100,
                                                    step = 5)
                            )
                        )
                    ),
                    shiny::hr(),

                    shiny::fluidRow(
                        #                   threshold_mz = 0.005,
                        shiny::column(4,
                            shiny::numericInput("threshold_mz",
                                                "Maximum mass difference",
                                                value = 0.005,
                                                min = 0,
                                                step = 0.001)
                        ),
                        #                   threshold_rt = 0.025,
                        shiny::column(4,
                            shiny::numericInput("threshold_rt",
                                                "Maximum retention time difference",
                                                value = 0.025,
                                                min = 0,
                                                step = 0.001)
                        )
                    ),
                    shiny::hr(),

                    shiny::fluidRow(
                        #                   compute_pearson_correlation = TRUE,
                        shiny::column(4,
                            shiny::checkboxInput("compute_pearson_correlation",
                                                 "Use Pearson correlation to compute clusters?",
                                                 value = TRUE),
                            shiny::br()
                        ),
                        #                   pearson_correlation_threshold = 0.8,
                        shiny::column(4,
                            shiny::conditionalPanel(
                                condition = "input.compute_pearson_correlation",  # JS expression
                                shiny::numericInput("pearson_correlation_threshold",
                                                    "Minimum correlation",
                                                    value = 0.8,
                                                    min = 0,
                                                    max = 1,
                                                    step = 0.05)
                            )
                        ),
                        #                   pearson_p_value = 0.05,
                        shiny::column(4,
                            shiny::conditionalPanel(
                                condition = "input.compute_pearson_correlation",
                                shiny::selectInput("pearson_p_value",
                                                   "Maximum p-value",
                                                   c("No filtering", 0.05, 0.01),
                                                   multiple = FALSE,
                                                   selected = 0.05)
                            )
                        )
                    ),
                    shiny::hr(),

                    "You can optionally import personal reference files for adducts and neutral losses.",
                    "By default, data displayed in the Datasets tab will be used.",
                    shiny::checkboxInput("user_refs",
                                         "Use personal reference files?",
                                         value = FALSE),
                    shiny::conditionalPanel(
                        condition = "input.user_refs",
                        shiny::tabsetPanel(
                            #                   user_pos_adducts_refs = NA,
                            referenceInput("user_pos_adducts_refs", "Positive adducts"),
                            #                   user_neg_adducts_refs = NA,
                            referenceInput("user_neg_adducts_refs", "Negative adducts"),
                            #                   user_pos_neutral_refs = NA,
                            referenceInput("user_pos_neutral_refs", "Neutral losses in positive mode"),
                            #                   user_neg_neutral_refs = NA,
                            referenceInput("user_neg_neutral_refs", "Neutral losses in negative mode")
                        )
                    ),

                    mainButton("button_clean", "Clean MS-DIAL data"),
                    shiny::verbatimTextOutput("clean")
                )
            ),


            shiny::tabPanel(
                "Keep top peaks by cluster",

                shiny::sidebarPanel(
                    shiny::h6("Keeps only the top peaks by cluster and by method."),
                    shiny::br(),
                    shiny::h5("Selection mode"),
                    "Peaks selected can be the most intense (intensity), most connected (degree), or both.",
                    "If there are ties, the number of peaks selected can be greater than the number requested by the user.",
                    shiny::h5("Exporting filtered peaks"),
                    "Copy MAT files corresponding to the selected peaks in a new folder for an faster analysis in MS-FINDER."
                ),

                shiny::mainPanel(
                    warning_wd(),
                    warning_perso("!(input.selection_criteria.includes('intensity') | input.selection_criteria.includes('degree'))",
                                  "Please select at least one selection criterion."),
                    shiny::fluidRow(
                        #                selection_criterion = "degree",
                        shiny::column(4,
                            shiny::checkboxGroupInput("selection_criteria",
                                                      "Selection criteria",
                                                      choiceNames = c("Intensity (most intense peaks)",
                                                                      "Degree (most connected peaks)"),
                                                      choiceValues = c("intensity", "degree"),
                                                      selected = "intensity")
                        ),
                        #                n = 1,
                        shiny::column(4,
                            shiny::numericInput("n",
                                                "Number of peaks to keep (by cluster and by method)",
                                                value = 1,
                                                min = 1,
                                                step = 1)
                        ),
                        #                export_filtered_peaks = TRUE,
                        shiny::column(4,
                            shiny::checkboxInput("export_filtered_peaks",
                                                 "Export final peaks in a new folder?",
                                                 value = TRUE)
                        )
                    ),

                    mainButton("button_keep", "Keep top peaks by cluster"),
                    shiny::verbatimTextOutput("keep")
                )
            ),


            shiny::tabPanel(
                "Launch MS-FINDER annotation",

                shiny::sidebarPanel(
                    shiny::h6("Annotates peaks based on files extracted from MSFinder."),
                    shiny::br(),
                    shiny::h5("Compound levels"),
                    "The list of compound levels to consider, in the given order (from more important to least important).",
                    shiny::h5("Biosource levels"),
                    "The list of biosource levels to consider, in the given order (from more important to least important).",
                    "They must correspond to the folders containing MS-FINDER files in your project directory.",
                    "The level 'generic' is always added as the last biosource level considered.",
                    shiny::h5("Levels scores"),
                    "A list of levels names and their corresponding multiplier to adapt final annotation scores."
                ),

                shiny::mainPanel(
                    warning_wd(),
                    # compound_levels = c("1a", "1b"),
                    levelsInput("compound_levels",
                                "compound levels in your annotation files",
                                "direct,extension,..."),
                    # biosoc_levels = c("genre", "family"),
                    levelsInput("biosoc_levels",
                                "biosource levels in your annotation process",
                                "genus,family,...",
                                blank_if_none = FALSE),
                    # levels_scores = list("1a" = 2, "1b" = 1.5, "genre" = 2, "family" = 1.5),
                    levelsInput("levels_scores",
                                "scores multipliers associated to your compound or biosource levels",
                                "direct:2,genus:2,family:1.5,..."),

                    mainButton("button_launch", "Launch MS-FINDER annotation"),
                    shiny::verbatimTextOutput("launch")
                )
            ),


            shiny::tabPanel(
                "Convert annotated peaks to MSP files",

                shiny::sidebarPanel(
                    shiny::h6("Convert the final CSV file post annotations to MSP format."),
                    shiny::br(),
                    shiny::h5("Minimum score"),
                    "Minimum annotation score needed to export peaks to the MSP files."
                ),

                shiny::mainPanel(
                    warning_wd(),
                    shiny::numericInput("min_score",
                                        "Minimum score to convert peak to MSP",
                                        value = 10,
                                        min = 0),

                    mainButton("button_convert", "Convert peaks to MSP files"),
                    shiny::verbatimTextOutput("convert")
                )
            ),


            shiny::tabPanel(
                "Datasets",

                shiny::sidebarPanel(
                    shiny::h6("Datasets used by default during the analysis if none are provided by the user."),
                    shiny::br(),
                    shiny::h5("Adducts"),
                    "Mass differences of adducts with the original compound M, in positive and negative mode.",
                    shiny::h5("Neutral losses"),
                    "Masses of neutral losses, in positive and negative mode.",
                    shiny::h5("Isotopes"),
                    "Masses of isotopes."
                ),

                shiny::mainPanel(
                    shiny::tabsetPanel(
                        dataset_panel("Positive adducts",                "adducts_pos"),
                        dataset_panel("Negative adducts",                "adducts_neg"),
                        dataset_panel("Neutral losses in positive mode", "nl_pos"),
                        dataset_panel("Neutral losses in negative mode", "nl_neg"),
                        shiny::tabPanel("Isotopes", shiny::dataTableOutput("isotopes"))
                    )
                )
            )
        )
    )



    server <- function(input, output, session) {
        # Project directory
        volumes <- c(Home = fs::path_home(), shinyFiles::getVolumes()())
        shinyFiles::shinyDirChoose(input,
                                   "project_directory",
                                   roots = volumes,
                                   session = session)

        output$project_directory_path <- shiny::renderPrint({
            if (is.integer(input$project_directory)) {
                cat("No project directory has been selected")
                shinyjs::disable("button_clean")
                shinyjs::disable("button_keep")
                shinyjs::disable("button_launch")
                shinyjs::disable("button_convert")
            } else {
                path <- shinyFiles::parseDirPath(volumes, input$project_directory)
                assign("analysis_directory", path, envir = mscleanrCache)
                cat(path)
                shinyjs::enable("button_clean")
                shinyjs::enable("button_keep")
                shinyjs::enable("button_launch")
                shinyjs::enable("button_convert")
            }
        })

        shiny::observeEvent(input$overwrite, { assign("overwrite", input$overwrite, envir = mscleanrCache) })


        # References
        user_ref_generic <- function(ref_file) {
            shiny::req(ref_file)
            data <- vroom::vroom(ref_file$datapath)
            shiny::validate(shiny::need(ncol(data) == 2, "Please upload a reference file with 2 columns (see documentation)."))
            return(data)
        }
        user_pos_adducts_refs <- shiny::reactive({ shiny::req(user_ref_generic(input$user_pos_adducts_refs))
                                                   user_ref_generic(input$user_pos_adducts_refs) })
        user_neg_adducts_refs <- shiny::reactive({ shiny::req(user_ref_generic(input$user_neg_adducts_refs))
                                                   user_ref_generic(input$user_neg_adducts_refs) })
        user_pos_neutral_refs <- shiny::reactive({ shiny::req(user_ref_generic(input$user_pos_neutral_refs))
                                                   user_ref_generic(input$user_pos_neutral_refs) })
        user_neg_neutral_refs <- shiny::reactive({ shiny::req(user_ref_generic(input$user_neg_neutral_refs))
                                                   user_ref_generic(input$user_neg_neutral_refs) })

        preview_ref <- function(ref_data) rbind(utils::head(ref_data, 2),
                                                stats::setNames(data.frame(a = "...", b = "..."), names(ref_data)))
        output$user_pos_adducts_refs <- shiny::renderTable(preview_ref(user_pos_adducts_refs()))
        output$user_neg_adducts_refs <- shiny::renderTable(preview_ref(user_neg_adducts_refs()))
        output$user_pos_neutral_refs <- shiny::renderTable(preview_ref(user_pos_neutral_refs()))
        output$user_neg_neutral_refs <- shiny::renderTable(preview_ref(user_neg_neutral_refs()))


        # clean_msdial_data
        clean_params <- shiny::eventReactive(input$button_clean, {
            check_ref <- function(ref) {
                tryCatch({ shiny::validate(shiny::need(ref, "not found")) ; return(ref) },
                         error = function(e) {})
            }
            get_ref   <- function(ref) if(is.null(check_ref(ref))) NA else as.data.frame(ref)
            list(
                filter_blk                    = "filter_blk" %in% input$filters,
                filter_blk_threshold          = input$filter_blk_threshold,
                filter_mz                     = "filter_mz" %in% input$filters,
                filter_rsd                    = "filter_rsd" %in% input$filters,
                filter_rsd_threshold          = input$filter_rsd_threshold,
                threshold_mz                  = input$threshold_mz,
                threshold_rt                  = input$threshold_rt,
                compute_pearson_correlation   = input$compute_pearson_correlation,
                pearson_correlation_threshold = input$pearson_correlation_threshold,
                pearson_p_value               = input$pearson_p_value,
                user_pos_adducts_refs         = get_ref(user_pos_adducts_refs()),
                user_neg_adducts_refs         = get_ref(user_neg_adducts_refs()),
                user_pos_neutral_refs         = get_ref(user_pos_neutral_refs()),
                user_neg_neutral_refs         = get_ref(user_neg_neutral_refs())
            )
        })

        output$clean <- shiny::renderPrint({
            shiny::req(input$button_clean > 0)
            shiny::req(exists("analysis_directory", envir = mscleanrCache))
            shiny::req(clean_params())
            waiter <- fun_waiter("Cleaning data...")
            waiter$show()
            tryCatch(
                do.call(clean_msdial_data, clean_params()),
                finally = waiter$hide()
            )
        })


        # keep_top_peaks
        keep_params <- shiny::eventReactive(input$button_keep, {
            shiny::req("degree" %in% input$selection_criteria | "intensity" %in% input$selection_criteria)
            if ("degree" %in% input$selection_criteria) {
                if ("intensity" %in% input$selection_criteria) mode <- "both"
                else                                           mode <- "degree"
            } else                                             mode <- "intensity"
            list(selection_criterion   = mode,
                 n                     = input$n,
                 export_filtered_peaks = input$export_filtered_peaks)
        })

        output$keep <- shiny::renderPrint({
            shiny::req(input$button_keep > 0)
            shiny::req(exists("analysis_directory", envir = mscleanrCache))
            shiny::req(keep_params())
            waiter <- fun_waiter("Keeping only selected peaks...")
            waiter$show()
            tryCatch(
                do.call(keep_top_peaks, keep_params()),
                finally = waiter$hide()
            )
        })


        # launch_msfinder_annotation
        clean <- function(l) {
            tmp <- unlist(lapply(strsplit(l, ",", fixed = TRUE), function(x) trimws(x)))
            tmp <- tmp[tmp != ""]
            if (length(tmp) == 0) return(NULL) else return(unique(tmp))
        }

        compounds <- shiny::reactive(clean(input$compound_levels))
        biosoc <- shiny::reactive({
            tmp <- clean(input$biosoc_levels)
            if (is.null(tmp)) c("generic") else tmp
        })
        scores <- shiny::reactive({
            tmp_scores <- lapply(clean(input$levels_scores), strsplit, split = ":", fixed = TRUE)
            if (length(tmp_scores) == 0) scores <- NULL
            else {
                scores <- list()
                for (x in tmp_scores) scores[x[[1]][1]] <- ifelse(is.na(x[[1]][2]) | x[[1]][2] == "",
                                                                  1,
                                                                  as.numeric(x[[1]][2]))
            }
            scores
        })

        output$compound_levels <- shiny::renderTable({
            shiny::req(compounds())
            data.frame("Rank" = c(1:length(compounds())), "Compound level" = compounds())
        })
        output$biosoc_levels   <- shiny::renderTable(
            data.frame("Rank" = c(1:length(unique(c(biosoc(), "generic")))),
                       "Biosource level" = unique(c(biosoc(), "generic")))
        )
        output$levels_scores   <- shiny::renderTable(
            data.frame(Level = names(scores()), Multiplier = unlist(scores(), use.names = FALSE))
        )

        launch_params <- shiny::eventReactive(input$button_launch, {
            list(compound_levels = compounds(),
                 biosoc_levels   = biosoc(),
                 levels_scores   = scores())
        })

        output$launch <- shiny::renderPrint({
            shiny::req(input$button_launch > 0)
            shiny::req(exists("analysis_directory", envir = mscleanrCache))
            shiny::req(launch_params())
            waiter <- fun_waiter("Annotating peaks with MS-FINDER data...")
            waiter$show()
            tryCatch(
                do.call(launch_msfinder_annotation, launch_params()),
                finally = waiter$hide()
            )
        })


        # convert_csv_to_msp
        convert_params <- shiny::eventReactive(input$button_convert, list(min_score = input$min_score))

        output$convert <- shiny::renderPrint({
            shiny::req(input$button_convert > 0)
            shiny::req(exists("analysis_directory", envir = mscleanrCache))
            shiny::req(convert_params())
            waiter <- fun_waiter("Converting peaks to MSP...")
            waiter$show()
            tryCatch(
                do.call(convert_csv_to_msp, convert_params()),
                finally = waiter$hide()
            )
        })

        # Datasets
        output$adducts_pos <- shiny::renderDataTable(mass_adducts_pos,      escape = FALSE, options = list(pageLength = 10))
        output$adducts_neg <- shiny::renderDataTable(mass_adducts_neg,      escape = FALSE, options = list(pageLength = 10))
        output$nl_pos      <- shiny::renderDataTable(mass_neutral_loss_pos, escape = FALSE, options = list(pageLength = 10))
        output$nl_neg      <- shiny::renderDataTable(mass_neutral_loss_neg, escape = FALSE, options = list(pageLength = 10))
        output$isotopes    <- shiny::renderDataTable(mass_isotopes,         escape = FALSE, options = list(pageLength = 10))

        dl <- function(data, filename) shiny::downloadHandler(filename = filename,
                                                              content = function(file) vroom::vroom_write(data, file, delim = ","))
        output$dl_adducts_pos <- dl(mass_adducts_pos,      "positive_adducts.csv")
        output$dl_adducts_neg <- dl(mass_adducts_neg,      "negative_adducts.csv")
        output$dl_nl_pos      <- dl(mass_neutral_loss_pos, "positive_neutral_losses.csv")
        output$dl_nl_neg      <- dl(mass_neutral_loss_neg, "negative_neutral_losses.csv")



        shiny::onStop(function() { assign("shiny_running", FALSE, envir = mscleanrCache) })
    }



    shiny::runApp(shiny::shinyApp(ui = ui,
                                  server = server,
                                  onStart = function() { assign("shiny_running", TRUE, envir = mscleanrCache) }))
}
