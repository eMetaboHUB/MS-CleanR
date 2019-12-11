
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


    ui <- shiny::tagList(
        shiny::navbarPage(
            theme = shinythemes::shinytheme("paper"),
            "mscleanr",


            shiny::tabPanel(
                "Global parameters",

                shiny::sidebarPanel(
                    shiny::h2("TODO"),
                    "Help message."
                ),

                shiny::mainPanel(
                    #                   project_directory,
                    # TODO
                    #fileInput("upload", NULL),
                    #                   overwrite = TRUE)
                    shiny::checkboxInput("overwrite",
                                         "Overwrite existing results?",
                                         value = TRUE)
                )
            ),


            shiny::tabPanel(
                "Clean MS-DIAL data",

                shiny::sidebarPanel(
                    shiny::h2("TODO"),
                    "Help message."
                ),

                shiny::mainPanel(
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

                    # TODO
                    shiny::helpText("You can optionally import personal reference files for adducts and neutral losses.",
                                    "By default, data displayed in the Datasets tab will be used."),
                    shiny::fluidRow(
                        #                   user_pos_adducts_refs = NA,
                        shiny::column(6,
                            shiny::fileInput("user_pos_adducts_refs",
                                             "(OPTIONAL) Upload a reference file for adducts in positive mode.",
                                             width = "100%")
                        ),
                        #                   user_neg_adducts_refs = NA,
                        shiny::column(6,
                            shiny::fileInput("user_neg_adducts_refs",
                                             "(OPTIONAL) Upload a reference file for adducts in negative mode.",
                                             width = "100%")
                        )
                    ),
                    shiny::fluidRow(
                        #                   user_pos_neutral_refs = NA,
                        shiny::column(6,
                            shiny::fileInput("user_pos_neutral_refs",
                                             "(OPTIONAL) Upload a reference file for neutral losses in positive mode.",
                                             width = "100%")
                        ),
                        #                   user_neg_neutral_refs = NA,
                        shiny::column(6,
                            shiny::fileInput("user_neg_neutral_refs",
                                             "(OPTIONAL) Upload a reference file for neutral losses in negative mode.",
                                             width = "100%")
                        )
                    ),

                    shiny::actionButton("button_clean",
                                        "Clean MS-DIAL data",
                                        width = "100%",
                                        class = "btn btn-primary"),
                    shiny::verbatimTextOutput("clean")
                )
            ),


            shiny::tabPanel(
                "Keep top peaks by cluster",

                shiny::sidebarPanel(
                    shiny::h2("TODO"),
                    "Help message."
                ),

                shiny::mainPanel(
                    shiny::fluidRow(
                        #                selection_criterion = "degree",
                        shiny::column(4,
                            shiny::checkboxGroupInput("selection_criteria",
                                                      "Selection criteria",
                                                      c("Degree (most connected peaks)", "Intensity (most intense peaks)"))
                        ),
                        #                n = 1,
                        shiny::column(4,
                            shiny::numericInput("n",
                                                "Number of peaks to keep (by cluster)",
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

                    shiny::actionButton("button_keep",
                                        "Keep top peaks by cluster",
                                        width = "100%",
                                        class = "btn btn-primary"),
                    shiny::verbatimTextOutput("keep")
                )
            ),


            shiny::tabPanel(
                "Launch MS-FINDER annotation",

                shiny::sidebarPanel(
                    shiny::h2("TODO"),
                    "Help message."
                ),

                shiny::mainPanel(
                    # TODO
                    # compound_levels = c("1a", "1b"),
                    shiny::textInput("compound_levels",
                                     "Indicate the compound levels in your annotation files, separated by commas (leave blank if none).",
                                     placeholder = "direct,extension,...",
                                     width = "90%"),
                    # biosoc_levels = c("genre", "family"),
                    shiny::textInput("biosoc_levels",
                                     "Indicate the biosource levels in your annotation process, separated by commas (leave blank if none).",
                                     placeholder = "genus,family,...",
                                     width = "90%"),
                    # levels_scores = list("1a" = 2, "1b" = 1.5, "genre" = 2, "family" = 1.5),
                    shiny::textInput("levels_scores",
                                     "Indicate the scores multipliers associated to your compound or biosource levels, separated by commas (leave blank if none).",
                                     placeholder = "direct:2,genus:2,family:1.5,...",
                                     width = "90%"),

                    shiny::actionButton("button_launch",
                                        "Launch MS-FINDER annotation",
                                        width = "100%",
                                        class = "btn btn-primary"),
                    shiny::verbatimTextOutput("launch")
                )
            ),


            shiny::tabPanel(
                "Convert annotated peaks to MSP files",

                shiny::sidebarPanel(
                    shiny::h2("TODO"),
                    "Help message."
                ),

                shiny::mainPanel(
                    shiny::numericInput("min_score",
                                        "Minimum score to convert peak to MSP",
                                        value = 10,
                                        min = 0),

                    shiny::actionButton("button_convert",
                                        "Convert peaks to MSP files",
                                        width = "100%",
                                        class = "btn btn-primary"),
                    shiny::verbatimTextOutput("convert")
                )
            ),


            shiny::tabPanel(
                "Datasets",

                shiny::sidebarPanel(
                    shiny::h2("TODO"),
                    "Help message."
                ),

                shiny::mainPanel(
                    shiny::tabsetPanel(
                        shiny::tabPanel("Positive adducts",                shiny::dataTableOutput("adducts_pos")),
                        shiny::tabPanel("Negative adducts",                shiny::dataTableOutput("adducts_neg")),
                        shiny::tabPanel("Neutral losses in positive mode", shiny::dataTableOutput("nn_pos")),
                        shiny::tabPanel("Neutral losses in negative mode", shiny::dataTableOutput("nn_neg")),
                        shiny::tabPanel("Isotopes",                        shiny::dataTableOutput("isotopes"))
                    )
                )
            )
        )
    )



    server <- function(input, output, session) {
        output$clean <- shiny::renderPrint({})
        output$keep <- shiny::renderPrint({})
        output$launch <- shiny::renderPrint({})
        output$convert <- shiny::renderPrint({})

        # Datasets
        output$adducts_pos <- shiny::renderDataTable(mass_adducts_pos, escape = FALSE, options = list(pageLength = 10))
        output$adducts_neg <- shiny::renderDataTable(mass_adducts_neg, escape = FALSE, options = list(pageLength = 10))
        output$nn_pos      <- shiny::renderDataTable(mass_neutral_loss_pos, escape = FALSE, options = list(pageLength = 10))
        output$nn_neg      <- shiny::renderDataTable(mass_neutral_loss_neg, escape = FALSE, options = list(pageLength = 10))
        output$isotopes    <- shiny::renderDataTable(mass_isotopes, escape = FALSE, options = list(pageLength = 10))
    }

    shiny::runApp(shiny::shinyApp(ui, server))
}
