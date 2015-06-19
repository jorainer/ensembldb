library(shiny)
## start with runApp("jo_test")

shinyUI(fluidPage(

            ## Application title
            titlePanel("Get gene/transcript/exon annotations"),

            fluidRow(
                shiny::column(3,
                              uiOutput("packages")
                              ),
                shiny::column(3,
                              ##                    div(
                              h4("EnsDb annotation:"),
                              textOutput("metadata_organism"),
                              textOutput("metadata_ensembl"),
                              textOutput("metadata_genome"),
                              " "
                              ##                    )
                              ),
                shiny::column(4,
                              h4("Hints:"),
                              tags$li("Enter comma or whitespace separated values so search for multiple e.g. genes."),
                              tags$li("Use a % sign and condition like for partial matching."))
            ),
            fluidRow(
                shiny::column(4,
                              " "
                              )
            ),
            fluidRow(
                shiny::column(2,
                              selectInput("type", NA,
                                          choices=c("Gene name", "Chrom name",
                                                    "Gene biotype", "Tx biotype"),
                                          selected="Gene name")
                              ),
                shiny::column(1,
                              selectInput("condition", NA,
                                          choices=c("=", "!=", "like", "in"),
                                          selected="=")
                              ),
                shiny::column(2,
                              textInput("geneName", NA, value="")
                              ),
                shiny::column(2,
                              submitButton("Go!")
                              )
            ),
            fluidRow(
                mainPanel(
                    tabsetPanel(
                        tabPanel('Genes',
                                 dataTableOutput("Genes")
                                 ),
                        tabPanel('Transcripts',
                                 dataTableOutput("Transcripts")
                                 ),
                        tabPanel('Exons',
                                 dataTableOutput("Exons")
                                 )
                    )
                    ##h4(textOutput("genename"))
                )
            )
        ))


