library(shiny)
library(gapminder)
library(plotly)
library(ggplot2)
library(readr)
library(NGLVieweR)
library(shinyjs)
library(shinythemes)
library(shinyWidgets)

data <- read_tsv("../data/part_a_data.tsv")
new_data <- data
new_data <- data.frame(new_data)
#headers <- c('GENE', 'ACC', 'Total_MUT', 'UniProt_MUT', 'ClinVar_MUT', 'COSMIC_MUT', 'Position', 'SignalP_TM', 'SignalP_noTM', 'N_terminus', 'Topology')
headers <- c('GENE', 'Total_MUT', 'UniProt_MUT', 'ClinVar_MUT', 'COSMIC_MUT', 'Position', 'SignalP_TM', 'SignalP_noTM', 'N_terminus', 'Topology')
#headers <- c('GENE', 'Total_MUT', 'UniProt_MUT', 'ClinVar_MUT', 'COSMIC_MUT', 'Position', 'SignalP_TM', 'SignalP_noTM', 'N-terminus', 'Topology')
new_data <- new_data[headers]
RV <- reactiveValues(data = data.frame(new_data))

data2 <- read_tsv("../data/part_c_data.tsv")
new_data2 <- data2
new_data2 <- data.frame(new_data2)
#print (new_data2['N-terminus'])
headers2 <- c('GENE', 'ACC', 'Total_MUT', 'UniProt_MUT', 'ClinVar_MUT', 'COSMIC_MUT', 'Position', 'TMD', 'SignalP_TM', 'SignalP_noTM', 'N_terminus')
new_data2 <- new_data2[headers2]
RV2 <- reactiveValues(data = data.frame(new_data2))

# Define UI for application that draws a histogram
shinyUI(navbarPage("SPC as a Quality Control Enzyme",
        tabPanel("N-terminal analysis",
            fluidPage(#titlePanel("PART A"),
                      tags$head(
                        #tags$style(type="text/css", "label{ display: table-cell; text-align: center; vertical-align: middle; } .form-group { display: table-row;}")
                      ),
                  fluidRow(
                    column(3,
                      wellPanel(
                          align = 'center',
                          useShinyjs(),
                          
                          ## ---------------------------------------------------
                          #div(style="display: inline-block; padding-left: 3%;", actionButton("inputhelp", "", icon = icon("question-circle"))),
                          #hr(style = "border: 1px solid;"),
                          div(style="padding-left: 3%; padding-bottom: 5%", strong('Process all proteins?', style="display: inline-block; padding-bottom: 5px; font-size: 18px"),
                            div(style="display: inline-block; padding-left: 3%", actionButton("processHelp", "", icon = icon("question-circle"))),
                          ),
                          div(style="display: inline-block; width: 100px;",
                            switchInput(inputId = "all", label="All", value = TRUE, size = "mini", onLabel = "Yes", offLabel = "No"),
                          ),
                          div(style="display: inline-block; vertical-align:center; padding-left: 3%; width: 25px", strong('OR')),
                          div(style="display: inline-block; vertical-align:top; padding-left: 3%; width: 200px;",
                            selectizeInput("gene", "Choose a protein(s):",
                                           choices = NULL, multiple = TRUE),
                          ),
                          hr(style = "border: 1px solid;"),
                          
                          ## ---------------------------------------------------
                          div(style="padding-left: 3%; padding-bottom: 5%", strong('SignalP cutoffs', style="display: inline-block; vertical-align:center; font-size: 18px"),
                            div(style="display: inline-block; padding-left: 3%;vertical-align:center", actionButton("signalpHelp", "", icon = icon("question-circle"))),
                          ),
                          div(style="display: inline-block;vertical-align:top;horizontal-align:center; width: 170px;",
                            sliderInput("SignalP_TM", label = "TM version", min = 0,
                                        max = 1.0, value = c(0.0, 0.6))
                          ),
                          div(style="display: inline-block;vertical-align:top; width: 25px;"),
                          div(style="display: inline-block;vertical-align:top; width: 170px;",
                            sliderInput("SignalP_noTM", label = "no TM version", min = 0, 
                                        max = 1.0, value = c(0.6, 1.0))
                          ),
                          hr(style = "border: 1px solid;"),
                          
                          ## ---------------------------------------------------
                          div(style="display: inline-block; padding-left: 3%; width: 150px; font-size: 18px", strong('# Mutations >')),
                          div(style="display: inline-block; padding-left: 3%; text-align: left; width: 100px", numericInput('Num_MUT', '', value = 0)),
                          div(style="display: inline-block; padding-left: 3%;", actionButton("mutationHelp", "", icon = icon("question-circle"))),
                          hr(style = "border: 0.5px solid;"),
                          
                          ## ---------------------------------------------------
                          div(style="display: inline-block; padding-left: 3%; width: 200px; font-size: 18px", strong('N-term localization?')),
                          div(style="display: inline-block;", actionButton("ntermHelp", "", icon = icon("question-circle"))),
                          checkboxGroupInput("local", label = "", 
                                             choices = levels(factor(data$`N_terminus`, exclude=c())),
                                             #selected = "Cytoplasmic", inline = TRUE,
                                             selected = levels(factor(data$`N_terminus`, exclude=c())),
                                             inline = TRUE
                                             ),
                          actionLink("selectall","Select All"),
                          hr(style = "border: 1px solid;"),
                          
                          ## ---------------------------------------------------
                          div(style="display: inline-block; padding-left: 3%; width: 175px; vertical-align: top; font-size: 18px", strong('Show only TM proteins?')),
                          div(style="display: inline-block; width: 100px;",
                          switchInput(inputId = "TM", label="Y/N", value = TRUE, onLabel = "Yes", offLabel = "No"),
                          ),
                          div(style="display: inline-block; padding-left: 3%; width:25px", actionButton("tmHelp", "", icon = icon("question-circle"))),
                          hr(style = "border: 1px solid;"),
                          
                          ## ---------------------------------------------------
                          div(style="display: inline-block; width: 175px; vertical-align: top; font-size: 18px", strong('Show only non-canonical cases?')),
                          div(style="display: inline-block; width: 100px; vertical-align: center",
                            switchInput(inputId = "NC", label="Y/N", value = TRUE, onLabel = "Yes", offLabel = "No"),
                          ),
                          div(style="display: inline-block; padding-left: 3%; width:25px; vertical-align: top", actionButton("canonicalHelp", "", icon = icon("question-circle"))),
                          #hr(style = "border: 1px solid;"),
                          
                          ## ---------------------------------------------------
                          #actionButton(inputId = 'submit', label = 'Update', icon = icon("sync"),)
                      )
                    ),
                    column(9,
                      wellPanel(
                      plotlyOutput(outputId = "scatter")
                      ),
                      wellPanel(
                          DT::dataTableOutput("table")
                      )  
                    )
                )
            )
        ),
        tabPanel("Internal TMDs analysis",
                 fluidPage(#titlePanel("PART A"),
                   tags$head(
                     #tags$style(type="text/css", "label{ display: table-cell; text-align: center; vertical-align: middle; } .form-group { display: table-row;}")
                   ),
                   fluidRow(
                     column(3,
                            wellPanel(
                              align = 'center',
                              useShinyjs(),
                              
                              ## ---------------------------------------------------
                              #div(style="display: inline-block;vertical-align:center;",
                              #    h1(id="big-heading", "Input panel"),
                              #    tags$style(HTML("#big-heading{color: black; font-size: 20px; text-align: center; text-decoration: underline;}"))
                              #),
                              
                              ## ---------------------------------------------------
                              div(style="padding-left: 3%; padding-bottom: 5%", strong('Process all proteins?', style="display: inline-block; padding-bottom: 5px; font-size: 18px"),
                                  div(style="display: inline-block; padding-left: 3%", actionButton("processHelp", "", icon = icon("question-circle"))),
                              ),
                              div(style="display: inline-block; width: 100px;",
                                  switchInput(inputId = "all2", label="All", value = TRUE, size = "mini", onLabel = "Yes", offLabel = "No"),
                              ),
                              div(style="display: inline-block; vertical-align:center; padding-left: 3%; width: 25px", strong('OR')),
                              div(style="display: inline-block; vertical-align:top; padding-left: 3%; width: 200px;",
                                  selectizeInput("gene2", "Choose a protein(s):",
                                                 choices = NULL, multiple = TRUE),
                              ),
                              hr(style = "border: 1px solid;"),
                              
                              ## ---------------------------------------------------
                              div(style="padding-left: 3%; padding-bottom: 5%", strong('SignalP cutoffs', style="display: inline-block; vertical-align:center; font-size: 18px"),
                                  div(style="display: inline-block; padding-left: 3%;vertical-align:center", actionButton("signalpHelp", "", icon = icon("question-circle"))),
                              ),
                              div(style="display: inline-block;vertical-align:top;horizontal-align:center; width: 170px;",
                                  sliderInput("SignalP_TM2", label = "TM version", min = 0,
                                              max = 1.0, value = c(0.0, 0.6))
                              ),
                              div(style="display: inline-block;vertical-align:top; width: 25px;"),
                              div(style="display: inline-block;vertical-align:top; width: 170px;",
                                  sliderInput("SignalP_noTM2", label = "no TM version", min = 0, 
                                              max = 1.0, value = c(0.6, 1.0))
                              ),
                              hr(style = "border: 1px solid;"),
                              
                              ## ---------------------------------------------------
                              div(style="display: inline-block; padding-left: 3%; width: 150px; font-size: 18px", strong('# Mutations >')),
                              div(style="display: inline-block; padding-left: 3%; text-align: left; width: 100px", numericInput('Num_MUT2', '', value = 0)),
                              div(style="display: inline-block; padding-left: 3%;", actionButton("mutationHelp", "", icon = icon("question-circle"))),
                              hr(style = "border: 0.5px solid;"),
                              
                              ## ---------------------------------------------------
                              div(style="display: inline-block; padding-left: 3%; width: 200px; font-size: 18px", strong('N-term localization?')),
                              div(style="display: inline-block;", actionButton("ntermHelp", "", icon = icon("question-circle"))),
                              checkboxGroupInput("local2", label = "", 
                                                 choices = levels(factor(data2$`N_terminus`, exclude=c())),
                                                 #selected = "Cytoplasmic", inline = TRUE,
                                                 selected = levels(factor(data2$`N_terminus`, exclude=c())),
                                                 inline = TRUE
                              ),
                              actionLink("selectall2","Select All"),
                              hr(style = "border: 1px solid;"),
                              
                              ## ---------------------------------------------------
                              div(style="display: inline-block; padding-left: 3%; width: 175px; vertical-align: top; font-size: 18px", strong('Show only TM proteins?')),
                              div(style="display: inline-block; width: 100px;",
                                  switchInput(inputId = "TM2", label="Y/N", value = TRUE, onLabel = "Yes", offLabel = "No"),
                              ),
                              div(style="display: inline-block; padding-left: 3%; width:25px", actionButton("tmHelp", "", icon = icon("question-circle"))),
                              hr(style = "border: 1px solid;"),
                              
                              ## ---------------------------------------------------
                              div(style="display: inline-block; width: 175px; vertical-align: top; font-size: 18px", strong('Show only non-canonical cases?')),
                              div(style="display: inline-block; width: 100px; vertical-align: center",
                                  switchInput(inputId = "NC2", label="Y/N", value = TRUE, onLabel = "Yes", offLabel = "No"),
                              ),
                              div(style="display: inline-block; padding-left: 3%; width:25px; vertical-align: top", actionButton("canonicalHelp", "", icon = icon("question-circle"))),
                              #hr(style = "border: 1px solid;"),
                              
                              ## ---------------------------------------------------
                              #actionButton(inputId = 'submit2', label = 'Update', icon = icon("sync"),)
                            )
                     ),
                     column(9,
                            wellPanel(
                              plotlyOutput(outputId = "scatter2", width = '50%')
                            ),
                            wellPanel(
                              DT::dataTableOutput("table2")
                            )  
                     )
                   )
                   )
        ),
        tabPanel("About",
                 fluidPage(theme = shinytheme("cerulean"),
                           fluidRow(
                             column(7,
                                    wellPanel(
                                      align = 'center',
                                      useShinyjs(),
                                      htmlOutput("about"),
                                      actionButton(inputId = 'github',
                                                   label = '',
                                                   icon = icon("github"),
                                                   onclick ="window.open('https://github.com/russelllab/spc', '_blank')",
                                                   style='padding:4px; font-size:250%'
                                                   )
                                    )
                             ),
                             column(5,
                                      imageOutput(outputId = "workflow")
                                    )
                            )
                           ),
                 tags$head(tags$style(HTML(".shiny-output-error-validation{color: red;}")))
        )
    )
)

