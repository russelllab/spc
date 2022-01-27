library(shiny)
library(gapminder)
library(plotly)
library(ggplot2)
library(readr)
library(NGLVieweR)
library(shinyjs)
library(shinythemes)
library(shinyWidgets)

data <- read_tsv("../data/data.tsv")
new_data <- data
new_data <- data.frame(new_data)
headers <- c('GENE', 'ACC', 'Total_MUT', 'UniProt.MUT', 'ClinVar.MUT', 'COSMIC.MUT', 'Position', 'SignalP_TM', 'SignalP_noTM', 'N.terminus', 'Topology')
new_data <- new_data[headers]
RV <- reactiveValues(data = data.frame(new_data))

# Define UI for application that draws a histogram
ui <- navbarPage("SPC as a Quality Control Enzyme",
        tabPanel("PART A",
            fluidPage(#titlePanel("PART A"),
                      #tags$footer(title = "Your footer here"),
                      #theme = shinytheme("cerulean"),
                      tags$head(
                        #tags$style(type="text/css", "label{ display: table-cell; text-align: center; vertical-align: middle; } .form-group { display: table-row;}")
                      ),
                  fluidRow(
                    column(2,
                      wellPanel(
                          align = 'center',
                          useShinyjs(),
                          
                          ## ---------------------------------------------------
                          div(style="display: inline-block;vertical-align:center;",
                              h1(id="big-heading", "Input panel"),
                              tags$style(HTML("#big-heading{color: black; font-size: 20px; text-align: center; text-decoration: underline;}"))
                          ),
                          
                          ## ---------------------------------------------------
                          div(style="display: inline-block; padding-left: 3%;", actionButton("inputhelp", "", icon = icon("question-circle"))),
                          hr(style = "border: 1px solid;"),
                          div(style="display: inline-block; padding-left: 3%;", strong('Process all proteins?', style="display: inline-block; padding-bottom: 5px; vertical-align:top;")),
                          switchInput(inputId = "val", label="YES", value = TRUE),
                          div(style="display: inline-block; padding-left: 3%;", strong('OR', style="display: inline-block; padding-bottom: 10px; vertical-align:top;")),
                          selectizeInput("gene", "Choose a protein(s):",
                                         choices = NULL, multiple = TRUE),
                          hr(style = "border: 1px solid;"),
                          
                          ## ---------------------------------------------------
                          sliderInput("SignalP_TM", label = "SignalP_TM", min = 0, 
                                      max = 1.0, value = c(0.4, 1.0)),
                          sliderInput("SignalP_noTM", label = "SignalP_noTM", min = 0, 
                                      max = 1.0, value = c(0.4, 1.0)),
                          hr(style = "border: 1px solid;"),
                          
                          ## ---------------------------------------------------
                          div(style="display: inline-block; padding-left: 3%; text-align: left;", numericInput('Num_MUT', 'Show proteins # mutations >', value = 0)),
                          hr(style = "border: 1px solid;"),
                          
                          ## ---------------------------------------------------
                          checkboxGroupInput("local", label = "N-term localization?", 
                                             choices = levels(factor(data$`N-terminus`, exclude=c())),
                                             selected = "Cytoplasmic", inline = TRUE),
                          actionLink("selectall","Select All"),
                          hr(style = "border: 1px solid;"),
                          
                          ## ---------------------------------------------------
                          div(style="display: inline-block; padding-left: 3%;", strong('Show only TM proteins?')),
                          switchInput(inputId = "val1", label="YES", value = TRUE),
                          hr(style = "border: 1px solid;"),
                          
                          ## ---------------------------------------------------
                          div(style="display: inline-block; padding-left: 3%;", strong('Show only non-canonical cases?')),
                          switchInput(inputId = "val2", label="YES", value = TRUE),
                          hr(style = "border: 1px solid;"),
                          
                          ## ---------------------------------------------------
                          actionButton(inputId = 'submit', label = 'Update', icon = icon("sync"),)
                      )
                    ),
                    column(10,
                      wellPanel(
                      plotlyOutput(outputId = "scatter", width = '50%')
                      ),
                      wellPanel(
                        DT::dataTableOutput("table")
                      )  
                    )
                )
            )
        ),
        tabPanel("PART C",
                 fluidPage(theme = shinytheme("cerulean")),
                 tags$head(tags$style(HTML(".shiny-output-error-validation{color: red;}"))),
                 mainPanel(tabsetPanel(
                     id = 'dataset2',
                     tabPanel("PART A", DT::dataTableOutput("table1")),
                     tabPanel("PART C", DT::dataTableOutput("table2"))
                        )
                 )
        ),
        tabPanel("ABOUT",
                 fluidPage(theme = shinytheme("cerulean")),
                 tags$head(tags$style(HTML(".shiny-output-error-validation{color: red;}")))
        ),
        tabPanel(icon = icon("github"),
                 fluidPage(theme = shinytheme("cerulean")),
                 tags$head(tags$style(HTML(".shiny-output-error-validation{color: red;}")))
        )
    )

spc_plot <- function(a, b, c, d, e, f, signal_no, signal_yes) {
    data <- read_tsv("../data/data.tsv")
    new_data <- data
    a <- levels(factor(a, exclude=c()))
    new_data <- new_data %>% filter(`N-terminus` %in% a)
    
    b <- levels(factor(b, exclude = c()))
    if (b == TRUE){
      b <- 'YES'
      new_data <- new_data %>% filter(`TM protein` %in% b)
    }
    
    if (c == FALSE) {
      new_data <- new_data %>% filter(`GENE` %in% d)
    }
    
    new_data <- new_data %>% filter(`Total_MUT` >= e)
    
    f <- levels(factor(f, exclude = c()))
    if (f == TRUE){
      f <- 0 #if show only non-canonical signalp proteins, then select 0
      new_data <- new_data %>% filter(`Signalp` %in% f)
    }
    
    new_data <- new_data %>% filter(`SignalP_noTM` >= signal_no[1] & `SignalP_noTM` <= signal_no[2])
    new_data <- new_data %>% filter(`SignalP_TM` >= signal_yes[1] & `SignalP_TM` <= signal_yes[2])
    #print (class(signal_no[1]))
    #print (signal_no[2])
    #print (new_data)
    
    p <- new_data %>%
        #ggplot(aes(text=new_data$GENE, SignalP_noTM, SignalP_TM, size = `Total_MUT`, color=`N-terminus`, shape=`TM protein`)) +
        #geom_point()
        ggplot(aes(text=new_data$GENE, SignalP_noTM, SignalP_TM, size = `Total_MUT`, color=`N-terminus`)) +
        geom_point()
    
    #result <- list('gg' = ggplotly(p, height =400, width = 1400), 'values' = new_data)
    #return (result)
    return(ggplotly(p, height =400, width = 1400))
    
}

spc_table <- function(a, b, c, d, e, f, signal_no, signal_yes) {
  data <- read_tsv("../data/data.tsv")
  new_data <- data
  a <- levels(factor(a, exclude=c()))
  new_data <- new_data %>% filter(`N-terminus` %in% a)

  b <- levels(factor(b, exclude = c()))
  if (b == TRUE){
    b <- 'YES'
    new_data <- new_data %>% filter(`TM protein` %in% b)
  }

  if (c == FALSE) {
    new_data <- new_data %>% filter(`GENE` %in% d)
  }
  
  new_data <- new_data %>% filter(`Total_MUT` >= e)
  
  f <- levels(factor(f, exclude = c()))
  if (f == TRUE){
    f <- 0 #if show only non-canonical signalp proteins, then select 0
    new_data <- new_data %>% filter(`Signalp` %in% f)
  }
  
  new_data <- new_data %>% filter(`SignalP_noTM` >= signal_no[1] & `SignalP_noTM` <= signal_no[2])
  new_data <- new_data %>% filter(`SignalP_TM` >= signal_yes[1] & `SignalP_TM` <= signal_yes[2])
  #print (class(signal_no[1]))
  #print (signal_no[2])
  #print (new_data)

  #result <- list('gg' = ggplotly(p, height =400, width = 1400), 'values' = new_data)
  #return (result)
  return(data.frame(new_data)[headers])
  #return (reactiveValues(data = data.frame(new_data)))
}
# Define server logic required to draw a histogram
server <- function(input, output, session) {
    observeEvent(input$selectall, {
      print('hello')
      if(input$selectall == 0) return(NULL) 
      else if (input$selectall%%2 == 0)
      {
        updateCheckboxGroupInput(session,"local", inline = TRUE, choices=levels(factor(data$`N-terminus`, exclude=c())))
        updateActionLink(session, 'selectall', label = 'Select All')
        print('hello here')
      }
      else
      {
        print(levels(factor(data$`N-terminus`, exclude=c())))
        updateCheckboxGroupInput(session,"local", inline = TRUE, choices=levels(factor(data$`N-terminus`, exclude=c())), selected=levels(factor(data$`N-terminus`, exclude=c())))
        updateActionLink(session, 'selectall', label = 'Select None')
        print('hello')
      }
    })
  
    observeEvent(input$val, {
        if (input$val == TRUE){
            shinyjs::disable("gene")
        }
        else {
            shinyjs::enable("gene")
        }
        })

    updateSelectizeInput(session, 'gene', choices = levels(factor(new_data$`GENE`)), server = TRUE)
    spc <- eventReactive(input$submit, {spc_plot(input$local, input$val1, input$val, input$gene, input$Num_MUT, input$val2, input$SignalP_noTM, input$SignalP_TM)})
    #result <- eventReactive(input$submit, {spc(input$local, input$tmp, input$val, input$gene, input$Num_MUT, input$SignalP_noTM, input$SignalP_TM)})
    #result
    #data <- result[2]
    #gg <- result$gg
    #spc <- ggplotly(p, height =400, width = 1400)
    output$scatter <- renderPlotly({spc()})
    #output$scatter <- renderPlotly({data()})
    observeEvent(input$submit, {RV$data <- spc_table(input$local, input$val1, input$val, input$gene, input$Num_MUT, input$val2, input$SignalP_noTM, input$SignalP_TM)})
    #data <- read_tsv("../data/data.tsv")
    #new_data <- data
    #new_data <- new_data %>% filter(`TM protein` %in% 'YES')
    output$table <- DT::renderDataTable({
      DT::datatable(RV$data,
                    filter = "top",
                    height = 1,
                    class = 'cell-border strip hover',
                    extensions = list("ColReorder" = NULL,
                                      "Buttons" = NULL,
                                      "FixedColumns" = list(leftColumns=1)),
                    options = list(
                      dom = 'lBRrftpi',
                      autoWidth=TRUE,
                      pageLength = 3,
                      lengthMenu = list(c(3, 5, 10, 15, 50, -1), c('3', '5', '10', '15','50', 'All')),
                      ColReorder = TRUE,
                      buttons =
                        list(
                          'copy',
                          'print',
                          list(
                            extend = 'collection',
                            buttons = c('csv', 'excel', 'pdf'),
                            text = 'Download'
                          )
                        )
                    )
      )
    }, server = TRUE)
    output$x5 = renderPrint({
      cat('\n\nSelected rows:\n\n')
      cat(input$table_rows_selected, sep = ', ')
    })
    
    ##observeEvent(input$table_cell_clicked, {
    #  info = input$table_cell_clicked
    #  cat('\n\nSelected rows:\n\n')
    #  cat(info, sep = ', ')
      # do nothing if not clicked yet, or the clicked cell is not in the 1st column
      #if (is.null(info$value) || info$col != 0) return()
      #updateTabsetPanel(session, 'x0', selected = 'Plot')
      #updateTextInput(session, 'x2', value = info$value)
    #})
}

# Run the application 
shinyApp(ui = ui, server = server)
