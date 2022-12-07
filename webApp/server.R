library(shiny)
library(gapminder)
library(plotly)
library(ggplot2)
library(readr)
library(NGLVieweR)
library(shinyjs)
library(shinythemes)
library(shinyWidgets)

#data <- read_tsv("../data/part_a_data.tsv")
data <- read_tsv("part_a_data.tsv")
new_data <- data
new_data <- data.frame(new_data)
#headers <- c('GENE', 'ACC', 'Total_MUT', 'UniProt_MUT', 'ClinVar_MUT', 'COSMIC_MUT', 'Position', 'SignalP_TM', 'SignalP_noTM', 'N_terminus', 'Topology')
headers <- c('GENE', 'Total_MUT', 'UniProt_MUT', 'ClinVar_MUT', 'COSMIC_MUT', 'Position', 'SignalP_TM', 'SignalP_noTM', 'N_terminus', 'Topology')
new_data <- new_data[headers]
RV <- reactiveValues(data = data.frame(new_data))

#data2 <- read_tsv("../data/part_c_data.tsv")
data2 <- read_tsv("part_c_data.tsv")
new_data2 <- data2
new_data2 <- data.frame(new_data2)
#print (new_data2['N-terminus'])
#headers2 <- c('GENE')
headers2 <- c('GENE', 'ACC', 'Total_MUT', 'UniProt_MUT', 'ClinVar_MUT', 'COSMIC_MUT', 'Position', 'TMD', 'SignalP_TM', 'SignalP_noTM', 'N_terminus')
new_data2 <- new_data2[headers2]
RV2 <- reactiveValues(data = data.frame(new_data2))


spc_plot <- function(x, a, b, c, d, e, f, signal_no, signal_yes) {
    if (x==1) {
      #datar <- read_tsv("../data/part_a_data.tsv")
      datar <- read_tsv("part_a_data.tsv")
    }
    else {
      #datar <- read_tsv("../data/part_c_data.tsv")
      datar <- read_tsv("part_c_data.tsv")
    }
    new_datar <- datar
    a <- levels(factor(a, exclude=c()))
    new_datar <- new_datar %>% filter(`N_terminus` %in% a)
    
    b <- levels(factor(b, exclude = c()))
    if (b == TRUE){
      b <- 'YES'
      new_datar <- new_datar %>% filter(`TM protein` %in% b)
    }
    
    if (c == FALSE) {
      new_datar <- new_datar %>% filter(`GENE` %in% d)
    }
    
    new_datar <- new_datar %>% filter(`Total_MUT` >= e)
    
    
    f <- levels(factor(f, exclude = c()))
    if (f == TRUE){
      f <- 0 #if show only non-canonical signalp proteins, then select 0
      new_datar <- new_datar %>% filter(`Signalp` %in% f)
    }
    
    
    new_datar <- new_datar %>% filter(`SignalP_noTM` >= signal_no[1] & `SignalP_noTM` <= signal_no[2])
    new_datar <- new_datar %>% filter(`SignalP_TM` >= signal_yes[1] & `SignalP_TM` <= signal_yes[2])
    
    p <- new_datar %>%
        #ggplot(aes(text=new_data$GENE, SignalP_noTM, SignalP_TM, size = `Total_MUT`, color=`N-terminus`, shape=`TM protein`)) +
        #geom_point()
        ggplot(aes(text=new_datar$GENE, SignalP_noTM, SignalP_TM, size = `Total_MUT`, color=`N_terminus`)) +
        geom_point()
    
    #result <- list('gg' = ggplotly(p, height =400, width = 1400), 'values' = new_data)
    #return (result)
    return(ggplotly(p, height =400, width = 1300))
    
}

spc_table <- function(x, a, b, c, d, e, f, signal_no, signal_yes) {
  if (x==1) {
    #datar <- read_tsv("../data/part_a_data.tsv")
    datar <- read_tsv("part_a_data.tsv")
  }
  else {
    #datar <- read_tsv("../data/part_c_data.tsv")
    datar <- read_tsv("part_c_data.tsv")
  }
  #data <- read_tsv("../data/part_a_data.tsv")
  new_datar <- datar
  a <- levels(factor(a, exclude=c()))
  new_datar <- new_datar %>% filter(`N_terminus` %in% a)

  b <- levels(factor(b, exclude = c()))
  if (b == TRUE){
    b <- 'YES'
    new_datar <- new_datar %>% filter(`TM protein` %in% b)
  }

  if (c == FALSE) {
    new_datar <- new_datar %>% filter(`GENE` %in% d)
  }
  
  new_datar <- new_datar %>% filter(`Total_MUT` >= e)
  #print (new_datar)
  #print (e)
  
  f <- levels(factor(f, exclude = c()))
  if (f == TRUE){
    f <- 0 #if show only non-canonical signalp proteins, then select 0
    new_datar <- new_datar %>% filter(`Signalp` %in% f)
  }
  
  new_datar <- new_datar %>% filter(`SignalP_noTM` >= signal_no[1] & `SignalP_noTM` <= signal_no[2])
  new_datar <- new_datar %>% filter(`SignalP_TM` >= signal_yes[1] & `SignalP_TM` <= signal_yes[2])

  if (x==1) {
    return(data.frame(new_datar)[headers])
  }
  else {
    return(data.frame(new_datar)[headers2])
  }
  #return (reactiveValues(data = data.frame(new_data)))
}
# Define server logic required to draw a histogram
shinyServer(
  function(input, output, session) {
    ### Help modals
    observeEvent(input$processHelp, {
      showModal(modalDialog(
        title = "Process all proteins?",
        HTML("<div class='text-justify' style='font-size:16px'><p>Choose to process either <strong>All</strong> the proteins (Yes) or only selected proteins (No) using the button on the left. 
             In case of latter, on the right an input field will be enabled where you can enter one or more proteins (gene names) of your choice.</p></div>"
            )
      ))
    })
    
    observeEvent(input$signalpHelp, {
      showModal(modalDialog(
        title = "SignalP Y-scores",
        HTML("<div class='text-justify' style='font-size:16px'><p>Select the ranges for Y-scores of the two version of <a href='https://services.healthtech.dtu.dk/service.php?SignalP-4.1'>SignalP</a>: <kbd>TM</kbd> and <kbd>no TM</kbd>.</p>
             <p><kbd>TM version</kbd>: Y-scores derived by running the SignalP program in <i><strong>transmembrane network</strong></i> mode</p>
             <p><kbd>no TM version</kbd>: Y-scores derived by running the SignalP program in <i><strong>no-transmembrane network</strong></i> mode</p>
             </div>"
        )
      ))
    })
    
    observeEvent(input$mutationHelp, {
      showModal(modalDialog(
        title = "Number of mutations",
        HTML("<div class='text-justify' style='font-size:16px'><p>Select proteins that harbor minimum number of total mutations i.e.
        unique number of mutations found in <a href='https://pubmed.ncbi.nlm.nih.gov/14681372/'>UniProt</a>,
        <a href='https://pubmed.ncbi.nlm.nih.gov/29165669/'>ClinVar</a>, and <a href='https://pubmed.ncbi.nlm.nih.gov/30371878/'>COSMIC</a>. 
             </p></div>"
        )
      ))
    })
    
    observeEvent(input$ntermHelp, {
      showModal(modalDialog(
        title = "N-terminus localization",
        HTML("<div class='text-justify' style='font-size:16px'><p>
              Select the cell compartments in which the protein's N-terminus must lie.
              </p>
             </div>"
        )
      ))
    })
    
    observeEvent(input$tmHelp, {
      showModal(modalDialog(
        title = "Transmembrane proteins",
        HTML("<div class='text-justify' style='font-size:16px'><p>Choose to display either only the transmembrane proteins (<strong>Yes</strong>) or
            otherwise (<strong>No</strong>). 
             </p></div>"
        )
      ))
    })
    
    observeEvent(input$canonicalHelp, {
      showModal(modalDialog(
        title = "Help",
        HTML("<div class='text-justify' style='font-size:16px'><p>Choose to display either only the proteins with non-canonical signal peptides (<strong>Yes</strong>) or
            otherwise (<strong>No</strong>). 
             </p></div>"
        )
      ))
    })
    
    ### Help modals
    observeEvent(input$processHelp2, {
      showModal(modalDialog(
        title = "Process all proteins?",
        HTML("<div class='text-justify' style='font-size:16px'><p>Choose to process either <strong>All</strong> the proteins (Yes) or only selected proteins (No) using the button on the left. 
             In case of latter, on the right an input field will be enabled where you can enter one or more proteins (gene names) of your choice.</p></div>"
        )
      ))
    })
    
    observeEvent(input$signalpHelp2, {
      showModal(modalDialog(
        title = "SignalP Y-scores",
        HTML("<div class='text-justify' style='font-size:16px'><p>Select the ranges for Y-scores of the two version of <a href='https://services.healthtech.dtu.dk/service.php?SignalP-4.1'>SignalP</a>: <kbd>TM</kbd> and <kbd>no TM</kbd>.</p>
             <p><kbd>TM version</kbd>: Y-scores derived by running the SignalP program in <i><strong>transmembrane network</strong></i> mode</p>
             <p><kbd>no TM version</kbd>: Y-scores derived by running the SignalP program in <i><strong>no-transmembrane network</strong></i> mode</p>
             </div>"
        )
      ))
    })
    
    observeEvent(input$mutationHelp2, {
      showModal(modalDialog(
        title = "Number of mutations",
        HTML("<div class='text-justify' style='font-size:16px'><p>Select proteins that harbor minimum number of total mutations i.e.
        unique number of mutations found in <a href='https://pubmed.ncbi.nlm.nih.gov/14681372/'>UniProt</a>,
        <a href='https://pubmed.ncbi.nlm.nih.gov/29165669/'>ClinVar</a>, and <a href='https://pubmed.ncbi.nlm.nih.gov/30371878/'>COSMIC</a>. 
             </p></div>"
        )
      ))
    })
    
    observeEvent(input$ntermHelp2, {
      showModal(modalDialog(
        title = "N-terminus localization",
        HTML("<div class='text-justify' style='font-size:16px'><p>
              Select the cell compartments in which the protein's N-terminus must lie.
              </p>
             </div>"
        )
      ))
    })
    
    observeEvent(input$tmHelp2, {
      showModal(modalDialog(
        title = "Transmembrane proteins",
        HTML("<div class='text-justify' style='font-size:16px'><p>Choose to display either only the transmembrane proteins (<strong>Yes</strong>) or
            otherwise (<strong>No</strong>). 
             </p></div>"
        )
      ))
    })
    
    observeEvent(input$canonicalHelp2, {
      showModal(modalDialog(
        title = "Non-canonical cases",
        HTML("<div class='text-justify' style='font-size:16px'><p>Choose to display either only the proteins with non-canonical signal peptides (<strong>Yes</strong>) or
            otherwise (<strong>No</strong>). 
             </p></div>"
        )
      ))
    })
    
    getPage<-function() {
      return(includeHTML("about.html"))
    }
    output$about<-renderUI({getPage()})
    
    output$workflow <- renderImage({
      list(src = "workflow.png", contentType = 'image/png',
           alt = "Workflow")
    }, deleteFile = FALSE)
    
    observeEvent(input$selectall, {
      if(input$selectall == 0) return(NULL) 
      else if (input$selectall%%2 == 0)
      {
        updateCheckboxGroupInput(session,"local", inline = TRUE, choices=levels(factor(data$`N_terminus`, exclude=c())))
        updateActionLink(session, 'selectall', label = 'Select All')
      }
      else
      {
        updateCheckboxGroupInput(session,"local", inline = TRUE, choices=levels(factor(data$`N_terminus`, exclude=c())), selected=levels(factor(data$`N_terminus`, exclude=c())))
        updateActionLink(session, 'selectall', label = 'Select None')
      }
    })
    
    observeEvent(input$selectall2, {
      if(input$selectall2 == 0) return(NULL) 
      else if (input$selectall2%%2 == 0)
      {
        updateCheckboxGroupInput(session,"local2", inline = TRUE, choices=levels(factor(data2$`N_terminus`, exclude=c())))
        updateActionLink(session, 'selectall2', label = 'Select All')
      }
      else
      {
        print(levels(factor(data$`N-terminus`, exclude=c())))
        updateCheckboxGroupInput(session,"local2", inline = TRUE, choices=levels(factor(data2$`N_terminus`, exclude=c())), selected=levels(factor(data2$`N_terminus`, exclude=c())))
        updateActionLink(session, 'selectall2', label = 'Select None')
      }
    })
  
    observeEvent(input$all, {
        if (input$all == TRUE){
            shinyjs::disable("gene")
        }
        else {
            shinyjs::enable("gene")
        }
        })
    
    observeEvent(input$all2, {
      if (input$all2 == TRUE){
        shinyjs::disable("gene2")
      }
      else {
        shinyjs::enable("gene2")
      }
    })

    updateSelectizeInput(session, 'gene', choices = levels(factor(new_data$`GENE`)), server = TRUE)
    updateSelectizeInput(session, 'gene2', choices = levels(factor(new_data2$`GENE`)), server = TRUE)
    

    observeEvent(c(input$all, input$SignalP_TM, input$SignalP_noTM, input$gene, input$Num_MUT, input$local, input$TM, input$NC), {RV$data <- spc_table(1, input$local, input$TM, input$all, input$gene, input$Num_MUT, input$NC, input$SignalP_noTM, input$SignalP_TM)})
    observeEvent(c(input$all, input$SignalP_TM, input$SignalP_noTM, input$gene, input$Num_MUT, input$local, input$TM, input$NC), {output$scatter <- renderPlotly({spc_plot(1, input$local, input$TM, input$all, input$gene, input$Num_MUT, input$NC, input$SignalP_noTM, input$SignalP_TM)})})
    
    observeEvent(c(input$all2, input$SignalP_TM2, input$SignalP_noTM2, input$gene2, input$Num_MUT2, input$local2, input$TM2, input$NC2), {RV2$data <- spc_table(2, input$local2, input$TM2, input$all2, input$gene2, input$Num_MUT2, input$NC2, input$SignalP_noTM2, input$SignalP_TM2)})
    observeEvent(c(input$all2, input$SignalP_TM2, input$SignalP_noTM2, input$gene2, input$Num_MUT2, input$local2, input$TM2, input$NC2), {output$scatter2 <- renderPlotly({spc_plot(2, input$local2, input$TM2, input$all2, input$gene2, input$Num_MUT2, input$NC2, input$SignalP_noTM2, input$SignalP_TM2)})})
    
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
    
    output$table2 <- DT::renderDataTable({
      DT::datatable(RV2$data,
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
}
)
