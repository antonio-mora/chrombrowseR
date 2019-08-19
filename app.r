library(shiny)
library(GenomicRanges)
library(R.utils)    	# includes gunzip
library(rtracklayer)	# includes liftOver
source("get_rao.R")
source("format_rao.R")

ui <- fluidPage(
  titlePanel("ChromBrowser"),
  sidebarLayout(
    sidebarPanel(
      helpText("Retrieving and formatting the chromatin interaction datasets from Rao et al., 2014"),
      selectInput("cell", 
                  label = "Cell type",
                  choices = c("GM12878", "HMEC", "HUVEC", "HeLa", "IMR90", "K562", "KBM7", "NHEK", "CH12-LX"),
                  selected = "GM12878")
    ),
    mainPanel(
      textOutput("selected_cell"),
	  htmlOutput("text"),
	  tableOutput("headTable"),
	  tableOutput("headTable2")
    )
  )
)

server <- function(input, output) {
	# Get Rao datasets from GEO:
	raoList <- reactive({
		get_rao(input$cell)
	})
	
	# Format Rao dataset ():
	raoLoops <- reactive({
		format_rao(raoList(), "loops", input$cell)
	})
	raoDomains <- reactive({
		format_rao(raoList(), "domains", input$cell)
	})

	# Output:
	output$text <- renderUI({
		str1 <- paste("You have selected: ", input$cell)
		str2 <- paste("Raw interaction data has been downloaded to: ", getwd())
		str3 <- paste("Formatted interaction data has been downloaded to: ", getwd())
		str4 <- paste("Loop dataset size: ", dim(raoLoops())[1])
		str5 <- paste("Domain dataset size: ", dim(raoDomains())[1])
		str6 <- "Tables (head only): "

		HTML(paste(str1, str2, str3, str4, str5, str6, sep = '<br/>'))
	})
    output$headTable <- renderTable({ 
		#head(raoList()$loops)
		head(raoLoops())
    })
    output$headTable2 <- renderTable({ 
		#head(raoList()$domains)
		head(raoDomains())
    })
}

shinyApp(ui = ui, server = server)
