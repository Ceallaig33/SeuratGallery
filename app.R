library(shiny)
library(Seurat)
library(dplyr)
library(Matrix)

# demoset downloaded from https://www.dropbox.com/s/63gnlw45jf7cje8/pbmc3k_final.rds?dl=0
# Full details at https://satijalab.org/seurat/vignettes.html


ui <- fluidPage(
  titlePanel("SeuratGallery"),
  fluidRow(
    column(3, 
           wellPanel(
             fileInput(inputId="file_name", label="Select Saved Seurat Object File"),
             radioButtons(inputId = "labelBoolean", label="Display labels?", choices = c("TRUE","FALSE"), selected="TRUE"),
             radioButtons(inputId = "legendBoolean", label="Hide legend?", choices = c("TRUE","FALSE"), selected="FALSE"),
             tags$hr(),
             sliderInput(inputId = "dotSize", label = "Set point size", value=1, min=0.01, max=10),
             sliderInput(inputId = "labelSize", label = "Set label size", value=4, min=0.5, max=10)
           ),
           wellPanel(
             textInput(inputId = "geneName", label="Enter Gene Name", value="ACTB"),
             actionButton("goButton", "Plot Gene")
           )
           
    ),
    column(9,
           verbatimTextOutput("console1"),
           plotOutput("UMAP"),
           tags$hr(),
           tabsetPanel(
             tabPanel("UMAP", plotOutput("featureUMAP")),
             tabPanel("Violins", plotOutput("vlnPlot"))
           )
           
    )
  )
)



server <- function(input, output) {
  
  # increase the max upload file size
  options(shiny.maxRequestSize=1000*1024^2)
  
  # set demo dataset file here
  demo = "pbmc3k_final.rds"
  
  SeuratObject <- reactive({
    if(is.null(input$file_name)){
      # Load a default demo dataset
      return(readRDS(demo))
    } else {
      # Load user-defined dataset
      return(readRDS(input$file_name$datapath))
    }
  })
  
  
  output$console1 <- renderPrint({
    sce <- SeuratObject()
    sce
    
  })
  
  # Display UMAP Plot as overview toggled options for display legends/labels and setting size. Default displays demo dataset
  output$UMAP <- renderPlot({
    DimPlot(SeuratObject(), reduction = 'umap', label = input$labelBoolean,
            pt.size = input$dotSize, label.size = input$labelSize, no.legend=input$legendBoolean)
  })
  
  fPlot <- eventReactive(input$goButton, {
    cols_use = viridis(n=100,option = "A", direction = -1)
    FeaturePlot(object = SeuratObject(),reduction = 'umap', features=input$geneName,
                pt.size = input$dotSize, cols = cols_use)
  })
  
  vPlot <- eventReactive(input$goButton, {
    VlnPlot(object = SeuratObject(),features =  input$geneName, size.x.use = NULL)
  })
  
  output$featureUMAP <- renderPlot(fPlot())
  output$vlnPlot <- renderPlot(vPlot())
  
}

shinyApp(ui = ui, server = server)

