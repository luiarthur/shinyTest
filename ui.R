shinyUI(fluidPage(
  titlePanel("Attraction IBP"),

  sidebarLayout(position="right",
    sidebarPanel(
      numericInput("distMatNum",label="Distance Matrix",value=1,min=1,max=3),
      uiOutput("distMat"),

      imageOutput("Image",height=700),
      br(),

      numericInput("matNum",label="Matrix #:",value=1,min=1),
      uiOutput("matrix")
    ),
    mainPanel( 
      dataTableOutput("matFreq")
    )
  )
))
