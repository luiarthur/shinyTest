shinyServer(function(input,output) {
  source("../ddibp.R",chdir=T)
  source("oneSim.R")
  library(xtable)
  
  output$distMat = renderUI({
    D <- eval(parse(text=paste0("D",input$distMatNum)))
    D <- print(xtable(D, align=rep("", ncol(D)+1), digits=0),
               floating=FALSE, tabular.environment="array", comment=FALSE, 
               print.results=FALSE,hline.after=NULL,
               include.rowname=FALSE,include.colname=FALSE)
    html <- paste0("$$", D, "$$")
    withMathJax(helpText(html))
  })
  #output$Image = renderImage({
  #  list(src=paste0("www/",D.name,".png"))
  #}, deleteFile=FALSE)

  #output$matFreq = renderDataTable({matfreq})

  #output$matrix = renderUI({
  #  M <- um[[input$matNum]]
  #  if (ncol(M)>0) {
  #    M <- print(xtable(M, align=rep("", ncol(M)+1)),
  #               floating=FALSE, tabular.environment="bmatrix", comment=FALSE, 
  #               print.results=FALSE,hline.after=NULL,
  #               include.rowname=FALSE,include.colname=FALSE)
  #    html <- paste0("$$", M, "$$")
  #    withMathJax(helpText(html))
  #  }
  #})
})
