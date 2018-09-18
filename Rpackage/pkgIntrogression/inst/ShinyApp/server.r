paramNames <- c("AB0", "Ab0", "aB0", "ab0" ,"bA", "ba", "dA", "da", "r")

plot_pop <- function(BDtestdata) {
  plot(BDtestdata$AB_mean ~ BDtestdata$t, ylim=c(0,100), type='l', col = "red", lwd = 5, ylab = "# Individuals", xlab = "Time", cex.axis = 1.5, cex.lab = 1.5)
  lines(BDtestdata$Ab_mean ~ BDtestdata$t, col = "red", lty = "dashed", lwd = 5)
  lines(BDtestdata$aB_mean ~ BDtestdata$t,  col = "blue", lwd = 5)
  lines(BDtestdata$ab_mean ~ BDtestdata$t, lty = "dashed",col = "blue", lwd = 5)
  legend(3600,100, legend=c("AB","Ab","aB","ab"), col = c("red","red","blue","blue"), lty = c("solid","dashed","solid","dashed"), lwd = 5, cex = 1.5)
}

plot_intro <- function(BDtestdata) {
  plot(BDtestdata$FA_mean ~ BDtestdata$t, ylim=c(0,1), type='l', col = "red", lwd = 5, ylab = "# fraction B-alleles", xlab = "Time", cex.axis = 1.5, cex.lab = 1.5)
  lines(BDtestdata$Fa_mean ~ BDtestdata$t, col = "blue", lwd = 5)
}

RunAllSimulation <- function(pars){
  pkgIntrogression::BDSim(100, 5000, pars, progressbar = FALSE)
}

function(input, output, session) {
    getParams <- function() {
      input[["recalc"]]

      params <- list()
      length(params) <- 10
      params <- lapply(paramNames, function(p) {input[[p]]})
      names(params) <- paramNames
      as.list(params)
  }

  navA <- eventReactive(input$recalc, RunAllSimulation(getParams())) 

  output$plot_popdynamic <- renderPlot({
    suppressWarnings(plot_pop(navA()))
  })

  output$plot_introgression <- renderPlot({
   suppressWarnings(plot_intro(navA()))
  })

  output$write_fixation <- renderText({ 
    paste("The fixation probability is: ", navA()$extinct)
  })

  #output$a_lociPlot <- renderPlotly({
  #allelefrequencies
  #dataallelemean <-tidyr::gather(as.data.frame(navA()[[3]]), locus, value, 2:ncol(as.data.frame(navA()[[3]])))
  #dataallelevar <-tidyr::gather(as.data.frame(navA()[[4]]), locus, value, 2:ncol(as.data.frame(navA()[[4]])))
  #dataallele <- dplyr::bind_cols(dataallelemean,dataallelevar[3])
  #dataallele <- dplyr::mutate(dataallele, col=factor(locus,locus))
  #allelefrequencyplot <- ggplot(dataallele, aes(col, value, ymin=value-sqrt(value1),ymax=value+sqrt(value1))) + 
  #geom_pointrange(aes(frame = Generation), shape=22) + 
  #theme_minimal() + 
  #xlab("Locus") +
  #ylab("Allelefrequency") + 
  #theme(axis.text.x = element_blank())

  #suppressWarnings(ggplotly(allelefrequencyplot))
  # })

}