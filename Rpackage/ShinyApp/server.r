paramNames <- c("ngen", "nloci", "nrep", "ninit1" ,"ninit0", "k", "b", "dA", "da", "rec")

plot_pop <- function(data) {
  ggplot(data[[2]]) + 
  geom_ribbon(aes(x=Generation, ymin=Popsize_avg-sqrt(Popsize_var), ymax=Popsize_avg+sqrt(Popsize_var)), fill = "grey", alpha = 0.5) +
  geom_ribbon(aes(x=Generation, ymin=Major0_avg-sqrt(Major0_avg), ymax=Major0_avg+sqrt(Major0_avg)), fill = "red", alpha=0.8) +
  geom_ribbon(aes(x=Generation, ymin=Major1_avg-sqrt(Major1_avg), ymax=Major1_avg+sqrt(Major1_avg)), fill = "green", alpha = 0.8) + 
  xlab("Time") + ylab("# individuals") + 
  geom_line(aes(x=Generation, y=Popsize_avg)) +
  geom_line(aes(x=Generation, y=Major0_avg)) +
  geom_line(aes(x=Generation, y=Major1_avg)) + 
  theme_minimal() + 
  theme(text=element_text(size=25))
}

plot_intro <- function(data) {
  ggplot(data[[2]]) + 
  geom_ribbon(aes(x=Generation, ymin=Introgressed0_avg-sqrt(Introgressed0_var), ymax=Introgressed0_avg+sqrt(Introgressed0_var)), fill = "red", alpha = 0.8) +
  geom_ribbon(aes(x=Generation, ymin=Introgressed1_avg-sqrt(Introgressed1_var), ymax=Introgressed1_avg+sqrt(Introgressed1_var)), fill = "green", alpha = 0.8) +
  geom_line(aes(x=Generation, y=Introgressed0_avg)) +
  geom_line(aes(x=Generation, y=Introgressed1_avg)) +
  theme_minimal() + 
  xlab("Time") + 
  ylab("Genetic rescue") + 
  theme(text=element_text(size=25))
}

RunAllSimulation <- function(pars){
  pkgIntrogression::ShinyInitializeSimulation(pars)
  
  withProgress(message="Running simulations",value = 0,{
    for(x in c(1:pars$nrep)){
      incProgress(1/pars$nrep, detail = paste("Doing replicate", x))
      pkgIntrogression::ShinyRunSimulation()
    }
  })

  pkgIntrogression::ShinyWriteOutputandCleanup()
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
    paste("The fixation probability is: ", navA()$fixation)
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