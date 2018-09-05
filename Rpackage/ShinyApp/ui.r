library(ggplot2)

renderInputs <- function() {
  wellPanel(
    fluidRow(
        column(3,
        sliderInput("ngen", "nr of generations:", min = 10, max = 100, value = 100,step = 20),
        sliderInput("nloci", "nr of loci:", min = 50, max = 200, value = 100, step = 25),
        sliderInput("nrep", "nr of repetitions:", min = 1, max = 100, value = 50, step = 10)
        ),
        column(3,
        sliderInput("ninit0", "nr of wildtype (a) individuals:", min = 10, max = 100, value = 20, step = 10),
        sliderInput("ninit1", "nr of rescuetype (A) individuals:", min = 1, max = 100, value = 2)
        ),
        column(3,
        sliderInput("k", "Carrying capacity of the population:", min = 10, max = 150, value = 20, step = 20),
        sliderInput("b", "Birthrate of individuals:", min = 0.0, max = 2.0, value = 1.1, step = 0.1)
        ),
        column(3,
        sliderInput("dA", "Deathrate A:", min = 0.0, max = 2.0, value = 1.0, step = 0.05),
        sliderInput("da", "Deathrate a:", min = 0.0, max = 2.0, value = 1.1, step = 0.01),
        sliderInput("rec", "Recombination rate:", min = 0.0, max = 0.5, value = 0.5, step = 0.05)
        )
    ),
    p(actionButton("recalc","Run simulation", icon("random")))
  )
}

fluidPage(theme="simplex.min.css",
  tags$style(type="text/css",
    "label {font-size: 12px;}",
    ".recalculating {opacity: 1.0;}"
  ),

  # Application title
  tags$h2("Genetic rescue due to introgressive hybridization"),
  p("Simulations implemented in c++"),
  hr(),

  fluidRow(
    column(6,
    renderInputs()
    )
    ,
    column(3,
    plotOutput("plot_popdynamic")
    ),
    column(3,
    plotOutput("plot_introgression")
    )
  ),
  fluidRow(
    textOutput("write_fixation")
  )
)