
## app.R ##
library(shinydashboard)

renderInputs <- function() {
  wellPanel(
    fluidRow(
        column(4,
        sliderInput("t", "nr of timesteps:", min = 10, max = 100, value = 100,step = 20),
        sliderInput("nrep", "nr of repetitions:", min = 1, max = 100, value = 50, step = 10),
        sliderInput("r", "Recombination rate:", min = 0.0, max = 0.5, value = 0.5, step = 0.05)
        ),
        column(4,
        sliderInput("AB0", "nr of AB individuals:", min = 0, max = 100, value = 0, step = 5),
        sliderInput("Ab0", "nr of Ab individuals:", min = 1, max = 10, value = 1, step = 1),
        sliderInput("aB0", "nr of aB individuals:", min = 1, max = 1000, value = 100, step = 10),
        sliderInput("ab0", "nr of ab individuals:", min = 0, max = 100, value = 0, step = 5)
        ),
        column(4,
        sliderInput("bA", "Birthrate of A", min = 0.0, max = 2.0, value = 1.1, step = 0.05),
        sliderInput("ba", "Birthrate of a:", min = 0.0, max = 2.0, value = 1.1, step = 0.05),
        sliderInput("dA", "Deathrate of A:", min = 0.0, max = 2.0, value = 1.05, step = 0.05),
        sliderInput("da", "Deathrate of a:", min = 0.0, max = 2.0, value = 1.15, step = 0.05)
        )
    ),
    p(actionButton("recalc","Run simulation", icon("random")))
  )
}

ui <- dashboardPage(
  dashboardHeader(title = "Introgression"),
  dashboardSidebar(),
  dashboardBody(
    # Boxes need to be put in a row (or column)
      fluidRow(
          column(4,
          renderInputs()
          )
          ,
          column(4,
          plotOutput("plot_popdynamic")
          ),
          column(4,
          plotOutput("plot_introgression")
          )
        ),
        fluidRow(
          textOutput("write_fixation")
        )
  )
)
