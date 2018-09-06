library(pkgIntrogression)

#test parameters:
IBStestpars()
BDtestpars()

#run simulations
BDtestdata <- BDSim(100000,1000, BDtestpars(), setthreads = 0)
IBStestdata <- IBSSim(50, IBStestpars(), setthreads = 4, progressbar = TRUE)

BDtestdata

testdata <- RcppIntrogressionSimulation(testpars,0)


tail(testdata$data$Introgressed1_avg,n=1)

## shiny stuff
IntrogressionShiny()
library(shiny)
runApp("/home/freek/Documents/GeneticRescue/Rpackage/ShinyApp")