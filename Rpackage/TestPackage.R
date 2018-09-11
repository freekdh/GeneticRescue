#F.J.H. de Haas
#testing R package pkgIntrogression
#Two models: BDSim and IBSSim

library(pkgIntrogression)
library(ggplot2)

#### BDSIM ####

#test parameters:
BDpars <- BDtestpars()
BDpars$aB0 <- 100

#run simulations
BDtestdata <- BDSim(10000,50000, BDpars, setthreads = 0)
#IBStestdata <- IBSSim(50, IBStestpars(), setthreads = 4, progressbar = TRUE)

#Genetic rescue:
plot(BDtestdata$FA_mean ~ BDtestdata$t, ylim=c(0,1), type='l', col = "red", lwd = 5)
lines(BDtestdata$Fa_mean ~ BDtestdata$t, col = "blue", lwd = 5)

#Population dynamics:
plot(BDtestdata$AB_mean ~ BDtestdata$t, ylim=c(0,100), type='l', col = "red", lwd = 5)
lines(BDtestdata$Ab_mean ~ BDtestdata$t, col = "red", lty = "dashed", lwd = 5)
lines(BDtestdata$aB_mean ~ BDtestdata$t,  col = "blue", lwd = 5)
lines(BDtestdata$ab_mean ~ BDtestdata$t, lty = "dashed",col = "blue", lwd = 5)

#Population dynamics:
plot(sqrt(BDtestdata$AB_var) ~ BDtestdata$t, ylim=c(0,100), type='l', col = "red", lwd = 5)
lines(sqrt(BDtestdata$Ab_var) ~ BDtestdata$t, col = "red", lty = "dashed", lwd = 5)
lines(sqrt(BDtestdata$aB_var) ~ BDtestdata$t,  col = "blue", lwd = 5)
lines(sqrt(BDtestdata$ab_var) ~ BDtestdata$t, lty = "dashed",col = "blue", lwd = 5)

BDtestdata$extinct

#### IBSSIM ####

IBStestdata <- IBSSim(generations = 100, parslist = IBStestpars())

View(IBStestdata$data)

## shiny stuff
IntrogressionShiny()
library(shiny)
runApp("/home/freek/Documents/GeneticRescue/Rpackage/ShinyApp")