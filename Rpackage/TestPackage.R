library(pkgIntrogression)
testpars <- list(
    b=1.1,
    dA=1,
    da=0.8,
    nloci=100,
    ninit0=20,
    ninit1=2,
    ngen=50,
    nrep=1000,
    rec=0.5,
    k=20)

ShinyInitializeSimulation(testpars)
ShinyRunSimulation()
ShinyWriteOutputandCleanup()

testdata <- RcppIntrogressionSimulation(testpars,0)
tail(testdata$data$Introgressed1_avg,n=1)

plot(testdata$data$Major0_avg)
plot(testdata$data$Major1_avg)
plot(testdata$data$Popsize_avg)
## shiny stuff
IntrogressionShiny()
library(shiny)
runApp("/home/freek/Documents/GeneticRescue/Rpackage/ShinyApp")

  ggplot(testdata[[2]]) + 
  geom_ribbon(aes(x=Generation, ymin=Popsize_avg-sqrt(Popsize_var), ymax=Popsize_avg+sqrt(Popsize_var)), fill = "grey", alpha = 0.8) +
  geom_ribbon(aes(x=Generation, ymin=Major0_avg-sqrt(Major0_avg), ymax=Major0_avg+sqrt(Major0_avg)), fill = "red", alpha=0.8) +
  geom_ribbon(aes(x=Generation, ymin=Major1_avg-sqrt(Major1_avg), ymax=Major1_avg+sqrt(Major1_avg)), fill = "green", alpha = 0.8) + 
  xlab("Time") + ylab("# individuals") + 
  geom_line(aes(x=Generation, y=Popsize_avg)) +
  geom_line(aes(x=Generation, y=Major0_avg)) +
  geom_line(aes(x=Generation, y=Major1_avg)) + 
  theme_minimal() + 
  theme(text=element_text(size=25))
