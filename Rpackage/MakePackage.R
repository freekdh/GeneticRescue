library(Rcpp)

Rcpp.package.skeleton(
    name = "pkgIntrogression", 
    code_files = c("Rscripts/ShinyExport.R"), 
    cpp_files = c(
    "./../Models/IBS_Simulations/IntrogressionSimulations.cpp",
    "./../Models/IBS_Simulations/IntrogressionSimulations.h",
    "./../Models/IBS_Simulations/Rcpp_output.cpp",
    "./../Models/IBS_Simulations/Rcpp_output.h",
    "./../Models/IBS_Simulations/Shiny_output.cpp", 
    "./../Models/IBS_Simulations/random.h", 
    "./../Models/IBS_Simulations/random.cpp", 
    "./../Models/IBS_Simulations/utils.h", 
    "./../Models/IBS_Simulations/utils.cpp",
    "./../Models/IBS_Simulations/Makevars"),
    example_code = FALSE,
    author = "F.J.H. de Haas",
    email = "dehaas@zoology.ubc.ca")
    
#Add to the DESCRIPTION file:
#Imports: Rcpp (>= 0.12.17), BH, RcppProgress
#LinkingTo: Rcpp, BH, RcppProgress