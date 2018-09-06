// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(BH)]]

#ifndef SHINYFUNCTION_H
#define SHINYFUNCTION_H

#include "IBSSimulations.h"
#include "Rcpp_output.h"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <Rcpp.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

static bool RUNSIMULATION_FUN = false, INITIALIZESIMULATION_FUN = false;
static Parameters Shiny_Pars;
static SimData Shiny_Data; 

// [[Rcpp::export]]
void ShinyIBS_init(const int &generations, const Rcpp::List &parslist){
    if(RUNSIMULATION_FUN == true) {std::cerr << "First write output or clear dataset : WriteOutputandCleanupt()" << std::endl;}
    else{
        // Prepare for simulations
        Parameters Shiny_Pars(generations,parslist);
        INITIALIZESIMULATION_FUN = true;
    }
}

// [[Rcpp::export]]
void ShinyIBS_run(){
    if(INITIALIZESIMULATION_FUN == true){
        while(RunSimulation(Shiny_Pars, Shiny_Data)==false);        
        RUNSIMULATION_FUN = true;
    }
    else{std::cerr << "First initialize simulation: InitializeSimulation()" << std::endl;}
}

// [[Rcpp::export]]
Rcpp::List ShinyIBS_write(){
    if(INITIALIZESIMULATION_FUN == true && RUNSIMULATION_FUN == true){
        return Rcpp_WriteOutput(Shiny_Pars,Shiny_Data);
    }
    else{
        std::cerr << "First initialize and run a simulation: InitializeSimulation() and RunSimulation()" << std::endl;
        return NULL;
    }
}


#endif