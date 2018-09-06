/*=============================================================================================================
                                                BDSimulations.cpp
===============================================================================================================

 Simulations of genetic rescue

 C++-code accompanying:

		(ms. in prep).

 Written by:
        F.J.H. de Haas
       	Theoretical Biology Group
        University of British Columbia
        the Netherlands

 Program version
		xx/xx/xxxx	:

=============================================================================================================*/

#include <assert.h>
#include <stdlib.h>
#include <iomanip>
#include <string>
#include <iostream>
#include "random.h"
#include "utils.h"
#include <progress.hpp>
#include <mutex>
#include <Rcpp.h>
#include <atomic>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#ifdef _OPENMP
    #include <omp.h>
#endif

enum typenames {AB,Ab,aB,ab};
using namespace boost::accumulators;

struct Parameters{
    Parameters(const int &in_AB0, const int &in_Ab0,const int &in_aB0, const int &in_ab0, const double &in_bA, const double &in_ba, double const& in_dA, double const& in_da, double const&in_r){
        AB0 = in_AB0;
        Ab0 = in_Ab0;
        aB0 = in_aB0;
        ab0 = in_ab0;
        bA = in_bA;
        ba = in_ba;
        dA = in_dA;
        da = in_da;
        r = in_r;
    }

    Parameters(Rcpp::List parslist){
        AB0 = parslist["AB0"];
        Ab0 = parslist["Ab0"];
        aB0 = parslist["aB0"];
        ab0 = parslist["ab0"];
        bA = parslist["bA"];
        ba = parslist["ba"];
        dA = parslist["dA"];
        da = parslist["da"];
        r = parslist["r"];
    }

    int AB0,Ab0,aB0,ab0;
    double bA,ba,dA,da;
    double r;
};

class BDPopulation{
    public:
    BDPopulation(const Parameters &pars) {
        type[AB].push_back(pars.AB0);
        type[Ab].push_back(pars.Ab0);
        type[aB].push_back(pars.aB0);
        type[ab].push_back(pars.ab0); 
    }
    bool Iterate(const Parameters &pars){
        /*parents become offspring*/
        type[AB].push_back(type[AB].back());
        type[Ab].push_back(type[Ab].back());
        type[aB].push_back(type[aB].back());
        type[ab].push_back(type[ab].back());
        /*offspring changes*/ 
        static double f[4];
        assert(pars.da>=0.0);
        assert(pars.dA>=0.0);
        if(CalculateFrequencies(f)==false) {extinct = true; return false;};
        const double sumdeathrate = (f[AB]+f[Ab])*pars.dA+(f[aB]+f[ab])*pars.da;
        const double sumbirthrate = (f[AB]+f[Ab])*pars.bA+(f[aB]+f[ab])*pars.ba;
        const double Pdeath = sumdeathrate/(sumbirthrate+sumdeathrate);
        const double Pbirth = 1-Pdeath;
        const double rD = pars.r*(f[AB]*f[ab]-f[Ab]*f[aB]);
        const double PAB_death = f[AB]*pars.dA/sumdeathrate;
        const double PAb_death = f[Ab]*pars.dA/sumdeathrate;
        const double PaB_death = f[aB]*pars.da/sumdeathrate;
        const double Pab_death = f[ab]*pars.da/sumdeathrate;
        const double PAB_birth = f[AB]-rD;
        const double PAb_birth = f[Ab]+rD;
        const double PaB_birth = f[aB]+rD;
        const double Pab_birth = f[ab]-rD;

        rnd::discrete_distribution birthdeathevent(8);
        birthdeathevent[birth_AB] = Pbirth * PAB_birth;
        birthdeathevent[birth_Ab] = Pbirth * PAb_birth;
        birthdeathevent[birth_aB] = Pbirth * PaB_birth;
        birthdeathevent[birth_ab] = Pbirth * Pab_birth;
        
        birthdeathevent[death_AB] = Pdeath * PAB_death;
        birthdeathevent[death_Ab] = Pdeath * PAb_death;
        birthdeathevent[death_aB] = Pdeath * PaB_death;
        birthdeathevent[death_ab] = Pdeath * Pab_death;

        const int sample = birthdeathevent.sample();
        switch(sample) {
            case birth_AB:      ++type[AB].back(); break;
            case birth_Ab:      ++type[Ab].back(); break;
            case birth_aB:      ++type[aB].back(); break;
            case birth_ab:      ++type[ab].back(); break;
            case death_AB:      --type[AB].back(); break;
            case death_Ab:      --type[Ab].back(); break;
            case death_aB:      --type[aB].back(); break;
            case death_ab:      --type[ab].back(); break;
        }
        return true;
    }
    inline int returnTypes(int index){
        return type[index].back();
    }
    bool getextinction(){return extinct;}

    private:
    bool extinct = false;
    enum distribution {birth_AB,birth_Ab,birth_aB,birth_ab,death_AB,death_Ab,death_aB,death_ab};        
    std::vector<int> type[4];
    bool CalculateFrequencies(double *f){
        const double sum=(double)(type[AB].back()+type[Ab].back()+type[aB].back()+type[ab].back());
        if(sum<=0){return false;}
        f[AB] = (double)type[AB].back()/sum;
        f[Ab] = (double)type[Ab].back()/sum;
        f[aB] = (double)type[aB].back()/sum;
        f[ab] = (double)type[ab].back()/sum;
        return true;
    }
};

struct RcppOutput{
    RcppOutput(const unsigned int &ngen){nofixcounter = 0;};
    void pushback_protect(BDPopulation* pop){
        if(pop->getextinction()){
            ++nofixcounter;
        }
        else{
            mu_acc.lock();
            ABvec(pop->returnTypes(AB));
            Abvec(pop->returnTypes(Ab));
            aBvec(pop->returnTypes(aB));
            abvec(pop->returnTypes(ab));
            mu_acc.unlock();
        }
    };

    Rcpp::List pushout(){
        return Rcpp::List::create(
            Rcpp::_["AB"] = mean(ABvec),
            Rcpp::_["Ab"] = mean(Abvec),
            Rcpp::_["aB"] = mean(aBvec),
            Rcpp::_["ab"] = mean(abvec),
            Rcpp::_["extinct"] = static_cast<int>(nofixcounter)
    );
    }

    private:
    std::mutex mu_acc;
    std::atomic<int> nofixcounter;
    accumulator_set<int, stats<tag::mean, tag::variance > > ABvec, Abvec, aBvec, abvec;
};

// [[Rcpp::export]]
Rcpp::List BDSim(const int &nrep, const int &tend, Rcpp::List parslist, int setthreads = 0, bool progressbar = true){
    rnd::set_seed();
    Parameters pars(parslist);
    #ifdef _OPENMP
        const static int maxthreads = omp_get_max_threads();
        if(setthreads>0) omp_set_num_threads(setthreads);
        else omp_set_num_threads(maxthreads);
        REprintf("Parallel activated : Number of threads=%i\n",omp_get_max_threads());   
    #endif
    Progress p(nrep, progressbar);

    /*collect data */    
    RcppOutput OUTPUT(tend);
    BDPopulation* arrayPopulation[nrep];

    #pragma omp parallel for
    for(int j = 0; j < nrep; ++j){
        /* run simulations in parallel*/
        /*#pragma omp task*/
        {
            arrayPopulation[j] = new BDPopulation(pars);
            for(int i = 0; i < tend; ++i){arrayPopulation[j]->Iterate(pars);}
        }
        /*#pragma omp taskwait*/
        /*#pragma omp task*/
        {
            OUTPUT.pushback_protect(arrayPopulation[j]);
            delete arrayPopulation[j];
            p.increment();
        }
    }
    
    return OUTPUT.pushout();
}