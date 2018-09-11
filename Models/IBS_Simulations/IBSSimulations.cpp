/*=============================================================================================================
                                                IntrogressionSimulations.cpp
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

#include <boost/dynamic_bitset.hpp>
#include <assert.h>
#include <stdlib.h>
#include <ctime>
#include <iomanip>
#include <chrono>
#include <string>
#include "random.h"
#include "IBSSimulations.h"

enum type {A,a};

class Individual{
  public:
    Individual(const type maketype, const int &ngenes) {
        maketype == A ? genome.set(ngenes, true) : genome.set(ngenes, false);
    }

    Individual(Individual* parent1, Individual* parent2, const int &ngenes){
        
        /*lots of things can be improved here, not my best coding*/

        genome.resize(ngenes);
        boost::dynamic_bitset<> temp1, temp2, temp3, r_localbit, m_localbit;
        r_localbit.resize(NLOCI);
        m_localbit.resize(NLOCI);

        // Recombination
        Make_localbit(r_localbit, pars.r_global);
        temp1 = parent1->genome & r_localbit;
        temp2 = parent2->genome & r_localbit.flip();
        temp3 = temp1 | temp2;
        
        // Mutation
        Make_localbit(m_localbit, pars.m_global);
        temp1 = temp3 & m_localbit;
        temp2 = temp3 | m_localbit;
        genome = temp2 & temp1.flip();
    }

    inline bool rescue(const Parameters &pars){if(genome[pars.index[0]]==1){return true;} else{return false;} ;}

    bool Genotype(const int &locus) {return genome[locus]; }

    int GenotypeCount() {return genome.count();}
    
    private:
    void Make_localbit(boost::dynamic_bitset<> &local, const boost::dynamic_bitset<> &global)
    {
        const int NLOCI = genome.size();
        const int global_size = global.size();

        assert(local.size() == NLOCI);
        assert(global_size > NLOCI);
        int start = rnd::integer(global_size - NLOCI);

        for (int i = 0; i < NLOCI; ++i)
        {
            local[i] = global[start + i];
        }
    };

    boost::dynamic_bitset<> genome; 
};

class Population{
    // constructor
    Population(const Parameters* Simpars) pars(Simpars) {
        
        // Create storage block for population metrics
        data = new Datablock;

        // Create initial population distribution
        for(int i = 0; i < Simpars->NINIT[0]; ++i)
            pop[i] = new Individual(SimPars.INIT_GENOME[0]);
        }
        for(int i = Simpars->NINIT[0]; i < Simpars->NINIT[0]+Simpars->NINIT[1]; ++i){
            pop[i] = new Individual(SimPars.INIT_GENOME[1]);
        }

        // Save initial population
        WriteToData();
    }

    // destructor 
    ~Population(){
        delete data;
    }

    // iterate population
    bool Iterate(){ 
        
        // Choose n offspring
        const int N = pop.size();
        const double r = 1.0 + (pars->BIRTHRATE-1.0) * (1.0 - ((double)N / (double)pars->K));
        assert(r >= 0.0);
        const int noffspring = rnd::poisson((double)N*r);
        std::vector<Individual*> offspring(noffspring);
        std::vector<Individual*> offspring2;

        // Create offspring by random mating of parents
        for(it = offspring.begin(); it != offspring.end(); ++it){
            *it = new Individual(pop[rnd::integer(nparents)],pop[rnd::integer(nparents)],pars);
        }

        // Kill parents
        for(Individual* ind : pop) delete ind;
        pop.clear();

        // Viability selection
        bool rescue = false;
        for(it = offspring.begin(); it != offspring.end(); ++it){
            if((*it)->rescue(pars)) {
                rescue = true;
                if(rnd::uniform() < pars->DEATHRATEA){
                    pop.push_back(*it);
                }
            }
            else{
                if(rnd::uniform() < pars->DEATHRATEa){
                    pop.push_back(*it);
                }
            }
        }

        // Save population
        WriteToData();
        return rescue;
    }

    // return pointer to Datablock
    inline Datablock* ptr2datablock(){
        return data;
    }

    private:
    std::vector<Individual*> pop;
    std::vector<Individual*>::iterator it;
    const Datablock* data;
    const Parameters* pars;
    
    void WriteToData(){
        //Save population size:
        const int populationsize = population.size();
        data->popsize.push_back(populationsize);

        //Save the rest
        std::vector<int> count[2];
        count[0].resize(pars->NLOCI,0);
        count[1].resize(pars->NLOCI,0);
        int temp0 = 0, temp1 = 0, ind0 = 0, ind1 = 0;
        for(Individual* ind : population){
            if(ind->rescue(pars)){
                ++ind1;
                temp1 += ind->GenotypeCount() - 1;
            }
            else{
                ++ind0;
                temp0 += ind->GenotypeCount(); 
            } 
            // Genome allele frequencies
            for(int j = 0; j < pars->NLOCI; ++j){
                ++count[ind->Genotype(j)][j];
            }
        }

        data->allele0.push_back(count[0]);
        data->allele1.push_back(count[1]);
        data->major0.push_back(ind0);
        data->major1.push_back(ind1);
        data->introgressed0.push_back((double)temp0/(double)((pars->NLOCI-1)*ind0));
        data->introgressed1.push_back((double)temp1/(double)((pars->NLOCI-1)*ind1));
    }
}

bool DoSimulation(const Parameters &SimPars, SimData &SimulationData){
    
    // Initialize testpop
    Population testpop(SimPars);
   
    // Iterate testpop (conditioned on A fixating)
    for (int i = 0; i < SimPars.NGEN; ++i){
        if(testpop.Iterate(SimPars)==false){
            delete testpop;
            ++SimulationData.nofixcounter;
            return false;
    }}    
  
    // Save successful testpop
    SimulationData.push_back_protect(testpop.ptr2datablock()); 
    delete testpop;
    return true;
}
