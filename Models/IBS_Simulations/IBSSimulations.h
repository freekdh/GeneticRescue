#ifndef GLOBALSTRUCT_H
#define GLOBALSTRUCT_H

#include <mutex>
#include <atomic>
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <Rcpp.h>
#include "random.h"

enum type {A,a};

struct Parameters{
    Parameters(){};
    Parameters(int argc, char *argv[]);    
    Parameters(const int &, const Rcpp::List &);

    void Initialize();
    
    double RECOMBINATIONRATE;
    double MUTATIONRATE = 0.0;
    int NGEN;
    int NLOCI;
    int NREP;
    int NINIT[2];
    int K;
    double BIRTHRATE;
    double DEATHRATEA;
    double DEATHRATEa;

    //Initialize:
    std::vector<int> index;
    boost::dynamic_bitset<> INIT_GENOME[2];
    boost::dynamic_bitset<> r_global, m_global;
};

struct DataBlock{
    public:
        std::vector<int> Nind;
        std::vector<int> NA;
        std::vector<int> Na;
        std::vector<double> FA;
        std::vector<double> Fa;
};
/*
class Individual{
        public:
        Individual(const type IndividualType);
        Individual(Individual* parent1, Individual* parent2);
        inline bool GetAllele(const int &locus);
        inline int GenotypeCount();
};
*/
struct SimData{  
    public:
    SimData(){nofixcounter = 0;}
    void push_back_protect(DataBlock* &datablock){
        mu_datablock.lock();
        DataSet.push_back(datablock);
        mu_datablock.unlock();
    };

    private:
    std::mutex mu_datablock;
    std::vector<DataBlock*> DataSet;
    std::atomic<int> nofixcounter;
};

bool DoSimulation(const Parameters &SimPars, SimData &SimulationData);


Parameters::Parameters(int argc, char *argv[]){
    MUTATIONRATE = 0.0;
    BIRTHRATE = std::atof(argv[1]);
    DEATHRATEA = std::atof(argv[2]);
    DEATHRATEa = std::atof(argv[3]);
    NLOCI = std::atoi(argv[4]);
    NINIT[0] = std::atoi(argv[5]);
    NINIT[1] = std::atoi(argv[6]);
    NGEN = std::atoi(argv[7]);
    NREP = std::atoi(argv[8]);
    RECOMBINATIONRATE = std::atof(argv[9]);
    K = std::atoi(argv[10]);

    Initialize();
}

Parameters::Parameters(const int &generations, const Rcpp::List &parslist){
    MUTATIONRATE = 0.0;
    BIRTHRATE = parslist["b"];
    DEATHRATEA = parslist["dA"];
    DEATHRATEa = parslist["da"];
    NLOCI = parslist["nloci"];
    NINIT[0] = parslist["ninit0"];
    NINIT[1] = parslist["ninit1"];
    NGEN = generations;
    NREP = parslist["nrep"];
    RECOMBINATIONRATE = parslist["rec"];
    K = parslist["k"];

    Initialize();
}

void Parameters::Initialize(){
    
    //Genetic architecture of selection
    index.clear();
    index.resize(NLOCI);
    index[0] = std::floor((double)NLOCI/2.0);

    std::vector<int> index;
    boost::dynamic_bitset<> INITGENOME0(NLOCI);
    boost::dynamic_bitset<> INITGENOME1(NLOCI);
    INITGENOME1.set();
    INIT_GENOME[0] = INITGENOME0;
    INIT_GENOME[1] = INITGENOME1;

    //Recombination and mutation
    const int GLOBALMAX = 100000;
    r_global.resize(GLOBALMAX);
    boost::dynamic_bitset<> a(1);
    a[0] = rnd::bernoulli(0.5);
    for (int i = 0; i < GLOBALMAX; ++i)
    {
        if (rnd::bernoulli(RECOMBINATIONRATE) == true)
        a.flip();
        r_global[i] = a[0];
    }
    m_global.resize(GLOBALMAX);
    for (int i = 0; i < GLOBALMAX; ++i)
    {
        m_global[i] = rnd::bernoulli(MUTATIONRATE);
    }
}

#endif
