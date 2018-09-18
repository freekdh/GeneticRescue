
#include <boost/dynamic_bitset.hpp>
#include <assert.h>
#include <stdlib.h>
#include <ctime>
#include <iomanip>
#include <chrono>
#include <string>
#include "random.h"
#include "IBSSimulations.h"

namespace PopGen{

    #ifndef random_h
        rnd::set_seed();
    #endif

    enum type {type_A,type_a};

    template <class _ploidytype>
    class Population{
        public:
        Population(const int n0[2], const Parameters* pop_pars, const Parameters* ind_pars) : Individual_pars(ind_pars), Population_pars(pop_pars) {
            N = n0[type_A]+n0[type_a];
            pop.reserve(N);
            for(int i = 0; i < n0[type_A]; ++i)
            pop.push_back(new _ploidytype(type_A, ind_pars));
            for(int i = 0; i < n0[type_a]; ++i)
            pop.push_back(new _ploidytype(type_a, ind_pars));
        }
        ~Population();
        /*
        Population Population::operator+ (const Population& pop2) const 
        {
            Population out;
            out.genome = this->genome + pop2.genome;
            return out;
        }
        */
        bool Iterate(const int &Selection_locus, const double prob_death[2]){             
            // Create offspring
            const int Nparents = pop.size();
            const int N = rnd::poisson(E_r(Nparents));
            offspring.clear();
            offspring.resize(N);

            // Random mating
            for(_ploidytype* ind : offspring){ ind = new _ploidytype(pop[rnd::integer(Nparents)],pop[rnd::integer(Nparents)]); }
            for(_ploidytype* ind : pop) delete ind;
            pop.clear();

            // Viability selection
            bool is_rescue = false;
            for(_ploidytype* ind : offspring){
                if(ind->GetRescue(Selection_locus)) {
                    is_rescue = true;
                    if(rnd::uniform() < prob_death[type_A]) pop.push_back(ind) ;
                }
                else{
                    if(rnd::uniform() < prob_death[type_a]) pop.push_back(ind) ;
                }
            }

            ++generation;
            return is_rescue;
        }

        private:
        typename std::vector<_ploidytype*>::iterator it;
        std::vector<_ploidytype*> offspring;
        int N;
        int generation = 0;
        inline double E_r(const int &Nparents){
            return (double)Nparents*(1.0 + Population_pars->r*(1.0 - ((double)Nparents / Population_pars->K)));
        }

        protected:
        Parameters* Individual_pars;
        Parameters* Population_pars;
        std::vector<_ploidytype*> pop;
    };

    class Individual{
    public:
        Individual(const Parameters* Simpars) : pars(Simpars){};
        virtual bool GetRescue(const int &locus);
        virtual int GenotypeCount();
        
        protected:
        const Parameters* pars;
        void Make_localbit(boost::dynamic_bitset<> &local, const boost::dynamic_bitset<> &global)
        {
            const int global_size = global.size();

            assert(local.size() == pars->NLOCI);
            assert(global_size > NLOCI);
            int start = rnd::integer(global_size - pars->NLOCI);

            for (int i = 0; i < pars->NLOCI; ++i)
            {
                local[i] = global[start + i];
            }
        };
    };

    class Haploid : public Individual{
        public:
        Haploid(const type &maketype, const Parameters* Simpars) : Individual(Simpars) {
            maketype == type_A ? haplotype.set(true) : haplotype.set(false);
        }
        Haploid(const Haploid* parent1, const Haploid* parent2, const Parameters*Simpars) : Individual(Simpars){
            boost::dynamic_bitset<> temp1, temp2, temp3, r_localbit(pars->NLOCI), m_localbit(pars->NLOCI);

            // Recombination
            Make_localbit(r_localbit, pars->r_global);
            temp1 = parent1->haplotype & r_localbit;
            temp2 = parent2->haplotype & r_localbit.flip();
            temp3 = temp1 | temp2;
            
            // Mutation
            Make_localbit(m_localbit, pars->m_global);
            temp1 = temp3 & m_localbit;
            temp2 = temp3 | m_localbit;
            haplotype = temp2 & temp1.flip();
        }
        int GenotypeCount(){
            return haplotype.count();
        }
        bool GetRescue(const int &locus){
            return haplotype[locus];
        } 

        private: 
        boost::dynamic_bitset<> haplotype; 
    };

    /*
    class Diploid : public Individual{
        
    };
    */
};

int main(){
    using namespace PopGen;
    Parameters* pars;
    const int n0[2] = {5,5};

    Population<Haploid> pop(n0, pars, pars);
    return 0;
}

/*
    class Diploid : public Individual{
        public:
        Polyploid(const Parameters* Simpars, const int &NPLOIDY) : public Individual(Simpars) {} 

        private: 
        boost::dynamic_bitset<> genome; 
    };
*/


/// Code graveyard: 

/*
static void RescueSimulation(Population* pop, const int Ngen,  const Parameters* Simpars, SimData* SimulationData){
    using namespace PopGen;  
    do{
        for (int i = 0; i < Ngen; ++i){
            if(pop->Iterate(Simpars)==false){
                delete pops;
                ++SimulationData.nofixcounter;
                break;
        }}
    }
    while(i < Ngen);

    SimulationData.push_back_protect(pop); 
}

        void SaveCurrentState(DataBlock* data, const int Nloci, const int &Selection_locus){
            data->Nind.push_back(N);
            std::vector<int> count[2];
            int F_A = 0, F_a = 0, N_A = 0, N_a = 0;
            for(_ploidytype* ind : pop){
                if(ind->GetAllele(Selection_locus)){
                    ++N_A;
                    F_A += ind->GenotypeCount() - 1;
                }
                else{
                    ++N_a;
                    F_a += ind->GenotypeCount(); 
                } 
            }
            data->NA.push_back(N_A);
            data->Na.push_back(N_a);
            data->FA.push_back((double)F_A/(double)((Nloci-1)*N_A));
            data->Fa.push_back((double)F_a/(double)((Nloci-1)*N_a));
        }
*/
