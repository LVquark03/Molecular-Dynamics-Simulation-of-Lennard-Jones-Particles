#ifndef _RANDNUMGEN_
#define _RANDNUMGEN_

#include<chrono>
#include<random>

// use by default mersenne-twister rng engine
template <class ntype=double, class engine=std::mt19937_64>

class randnumgen {
  std::random_device rd;
  using myclock = std::chrono::high_resolution_clock;
  engine rng;
  using ud = std::uniform_real_distribution<ntype>;
  ud* unidst; // puntatore a uniform_real_distribution, così da poter cambiare il range di generazione di numeri casuali
public:
  void seed(int s); // seeding fisso
  void rseed(void); // seeding casuale

  ntype ranf() {return (*unidst)(rng); }// generazione di numeri casuali in [0.0,1.0)

  randnumgen()
  {
      unidst = new ud(0.0, 1.0);
      // equivalent way to set the range of random numbers:
      // unidst->param(std::uniform_real_distribution<double>::param_type(0.0, 1.0));
  }
  ~randnumgen(){delete unidst;}
};



template <class ntype, class engine>
void randnumgen<ntype, engine>::seed(int s)
{
    std::seed_seq seq{s + 1, s + 2, s + 2}; // passi una sequenza di valori di seed al motore
    rng.seed(seq);// il motore utilizza questi valori per impostare il suo stato interno.
}

template <class ntype, class engine>
void randnumgen<ntype, engine>::rseed(void)
{
    unsigned int t = myclock::now().time_since_epoch().count();
    std::seed_seq seq{rd() ^ t, rd() ^ t, rd() ^ t}; // ^ is bitwise XOR
    rng.seed(seq);
}

randnumgen rng; // inizializzazione di un oggetto di tipo randnumgen

#endif