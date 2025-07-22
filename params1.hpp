#ifndef _PARAMS_
#define _PARAMS_
#include "/Users/leopoldo/Desktop/Computing/CSM/MD/definitivo/libraries/pvector.hpp"
#include <fstream>

#define species 2

class simpars
{
  using ntype=double;
public:
  int nx, ny, nz; /* nx*ny*nz particelle */
  double T, P; // temperature and pressure
  int Np,NpL[species] ; // numero di particelle
  long int maxadjstps, eqstps, adjstps, save_mgl_snapshot;
  long int savemeasure, outstps, totsteps; // Nsteps = simulations steps, outstps steps print something on current simulation status
  float DT; // simulation time
  double rho, rc; // density
  int simtype; // simulation type (see below)
  int seed; // -1 means random
  pvector<double, 3> L; // box
  double sigma, epsilon, mass; // Lennard-Jones parameters
  double masses[species],rate[species];
  double dt;
  simpars()
    {
      std::ifstream ri;
      
      simtype = 0; // 0 NTV, 1 NPT
      nx = 8;  // number of particles along each direction
      ny = 8;
      nz = 8;
      sigma=1.0;
      epsilon=1.0;
      rho = 0.5; 
      rc = 2.5;
      seed=0;
      mass=1.0;
      masses[0]=1.;
      masses[1]=100.;
      rate[0]=0.1;
      rate[1]=.9;
      ///
      adjstps = -200;
      maxadjstps = 2000;
      eqstps=-1;
      save_mgl_snapshot = 1000; 
      savemeasure=5; 
      outstps=1000;
      DT=2;
      dt = 0.01;
      totsteps = static_cast<long int>(DT/dt);
      T = 2.0;
      P = 3.838; //se P*=\beta*P*v0, 1 < P* < 10 dove v0 è il volume di una particella
    }
  void set_rate(double r);
};
void simpars::set_rate(double r)
{
  rate[0] = r;
  rate[1] = 1 - r;
}
#endif

