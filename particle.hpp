#ifndef _PARTICLE_
#define _PARTICLE_

#include "./libraries/pvector.hpp"
#include "./params1.hpp"
using ntype = double;
class particle
{
  simpars pars;
protected:
  ntype vcut;
public:
  ntype sigma, epsilon, rc, m;
  pvector<ntype,3> r, v, f; // particle's position and velocity

  void set_vcut(void){vcut = 4.0*epsilon*(pow(sigma/rc,12.0)-pow(sigma/rc,6));};
  void set_sigma(ntype sig){sigma = sig;}
  void set_epsilon(ntype eps){epsilon = eps;}
  void set_rcut(ntype rcut){rc = rcut;}
  void set_mass(ntype mass){m = mass;}

  particle()
    {
      sigma=pars.sigma;
      epsilon=pars.epsilon;
      m=pars.mass;
      rc=pars.rc;
      set_vcut(); // oppure metto a zero vcut = 0.0;
    }

  // methods for MD
  void expiLp(ntype dt){v = v + f*dt/m;}// implementation of action of operator exp(iLp*dt)
  void expiLq(ntype dt){r = r + v*dt;}// implementation of action of operator exp(iLq*dt)
};

class particleLJ: public particle
{
public:
  ntype vij(particleLJ P, pvector<ntype,3> L);
  pvector<ntype,3> fij(particleLJ P, pvector<ntype,3> L, ntype &vij, ntype& vijs, ntype &wij);
};

ntype particleLJ::vij(particleLJ P, pvector<ntype, 3> L)
{
  ntype ene;
  pvector<ntype, 3> Dr;

  Dr = r - P.r;
  // MINIMUM IMAGE CONVENTION
  Dr = Dr - L.mulcw(rint(Dr.divcw(L))); // Dr - L*rint(Dr/L)
  ntype rsq, rn = Dr.norm();
  rsq = rn * rn;
  if (rsq < rc * rc) // interaction potential cut-off
    ene = 4.0 * epsilon * (pow(sigma / rn, 12.0) - pow(sigma / rn, 6));
  else
    ene = 0.0;
  return ene;
}

pvector<ntype, 3> particleLJ::fij(particleLJ P, pvector<ntype, 3> L, ntype &vij, ntype &vijs, ntype &wij)
{
  ntype a,b;
  pvector<ntype, 3> fijv,Dr;
  Dr = r - P.r;
  // MINIMUM IMAGE CONVENTION
  Dr = Dr - L.mulcw(rint(Dr.divcw(L))); // Dr - L*rint(Dr/L)
  ntype rsq, rn = Dr.norm();
  rsq = rn * rn;
  if (rsq < rc * rc)
  {
    b = pow(sigma / rn, 6.0);
    a = b*b;
    fijv = Dr *24.0 * epsilon * (2.0 * a - b) / (rn*rn);
    vij = 4.0 * epsilon * (a - b);
    vijs = vij - vcut;
    wij = fijv * Dr;    
  }
  else
  {
    fijv = {0.0, 0.0, 0.0};
    vij = 0.0;
    vijs = 0.0;
    wij = 0.0;
  }

  return fijv;
}

#endif 
