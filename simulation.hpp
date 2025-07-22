#ifndef _SIMCLASS_
#define _SIMCLASS_

#include <vector>
#include <string>
#include <fstream>
#include "./params1.hpp"
#include "./particle.hpp"
#include "./libraries/randnumgen.hpp"
#include <iomanip> // for setprecision()
#include <chrono>

using ntype=double;

template<typename particle_type>
class sim
{
  using simp = simpars;
protected:
  simp params;  // parametri per composizione
  std::vector<particle_type> p;
  ///
  ntype dt; 
  ntype calcenergyi(int i, int opt=0);// if opt=1 calculate energies only for i < j
  ntype totenergy(void);
  void pbc(int i);
  void save_mgl_snapshot(long int t);

public:
  void prepare_initial_conf(void);

  void set_sim_type(int type){params.simtype = type;}
  void set_dt(ntype timestep) { dt = timestep; } // New method to set dt
  void set(double dt, double r)
  {
    this->dt = dt;
    params.dt = dt;
    params.totsteps = static_cast<long int>(params.DT/dt);
    params.set_rate(r);
  }

  void init_rng(int n) 
    {
      if (n < 0)
        rng.rseed();
      else
        rng.seed(n);
    }

  void run(void) {std::cout << "Running simulation\n";}; 
};

// MD class
template <typename particle_type>
class mdsim : public sim<particle_type>
{
  using ntype = double;
  ntype gauss(void); // Generate gaussian number by Box Muller algorithm
protected:
  using bc = sim<particle_type>;
  using bc::params, bc::p, bc::pbc, bc::dt;
  ntype Us, U, E0;
  void calc_forces(ntype &Us, ntype &U);
  ntype calcK(void);
  void init_measures(void);
  void save_measures(long int t);

public:
  void prepare_initial_conf(void);
  void algo(void);
  void print_step_results(int t);
  void print_final_results(void);
  void simulation_result(std::chrono::duration<double> d);
  void run(void);
};

// Multi-step class
template <typename particle_type>
class multistep : public mdsim<particle_type>
{
  using ntype = double;
  using bc = mdsim<particle_type>;
  using bc::params, bc::pbc, bc::p, bc::Us, bc::U, bc::E0, bc::dt;
  int M,N1,N2;
  ntype eps = 1.; // epsilon è usato per variare M, di default è 1

  void calc_forces(int flag);  

public:
  void set_eps(ntype eps){this->eps=eps;}
  int get_M(void) {return M;}

  int get_N1(void) {return N1;}
  int get_N2(void) {return N2;}
  void prepare_initial_conf(void);
  void sort_particles(void);  // Metodo per ordinare le particelle in base alla massa all'interno del vettore p
  void algo(double dt, int Ni, int Nf, int flag);
  void simulation_result(std::chrono::duration<double> d);
  void run(void);
};

//################################
// METODI CLASSE MULTISTEP
template <typename particle_type>
void multistep<particle_type>::prepare_initial_conf(void)
{
  bc::prepare_initial_conf();
  sort_particles();
  M = int(sqrt(params.masses[1] / params.masses[0]) * eps); // Definisco M
  if(M%2==1) M++;
  N1 = params.NpL[0];
  N2 = params.NpL[1];
}

template <typename particle_type>
void multistep<particle_type>::sort_particles(void)
{
  std::vector<particle_type> pvec = p;
  int c1 = 0, c2 = params.Np - 1;
  for (int i = 0; i < params.Np; i++)
  {
    if (pvec[i].m == params.masses[0])
      p[c1++] = pvec[i]; // c1 viene prima usato e poi incrementato
    else
      p[c2--] = pvec[i];
  }
}

template <typename particle_type>
void multistep<particle_type>::calc_forces(int flag)// flag=1 per il primo gruppo di particelle, flag=2 per il secondo
{
  int starti, stopi, startj1, stopj1, stopj2; // starti e stopi sono gli indici di partenza e fine per il ciclo sulle particelle

  if (flag == 1)
  {
    starti = 0, stopi = N1, startj1 = N1, stopj1 = params.Np, stopj2 = N1;
  }
  if (flag == 2)
  {
    starti = N1, stopi = params.Np, startj1 = 0, stopj1 = N1, stopj2 = params.Np;
  }


  ntype vij, vijs, wij;
  for (int i = starti; i < stopi; i++){p[i].f = {0.0, 0.0, 0.0};}
  for (int i = starti; i < stopi; i++)
  {
    for (int j = startj1; j < stopj1; j++) // Calcolo forze tra particelle di gruppi diversi
    {
      pvector fijv = p[i].fij(p[j], params.L, vij, vijs, wij);
      p[i].f = p[i].f + fijv;
    }
    for (int j = i + 1; j < stopj2; j++) // Calcolo forze tra particelle dello stesso gruppo
    {
      pvector fijv = p[i].fij(p[j], params.L, vij, vijs, wij);
      p[i].f = p[i].f + fijv;
      p[j].f = p[j].f - fijv;
    }
  }
}

// Velocity Verlet algorithm for multistep ( evolve solo le particelle del gruppo flag)
template <typename particle_type>
void multistep<particle_type>::algo(double dt, int Ni, int Nf, int flag) 
{
  for (int i = Ni; i < Nf; i++) 
  {
    p[i].expiLp(dt / 2);
    p[i].expiLq(dt);
    pbc(i);
  }
  calc_forces(flag);
  for (int i = Ni; i < Nf; i++)
  {
    p[i].expiLp(dt / 2);
  }
}

template <typename particle_type>
void multistep<particle_type>::run(void)
{
  bc::calc_forces(Us, U);
  E0 = bc::calcK() + Us;
  std::fstream fo;
  bc::init_measures();

  auto start = std::chrono::high_resolution_clock::now();
  for (long int t = 0; t < params.totsteps ; t++)
  {
    for (int i = 0; i < M/2; i++) {algo(dt/M, 0, N1, 1);} // M/2 steps for the first group
    calc_forces(2);
    algo(dt, N1, params.Np, 2); // One steps for the second group
    calc_forces(1);
    for (int i = 0; i < M/2; i++) {algo(dt/M, 0, N1, 1);} // M/2 steps for the first group
    bc::calc_forces(Us, U); // Usare calc_forces(1) se non c'é bisogno di calcolare l'energia step per step
    bc::print_step_results(t);
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end - start;

  //bc::print_final_results();
  simulation_result(duration);
}

// Metodo per salvare i risultati della simulazione
template <typename particle_type>
void multistep<particle_type>::simulation_result(std::chrono::duration<double> d)
{
  std::fstream f;
  f.open("./results/MS.dat", std::ios::out | std::ios::app);
  f << dt << " "
    << dt / M << " "
    << E0 << " "
    << M << " "
    << N1 << " "
    << N2 << " "
    << d.count()<<"\n";
  f.sync();
  f.close();
}

//################################
//################################
//SIM class methods
template<typename particle_type>
ntype sim<particle_type>::calcenergyi(int i, int opt)
{
  int j;
  ntype enei = 0.0;
  for (j = 0; j < params.Np; j++) // pars.Np è il numero totale di particelle
  {
    if (opt == 1 && i >= j)
      continue;
    if (i == j)
      continue;
    enei += p[i].vij(p[j], params.L);
    // pars.L è un vettore con i lati del box
  }
  return enei;
}

template <typename particle_type>
ntype sim<particle_type>::totenergy(void)
{
  ntype ene = 0.0;
  for (auto i = 0; i < p.size(); i++)
  {
    ene += calcenergyi(i, 1);
  }
  return ene;
}

template <typename particle_type>
void sim<particle_type>::pbc(int i)
{
  auto Dr = p[i].r;
  Dr = params.L.mulcw(rint(Dr.divcw(params.L))); // L*rint(Dr/L)
  p[i].r = p[i].r - Dr;
}

template <typename particle_type>
void sim<particle_type>::save_mgl_snapshot(long int t)
{
  std::fstream f;
  f.open("snapshot.dat", std::ios::out|std::ios::trunc);
  f << "# snapshot at time prova" << t << "\n";
  for (int i=0; i < params.Np; i++)
    {
      f << p[i].r(0) << " " << p[i].r(1) << " " << p[i].r(2) << " " << p[i].m <<"\n";
    }
  f << "e\n";
  f.sync();
  f.close();
}

template <typename particle_type>
void sim<particle_type>::prepare_initial_conf(void)
{
  // SC
  int ix, iy, iz;
  int cc = 0,count=0;
  params.Np = params.nx * params.ny * params.nz;
  p.resize(params.Np);
  ntype vcell = pow(params.sigma, 3.0);
  ntype rhomax = 1.0 / vcell;
  ntype sf;
  sf = cbrt(rhomax / params.rho);
  params.L = {ntype(params.nx), ntype(params.ny), ntype(params.nz)};
  ntype clen = sf * params.sigma; //dimensione della cella unitaria a=1/(rho)^1/3 = clen
  params.L *= clen; //quindi 1D: L=a*n dove n è il numero di particelle
  for (ix = 0; ix < params.nx; ix++)
    for (iy = 0; iy < params.ny; iy++)
      for (iz = 0; iz < params.nz; iz++)
      {
        p[cc].r = {ix * clen, iy * clen, iz * clen};
        p[cc].r -= params.L * 0.5; //per centrarla in (-L/2,L/2)
        // Settare la massa secondo il rate impostato
        if (rng.ranf() < params.rate[0]) 
        {
          p[cc].set_mass(params.masses[0]);
          count++;
        }
        else p[cc].set_mass(params.masses[1]);
        cc++;
      }
      // Calcolo del numero di particelle per ogni gruppo
      params.NpL[0]=count;
      params.NpL[1]=params.Np-count;
}

//################################
// MD METHODS

template <typename particle_type>
void mdsim<particle_type>::init_measures(void)
{
  std::fstream f;
  f.open("totenergy.dat", std::ios::out | std::ios::trunc);
  f.close();
}

template <typename particle_type>
ntype mdsim<particle_type>::calcK(void)
{
  ntype K = 0.0;
  for (int i = 0; i < params.Np; i++){K += 0.5 * p[i].m * p[i].v * p[i].v;}
  return K;
}

template <typename particle_type>
void mdsim<particle_type>::save_measures(long int t)
{
  std::fstream f;
  f.open("totenergy.dat", std::ios::out | std::ios::app);
  ntype K = calcK();
  // K+Us is the total conserved energy, where Us is the shifted potential energy
  f << t << " " << std::setprecision(15) << K + Us << " ";
  f << std::setprecision(8)<< K << " " << Us << " ";
  f << std::setprecision(15) << K + Us - E0 << "\n";
  f.sync();
  f.close();
}

template <typename particle_type>
void mdsim<particle_type>::calc_forces(ntype& Us, ntype& U)
{
  U = 0.0;
  Us = 0.0;
  ntype vij, vijs, wij;
  for (int i = 0; i < params.Np; i++){p[i].f = {0.0, 0.0, 0.0};}
  for (int i = 0; i < params.Np; i++)
    {
      for (int j = i + 1; j < params.Np; j++)
      {
        pvector fijv = p[i].fij(p[j], params.L, vij, vijs, wij);
        p[i].f = p[i].f + fijv;
        p[j].f = p[j].f - fijv;
        U += vij;
        Us += vijs;
      }
    }
}

template <typename particle_type>
ntype mdsim<particle_type>::gauss(void)
{
  // gaussian of variance 1 and 0 mean
  double x1, x2;
  do
  {
    x1 = rng.ranf();
    x2 = rng.ranf();
  } while (x1 == 0 || x2 == 0);
  return cos(2 * M_PI * x2) * sqrt(-2.0 * log(x1));
}

template <typename particle_type>
void mdsim<particle_type>::prepare_initial_conf(void)
{
  bc::prepare_initial_conf(); // place particles on a SC lattice

  pvector<ntype, 3> vcm = {0.0, 0.0, 0.0};
  ntype M = 0.0;
  // set initial velocities
  for (int i = 0; i < params.Np; i++)
  {
    ntype sigma = sqrt(params.T / p[i].m);
    p[i].v = {gauss()*sigma, gauss()*sigma, gauss()*sigma};
    p[i].set_vcut();

    vcm = vcm + p[i].v*p[i].m;
    M += p[i].m;    
  }
  // Remove center of mass velocity
  vcm = vcm / M;
  for (int i = 0; i < params.Np; i++) p[i].v = p[i].v - vcm; 
}

// Velocity Verlet algorithm
template <typename particle_type>
void mdsim<particle_type>::algo()
{
  for (int i = 0; i < params.Np; i++)
  {
    p[i].expiLp(dt / 2);
    p[i].expiLq(dt);
    pbc(i);
  }
  calc_forces(Us, U);
  for (int i = 0; i < params.Np; i++)
  {
    p[i].expiLp(dt / 2);
  }
}

template <typename particle_type>
inline void mdsim<particle_type>::print_step_results(int t)
{
  if (params.outstps > 0 && t % params.outstps == 0)
  {
    ntype K = calcK();
  }
  if (t > params.eqstps && params.savemeasure > 0 && (t % params.savemeasure == 0))
    save_measures(t);

  if (t > 0 && params.save_mgl_snapshot > 0 && t % params.save_mgl_snapshot == 0)
  {
    bc::save_mgl_snapshot(t);
  }
}

template <typename particle_type>
void mdsim<particle_type>::print_final_results(void)
{
  std::cout << "#################\n";
  std::cout << "dt\t E \t\t\t sigma(E)\t totalsteps  " << "\n";
  std::cout << dt << "\t" << std::setprecision(15) << calcK() + Us <<"\t"<< abs(calcK() + Us - E0) << "\t" << params.totsteps << "\n";
}

template <typename particle_type>
void mdsim<particle_type>::run(void)
{
  calc_forces(Us, U);
  E0 = calcK() + Us;
  std::fstream fo;
  init_measures();
  auto start = std::chrono::high_resolution_clock::now();
  for (long int t = 0; t < params.totsteps; t++)
  {
    algo();
    print_step_results(t);
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end - start;

  //print_final_results();
  simulation_result(duration);
}

template <typename particle_type>
void mdsim<particle_type>::simulation_result(std::chrono::duration<double> d)
{
  std::fstream f;
  f.open("./results/MD.dat", std::ios::out | std::ios::app);
  f << dt << " "
    << E0 << " "
    << d.count() << "\n";
  f.sync();
  f.close();
}

#endif