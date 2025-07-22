#ifndef _PVECTOR_
#define _PVECTOR_

#include<initializer_list> // inizializzazione vettori
#include<iostream> // input/output
#include<cmath>
#include<string> // strings
#include "./randnumgen.hpp"
template <typename ntype=double, int NT=3> // template class

class pvector
{
  ntype v[NT]; 
public:
  pvector(); // constructor
  pvector(std::initializer_list<ntype> list); 
  ~pvector(){};

  // get,set,show
  ntype get(int i) const{return v[i];}
  ntype operator()(int idx) const{return v[idx];}//v(i)
  ntype &operator()(int idx){return v[idx];} //v(i)
  ntype set(int i, ntype val){return v[i] = val;}
  void show(std::string s = "") const;

  //functions
  pvector sum(const pvector& v2) const {return (*this)+v2;}
  ntype norm(void) const {return sqrt((*this)*(*this));} //per gli int non va bene
  pvector rint();
  pvector mulcw(const pvector& v2) const;
  pvector divcw(const pvector& v2) const;

  // operator overloading
  pvector operator+(const pvector& v2) const;
  pvector operator-(const pvector& v2) const;
  pvector operator*(ntype x) const;
  pvector operator/(ntype x) const;  
  pvector operator+=(const pvector& v2) {return (*this)=(*this)+v2;} //+=
  pvector operator-=(const pvector &v2) { return (*this)=(*this)-v2; } //-=
  pvector operator*=(ntype x) {return (*this)=(*this)*x;} // *=
  pvector operator/=(ntype x) {return (*this)=(*this)/x;} // /=
  ntype operator*(const pvector& v2) const; // scalar product
  bool operator==(const pvector& v2) const; // confronto
  pvector operator^(const pvector& v2) const; // cross product

  // random generator
  pvector& random(const ntype& L);
  void random_orient(void);

  // friend functions
  template <typename ntype1, int NT1>
  friend std::ostream& operator<<(std::ostream& os, const pvector<ntype1,NT1>& vec);

  friend pvector operator*(ntype x, const pvector &vec) {return vec*x;}// scalar times vector
  };
//########################################
//########################################
template <typename ntype, int NT>
pvector<ntype,NT>::pvector()
{
  int i;
  for(i=0; i < NT; i++){v[i]=0;}
}

template <typename ntype, int NT>
pvector<ntype,NT>::pvector(std::initializer_list<ntype> list)
{
  int c=0;
  for (ntype el: list)
    {
      if (c < NT){v[c] = el;}
      c++;
    }
  for (;c < NT; c++){v[c]=0.0;}
}

template <typename ntype, int NT>
void pvector<ntype,NT>::show(std::string s) const
{
  std::cout << s << "(";
  for (int i=0; i < NT; i++)
    {
      std::cout << v[i];
      if (i < NT-1) 
        std::cout << ",";
    }
  std::cout << ")\n";
}

template <typename ntype, int NT>
pvector<ntype,NT> pvector<ntype,NT>::rint()
{
  pvector<ntype,NT> vt;
  for (auto i=0; i < NT; i++)
    {
      vt.v[i] = rint(v[i]);
    }
  return vt;
}

template <typename ntype, int NT>
pvector<ntype, NT> pvector<ntype, NT>::mulcw(const pvector &v2) const
{
  pvector<ntype, NT> vs;
  for (int i = 0; i < NT; i++){vs.v[i] = v[i] * v2.v[i];}
  return vs;
}

template <typename ntype, int NT>
pvector<ntype, NT> pvector<ntype, NT>::divcw(const pvector &v2) const
{
  pvector<ntype, NT> vs;
  for (int i = 0; i < NT; i++){vs.v[i] = v[i] / v2.v[i];}
  return vs;
}

// operator overloading
template <typename ntype, int NT>
pvector<ntype,NT> pvector<ntype,NT>::operator+(const pvector& v2) const
{
  pvector vs;
  for (int i=0; i < NT; i++){vs.v[i] = v[i]+v2.v[i];}
  return vs;
}

template <typename ntype, int NT>
pvector<ntype,NT> pvector<ntype,NT>::operator-(const pvector& v2) const
{
  pvector vs;
  for (int i=0; i < NT; i++){vs.v[i] = v[i]-v2.v[i];}
  return vs;
}

template <typename ntype, int NT>
pvector<ntype,NT> pvector<ntype,NT>::operator*(ntype x) const
{
  pvector vs;
  for (int i=0; i < NT; i++){vs.v[i] = v[i]*x;}
  return vs;
}

template <typename ntype, int NT>
pvector<ntype,NT> pvector<ntype,NT>::operator/(ntype x) const
{
  pvector vs;
  for (int i=0; i < NT; i++){vs.v[i] = v[i]/x;}
  return vs;
}

template <typename ntype, int NT>
ntype pvector<ntype,NT>::operator*(const pvector& v2) const
{
  ntype sp=0;
  for (int i=0; i < NT; i++){sp += v[i]*v2.v[i];}
  return sp;
}

template <typename ntype, int NT>
bool pvector<ntype,NT>::operator==(const pvector& v2) const
{
  for (int i=0; i < NT; i++)
    {
      if (v[i] != v2.v[i]) return 0;
    }
  return 1;
}

template <typename ntype, int NT>
pvector<ntype,NT> pvector<ntype,NT>::operator^(const pvector& v2) const
{
  if (NT==3)
    {
      pvector vt;
      vt.v[0] = v[1]*v2.v[2]-v[2]*v2.v[1];
      vt.v[1] = v[2]*v2.v[0]-v[0]*v2.v[2];
      vt.v[2] = v[0]*v2.v[1]-v[1]*v2.v[0];
      return vt;
    }
  else
    {
      std::cout << "Cross product not defined\n";
      exit(1);
    }
}

// friend functions
template <typename ntype, int NT>
std::ostream& operator<<(std::ostream& os, const pvector<ntype,NT>& vec)
{
  os << "(";
  for (int i=0; i < NT; i++)
    {
      os << vec.v[i];
      if (i < NT-1)
        os << ",";
    }
  os << ")";
  return os;
}

// random generator
template <typename ntype, int NT>
pvector<ntype,NT>& pvector<ntype,NT>::random(const ntype& L)
{
  for (int i = 0; i < NT; i++)
  {
    v[i] = (rng.ranf() - 0.5) * L; // assign a random value in [-L/2,L/2]
  }
  return (*this);
}

template <typename ntype, int NT>
void pvector<ntype,NT>::random_orient(void)
{
  ntype rS, S, V1, V2;
  if (NT==3)
    {
      do
        {
          V1 = 2.0*rng.ranf()-1.0;
          V2 = 2.0*rng.ranf()-1.0;
          S = V1*V1+V2*V2;
        }
      while (S >= 1.0);
      rS = sqrt(1.0-S);
      (*this) = {2.0*rS*V1, 2.0*rS*V2, 1.0-2.0*S};
    }
  else std::cout << "[random_orient] Only 3D vectors are supported\n";
}

// rint function outside the class
template<typename ntype, int NT>
pvector<ntype,NT> rint(const pvector<ntype,NT>& vec)
{
  pvector<ntype,NT> vt;

  for (int i=0; i < NT; i++)
    {
      vt(i) = rint(vec(i));
    }
  return vt;
}


using pvec3d=pvector<double,3>;
#endif
