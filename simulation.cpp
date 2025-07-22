#include "./simulation.hpp"
#include <iostream>


int main(int argc, char **argv)
{
  // Pass to the program: dt, ratio of I kind, kind of simultaion(1 for SVV, 2 for MTS), number of simulations, epsilon
  if (argc < 6)
  {
    std::cerr << "Usage: " << argv[0] << " <dt>" << std::endl;
    return 1;
  }

  double dt = std::stod(argv[1]);
  double r = std::stod(argv[2]);
  double flag = std::stod(argv[3]);
  int count = std::stoi(argv[4]);
  double eps = std::stod(argv[5]);
  int M,N1,N2;
  M=-1;

// PER FARE SIMULAZIONI CON MULTI-TIME STEP METHOD
if (flag==2)
{
  std::cout <<"MS : " << r <<" :\t";
  for (size_t i = 0; i < count; i++)
  {
    multistep<particleLJ> md;
    md.set(dt, r);
    md.set_eps(eps);
    md.init_rng(0);            
    md.prepare_initial_conf(); 
    md.run();                  
    // std::cout<<i+1<<" ";
    // if (i+1==count) M=md.get_M();
  }
  //std::cout << "\n";
}
//std::cout<<M<<"\n";

// PER FARE SIMULAZIONI CON STANDARD VV

if (flag == 1)
{
  std::cout << "MD : " << r << " :\t";
  for (size_t i = 0; i < count; i++)
  {
    mdsim<particleLJ> md;
    md.set(dt, r);
    md.init_rng(0);            
    md.prepare_initial_conf(); 
    md.run();                 
    std::cout << i+1<<" ";
  }
}
std::cout << "\n";
return 0;
}