/** \file
  * Here is a very basic test of the integrators
  * the functions are structs which contain a field
  * that points to the actual function that is to be integrated
  * the use of a struct is convenient if one wants to add spinors
  * or four vectors or other shizzle
  * \author Camille
  */

#include <iostream>
#include <numint/numint.hpp>
#include <numint/typedef.hpp>
#include <cstdlib>
using std::cout;
using std::endl;

void my_1D_func(double, void* param, double& ret);
void my_1D_quad_romberg_integration();

void my_MD_func(const numint::array<double,3> &x, void* param, double& ret);
void my_MD_cube_romberg_integration();


int main(int argc, char* argv[])
{
  //my_1D_quad_romberg_integration();
  my_MD_cube_romberg_integration();
  return 0;
}


void my_1D_func(double x, void* param, double& ret)
{
  ret = sin(x)/x;
}

void my_1D_quad_romberg_integration()
{
  numint::function<double> mystructfunction;
  mystructfunction.func = &my_1D_func;
  double prec = 1.E-04;
  double ret;
  unsigned int count=0;
  numint::quad_romb(mystructfunction,1e-8,1e4,1.E-08,prec,ret,count);
  cout << " returned value is " << ret << endl;
}

void my_MD_func(const numint::array<double,3>& x, void* param, double& ret)
{
  ret = pow(x[0],2)*sin(x[1]);
}
void my_MD_cube_romberg_integration()
{
  numint::mdfunction<double ,3 > mystructmdfunction;
  mystructmdfunction.func = &my_MD_func;
  numint::array<double,3> start = {{0.,0.,0.}};
  numint::array<double,3> stop = {{1.,PI,2.*PI}};
  double prec = 1.E-06; double ret;
  unsigned int count=0;
  
  numint::cube_romb(mystructmdfunction,start,stop,1.E-06,prec,ret,count);
  cout << " returned value is " << ret << endl;
}
