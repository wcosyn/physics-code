#include <iostream>
#include <cstdlib>

using namespace std;

#include <MeanFieldNucleusThick.hpp>
#include <Utilfunctions.hpp>
#include <constants.hpp>
#include <vector>
#include <numint/numint.hpp>
#include <gsl/gsl_sf_bessel.h>


void adap_intr(numint::vector_z &, double r, MeanFieldNucleusThick *pNucleus, int shell, double p);

int main(int argc, char *argv[])
{
  
  string homedir=argv[2];
  int nucleus=atoi(argv[1]);
  
  MeanFieldNucleusThick Nucleus(nucleus,homedir);
 //   cout << max << endl;
  struct Ftor {

    static void exec(const numint::array<double,1> &x, void *param, numint::vector_z &ret) {
      Ftor &p = * (Ftor *) param;
      p.f(ret,x[0],p.pNucleus,p.shell,p.p);
    }
    MeanFieldNucleusThick *pNucleus;
    int shell;
    double p;
    void (*f)(numint::vector_z &, double r, MeanFieldNucleusThick *pNucleus, int shell, double p);

  };
  vector<complex<double> > FT;
  double k=0.;
  for(int i=0;i<40;i++){
    double tot=0.;
    int totsh=0;
    double p=0.05*i*HBARC;
    cout << 0.05*i << " ";
    for(int shell=0;shell < Nucleus.getTotalLevels() ; shell++){
	Ftor F;
	F.pNucleus = &Nucleus;
	F.shell=shell;
	F.p=p;
	numint::mdfunction<numint::vector_z,1> mdf;
	mdf.func = &Ftor::exec;
	mdf.param = &F;

	unsigned neval = 0;
	numint::array<double,1> lower = {{0.}};
	numint::array<double,1> upper = {{Nucleus.getRange()}};
	
	
	F.f=adap_intr;
	
	unsigned count=0;
	numint::cube_romb(mdf,lower,upper,1.E-20,1.E-05,FT,count,0);
	tot+=(norm(FT[0])+norm(FT[1]))*2./PI*2.*abs(Nucleus.getKappas()[shell]);
// 	cout << (norm(FT[0])+norm(FT[1]))*2./PI*2.*abs(Nucleus.getKappas()[shell]) << " ";
	totsh+=2.*abs(Nucleus.getKappas()[shell]);
      //    numint::cube_adaptive(mdf,lower,upper,1.E-20,1.E-03,2.E06,totalcross,count,0);
	
    //       cout << Q2/1.E06 << " " << Ein << " " << totalcross[0]/totalcross[4] << " " << 
    // 	      totalcross[1]/totalcross[4] << " " << totalcross[2]/totalcross[4] << " " << totalcross[3]/totalcross[4] << " " << totalcross[4] << " " << count << endl;
	
      }    
      tot/=totsh;
      cout << tot << endl;
      k+=p*p*tot;
  }
   //cout << k*pow(INVHBARC,3)*3. << endl;
}

void adap_intr(numint::vector_z & results, double r, MeanFieldNucleusThick *pNucleus, int shell, double p){
   results = vector<complex<double> >(2,0.);
   results[1]=results[0]= r*gsl_sf_bessel_jl(pNucleus->getL_array()[shell],p*r*INVHBARC);
   results[0]*=pNucleus->getWave_G(shell,r);
   results[1]*=pNucleus->getWave_F(shell,r);
   
}
