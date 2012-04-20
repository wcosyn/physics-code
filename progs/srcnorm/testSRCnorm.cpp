#include <iostream>
#include <cstdlib>

using namespace std;

#include <MeanFieldNucleusThick.hpp>
#include <FsiCorrelator.hpp>
#include <FastParticle.hpp>
#include <constants.hpp>
#include <GlauberGrid.hpp>
#include <OneGlauberGrid.hpp>
#include <GlauberGridThick.hpp>
#include <Utilfunctions.hpp>

void intr(const double r, double *result, va_list ap);
void intcostheta(const double r, double *result, va_list ap);
void intphi(const double r, double *result, va_list ap);

int main(int argc, char *argv[])
{
  
  string homedir="/home/wim/Code/share";
  MeanFieldNucleusThick CarbonThick(atoi(argv[1]),HOMEDIR);
  FsiCorrelator CarbonCorr(&CarbonThick,90,90,HOMEDIR);
  CarbonCorr.printCorrGridAll();
  //CarbonCorr.printCorrGridProton();
  
//   for(int i=0;i<=15;i++){
//     double r=CarbonThick.getRange()/15.*i;
//     for(int j=0;j<=15;j++){
//       double theta=PI/15.*j;
//       double costheta,sintheta;
//       sincos(theta,&sintheta,&costheta);
//       for(int k=0;k<=6;k++){
// 	double phi=2.*PI/6.*k;
// 	double cosphi,sinphi;
// 	sincos(phi,&sinphi,&cosphi);
// 	double result;
// 	double restimate=0.,thetaestimate=0.,phiestimate=0.;
// 	rombergerN(intr,0.,CarbonThick.getRange(),1,&result,1.E-05,3,8,&restimate,r,costheta,sintheta,cosphi,sinphi,&CarbonThick,&CarbonCorr,&thetaestimate,&phiestimate);
// 	cout << i << " " << j << " " << k << " " << CarbonCorr.getCorrGridProton_interp(r,costheta)*result << endl;
//       }
//     }
//   }
 
  
  
  return 0;
}

void intr(const double r, double *result, va_list ap){
  double rhit = va_arg(ap,double);
  double costhetahit = va_arg(ap,double);
  double sinthetahit = va_arg(ap,double);
  double cosphihit = va_arg(ap,double);
  double sinphihit = va_arg(ap,double);
  MeanFieldNucleusThick* nucleus = va_arg(ap,MeanFieldNucleusThick*);
  FsiCorrelator* corr = va_arg(ap,FsiCorrelator*);
  double *pthetaestimate = va_arg(ap,double*);
  double *pphiestimate = va_arg(ap,double*);
  corr->setRinterp(r);
  rombergerN(intcostheta,-1.,1.,1,result,PREC,3,8,pthetaestimate,r,rhit,costhetahit,sinthetahit
	  ,cosphihit,sinphihit,corr,pphiestimate);
  *result*=nucleus->getProtonDensity(r);
  return;
}

void intcostheta(const double costheta, double *result, va_list ap){
  double r = va_arg(ap,double);
  double rhit = va_arg(ap,double);
  double costhetahit = va_arg(ap,double);
  double sinthetahit = va_arg(ap,double);
  double cosphihit = va_arg(ap,double);
  double sinphihit = va_arg(ap,double);
  FsiCorrelator* corr = va_arg(ap,FsiCorrelator*);
  double *pphiestimate = va_arg(ap,double*);
 
  double sintheta=sqrt(1.-costheta*costheta);
  rombergerN(intphi,0.,2.*PI,1,result,PREC,3,8,pphiestimate,r,costheta,sintheta,rhit,costhetahit,sinthetahit
	  ,cosphihit,sinphihit,corr);
  *result*=corr->getCorrGridProton_interp(costheta);
}

void intphi(const double phi, double *result, va_list ap){
  double r = va_arg(ap,double);
  double costheta = va_arg(ap,double);
  double sintheta = va_arg(ap,double);
  double rhit = va_arg(ap,double);
  double costhetahit = va_arg(ap,double);
  double sinthetahit = va_arg(ap,double);
  double cosphihit = va_arg(ap,double);
  double sinphihit = va_arg(ap,double);
  FsiCorrelator* corr = va_arg(ap,FsiCorrelator*);
  
  double sinphi,cosphi;
  sincos(phi,&sinphi,&cosphi);
  *result= corr->correlation(normr(rhit,costhetahit,sinthetahit,cosphihit,sinphihit,r,costheta,sintheta,cosphi,sinphi));
}