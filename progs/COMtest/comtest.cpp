#include <iostream>
#include <cstdlib>
#include <cmath>


using namespace std;

#include <constants.hpp>
#include <Utilfunctions.hpp>
#include <TDeuteron.h>
#include <TVector3.h>

void int_pcostheta(double costheta, double *result, va_list ap);
void int_pphi(double phi, double *result, va_list ap);
void int_pr(double pr, double *result, va_list ap);
void int_costheta(double costheta, double *result, va_list ap);
void int_phi(double phi, double *result, va_list ap);
double sigma=0.;

int main(int argc, char *argv[])
{
  TDeuteron::Wavefunction *wf = TDeuteron::Wavefunction::CreateWavefunction(argv[1]);
  double pestimate=0.,cthestimate=0.,phiestimate=0.;
  double total=0.;
//   for(int i=0;i<=20;i++){
//     double p = 50.*i;
//     for(int j=0;j<=10;j++){
//       double costheta =1.-0.2*j;
//       double sintheta= sqrt(1.-costheta*costheta);
//       for(int k=0;k<=5;k++){
// 	double phi = 2.*PI/5.*k;
// 	double cosphi = cos(phi);
// 	double sinphi=sin(phi);
// 	TVector3 pvec(p*sintheta*cosphi,p*sintheta*sinphi,p*costheta);
// 	double result;
// 	rombergerN(int_pr,0.,1000.,1,&result,PREC,3,8,&pestimate,&pvec,wf,&cthestimate,&phiestimate);
// 	cout << p << " " << costheta << " " << phi << " " << result << " " << (pow(wf->GetUp(p),2.)+pow(wf->GetWp(p),2.))/(4.*PI) << endl;
// 	total+=p*p*result/((k==0||k==5||j==0||j==10)?1.:2.);
//       }
//     }
//   }
//   cout << total*50*0.2*2.*PI/5. << endl;
  
  double cosestimate=0,phiphiestimate=0.;
  for(int i=0;i<=100;i++){
    double p = 10.*i;
    double result1, result2, result3;
//     sigma=160.;
//     rombergerN(int_pcostheta,-1.,1.,1,&result1,PREC,3,8,&cosestimate,p,wf,&phiphiestimate,&pestimate,&cthestimate,&phiestimate);
//     sigma=180.;
//     rombergerN(int_pcostheta,-1.,1.,1,&result2,PREC,3,8,&cosestimate,p,wf,&phiphiestimate,&pestimate,&cthestimate,&phiestimate);
//     sigma=200.;
//     rombergerN(int_pcostheta,-1.,1.,1,&result3,PREC,3,8,&cosestimate,p,wf,&phiphiestimate,&pestimate,&cthestimate,&phiestimate);    
//     cout << p << " " << (pow(wf->GetUp(p),2.)+pow(wf->GetWp(p),2.)) << " " << result1 << " " << result2 << " " << result3 << endl;
    sigma=164.;
    rombergerN(int_pcostheta,-1.,1.,1,&result1,PREC,3,8,&cosestimate,p,wf,&phiphiestimate,&pestimate,&cthestimate,&phiestimate);
    cout << p << " " << (pow(wf->GetUp(p),2.)+pow(wf->GetWp(p),2.)) << " " << result1 << endl;
  
    //     total+=p*p*result1;
  }
  
//   cout << total*50. << endl;
  delete wf;
}
 
void int_pcostheta(double costheta, double *result, va_list ap){
  double p = va_arg(ap,double);
  TDeuteron::Wavefunction *wf = va_arg(ap,TDeuteron::Wavefunction *);
  double *pphiphiestimate = va_arg(ap,double*);
  double *ppestimate = va_arg(ap,double*);
  double *pcosestimate = va_arg(ap,double*);
  double *pphiestimate = va_arg(ap,double*);
  
  double sintheta=sqrt(1.-costheta*costheta);
  rombergerN(int_pphi,0.,2.*PI,1,result,PREC,3,8,pphiphiestimate,p,costheta,sintheta,wf,ppestimate,pcosestimate,pphiestimate);
  return;
}


void int_pphi(double phi, double *result, va_list ap){
  double p = va_arg(ap,double);
  double costheta = va_arg(ap,double);
  double sintheta = va_arg(ap,double);
  TDeuteron::Wavefunction *wf = va_arg(ap,TDeuteron::Wavefunction *);
  double *ppestimate = va_arg(ap,double*);
  double *pcosestimate = va_arg(ap,double*);
  double *pphiestimate = va_arg(ap,double*);
  
  TVector3 pvec(p*sintheta*cos(phi),p*sintheta*sin(phi),p*costheta);
  rombergerN(int_pr,0.,1000.,1,result,PREC,3,8,ppestimate,&pvec,wf,pcosestimate,pphiestimate);
  return;
}
 
 
void int_pr(double pr, double *result, va_list ap){
  TVector3 *pvec = va_arg(ap,TVector3*);
  TDeuteron::Wavefunction *wf = va_arg(ap,TDeuteron::Wavefunction *);
  double *pcosestimate = va_arg(ap,double*);
  double *pphiestimate = va_arg(ap,double*);
 
  rombergerN(int_costheta,-1.,1.,1,result,PREC,3,8,pcosestimate,pr,pvec,wf,pphiestimate);
  //double sigma=180.;
  *result*=pr*pr/(pow(sigma*sqrt(2.*PI),3.))*exp(-4.*pr*pr/(2.*sigma*sigma))*8.;
  return;
}


void int_costheta(double costheta, double *result, va_list ap){
  double pr = va_arg(ap,double);
  TVector3 *pvec = va_arg(ap,TVector3*);
  TDeuteron::Wavefunction *wf = va_arg(ap,TDeuteron::Wavefunction *);
  double *pphiestimate = va_arg(ap,double*);
  double sintheta=sqrt(1.-costheta*costheta);
  
  rombergerN(int_phi,0.,2.*PI,1,result,PREC,3,8,pphiestimate,pr,costheta,sintheta,pvec,wf);
  return;
}


void int_phi(double phi, double *result, va_list ap){
  double pr = va_arg(ap,double);
  double costheta = va_arg(ap,double);
  double sintheta = va_arg(ap,double);
  TVector3 *pvec = va_arg(ap,TVector3*);
  TDeuteron::Wavefunction *wf = va_arg(ap,TDeuteron::Wavefunction *);
  double *pphiestimate = va_arg(ap,double*);
  
  TVector3 prel = *pvec - TVector3(pr*sintheta*cos(phi),pr*sintheta*sin(phi),pr*costheta);
//   cout << prel.Mag() << endl;
  *result=(pow(wf->GetUp(prel.Mag()),2.)+pow(wf->GetWp(prel.Mag()),2.))/(4.*PI); 
  return;
}
