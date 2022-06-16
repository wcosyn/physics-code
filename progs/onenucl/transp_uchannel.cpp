//program used to do some transparency calculations for bill Lee et al. for backward proton production
// this is A(e,e'p pi)A-1 but the pion is unobserved

//cross section is factorized, see notes writings/Notes/QE/electro-photo-transp.xopp
// we integrate over solid angles of pion in (pion,A-1) center of mass system

//elementary process is parametrized using input from Lee & Huber
// see python script  ~/Calculations/A\(ee\'p\)/backwards_T/dsigmaT_fit.py

//run ./transp_uchannel  [Q2,GeV^2]  [Eout, GeV]  [theta_e, degr] [pNout, GeV] [theta_p, rad] [precision in integration]



#include <iostream>
#include <cstdlib>

using namespace std;

#include <MeanFieldNucleusThick.hpp>
#include <constants.hpp>
#include <TKinematics2to3.h>
#include <TElectronKinematics.h>
#include <Utilfunctions.hpp>
#include <vector>
#include <numint/numint.hpp>
#include <OneGlauberGrid.hpp>
#include <TMFSpinor.hpp>

void getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, double Q2, double omega, int shell);
void adap_intpiCM(numint::vector_d &, double costhetapiCM, double phipiCM, MeanFieldNucleusThick *pNucleus, TElectronKinematics *elec,
		OneGlauberGrid *pgrid, double Q2, double omega, double EN, double thetaN);

//used in CM momenta calculation, for instance p_iCM = sqrt(lambda(s,m1^2,m2^2))/(2 sqrt{s})
double lambda(double s, double m1sq, double m2sq){ return (s*s-2*(m1sq+m2sq)*s+pow(m1sq-m2sq,2.));}


// function to compute nuclear cross section
void  getDiffCross(std::vector<double> &cross, MeanFieldNucleusThick * pnucleus, TKinematics2to3 * pkin, TElectronKinematics *pelec, OneGlauberGrid *pgrid,
			     int shellindex, int maxEval);

//u and Q2 in GeV-2
//fit from plot in proposal 2008.10768 and python script mentioned at the start
double getSigmaT(double u, double Q2);

//calculates distorted mom distribution
void rhoD(std::vector<double> &cross, MeanFieldNucleusThick * pnucleus, TKinematics2to3 * pkin, OneGlauberGrid * pgrid, int shellindex, int maxEval);

//for radial integral in distorted momentum distribution
struct Ftor_one {

    /*! integrandum function */
    static void exec(const numint::array<double,3> &x, void *param, numint::vector_z &ret) {
      Ftor_one &p = * (Ftor_one *) param;
      p.f(ret,x[0],x[1],x[2],p.pnucleus,p.pgrid,p.pmiss, p.pspinoru, p.shellindex, p.mindex);
    }
    MeanFieldNucleusThick *pnucleus;
    OneGlauberGrid *pgrid;
    TVector3 pmiss;
    TSpinor *pspinoru;
    int shellindex;
    int mindex;
    void (*f)(numint::vector_z & res, double r, double costheta, double phi, MeanFieldNucleusThick * pnucleus, OneGlauberGrid *pgrid,
        TVector3 &pmiss, TSpinor * pspinoru, int shellindex, int mindex);
  };

void int_r3(numint::vector_z & res, double r, double costheta, double phi, MeanFieldNucleusThick * pnucleus, OneGlauberGrid *grid,
        TVector3 &pmiss, TSpinor *spinoru, int shellindex, int mindex);


int main(int argc, char *argv[])
{
  
  string homedir=HOMEDIR;
  //kinematics input
  double Q2=atof(argv[1])*1.E06; // input in GeV^2 [8.0,9.4,11.4,14.2]
  double Ein=10.6*1.E03;  //input in GeV 10.6 fixed
  double Eout = atof(argv[2])*1.E03; // scattered electron momentum [MeV]
  double thetae = atof(argv[3])*DEGRTORAD; // electron scatt angle
  double pN = atof(argv[4])*1.E03; // outgoing proton momentum [MeV]
  double EN = sqrt(MASSP*MASSP+pN*pN);
  double thetaN = atof(argv[5])*DEGRTORAD; // final proton scatt angle (hall lab frame)
  int nucl = atoi(argv[7]);

  //calculate kinematics
  double omega=Ein-Eout;
  double x=Q2/(2.*MASSP*omega);
  //double q=sqrt(Q2+omega*omega);

  double qx = -Eout*sin(thetae);
  double qz = Ein-Eout*cos(thetae);
  double thetaq = atan2(qx,qz);
  double q = sqrt(qx*qx+qz*qz);
  Q2 = q*q-omega*omega;

  cout << thetaN << " " << thetaq << endl;
  cout << q*cos(thetaq) + Eout*cos(thetae) << endl;
  cout << q*sin(thetaq) + Eout*sin(thetae) << endl;
  thetaN+=thetaq;
  //cout << Q2 << " " << omega << " " << q << " " << Q2/4./Ein/(Ein-omega) << endl;
  
  //glauber parameters
  double prec=atof(argv[6]);
  bool userset=0;//atoi(argv[8]);
  double screening=0.;//atof(argv[9]);
  double lc_mod = 1.;
  double nkt_mod = 1.;
  
  //create objects for nucleus, kinematics and cross section
  MeanFieldNucleusThick Nucleus(nucl,homedir);  //1,6,5,7 for Carbon,Al,Cu,Au
  TElectronKinematics *elec = TElectronKinematics::CreateWithBeamEnergy(Ein);
 

  // for(int i=0;i<=20;i++){
  //   for(int j=0;j<=20;j++){
  //     double costheta = -1+2.*i/20.;
  //     double phi = 2.*PI*j/20.;
  //     TKinematics2to3 kin("","",7,TKinematics2to3::kN,"qsquared:wlab:pn:thn:kycosthkcm:kyphikcm",Q2,omega,EN,thetaN,costheta,phi,
  //               Nucleus.getMassA(), MASSP, M_PI0, Nucleus.getMassA_min_proton()+Nucleus.getExcitation()[1]);
  //     if(kin.IsPhysical()){
  //       double pm=kin.GetPy();
  //       double u = kin.GetTng();
  //       double t = kin.GetTkg();
  //     cout << costheta << " " << phi << " ";
  //       cout << pm << " " << kin.GetPk() << " " << u*1.E-06 << " " << t*1.E-06 << " " << kin.GetY3()[0] << " " << kin.GetY3()[1] << " " << kin.GetY3()[2] << endl;
  //     cout << endl;
  //     }
  //     //else cout << "not" << endl;

  //   }
  // }
  // exit(1);




  //proton momentum is input
  OneGlauberGrid grid(120,36,&Nucleus,prec,2,HOMEDIR);
  FastParticle proton(0, 0, pN,0.,0.,Q2*1.E-06,0.,lc_mod, nkt_mod, HOMEDIR);
  grid.clearParticles();
  grid.addParticle(proton);
  grid.updateGrids();


  //arrays to collect results
  //all p levels  + total proton, times 3 (glauber, CT, pw)
  vector<double> totalcross(3*Nucleus.getPLevels()+3,0.); 
//   double cthmax[Nucleus.getTotalLevels()];
//   double cthmin[Nucleus.getTotalLevels()];
//   double max=-1.;

//   //determine integration limits so pm<300 MeV
//   for(int shell=0;shell<Nucleus.getTotalLevels();shell++) {
//     cthmax[shell]=1.; cthmin[shell]=-1.;
//     getBound(cthmax[shell],cthmin[shell],Nucleus,Q2,omega,shell);
// //     cout << cthmax[shell] << endl;
//     if(max<cthmax[shell]) max=cthmax[shell];
//   }
//   cout << max << endl;

  //integration structure.
  struct Ftor {

    static void exec(const numint::array<double,2> &x, void *param, numint::vector_d &ret) {
      Ftor &p = * (Ftor *) param;
      p.f(ret,x[0],x[1],p.pNucleus,p.elec,p.pgrid,p.Q2,p.omega,p.EN,p.thetaN);
    }
    MeanFieldNucleusThick *pNucleus;
    TElectronKinematics *elec;
    OneGlauberGrid *pgrid;
    double Q2;
    double omega;
    double EN;
    double thetaN;
    void (*f)(numint::vector_d &, double cospicm, double phiphicm, MeanFieldNucleusThick *pNucleus, TElectronKinematics *elec,
	            OneGlauberGrid *pgrid, double Q2, double omega, double EN, double thetaN);

  };

  //initialize integration object + parameters needed to calculate integrals
  Ftor F;
  F.pNucleus = &Nucleus;
  F.elec=elec;
  F.pgrid=&grid;
  F.Q2=Q2;
  F.omega=omega;
  F.EN=EN;
  F.thetaN=thetaN;

  numint::mdfunction<numint::vector_d,2> mdf;
  mdf.func = &Ftor::exec;
  mdf.param = &F;

  unsigned neval = 0;
  numint::array<double,2> lower = {{-1.,0.}};
  numint::array<double,2> upper = {{1.,2.*PI}};
  


  F.f=adap_intpiCM;
  unsigned count=0;
   numint::cube_adaptive(mdf,lower,upper,1.E-20,1.E-02,1.E02,1.E04,totalcross,count,0);
  
  cout << Q2/1.E06 << " " << omega << " " << q  << " ";
  //proton shells transparencies
  //cout << "proton shells T s1/2, p3/2 no CT"
  for(int shell=0;shell<Nucleus.getPLevels();shell++) {
    cout << totalcross[3*shell]/totalcross[3*shell+2] << " ";
  }
  //proton shells T s1/2, p3/2 with CT"
  for(int shell=0;shell<Nucleus.getPLevels();shell++) {
    cout << totalcross[3*shell+1]/totalcross[3*shell+2] << " ";
  }
  //cout << "proton shells total cross no CT, CT, plane-wave"
  cout << totalcross[3*Nucleus.getPLevels()] << " " << totalcross[3*Nucleus.getPLevels()+1] << " " << totalcross[3*Nucleus.getPLevels()+2] << endl;

  //total proton and neutron transparency, proton no CT, with CT; neutron no CT, with CT
  cout << totalcross[3*Nucleus.getPLevels()]/totalcross[3*Nucleus.getPLevels()+2] 
        << " " << totalcross[3*Nucleus.getPLevels()+1]/totalcross[3*Nucleus.getPLevels()+2];
 
  
  delete elec;
  return 0;
}


//integrandum
void adap_intpiCM(numint::vector_d & results, double costhetapiCM, double phipiCM, MeanFieldNucleusThick *pNucleus, TElectronKinematics *elec,
		OneGlauberGrid *pgrid, double Q2, double omega, double EN, double thetaN){
		  
  results=numint::vector_d(3*pNucleus->getPLevels()+3,0.);
  for(int shell=0;shell<pNucleus->getPLevels();shell++) {
    //7 for user input masses, otherwise fixed by pieter's processes
    // in Pieter's language
    // N is the outgoing nucleon
    // K is the pion
    // Y is the A-1 residual nucleus
    TKinematics2to3 kin("","",7,TKinematics2to3::kN,"qsquared:wlab:pn:thn:kycosthkcm:kyphikcm",Q2,omega,EN,thetaN,costhetapiCM,phipiCM,
                pNucleus->getMassA(), MASSP, M_PI0, pNucleus->getMassA_min_proton()+pNucleus->getExcitation()[shell]);

    if(!kin.IsPhysical()){
        double pm=kin.GetPy();
        cout << "bla " << pm << endl;
        }
    else{
      numint::vector_d cross=numint::vector_d(3,0.);
      getDiffCross(cross, pNucleus, &kin, elec, pgrid, shell,20000);
      results[3*shell]+=cross[0];
      results[3*shell+1]+=cross[1];
      results[3*shell+2]+=cross[2];
      results[3*pNucleus->getPLevels()]+=cross[0];
      results[3*pNucleus->getPLevels()+1]+=cross[1];
      results[3*pNucleus->getPLevels()+2]+=cross[2];
      // cout << 2*shell << " " << 2*shell+1 << " " << 2*pNucleus->getTotalLevels() << " " << 2*pNucleus->getTotalLevels()+1 << endl;
      // cout << "0 " << shell << " " << costhetacm << " " << pm << " "  << acos(kin.GetCosthklab())*RADTODEGR << " " 
      // << acos(kin.GetCosthYlab())*RADTODEGR << " " << kin.GetPklab() << " " << kin.GetPYlab() 
      // << " " << kin.GetKlab() << " " << kin.GetWlab() <<  " " << results[3*shell] << " " << results[3*shell+2] << endl;
      cout << shell << " " << costhetapiCM << " " << phipiCM << " " << kin.GetPy() << " " << kin.GetPk() << " " << results[3*shell] << " " << results[3*shell+1] << " " << results[3*shell+2] << endl;
    }
  }

    
}

void  getDiffCross(std::vector<double> &cross, MeanFieldNucleusThick * pnucleus, TKinematics2to3 * pkin, TElectronKinematics *pelec, OneGlauberGrid *pgrid,
			     int shellindex, int maxEval){
  
  cross=vector<double>(3,0.);
  double front = pkin->GetMy()/(4.*pkin->GetSky())*sqrt(lambda(pkin->GetSky(),pow(pkin->GetMk(),2.),pow(pkin->GetMy(),2.))
                                                *lambda(pkin->GetSkn(),pow(pkin->GetMn(),2.), -pkin->GetQsquared())
                                                /lambda(pkin->GetStot(),pow(pkin->GetMd(),2.),-pkin->GetQsquared()));
  //extra 1/2pi because of dsigmaT/2pi
  double gamma = ALPHA/(4.*PI*PI*PI)*(pelec->GetBeamEnergy(pkin->GetQsquared(),pkin->GetWlab())-pkin->GetWlab())
                  /pelec->GetBeamEnergy(pkin->GetQsquared(),pkin->GetWlab())*(pkin->GetSkn()-MASSP*MASSP)/(2.*MASSP)/pkin->GetQsquared()
                  /(1-pelec->GetEpsilon(pkin->GetQsquared(),pkin->GetWlab()));
  double sigmaT = getSigmaT(pkin->GetTng()*1.E-06,pkin->GetQsquared()*1.E-06);

  rhoD(cross,pnucleus,pkin,pgrid,shellindex,maxEval);
  for(int i=0;i<3;i++){
    cout << cross[i] << " ";
    cross[i]*=front*gamma*sigmaT;
  }
  return;
}



double getSigmaT(double u, double Q2){
  double a=1.07492722, b=-0.92753151, c=-1.58234828, d=2.84226756;

  return exp(a+b*Q2+c/Q2+d*u);

}


void rhoD(std::vector<double> &cross, MeanFieldNucleusThick * pnucleus, TKinematics2to3 * pkin, OneGlauberGrid * pgrid, int shellindex, int maxEval){
  //careful!! we compute in a frame with p_f along z-axis, so have to rotate the pmiss vector!
  cross=vector<double>(3,0.);
  double thetaN=pkin->GetThn();
  TVector3 pN = pkin->GetN3();
  //cout << "pre " << pN[0] << " " << pN[1] << " " << pN[2] << endl;
  pN.RotateY(-thetaN);
  //cout << "post " << pN[0] << " " << pN[1] << " " << pN[2] << endl;
  TVector3 pmiss = -pkin->GetY3();
  pmiss.RotateY(-thetaN);
  FourVector<double> pmiss4(sqrt(MASSP*MASSP+pmiss.Mag2()),pmiss[0],pmiss[1],pmiss[2]);
  TSpinor spinoru(pmiss4,MASSP,TSpinor::Polarization(0.,0.,TSpinor::Polarization::kUp),TSpinor::kUnity);

  for(int m=-pnucleus->getJ_array()[shellindex];m<=pnucleus->getJ_array()[shellindex];m+=2){
     //only half of the m values due to parity symmetry, careful!!! 
    //cout << shellindex << " " << m << endl;
    pgrid->clearKnockout();
    pgrid->addKnockout(shellindex,m);    
    for(int i=0;i<1;i++){
      
      numint::array<double,3> lower = {{0.,-1.,0.}};
      numint::array<double,3> upper = {{pnucleus->getRange(),1.,2.*PI}};
      Ftor_one F;
      F.pnucleus = pnucleus;
      F.pgrid = pgrid;
      F.pmiss = pmiss;
      F.pspinoru = &spinoru;
      F.shellindex=shellindex;
      F.mindex = m;
      numint::mdfunction<numint::vector_z,3> mdf;
      mdf.func = &Ftor_one::exec;
      mdf.param = &F;
      numint::vector_z ret(3,0.);
      F.f=int_r3;
      unsigned count=0;
      numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-04,1E03,maxEval,ret,count,0);
      
      //factor of 2 because of parity symmetry
      for(int jj=0;jj<3;jj++) cross[jj]+= 2.*norm(ret[jj]);

    }
  }

}

void int_r3(numint::vector_z & results, double r, double costheta, double phi, MeanFieldNucleusThick * pnucleus, OneGlauberGrid * pgrid,
        TVector3 &pmiss, TSpinor * pspinoru, int shellindex, int mindex){

  results=numint::vector_z(3,0.);
  double sintheta=sqrt(1.-costheta*costheta);
  
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  TMFSpinor wave(*pnucleus,shellindex,mindex,r,costheta,phi);
  complex<double> exp_pr=r*exp(-I_UNIT*INVHBARC*(pmiss*TVector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)));
  complex<double> temp = TSpinor::Bar(*pspinoru)*wave;
  for(int i=0; i<3; i++) results[i] = temp*exp_pr;
  results[0]*=pgrid->getFsiGridFull_interp3(r,costheta,phi);
  results[1]*=pgrid->getFsiCtGridFull_interp3(r,costheta,phi);

  return;
}


// // find cm scatt angles so that missing momentum stays below 300 MeV 
// void getBound(double &high, double &low, MeanFieldNucleusThick &nucleus, double Q2, double omega, int shell){
//   TKinematics2to2 kin("","",nucleus.getMassA(),
// 		      (shell<nucleus.getPLevels()? nucleus.getMassA_min_proton(): nucleus.getMassA_min_neutron())+nucleus.getExcitation()[shell],
// 		      MASSP,"qsquared:wlab:costhkcm",Q2,omega,(high+low)/2.);
//   double pm=kin.GetPklab();
//   if(pm<300.) low=(high+low)/2.;
//   else high=(high+low)/2.;
//   if((high-low)<1.E-04) { return;}
//   else getBound(high,low,nucleus,Q2,omega,shell);  
// }

