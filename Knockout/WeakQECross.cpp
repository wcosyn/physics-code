#include "WeakQECross.hpp"

#include <cmath>
#include <TSpinor.h>
#include <OneGlauberGrid.hpp>
#include <GlauberGridThick.hpp>
#include <TMFSpinor.hpp>
#include <cassert>
#include <constants.hpp>
#include "GammaStructure.h"

using namespace std;

WeakQECross::WeakQECross(TElectronKinematics *elec, MeanFieldNucleusThick *pnucleus, double precision, int integr,
	     string dir, bool cc, double M_A_in, bool user_sigma, double gA_s_in, double r_s2_in, 
	      double mu_s_in, double sigma_screening)
:homedir(dir),electron(elec),pnucl(pnucleus),reacmodel(NULL), prec(precision),integrator(integr),
charged(cc), M_A(M_A_in), usersigma(user_sigma), gA_s(gA_s_in), r_s2(r_s2_in), mu_s(mu_s_in), sigmascreening(sigma_screening){
  if(cc){
     cerr << "ERROR in WeakQECross::constructor "
	<< "invalid, use TLeptonKinematics constructor for a CC process.\n";
      assert(1==0);
  }
}

WeakQECross::WeakQECross(TLeptonKinematics *lep, 
			 MeanFieldNucleusThick *pnucleus, double precision, int integr,
	     string dir, bool cc, double M_A_in, bool user_sigma, double gA_s_in, double r_s2_in, 
	      double mu_s_in, double sigma_screening)
:homedir(dir),lepton(lep),pnucl(pnucleus),reacmodel(NULL), prec(precision),integrator(integr),
charged(cc), M_A(M_A_in), usersigma(user_sigma), gA_s(gA_s_in), r_s2(r_s2_in), mu_s(mu_s_in), sigmascreening(sigma_screening){
  if(!cc){
     cerr << "ERROR in WeakQECross::constructor "
	<< "invalid, use TElectronKinematics constructor for a NC process.\n";
      assert(1==0);
  }
}


WeakQECross::~WeakQECross(){
  
}

double WeakQECross::getElWeakQECross(double Q2, double E_in, int current, double phi, int maxEval, bool phi_int,
 				    bool neutrino, bool proton){

   double omega=0.5*(MASSN*MASSN-MASSP*MASSP+Q2)/MASSP;
   double qvec=sqrt(omega*omega+Q2);              // spacelike Q2???
   //double Q2overkk=Q2/qvec/qvec;
   //double tan2=electron.GetTan2HalfAngle(kin);
   reacmodel=new WeakQEHadronCurrent(pnucl,prec,integrator,homedir,maxEval,charged, 
 				    M_A, getUsersigma(), gA_s, r_s2, mu_s,getSigmascreening());
 
   int numb_of_resp=5;
   double kinfactors[numb_of_resp];
   double response[numb_of_resp];
   for(int i=0;i<numb_of_resp;i++) {response[i]=0; kinfactors[i]=0;}
   
   //compute response[0] functions
   for(int i=0;i<numb_of_resp;i++) response[i]=0.;
   for(int spinin=-1;spinin<=1;spinin+=2){
     for(int spinout=-1;spinout<=1;spinout+=2){
       complex<double> j0=reacmodel->getFreeMatrixEl(Q2,proton,current, spinin,spinout,0);
       complex<double> j1=reacmodel->getFreeMatrixEl(Q2,proton,current, spinin,spinout,1);
       complex<double> j2=reacmodel->getFreeMatrixEl(Q2,proton,current, spinin,spinout,2);
       complex<double> j3=reacmodel->getFreeMatrixEl(Q2,proton,current, spinin,spinout,3);
//     cout << "current " << spinout << " " << j0 << " " << j1 << " " << j2 << " " << j3 << endl;
       response[0]+=norm(j0);
       response[1]+=-2.*real(conj(j0)*j3);    // *.2 because double counting in crosssection
       response[2]+=norm(j3);
       response[3]+=norm(j1)+norm(j2);
       if(neutrino)  response[4]+=2.*imag(conj(j1)*j2);  // *.2 because double counting in cs
       if(!neutrino) response[4]-=2.*imag(conj(j1)*j2);
     }
   }
// cout<<"Response "<<response[0]<<" "<<response[1]<<" "<<response[2]<<" "<<response[3]<<" "<<response[4]<<endl;
   double mmu=lepton->GetLeptonMass();   
   double eps=E_in;       
   double eeps=eps-omega;
   double kk=sqrt(eeps*eeps-mmu*mmu);
   double delta2=mmu*mmu/abs(Q2);
   double rho=abs(Q2)/(qvec*qvec);  
   double rrho=qvec/(2.*eps-omega);                
   double tan2theta=abs(Q2)/((eeps+eps)*(eps+eeps)-qvec*qvec);
   //kinematical factors
   kinfactors[0]=1.-delta2*tan2theta;
   kinfactors[1]=omega/qvec+delta2/rrho*tan2theta;     
   kinfactors[2]=omega*omega/(qvec*qvec)+(1.+2.*omega/(qvec*rrho)+rho*delta2)*delta2*tan2theta;
   kinfactors[3]=tan2theta+0.5*rho-delta2*(omega/qvec+0.5*rho*rrho*delta2)*tan2theta/rrho;  
   kinfactors[4]=(1.-omega*rrho*delta2/qvec)*tan2theta/rrho;      
// cout<<"Kinfactors "<<kinfactors[0]<<" "<<kinfactors[1]<<" "<<kinfactors[2]<<" "<<kinfactors[3]<<" "<<kinfactors[4]<<endl;
   double cross=0;
   for(int i=0;i<=4;i++) cross+=kinfactors[i]*response[i];
   cross*=G_FERMI*G_FERMI*0.975*0.975/(2*PI*PI)*kk*eeps*cos(atan(sqrt(tan2theta)))*cos(atan(sqrt(tan2theta)));
   cross*=HBARC*HBARC;
// cout<<"Factors "<<kk<<" "<<eeps<<" "<<cos(atan(sqrt(tan2theta)))<<" "<<cos(atan(sqrt(tan2theta)))<<endl;
   if (phi_int) cross*=phi;
   
   delete reacmodel;
   
   return cross;
   
 }

/*
double WeakQECross::getElWeakQECross(TKinematics2to2 &kin, int current, double phi, int maxEval, bool phi_int,
 				    bool neutrino, bool proton){
   double Q2=kin.GetQsquared();
   double qvec=kin.GetKlab();
   double Q2overkk=Q2/qvec/qvec;
   //double tan2=electron.GetTan2HalfAngle(kin);
   reacmodel=new WeakQEHadronCurrent(pnucl,prec,integrator,homedir,maxEval,charged, 
 				    M_A, getUsersigma(), gA_s, r_s2, mu_s,getSigmascreening());
 
   int numb_of_resp;
   double kinfactors[numb_of_resp];
   double response[numb_of_resp];
   
   //compute response[0] functions
   for(int i=0;i<numb_of_resp;i++) response[i]=0.;
   for(int spinin=-1;spinin<=1;spinin+=2){
     for(int spinout=-1;spinout<=1;spinout+=2){
       complex<double> jmin=reacmodel->getFreeMatrixEl(kin,proton, current, spinin,spinout,-1);
       complex<double> j0=reacmodel->getFreeMatrixEl(kin,proton,current, spinin,spinout,0);
       complex<double> jplus=reacmodel->getFreeMatrixEl(kin,proton,current, spinin,spinout,1);
       complex<double> jz=reacmodel->getFreeMatrixEl(kin,proton,current, spinin,spinout,3);
 //       cout << jmin << " " << j0 << " " << jplus << endl;
       response[0]+=norm(j0);
       response[1]+=norm(jmin)+norm(jplus);
       response[2]+=2.*real(conj(jplus)*jmin);
       response[3]+=2.*real(conj(j0)*(jplus-jmin));
     }
   }
   
   double Einon=sqrt(kin.GetHyperonMass()*kin.GetHyperonMass()+kin.GetKlab()*kin.GetKlab()
    +kin.GetPYlab()*kin.GetPYlab()-2.*kin.GetPYlab()*kin.GetKlab()*kin.GetCosthYlab());
   //average incoming proton spin
 //   response[0]/=2.;
 //   response[1]/=2.;
 //   response[2]/=2.;
 //   response[3]/=2.;
   delete reacmodel;
   //combine everything
   return 1.;
//   kin.GetPYlab()*mott*kin.GetHyperonMass()*kin.GetHyperonMass()/(Einon)*(kinfactors[0]*response[0][0]
// 	  +kinfactors[1]*response[0][1]+kinfactors[2]*response[0][2]*cos(2.*phi)+kinfactors[3]*response[0][3]*cos(phi));
//   
 } */

double WeakQECross::getDiffWeakQECross(TKinematics2to2 &kin, int current, int thick, int SRC, int CT, int pw, int shellindex, 
			   double phi, int maxEval, bool lab, bool phi_int){
  
  //electron kinematics
  double Q2=kin.GetQsquared();
  double qvec=kin.GetKlab();
  double Q2overkk=Q2/qvec/qvec;
  //kinematical front factor
  //to translate the 2to2kinematics language to our particles:
  //p is A
  //hyperon Y is fast final nucleon
  //kaon is residual A-1 nucleus
  
  //m_n*m_{A-1}*p_f/E_{A-1}/f_rec (for lab)
  double frontfactor=lab? (kin.GetHyperonMass()*kin.GetMesonMass()*kin.GetPYlab()/(8*pow(PI,3.))
      /abs(kin.GetEklab()+kin.GetEYlab()*(1-qvec*kin.GetCosthYlab()/kin.GetPYlab()))) :
	    kin.GetHyperonMass()*kin.GetMesonMass()*kin.GetPkcm()/(8*pow(PI,3.)*kin.GetW());
  reacmodel=new WeakQEHadronCurrent(pnucl,prec,integrator,homedir,maxEval,charged, 
				    M_A, getUsersigma(), gA_s, r_s2, mu_s,getSigmascreening());
  
  //CC
  if(charged){
    double Ebeam=lepton->GetBeamEnergy(kin);
    double Eout=Ebeam-kin.GetWlab();
    double leptonmass2=lepton->GetLeptonMass()*lepton->GetLeptonMass();
    double mott=sqrt(1.-leptonmass2/Eout/Eout)*pow(G_FERMI*COS_CAB*Eout/(Q2/M_W/M_W+1.)/PI/2.,2.); //[MeV]^2
    
    double kinfactors[9];
    double response[9];

    FourVector<double> k_in,k_out;
    lepton->GetLeptonVectors(kin,k_in,k_out);
    double massfactor=sqrt(1.-leptonmass2/k_out[0]/k_out[0]);
    
    
    
    double qvec = kin.GetKlab();
    double nu = kin.GetWlab();
    double Q2=kin.GetQsquared();
    double Q=sqrt(Q2);
    double costhl=lepton->GetCosScatterAngle(kin);
    double sinthl=sqrt(1.-costhl*costhl);
    
    kinfactors[0]=1.+massfactor*costhl-2.*leptonmass2*nu/k_out[0]/M_W/M_W
	+leptonmass2*pow(nu/M_W/M_W,2.)*(1.-massfactor*costhl); //v_L (part corresponding with unbroken current)
    kinfactors[1] = -leptonmass2/qvec/k_out[0] + leptonmass2/M_W/M_W/qvec*(nu*(1.-massfactor*costhl)+Q2/k_out[0])
      -nu/qvec*leptonmass2*pow(Q/M_W/M_W,2.)*(1-massfactor*costhl); //second contrib to v_L, mix z,0
    kinfactors[2] = leptonmass2/qvec/qvec*pow(1.-Q2/M_W/M_W,2.)*(1-massfactor*costhl); //third contrib to v_L, z*z
    kinfactors[3]=1.-massfactor*costhl+k_in[0]*k_out[0]/qvec/qvec*massfactor*massfactor*sinthl*sinthl; //v_T
    kinfactors[4]=-k_in[0]*k_out[0]/qvec/qvec*massfactor*massfactor*sinthl*sinthl; //v_TT
    kinfactors[5]=-sinthl/sqrt(2.)/qvec*(k_in[0]+k_out[0]-nu*leptonmass2/M_W/M_W); //v_TL
    kinfactors[6] =-sinthl*leptonmass2/sqrt(2.)/qvec/qvec*(1.-Q2/M_W/M_W); //v_TL contrib from mixed current
    kinfactors[7]=(k_in[0]+k_out[0])/qvec*(1.-massfactor*costhl)-leptonmass2/qvec/k_out[0]; //v_T'
    kinfactors[8]=-massfactor*sinthl/sqrt(2.); //v_TL'
    
// testing kin factors
//     const FourVector<GammaStructure> gamma_mu=FourVector<GammaStructure>(GammaStructure(0.,0.,1.),
// 											GammaStructure(0.,0.,0.,1.),
// 				      GammaStructure(0.,0.,0.,0.,1.),GammaStructure(0.,0.,0.,0.,0.,1.));
// 
// 
//     const GammaStructure gamma_5(0.,1.);
// 
//     FourVector<complex<double> > polVectorPlus(0.,
// 						      -1./sqrt(2.),
// 						      complex<double>(0.,-1./sqrt(2.)),
// 						  0.);
//     FourVector<complex<double> > polVectorMin(0.,
// 						    1./sqrt(2.),
// 						    complex<double>(0.,-1./sqrt(2.)),
// 						    0.);
//     FourVector<complex<double> > polVector0(qvec/Q,0.,0.,nu/Q);
//     FourVector<complex<double> > polVectorZ(nu/Q,0.,0.,qvec/Q);
//     polVectorZ*=1.-Q2/M_W/M_W;
    

//     cout << -Trace(((polVectorMin*gamma_mu)*(k_in*gamma_mu)*(polVectorPlus*gamma_mu)*(k_out*gamma_mu)).value())/(4.*k_in[0]*k_out[0]) 
//     << " " << 1-massfactor*costhl+k_in[0]*k_out[0]/qvec/qvec*massfactor*massfactor*sinthl*sinthl << endl;
// 
//     cout << -Trace(((polVectorMin*gamma_mu)*(k_in*gamma_mu)*(polVectorMin*gamma_mu)*(k_out*gamma_mu)).value())/(4.*k_in[0]*k_out[0]) 
//     << " " << -k_in[0]*k_out[0]/qvec/qvec*massfactor*massfactor*sinthl*sinthl << endl;
//     
//     cout << -Trace(((polVectorMin*gamma_mu)*gamma_5*(k_in*gamma_mu)*(polVectorPlus*gamma_mu)*(k_out*gamma_mu)).value())/(4.*k_in[0]*k_out[0]) << " " 
//     << (k_in[0]+k_out[0])/qvec*(1.-massfactor*costhl)-leptonmass2/qvec/k_out[0] << endl;
//     
//     cout << Trace(((polVector0*gamma_mu)*gamma_5*(k_in*gamma_mu)*(polVectorPlus*gamma_mu)*(k_out*gamma_mu)).value())/(4.*k_in[0]*k_out[0]) << " " 
//     << massfactor*sinthl*Q/qvec/sqrt(2.) << endl;
//     
//     cout << Trace(((polVectorZ*gamma_mu)*gamma_5*(k_in*gamma_mu)*(polVectorPlus*gamma_mu)*(k_out*gamma_mu)).value())/(4.*k_in[0]*k_out[0]) << endl;
//     cout << Trace(((polVector0*gamma_mu)*gamma_5*(k_in*gamma_mu)*(polVectorZ*gamma_mu)*(k_out*gamma_mu)).value())/(4.*k_in[0]*k_out[0]) << endl;
// 
//     cout << Trace(((polVector0*gamma_mu)*(k_in*gamma_mu)*(polVectorPlus*gamma_mu)*(k_out*gamma_mu)).value())/(4.*k_in[0]*k_out[0]) 
//     << " " << sinthl/qvec/qvec*massfactor/sqrt(2.)*(Q*(k_in[0]+k_out[0])-nu*leptonmass2/Q) << endl;
// 
//     cout << Trace(((polVectorZ*gamma_mu)*(k_in*gamma_mu)*(polVectorPlus*gamma_mu)*(k_out*gamma_mu)).value())/(4.*k_in[0]*k_out[0]) 
//     << " " << -sinthl/qvec*massfactor/sqrt(2.)*leptonmass2/Q*(1.-Q2/M_W/M_W) << endl;
//     
//     cout << std::setprecision(9) <<Trace(((polVectorZ*gamma_mu)*(k_in*gamma_mu)*(polVectorZ*gamma_mu)*(k_out*gamma_mu)).value())/2. 
//     << " " << leptonmass2*(1.+leptonmass2/Q2)*pow(1.-Q2/M_W/M_W,2.) << endl;
//     
//     cout << std::setprecision(9) <<Trace(((polVector0*gamma_mu)*(k_in*gamma_mu)*(polVectorZ*gamma_mu)*(k_out*gamma_mu)).value())/2. 
//     << " " << -leptonmass2/qvec*((k_in[0]+k_out[0])-leptonmass2*nu/Q2)*(1.-Q2/M_W/M_W) << endl;
// 
//     cout << std::setprecision(9) <<Trace(((polVector0*gamma_mu)*(k_in*gamma_mu)*(polVector0*gamma_mu)*(k_out*gamma_mu)).value())/2. 
//     << " " << (pow(k_in[0]+k_out[0],2.)/qvec/qvec-1.)*Q2-leptonmass2*(1.+2.*(k_in[0]*k_in[0]-k_out[0]*k_out[0])/qvec/qvec)
//     +pow(leptonmass2,2.)*nu*nu/qvec/qvec/Q2 << " " << -Q2-leptonmass2+pow(Q2*(k_in[0]+k_out[0])-leptonmass2*nu,2.)/qvec/qvec/Q2 << endl;
// 
//     cout << std::setprecision(9) <<(qvec*qvec*Trace(((polVector0*gamma_mu)*(k_in*gamma_mu)*(polVector0*gamma_mu)*(k_out*gamma_mu)).value()) 
//       - 2.*qvec*nu*Trace(((polVector0*gamma_mu)*(k_in*gamma_mu)*(polVectorZ*gamma_mu)*(k_out*gamma_mu)).value())
//       +nu*nu*Trace(((polVectorZ*gamma_mu)*(k_in*gamma_mu)*(polVectorZ*gamma_mu)*(k_out*gamma_mu)).value()))/(4.*k_in[0]*k_out[0]*Q2) << 
//       " " << 1.+massfactor*costhl -2.*leptonmass2*nu/k_out[0]/M_W/M_W+leptonmass2*pow(nu/M_W/M_W,2.)*(1.-massfactor*costhl)<< endl;
//     
//     cout << std::setprecision(9) <<
//       (Trace(((polVector0*gamma_mu)*(k_in*gamma_mu)*(polVectorZ*gamma_mu)*(k_out*gamma_mu)).value())
//       -nu/qvec*Trace(((polVectorZ*gamma_mu)*(k_in*gamma_mu)*(polVectorZ*gamma_mu)*(k_out*gamma_mu)).value()))/(4.*k_in[0]*k_out[0]) << 
//       " " << -leptonmass2/qvec/k_out[0] + leptonmass2/M_W/M_W/qvec*(nu*(1.-massfactor*costhl)+Q2/k_out[0])
//       -nu/qvec*leptonmass2*pow(Q/M_W/M_W,2.)*(1-massfactor*costhl) << endl;
// 
//     cout << std::setprecision(9) <<
//       Q2/qvec/qvec*Trace(((polVectorZ*gamma_mu)*(k_in*gamma_mu)*(polVectorZ*gamma_mu)*(k_out*gamma_mu)).value())/(4.*k_in[0]*k_out[0]) << 
//       " " << leptonmass2/qvec/qvec*pow(1.-Q2/M_W/M_W,2.)*(1-massfactor*costhl) << endl;
// //    z,z comparison pascal me   
//     cout << std::setprecision(9) << nu*nu/qvec/qvec*(1.+massfactor*costhl)+2.*nu/qvec*leptonmass2/qvec/k_out[0]+leptonmass2/qvec/qvec*(1.-massfactor*costhl) <<
//     " " << 1.+massfactor*costhl-2.*k_in[0]*k_out[0]*massfactor*massfactor*sinthl*sinthl/qvec/qvec << endl;


    double extraresponse[4];
    for(int i=0;i<4;i++) extraresponse[i]=0.;
    
    
    for(int i=0;i<9;i++) response[i]=0.;
    for(int m=-pnucl->getJ_array()[shellindex];m<=pnucl->getJ_array()[shellindex];m+=2){
	Matrix<2,4> J;
	Matrix<2,4> Vector;
	Matrix<2,4> Axial;
	reacmodel->getMatrixEl(kin,Vector,shellindex,m,CT,pw, current, SRC, thick,1);
	reacmodel->getMatrixEl(kin,Axial,shellindex,m,CT,pw, current, SRC, thick,0);
	J=Vector+Axial;
	for(int i=0;i<1;i++){ //exploit the parity symmetry of vector and axial parts to simplify things a bit
	  
//  	  cout << i << " " << m << " " << J(i,0) << " " << J(i,1) << " " <<  J(i,2) << " " << J(i,3) << endl;
	  response[0]+=2.*(norm(Vector(i,0)-kin.GetWlab()/qvec*Vector(i,3))+
		  norm(Axial(i,0)-kin.GetWlab()/qvec*Axial(i,3))); //W_L1
	  response[1]+=4.*real(Vector(i,3)*conj(Vector(i,0)-kin.GetWlab()/qvec*Vector(i,3))
	    +Axial(i,3)*conj(Axial(i,0)-kin.GetWlab()/qvec*Axial(i,3)));//W_L2
	  response[2]+=2.*(norm(Vector(i,3))+norm(Axial(i,3))); //W_L3
	  response[3]+=2.*(norm(Vector(i,1))+norm(Vector(i,2))+norm(Axial(i,1))+norm(Axial(i,2))); //W_T
	  response[7]+=4.*real(conj(Vector(i,2))*Axial(i,2)-conj(Vector(i,1))*Axial(i,1)); //W_T'
	  if(!phi_int){
	    response[4]+=4.*real(Vector(i,2)*conj(Vector(i,1))+Axial(i,2)*conj(Axial(i,1))); //W_TT
	    response[5]+=4.*real((Vector(i,0)-kin.GetWlab()/qvec*Vector(i,3))*conj(Vector(i,2)-Vector(i,1))+
		(Axial(i,0)-kin.GetWlab()/qvec*Axial(i,3))*conj(Axial(i,2)-Axial(i,1))); //W_LT1
	    response[6]+=4.*real((Vector(i,3))*conj(Vector(i,2)-Vector(i,1))
	      +(Axial(i,3))*conj(Axial(i,2)-Axial(i,1))); //W_LT2
	    response[8]+=4.*imag((Vector(i,0)-kin.GetWlab()/qvec*Vector(i,3))*conj(Vector(i,2)-Vector(i,1))
	      +(Axial(i,0)-kin.GetWlab()/qvec*Axial(i,3))*conj(Axial(i,2)-Axial(i,1))); //W_LT'
	    
	    //these are zero in electron scattering but not for CC due to parity rules!!!!
	    extraresponse[0]+=-4.*imag(Vector(i,2)*conj(Axial(i,1))+Axial(i,2)*conj(Vector(i,1)));
	    extraresponse[1]+=4.*imag((Vector(i,0)-kin.GetWlab()/qvec*Vector(i,3))*conj(Axial(i,2)+Axial(i,1))
	      +(Axial(i,0)-kin.GetWlab()/qvec*Axial(i,3))*conj(Vector(i,2)+Vector(i,1)));
	    extraresponse[2]+=4.*imag((Vector(i,3))*conj(Axial(i,2)+Axial(i,1))+
	      (Axial(i,3))*conj(Vector(i,2)+Vector(i,1)));
	    extraresponse[3]+=4.*real((Vector(i,0)-kin.GetWlab()/qvec*Vector(i,3))*conj(Axial(i,2)+Axial(i,1))
	      +(Axial(i,0)-kin.GetWlab()/qvec*Axial(i,3))*conj(Vector(i,2)+Vector(i,1)));
	  }
	}
    }
//     for(int i=0;i<9;i++) cout << kinfactors[i] << " " << response[i] << endl;
//     for(int i=0;i<4;i++) cout << extraresponse[i] << endl;
//     cout << kinfactors[0]*response[0]+kinfactors[1]*response[1]+kinfactors[2]*response[2] << " "
//      << response[3] << " " << response[7] << endl;
    double result=0.;
    if(!phi_int) result=kinfactors[0]*response[0]+kinfactors[1]*response[1]+kinfactors[2]*response[2]
	      +kinfactors[3]*response[3]+kinfactors[4]*(response[4]*cos(2.*phi)+extraresponse[0]*sin(2.*phi))
	      +(kinfactors[5]*response[5]+kinfactors[6]*response[6])*cos(phi)
	      +(kinfactors[5]*extraresponse[1]+kinfactors[6]*extraresponse[2])
	      +(shellindex<pnucl->getPLevels()?1.:-1.)*
	      (kinfactors[7]*response[7]+kinfactors[8]*(response[8]*sin(phi)+extraresponse[3]*cos(phi)));
    else result=2.*PI*(kinfactors[0]*response[0]+kinfactors[1]*response[1]+kinfactors[2]*response[2]
	      +kinfactors[3]*response[3]+(shellindex<pnucl->getPLevels()?1.:-1.)*kinfactors[7]*response[7]);
//  else result=2.*PI*(kinfactors[7]*response[7]);
    delete reacmodel;
//     cout <<"mott " << mott << " " << frontfactor << endl;
//     cout << kinfactors[3] << " " << kinfactors[7] << " " << kinfactors[0] << " " << response[3] << " " << response[7] << " " << response[0] << " " 
//     << massfactor << " " << mott << " " << frontfactor << " " << result << endl;
    return mott*frontfactor*result/HBARC;
  }  
  //NC
  else{
     double mott=(1.+electron->GetCosScatterAngle(kin))*
	pow(electron->GetBeamEnergy(kin)-kin.GetWlab()*G_FERMI/(Q2/M_Z/M_Z+1.)/PI,2.)/2.;
    double tan2=pow(tan(acos(electron->GetCosScatterAngle(kin))/2.),2.);
    double kinfactors[6];
    kinfactors[0]=1.; //v_L
    kinfactors[1]=tan2+Q2overkk/2.; //v_T
    kinfactors[2]=-Q2overkk/2.; //v_TT
    kinfactors[3]=-sqrt((tan2+Q2overkk)/2.); //v_TL
    kinfactors[4]=sqrt(tan2*(tan2+Q2overkk)); //v'_T
    kinfactors[5]=-1./sqrt(2)*sqrt(tan2); //v'_TL
    
    //compute response functions
    double response[6];
    double extraresponse[3];
    for(int i=0;i<6;i++) response[i]=0.;
    for(int i=0;i<3;i++) extraresponse[i]=0.;
    
    for(int m=-pnucl->getJ_array()[shellindex];m<=pnucl->getJ_array()[shellindex];m+=2){
	Matrix<2,4> J;
	Matrix<2,4> Vector;
	Matrix<2,4> Axial;
	reacmodel->getMatrixEl(kin,Vector,shellindex,m,CT,pw, current, SRC, thick,1);
	reacmodel->getMatrixEl(kin,Axial,shellindex,m,CT,pw, current, SRC, thick,0);
	J=Vector+Axial;
	for(int i=0;i<1;i++){ //exploit parity symmetries, only one polarization needed
//  	  cout << i << " " << m << " " << J(i,0) << " " << J(i,1) << " " <<  J(i,2) << " " << J(i,3) << endl;
	  response[0]+=2.*(norm(Vector(i,0)-kin.GetWlab()/qvec*Vector(i,3))+
		  norm(Axial(i,0)-kin.GetWlab()/qvec*Axial(i,3))); //W_L
	  response[1]+=2.*(norm(Vector(i,1))+norm(Vector(i,2))+norm(Axial(i,1))+norm(Axial(i,2))); //W_T
	  response[4]+=4.*real(conj(Vector(i,2))*Axial(i,2)-conj(Vector(i,1))*Axial(i,1)); //W_T'
	  if(!phi_int){
	    response[2]+=4.*real(Vector(i,2)*conj(Vector(i,1))+Axial(i,2)*conj(Axial(i,1))); //W_TT
	    response[3]+=4.*real((Vector(i,0)-kin.GetWlab()/qvec*Vector(i,3))*conj(Vector(i,2)-Vector(i,1))+
		(Axial(i,0)-kin.GetWlab()/qvec*Axial(i,3))*conj(Axial(i,2)-Axial(i,1))); //W_LT
	    response[5]+=4.*imag((Vector(i,0)-kin.GetWlab()/qvec*Vector(i,3))*conj(Vector(i,2)-Vector(i,1))
	      +(Axial(i,0)-kin.GetWlab()/qvec*Axial(i,3))*conj(Axial(i,2)-Axial(i,1))); //W_LT'
	    extraresponse[0]+=-4.*imag(Vector(i,2)*conj(Axial(i,1))+Axial(i,2)*conj(Vector(i,1)));
	    extraresponse[1]+=4.*imag((Vector(i,0)-kin.GetWlab()/qvec*Vector(i,3))*conj(Axial(i,2)+Axial(i,1))
	      +(Axial(i,0)-kin.GetWlab()/qvec*Axial(i,3))*conj(Vector(i,2)+Vector(i,1)));
	    extraresponse[2]+=4.*real((Vector(i,0)-kin.GetWlab()/qvec*Vector(i,3))*conj(Axial(i,2)+Axial(i,1))
	      +(Axial(i,0)-kin.GetWlab()/qvec*Axial(i,3))*conj(Vector(i,2)+Vector(i,1)));
	  }
	}
    }
//     for(int i=0;i<6;i++) cout << kinfactors[i] << " " << response[i] << endl;
//     for(int i=0;i<3;i++) cout << extraresponse[i] << endl;
    double result=0.;
    //combine everything
   if(!phi_int) result=kinfactors[0]*response[0]+kinfactors[1]*response[1]
			+kinfactors[2]*(response[2]*cos(2.*phi)+extraresponse[0]*sin(2.*phi))
			+kinfactors[3]*(response[3]*cos(phi)+extraresponse[1]*sin(phi))
	      +(shellindex<pnucl->getPLevels()?1.:-1.)*(kinfactors[4]*response[4]
							+kinfactors[5]*(response[5]*sin(phi)+extraresponse[2]*cos(phi)));
    else result=2.*PI*(kinfactors[0]*response[0]+kinfactors[1]*response[1]
      +(shellindex<pnucl->getPLevels()?1.:-1.)*kinfactors[4]*response[4]);
    delete reacmodel;
    return mott*frontfactor*result/HBARC;
  }  
}

// void  WeakQECross::getAllDiffWeakQECross(std::vector<double> &cross, TKinematics2to2 &kin, int current, 
// 			     int shellindex, int thick, double phi, int maxEval, bool lab, bool phi_int, bool neutrino){
//   
//   cross=vector<double>(5,0.);
//   
// }


