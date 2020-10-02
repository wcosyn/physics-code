#include <ElementaryUtils/logger/CustomException.h>
#include <ElementaryUtils/logger/LoggerManager.h>
#include <partons/Partons.h>
#include <partons/services/automation/AutomationService.h>
#include <partons/modules/running_alpha_strong/RunningAlphaStrongStandard.h>
#include <partons/ServiceObjectRegistry.h>
// #include <QtCore/qcoreapplication.h>
#include <string>
#include <vector>


#include <partons/ModuleObjectFactory.h>
#include <partons/services/GPDService.h>
#include <partons/services/ObservableService.h>
#include <partons/ServiceObjectRegistry.h>

#include <partons/modules/gpd/GPDMMS13.h>
#include <partons/modules/gpd/GPDVGG99.h>
#include <partons/modules/gpd/GPDGK16.h>
#include <partons/modules/gpd/GPDGK19.h>
#include <partons/modules/gpd/GPDGK16Numerical.h>


#include <DiDVCS.hpp>
#include <constants.hpp>


using namespace std;



#include <TVector3.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <TF2.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>

double me=511e-6;
double Mp=0.938272;
double mpi=0.13957;
double MassRho=0.77545;
double GammaRho=0.145;

double s1_min;
double s1_max;
double th_gen;
double mrho;

TLorentzVector BeamL;
TLorentzVector BeamH;
TLorentzVector Rho0,eplus,eminus,piplus,piminus;
TLorentzVector QPrime;
TLorentzVector Proton;
TLorentzVector ScattL;
TLorentzVector q;
TLorentzVector cms;
TLorentzVector S2;


double s;

TF1* virt_long=new TF1("virt long","TMath::Power(TMath::Sin(x),2)",0,2*TMath::Pi());//Angular dependence for the qprime decay
TF1* rho_long=new TF1("rho long","TMath::Power(TMath::Cos(x),2)",0,2*TMath::Pi());//Angular dependence for the rho0 decay

double theta_e,Ep, t_min;
TRandom3 *rando=new TRandom3();

double xb;

double Ecm,pcm; //Energy/momentum in center-of-mass of (rho+q')
double Ecm_pr,pcm_pr;  //Energy/momentum in center-of-mass of (p')
double Ecm_rho,pcm_rho; //Energy/momentum in center-of-mass of (rho)

TVector3 dir_pr;
TVector3 dir_qp;
TVector3 rot_dir;

//Method to compute the lower and upper bound of the squared momentum transfer to the recoil baryon (Proton, Nucleus, Delta...). It depends on q2, xb, the mass of the recoil baryon mrecoil and the mass of the produced particle mpp 
Double_t tmin(Double_t q, Double_t x, double mrecoil, double mpp){
 Double_t W=TMath::Sqrt(q*(1/x-1)+Mp*Mp);
 Double_t E1cm=(W*W+q+Mp*Mp)/2./W;
 Double_t P1cm=TMath::Sqrt(E1cm*E1cm-Mp*Mp);
 Double_t E3cm=(W*W-mpp*mpp+mrecoil*mrecoil)/2./W;
 Double_t P3cm=TMath::Sqrt(E3cm*E3cm-mrecoil*mrecoil);
 Double_t result=(TMath::Power(q+mpp*mpp+Mp*Mp-mrecoil*mrecoil,2)/4./W/W-TMath::Power(P1cm-P3cm,2));
 return result;
}

Double_t tmax(Double_t q, Double_t x, double mrecoil, double mpp){
 Double_t W=TMath::Sqrt(q*(1/x-1)+Mp*Mp);
 Double_t E1cm=(W*W+q+Mp*Mp)/2./W;
 Double_t P1cm=TMath::Sqrt(E1cm*E1cm-Mp*Mp);
 Double_t E3cm=(W*W-mpp*mpp+mrecoil*mrecoil)/2./W;
 Double_t P3cm=TMath::Sqrt(E3cm*E3cm-mrecoil*mrecoil);
 Double_t result=(TMath::Power(q+mpp*mpp+Mp*Mp-mrecoil*mrecoil,2)/4./W/W-TMath::Power(P1cm+P3cm,2));
 return result;
}


//THE method which computes all LorentzVector. Q2 and xb has already been set and are kept in memory. Then remains t and phi to define the entire final set
//In truth... it must be validated on the recoil baryon
void ComputeLeptScatt(double q2, double y){
  theta_e=2*TMath::ATan(TMath::Sqrt(q2/(1-y)/4./BeamL.E()/BeamL.E()));
  Ep=q2/4./BeamL.E()/TMath::Sin(theta_e/2.)/TMath::Sin(theta_e/2.);
 ScattL.SetPxPyPzE(Ep*TMath::Sin(TMath::Pi()+theta_e),0,Ep*TMath::Cos(TMath::Pi()+theta_e),Ep);
 
  q=BeamL-ScattL;
  cms=q+BeamH;
  //q2=sxy is approximative... need to readjust the value of xb for consistency
  xb=-q.Mag2()/2./(BeamH*q);

   //The collision is easier to understand in the center-of-mass frame
   ScattL.Boost(-cms.BoostVector());
   BeamH.Boost(-cms.BoostVector());
   BeamL.Boost(-cms.BoostVector());
  q.Boost(-cms.BoostVector());
}


void ComputePrQprime(double s2, double q2prime, double t_N, double phi_trento){
  double Es2_pr=(s2-q2prime+Mp*Mp)/2./TMath::Sqrt(s2);//Energy in cms of s2 for the proton
  double Es2_qp=(s2+q2prime-Mp*Mp)/2./TMath::Sqrt(s2);//Energy in cms of s2 for the q2'
  double ps2=TMath::Sqrt(((s2-pow(Mp+sqrt(q2prime),2))*(s2-pow(Mp-sqrt(q2prime),2))))/2./TMath::Sqrt(s2);//momentum of both particle in particle s2
  
  //Let's go to S2 cms
  BeamH.Boost(-S2.BoostVector());

  double theta=acos(((t_N-2*Mp*Mp)/2.+BeamH.E()*Es2_pr)/BeamH.Vect().Mag()/ps2);
 
  TVector3 unit=BeamH.Vect().Unit();
  unit.Rotate(theta,unit.Orthogonal());
  unit.Rotate(phi_trento,BeamH.Vect());
  Proton.SetPxPyPzE(ps2*unit.X(),ps2*unit.Y(),ps2*unit.Z(),Es2_pr);
  QPrime.SetPxPyPzE(-Proton.Px(),-Proton.Py(),-Proton.Pz(),Es2_qp);

  //Back to gamma*p cms frame
  BeamH.Boost(S2.BoostVector());
  Proton.Boost(S2.BoostVector());
  QPrime.Boost(S2.BoostVector());
  
  //Back to Lab frame
  BeamH.Boost(cms.BoostVector());
  Proton.Boost(cms.BoostVector());
  QPrime.Boost(cms.BoostVector());
  Rho0.Boost(cms.BoostVector());
  ScattL.Boost(cms.BoostVector());
  BeamL.Boost(cms.BoostVector());
  q.Boost(cms.BoostVector());
 }

void ComputeRho(double mpp, double mrecoil, double t_rho, double phi_trento){
  Ecm=(cms.Mag2()-mrecoil*mrecoil+mpp*mpp)/2./cms.Mag();
  pcm=TMath::Sqrt(Ecm*Ecm-mpp*mpp);
  
  Ecm_pr=(cms.Mag2()+mrecoil*mrecoil-mpp*mpp)/2./cms.Mag();
  pcm_pr=TMath::Sqrt(Ecm_pr*Ecm_pr-mrecoil*mrecoil);


   //t=2q_cm q'_cm(1-cos (theta_gamma*/meson))
  
  double theta_cm=acos(1-(tmin(-q.Mag2(), xb, mrecoil, mpp)-t_rho)/2./pcm/q.Vect().Mag());
   //Phi is a rotation of all four momentum around the axis defined by the virtual photon.
  TVector3 dir_pr=q.Vect().Unit();
  TVector3 dir_qp=q.Vect().Unit();
  TVector3 rot_dir=q.Vect().Orthogonal().Unit();
   dir_pr.Rotate(theta_cm,rot_dir);
   dir_pr.Rotate(phi_trento,q.Vect().Unit());
    Rho0.SetPxPyPzE(pcm*dir_pr.X(),pcm*dir_pr.Y(),pcm*dir_pr.Z(),Ecm);

    S2=q+BeamH-Rho0;
}

void GammaPrimeDecay(){
  double th_gen=virt_long->GetRandom();
  double pel=TMath::Sqrt(QPrime.Mag2()/4.-me*me);
  TVector3 unit_dir=QPrime.Vect().Unit();
  unit_dir.Rotate(th_gen,unit_dir.Orthogonal()); //theta angle with photon direction
  unit_dir.Rotate(rando->Uniform(0,2*TMath::Pi()),QPrime.Vect()); //rotation around the photon direction
  eminus.SetPxPyPzE(pel*unit_dir.X(),pel*unit_dir.Y(),pel*unit_dir.Z(),QPrime.Mag()/2.);
  eplus.SetPxPyPzE(-pel*unit_dir.X(),-pel*unit_dir.Y(),-pel*unit_dir.Z(),QPrime.Mag()/2.);

  //BackToLabFrame
  eminus.Boost(QPrime.BoostVector());
   eplus.Boost(QPrime.BoostVector());
  
}

void RhoDecay(){
  double th_gen_rho=rho_long->GetRandom();
  double ppi=TMath::Sqrt(mrho*mrho/4.-mpi*mpi);
  TVector3 unit_dir=Rho0.Vect().Unit();
  unit_dir.Rotate(th_gen_rho,unit_dir.Orthogonal()); //theta angle with rho direction
  unit_dir.Rotate(rando->Uniform(0,2*TMath::Pi()),Rho0.Vect()); //rotation around the rho direction
  piminus.SetPxPyPzE(ppi*unit_dir.X(),ppi*unit_dir.Y(),ppi*unit_dir.Z(),mrho/2.);
  piplus.SetPxPyPzE(-ppi*unit_dir.X(),-ppi*unit_dir.Y(),-ppi*unit_dir.Z(),mrho/2.);

  //BackToLabFrame
  piminus.Boost(Rho0.BoostVector());
   piplus.Boost(Rho0.BoostVector());
}

//This function is the main function...
//ScattL is the scattered lepton
//Proton is the recoiled proton
//eplus and eminus are the decay product of q'
//piplus piminus are the decay product of rho0
void GetParticles(double EneL, double EneH, double q2, double y, double q2prime, double s2, double t_N, double t_rho){
  //Set Beam particle
  BeamL.SetPxPyPzE(0,0,-EneL,EneL);
  BeamH.SetPxPyPzE(0,0,EneH,TMath::Sqrt(Mp*Mp+EneH*EneH));
  s=(BeamL+BeamH).Mag2();
  ComputeLeptScatt(q2,y);//Compute ScattL and xb... all 4 momenta are exprressed in gamma*p cms.
  
  while (mrho<0.4||mrho>1.2) mrho=rando->BreitWigner(MassRho,GammaRho); 
  
  ComputeRho(mrho,TMath::Sqrt(s2),t_rho,0); //Give rho 4-momentum and last argument is the angle around the virtual photon in gamma*p cms.
  
  ComputePrQprime(s2, q2prime, t_N, 3.141); //Give rho 4-momenta of p' and gamma'. last argument is the angle of p' gamma' around the incoming hadron in gamma*p cms. ALL THE 4-MOMENTA ARE NOW EXPRESSED in LAB FRAME.

  GammaPrimeDecay();
  RhoDecay();

  /*cout<<t_N<<" "<<(Proton-BeamH).Mag2()<<endl;
  cout<<t_rho<<" "<<(Rho0-q).Mag2()<<endl;
  cout<<q2prime<<" "<<QPrime.Mag2()<<" "<<(eplus+eminus).Mag2()<<endl;
  cout<<s2<<" "<<(QPrime+Proton).Mag2()<<endl;
  cout<<q2<<" "<<-(BeamL-ScattL).Mag2()<<endl;
  ScattL.Print();
  eplus.Print();
  eminus.Print();
  piplus.Print();
  piminus.Print();
  Proton.Print();*/
 

  return;
}



   /**
     * @brief list of possible GPD models
     * 
     */
    enum GPD_model_type { GK16Numerical,GK16,GK19,MMS13,MPSSW13,VGG99,Vinnikov06 };
    /**
     * @brief a string to type map so you can lookup int values using strings like "VGG99" and so on... 
     * use in constructor as Twovector_Nucl(Twovector_Nucl::TypeNames.at("VGG99"); (don't use [] acces op.)
     * 
     */
    static  std::map<std::string,GPD_model_type> TypeNames; 
    /**
     * @brief initialise typename lookup map, use in constructor as Twovector_Nucl(Twovector_Nucl::TypeNames.at("VGG99"); (don't use [] acces op.) 
     * 
     * @return std::map<std::string,GPD_model_type> 
     */
    static std::map<std::string,GPD_model_type> initTypeNames(){ std::map<std::string,GPD_model_type> m; 
                                                                    m["GK16Numerical"]=GK16Numerical;
                                                                    m["GK16"]=GK16;
                                                                    m["GK19"]=GK19;
                                                                    m["MMS13"]=MMS13;
                                                                    m["MPSSW13"]=MPSSW13;
                                                                    m["VGG99"]=VGG99;
                                                                    m["Vinnikov06"]=Vinnikov06;return m;} 




/*
 * Main function.
 */
int main(int argc, char** argv) {


    //double xi=atof(argv[1]);
    TypeNames=initTypeNames();
    //cout << argv[1] << endl;
    GPD_model_type gpdmodel = TypeNames.at(argv[1]);
    double Q2in = atof(argv[2]); // virtual photon scale squared [GeV^2]
    double Q2out = atof(argv[3]); // virtual photon scale squared [GeV^2]
    double t_N = atof(argv[4]);
    double t_rho = atof(argv[5]);
    double El_beam = atof(argv[6]);
    double p_beam = atof(argv[7]);    
    double s2 = atof(argv[8]);
    double y = atof(argv[9]);
        

    // Init Qt4
    //QCoreApplication a(argc, argv);
    PARTONS::Partons* pPartons = 0;
    
 
    try {

        // Init PARTONS application
        pPartons = PARTONS::Partons::getInstance();
        pPartons->init(argc, argv);

        // Retrieve GPD service
        PARTONS::GPDService* pGPDService =
                PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getGPDService();

        // Create GPD module with the BaseModuleFactory
        PARTONS::GPDModule* pGPDModel;
        
        switch(gpdmodel){
            case(GK19):
                pGPDModel=PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                        PARTONS::GPDGK19::classId);
                break;
            case(GK16):
                pGPDModel=PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                        PARTONS::GPDGK16::classId);
                break;
            case(GK16Numerical):
                pGPDModel=PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                        PARTONS::GPDGK16Numerical::classId);
                break;
            case(MMS13):
                pGPDModel=PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                        PARTONS::GPDMMS13::classId);
                break;
            case(VGG99):
                pGPDModel=PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                        PARTONS::GPDVGG99::classId);
                break;
            default:
                cerr << "Invalid GPD model name chosen " << endl;
                exit(1);
        } 
                
        PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule = PARTONS::Partons::getInstance()->getModuleObjectFactory()->newRunningAlphaStrongModule(
                    PARTONS::RunningAlphaStrongStandard::classId);

        //cout << "alpha_s " << alpha_s << endl;
        DiDVCS gimme_xs(pGPDService, pGPDModel, pRunningAlphaStrongModule);
        vector< double > results(2,0.);

        double Ep=sqrt(MASSP*MASSP+p_beam*p_beam);
        double s_eN = MASSP*MASSP + 2.*El_beam*(Ep+p_beam);

        gimme_xs.getCross_DiDVCS(results, 2., t_rho, t_N, Q2in, Q2out, s2, s_eN, y, 1.E04);
        GetParticles(El_beam,p_beam,Q2in,y,Q2out,s2,t_N,t_rho);

        cout << results[0] << " " << results[1] << endl;
        Rho0.Print();
        QPrime.Print();
        ScattL.Print();
        Proton.Print();
        


        // Remove pointer references
        // Module pointers are managed by PARTONS
        PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
                pGPDModel, 0);
        pGPDModel = 0;
        

    }
    // Appropriate catching of exceptions is crucial for working of PARTONS.
    // PARTONS defines its own type of exception, which allows to display class name and function name
    // where the exception has occurred, but also a human readable explanation.
    catch (const ElemUtils::CustomException &e) {

        // Display what happened
        pPartons->getLoggerManager()->error(e);

        // Close PARTONS application properly
        if (pPartons) {
            pPartons->close();
        }
    }
    // In a case of standard exception.
    catch (const std::exception &e) {

        // Display what happened
        pPartons->getLoggerManager()->error("main", __func__, e.what());

        // Close PARTONS application properly
        if (pPartons) {
            pPartons->close();
        }
    }

    // Close PARTONS application properly
    if (pPartons) {
        pPartons->close();
    }

    return 0;
}
