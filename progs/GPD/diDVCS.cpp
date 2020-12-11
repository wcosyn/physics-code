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
#include <TwoVector_Nucl.hpp>
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
    //double Q2out = atof(argv[3]); // virtual photon scale squared [GeV^2]
    double t_N = atof(argv[3]);  //NEGATIVE!!!
    double t_rho = atof(argv[4]); //NEGATIVE!!!
    double El_beam = atof(argv[5]);
    double p_beam = atof(argv[6]);    

    if(t_N >0.|| t_rho > 0.) {cerr << "Momentum transfer should be negative!  Exiting.." << endl; exit(1);}

    //we loop over s2 below
    //double s2 = atof(argv[7]);

    //we looped over y below when interested in electroproduction rates
    //double y = atof(argv[9]);
        

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
        TwoVector_Nucl gimme_xs_2meson(pGPDService, pGPDModel, pRunningAlphaStrongModule);
        vector< double > results(2,0.);

        double Ep=sqrt(MASSP_G*MASSP_G+p_beam*p_beam);
        double s_eN = MASSP_G*MASSP_G + 2.*El_beam*(Ep+p_beam);

        // cout << "sq_s_eN " << sqrt(s_eN) << endl;
        // cout << "sq_s_gN " << sqrt(s_gN) << endl;
        // cout << "x " << Q2in/s_eN/y << endl;
        double ximax = sqrt(1./(1.-4.*MASSP_G*MASSP_G/t_N));
        //cout << 2.*ximax/(1+ximax)*s2 << " " << ximax << " " << s2 << endl;
        
        double s2array[5] = {6.,10.,25.,50.,100.};
        
        for(int s2index = 0; s2index <5 ; s2index++){
            double s2 = s2array[s2index];

            for (int j=0;j<1;j++){
                double y=0.7/pow(2.,j);
                double s_gN = s_eN*y;
                for (int i=0;i<10;i++){
                    double Q2out = 2.+(2.*ximax/(1+ximax)*s2-2.)/10.*i;
                    double xi_didvcs = (Q2in +s_gN*(Q2out/s2))/(2.*s_gN+Q2in -s_gN*(Q2out/s2));

                    double trho_2meson = -(2.+((2.*ximax)/(1.-ximax)*s2-2.)/10.*i);
                    double xi_2meson= (1-s2/(s2-trho_2meson))/(1+s2/(s2-trho_2meson));

                    //cout << t_N << " " << -4.*MASSP_G*MASSP_G*xi*xi/(1-xi*xi) << endl;
                    double s1_dvcs=s_gN*Q2out/s2;
                    double s1_2meson = 2.*xi_2meson*(s_gN+Q2in)/(1.+xi_2meson)-Q2in;
                    //t_N = -4.*MASSP_G*MASSP_G*xi*xi/(1-xi*xi);  //put t_n equal to tmin
                    
                    cout << y << " " << Q2in << " " << Q2in/s_eN/y /* bjorken x */ << " " << s2 << " " << t_N << " ";
                    cout << Q2out << " " << s1_dvcs/s_gN << " " << s2/s1_dvcs << " " << xi_didvcs << " ";
                    if ((t_N>-4.*MASSP_G*MASSP_G*xi_didvcs*xi_didvcs/(1-xi_didvcs*xi_didvcs)) || s2>0.5*s1_dvcs){
                        cout << "0. 0. 0." << endl;
                    }
                    else{
                        gimme_xs.getCross_DiDVCS(results, 2., t_rho, t_N, Q2in, Q2out, s2, s_eN, y, 1.E04,0);
                        // dsigma/ds2 dQout^2  dtrho dtN   | dsigma/dxi dQout^2  dtrho dtN  | dsigma/ds2 dQout^2  dtrho dtN dy dQ^2 
                        cout << results[0] << " " << results[0]*Q2out/2/xi_didvcs/xi_didvcs << " " << results[1] << endl;
                    }

                    //comparison code for twomeson production

                    // cout << trho_2meson << " " << s1_2meson/s_gN << " " << s2/s1_dvcs << " " << xi_2meson << " ";
                    // if((t_N>-4.*MASSP_G*MASSP_G*xi_2meson*xi_2meson/(1-xi_2meson*xi_2meson)) || s2>0.5*s1_2meson){
                    //     cout << "0. 0. " << endl;
                    // }
                    // else{
                    //     vector< double > resultsL_2meson(12,0.);
                    //     vector< double > resultsT_2meson(12,0.);
                    //     vector< double > resultse_2meson(12,0.);

                    //     gimme_xs_2meson.getCross_twovector_ds2(resultsL_2meson, 2., s_eN,y,Q2in,trho_2meson,t_N,s2,TwoVector_Nucl::kgammaL,TwoVector_Nucl::krhoL,1E04);
                    //     gimme_xs_2meson.getCross_twovector_ds2(resultsT_2meson, 2., s_eN,y,Q2in,trho_2meson,t_N,s2,TwoVector_Nucl::kgammaT,TwoVector_Nucl::krhoL,1E04);

                    //     gimme_xs_2meson.getElectro_Cross_twovector(resultse_2meson, resultsL_2meson, resultsT_2meson, y, Q2in);
                    //     cout << resultsL_2meson[0] << " " << resultse_2meson[0] << endl;
                    // }
                }

            }
            cout << endl << endl;
        }

        // GetParticles(El_beam,p_beam,Q2in,y,Q2out,s2,t_N,t_rho);

        // cout << results[0] << " " << results[1] << endl;
        // cout << "Rho "; Rho0.Print();
        // cout << "Qprime "; QPrime.Print();
        // cout << "eprime "; ScattL.Print();
        // cout << "proton "; Proton.Print();
        


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
