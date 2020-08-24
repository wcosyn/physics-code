#include <ElementaryUtils/logger/CustomException.h>
#include <ElementaryUtils/logger/LoggerManager.h>
#include <partons/Partons.h>
#include <partons/services/automation/AutomationService.h>
#include <partons/modules/running_alpha_strong/RunningAlphaStrongStandard.h>
#include <partons/ServiceObjectRegistry.h>
#include <QtCore/qcoreapplication.h>
#include <string>
#include <vector>


#include <partons/ModuleObjectFactory.h>
#include <partons/services/GPDService.h>
#include <partons/services/ObservableService.h>
#include <partons/ServiceObjectRegistry.h>

#include <partons/modules/gpd/GPDMMS13.h>
#include <partons/modules/gpd/GPDVGG99.h>
#include <partons/modules/gpd/GPDGK16.h>
#include <partons/modules/gpd/GPDGK16Numerical.h>


#include <TwoVector_Nucl.hpp>
#include <constants.hpp>


using namespace std;


   /**
     * @brief list of possible GPD models
     * 
     */
    enum GPD_model_type { GK16Numerical,MMS13,MPSSW13,VGG99,Vinnikov06 };
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
    int meson_type = atoi(argv[2]);
    TwoVector_Nucl::Rho_pol kmeson;
    switch(meson_type){
        case 0:
            kmeson=TwoVector_Nucl::krhoT;
            break;
        case 1:
            kmeson=TwoVector_Nucl::krhoL;
            break;
        case 2:
            kmeson=TwoVector_Nucl::kaxial;
            break;
        default:
        ;
    }

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
        TwoVector_Nucl gimme_xs(pGPDService, pGPDModel, pRunningAlphaStrongModule);

        double masse2=pow(511.E-06,2.);
        double Q2max=atof(argv[3]);
        double sum=0.;
        for(int i=100;i<1000;i++){
            double y=i/1000.;
            sum+=ALPHA/2./PI*(2.*masse2*y*(1./Q2max-(1-y)/masse2/y/y)+(((1-y)*(1-y)+1)*log(Q2max*(1-y)/masse2/y/y)/y));
        }
        cout << sum*1.04692/2./1000. << endl;
        exit(1);


        vector< double > resultsL(12,0.), resultsT(12,0.), results_total(12,0.);
        double xi=0.15;
        double psq=4.;
        for(int j=-24; j<=-23; j++){
            double Q2=pow(10.,j/20.*5.+0.125);
            gimme_xs.getCross_twovector(resultsL, 1.,xi,Q2,psq, TwoVector_Nucl::kgammaL, TwoVector_Nucl::krhoL, 2E05);
            gimme_xs.getCross_twovector(resultsT, 1.,xi,Q2,psq, TwoVector_Nucl::kgammaT, TwoVector_Nucl::krhoL, 1E04);
            cout << Q2 << " " << resultsL[0] << " " << resultsT[0] << endl;
            // for(int i=0;i<20;i++){
            //     double y = 0.1+0.9/20.*i;
            //     gimme_xs.getElectro_Cross_twovector(results_total, resultsL, resultsT, y,Q2);
            //     cout << Q2 << " " << y << " ";
            //     for(int kk=0;kk<12;kk++){
            //         cout << results_total[kk] << " " << resultsL[kk] << " " << resultsT[kk] << " ";
            //     }
            //     cout << endl;
            // }
        }
        // cout << endl << endl;
        // for(int index=0;index<12;index++) cout << results[index] << " ";
 


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
