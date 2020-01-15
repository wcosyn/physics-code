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


#include <TwoVector_Deut.hpp>

#include <locale.h>

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
    static const std::map<std::string,GPD_model_type> TypeNames; 
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


    double xi=atof(argv[1]);
    double Qsq = atof(argv[2]); // virtual photon scale squared [GeV^2]
    bool longitudinal = atoi(argv[3]);
    std::string wfname = argv[4];
    std::string pdfname = argv[5];
    
    // Init Qt4
    //QCoreApplication a(argc, argv);
    //std::locale::global( std::locale( "" ) );
    //setlocale(LC_NUMERIC,"en_US");
    PARTONS::Partons* pPartons = 0;
     
    try {

        // Init PARTONS application
        pPartons = PARTONS::Partons::getInstance();
        pPartons->init(argc, argv);

        // Retrieve GPD service
        PARTONS::GPDService* pGPDService =
                PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getGPDService();

        // Create GPD module with the BaseModuleFactory
        PARTONS::GPDModule* pGPDModel =
                PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                        PARTONS::GPDMMS13::classId);
        PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule = PARTONS::Partons::getInstance()->getModuleObjectFactory()->newRunningAlphaStrongModule(
                    PARTONS::RunningAlphaStrongStandard::classId);

        //cout << "alpha_s " << alpha_s << endl;
        TwoVector_Deut gimme_xs(pGPDService, pGPDModel, pRunningAlphaStrongModule, wfname, pdfname);

        for(int i=0;i<=0;i++){
        
            double psq = 2.+i*0.5; //pomeron scale squared [GeV^2]
            double result_L=0.,result_T=0.;
            if(longitudinal){
                result_L=gimme_xs.getCross_gammaL_rhoL(1.,xi,Qsq,psq);
                result_T=gimme_xs.getCross_gammaT_rhoL(1.,xi,Qsq,psq);
            }
            else{
                result_L=gimme_xs.getCross_gammaL_rhoT(1.,xi,Qsq,psq);
                result_T=gimme_xs.getCross_gammaT_rhoT(1.,xi,Qsq,psq);
            }

 
            cout << Qsq << " " << xi << " " << psq << " " << result_T << " " << result_L << endl;

        }
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
