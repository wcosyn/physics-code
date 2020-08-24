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
    double Qsq = atof(argv[2]); // virtual photon scale squared [GeV^2]
    int meson_type = atoi(argv[3]);
    bool gammaT = atoi(argv[4]);
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

        //Testing chiral even GPDs

        // double mandelstam_t= -4.*MASSP_G*MASSP_G*xi*xi/(1-xi*xi);//tmin [GeV^2]
        // double scale=1.;
        // for(int i =-50;i<=50;i++){
        //     double x=i*0.01;
        //     PARTONS::GPDKinematic gpdKinematic(x,xi,mandelstam_t, scale*scale, scale*scale);

        //     PARTONS::GPDResult gpdResult = pGPDService->computeGPDModel(gpdKinematic,pGPDModel);
        // // (Hq^u-Hq^d) * z^2 * barz^2 / u / baru * P [z*barz,u*baru from 2 DA taken into account, factors 6 afterwards]
        //     cout << x << " " << xi << " " << gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
        //     << " " << gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution() <<
        //     " " << gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()-
        //     gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution() << " " << 
        //     " " << gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
        //     << " " << gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution() <<
        //     " " << gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()-
        //     gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution() << " " << endl;

        //     }
        // exit(1);

        //xi dependence, model comparisons
        // for(int j=0;j<=50;j++){
        //     double xi=0.1+0.01*j;
        //     for(int i=0;i<=3;i++){
            
        //         double psq = 2.+i*4; //pomeron scale squared [GeV^2]
        //         cout << Qsq << " " << xi << " " << psq << " "; 
        //         // double result_L=0.,result_T=0.;
        //         // if(meson_type){
        //         //     result_L=gimme_xs.getCross_gammaL_rhoL(1.,xi,Qsq,psq);
        //         //     result_T=gimme_xs.getCross_gammaT_rhoL(1.,xi,Qsq,psq);
        //         // }
        //         // else{
        //         //     result_L=gimme_xs.getCross_gammaL_rhoT(1.,xi,Qsq,psq,1);
        //         //     result_T=gimme_xs.getCross_gammaT_rhoT(1.,xi,Qsq,psq,1);
        //         // }
        //         // cout << result_T << " " << result_L << endl;


        //         vector< double > results(12,0.);
        //         gimme_xs.getCross_twovector(results, 1.,xi,Qsq,psq, gammaT? TwoVector_Nucl::kgammaT : TwoVector_Nucl::kgammaL, kmeson, 1E04);
        //         for(int index=0;index<12;index++) cout << results[index] << " ";
        //         cout << endl;
        //     }
        //     cout << endl << endl;
        // }


        //pt2 dependence
        double xi=0.3;
        for(int i=0;i<=18;i++){
        
            double psq = 2.+i*0.5; //pomeron scale squared [GeV^2]
            cout << Qsq << " " << xi << " " << psq << " "; 
            // double result_L=0.,result_T=0.;
            // if(meson_type){
            //     result_L=gimme_xs.getCross_gammaL_rhoL(1.,xi,Qsq,psq);
            //     result_T=gimme_xs.getCross_gammaT_rhoL(1.,xi,Qsq,psq);
            // }
            // else{
            //     result_L=gimme_xs.getCross_gammaL_rhoT(1.,xi,Qsq,psq,1);
            //     result_T=gimme_xs.getCross_gammaT_rhoT(1.,xi,Qsq,psq,1);
            // }
            // cout << result_T << " " << result_L << endl;


            vector< double > results(12,0.);
            gimme_xs.getCross_twovector(results, 1.,xi,Qsq,psq, gammaT? TwoVector_Nucl::kgammaT : TwoVector_Nucl::kgammaL, kmeson, 1E04);
            for(int index=0;index<12;index++) cout << results[index] << " ";
            cout << endl;
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
