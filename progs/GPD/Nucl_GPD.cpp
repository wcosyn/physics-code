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
#include <Deut_Conv_GPD_T.hpp>
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
        PARTONS::GPDModule *pGK16, *pMMS13, *pVGG99;
        
        pGK16=PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                        PARTONS::GPDGK16Numerical::classId);
        pMMS13=PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                        PARTONS::GPDMMS13::classId);
        pVGG99=PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(
                        PARTONS::GPDVGG99::classId);
                
        PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule = PARTONS::Partons::getInstance()->getModuleObjectFactory()->newRunningAlphaStrongModule(
                    PARTONS::RunningAlphaStrongStandard::classId);


//        Testing chiral even GPDs
        for(int i=0;i<=30;i++){
            double x=0.1;
            double xi=0.1+i*0.01;
            double mandelstam_t = -0.45;//-4.*MASSP_G*MASSP_G*xi*xi/(1-xi*xi);
            double scale = 1;
            cout << xi << " " << mandelstam_t << " ";
            PARTONS::GPDKinematic gpdKinematic(x,xi,mandelstam_t, scale*scale, scale*scale);
            PARTONS::GPDKinematic gpdKinematic_min(-x,xi,mandelstam_t, scale*scale, scale*scale);
            PARTONS::GPDResult gpdResult_GK16 = pGPDService->computeGPDModel(gpdKinematic,pGK16);
            PARTONS::GPDResult gpdResult_min_GK16 = pGPDService->computeGPDModel(gpdKinematic_min,pGK16);
            double Etu_GK16 = gpdResult_GK16.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()+
                        gpdResult_min_GK16.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
            double Etd_GK16 = gpdResult_GK16.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution()+
                        gpdResult_min_GK16.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();

            cout  << Etu_GK16 << " " << Etd_GK16 << " " <<  Etu_GK16-Etd_GK16 << endl; 

        }
        exit(1);


        for(int i=0;i<=2;i++){
            double xi=0.1+i*0.2;
            double mandelstam_t= -4.*MASSP_G*MASSP_G*xi*xi/(1-xi*xi);//tmin [GeV^2]
            double scale=1.;
            for(int i =0;i<=50;i++){
                double x=i*xi/50.;
                PARTONS::GPDKinematic gpdKinematic(x,xi,mandelstam_t, scale*scale, scale*scale);
                PARTONS::GPDKinematic gpdKinematic_min(-x,xi,mandelstam_t, scale*scale, scale*scale);

                PARTONS::GPDResult gpdResult_GK16 = pGPDService->computeGPDModel(gpdKinematic,pGK16);
                PARTONS::GPDResult gpdResult_min_GK16 = pGPDService->computeGPDModel(gpdKinematic_min,pGK16);
                PARTONS::GPDResult gpdResult_MMS13 = pGPDService->computeGPDModel(gpdKinematic,pMMS13);
                PARTONS::GPDResult gpdResult_min_MMS13 = pGPDService->computeGPDModel(gpdKinematic_min,pMMS13);
                PARTONS::GPDResult gpdResult_VGG99 = pGPDService->computeGPDModel(gpdKinematic,pVGG99);
                PARTONS::GPDResult gpdResult_min_VGG99 = pGPDService->computeGPDModel(gpdKinematic_min,pVGG99);

                GPD_T_Nucl_grid gpdTgrid("MSTW");
                TransGPD_set gpdT = gpdTgrid.getGK_param(x,xi,mandelstam_t*1.E06, scale);
                TransGPD_set gpdT_min = gpdTgrid.getGK_param(-x,xi,mandelstam_t*1.E06, scale);

                double Hu_GK16 = gpdResult_GK16.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()+
                            gpdResult_min_GK16.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
                double Hd_GK16 = gpdResult_GK16.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution()+
                            gpdResult_min_GK16.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
                double Eu_GK16 = gpdResult_GK16.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()+
                            gpdResult_min_GK16.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
                double Ed_GK16 = gpdResult_GK16.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution()+
                            gpdResult_min_GK16.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
                double Htu_GK16 = gpdResult_GK16.getPartonDistribution(PARTONS::GPDType::Ht).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()+
                            gpdResult_min_GK16.getPartonDistribution(PARTONS::GPDType::Ht).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
                double Htd_GK16 = gpdResult_GK16.getPartonDistribution(PARTONS::GPDType::Ht).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution()+
                            gpdResult_min_GK16.getPartonDistribution(PARTONS::GPDType::Ht).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
                double Etu_GK16 = gpdResult_GK16.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()+
                            gpdResult_min_GK16.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
                double Etd_GK16 = gpdResult_GK16.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution()+
                            gpdResult_min_GK16.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();

                double Hu_VGG99 = gpdResult_VGG99.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()+
                            gpdResult_min_VGG99.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
                double Hd_VGG99 = gpdResult_VGG99.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution()+
                            gpdResult_min_VGG99.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
                double Eu_VGG99 = gpdResult_VGG99.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()+
                            gpdResult_min_VGG99.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
                double Ed_VGG99 = gpdResult_VGG99.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution()+
                            gpdResult_min_VGG99.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
                double Htu_VGG99 = gpdResult_VGG99.getPartonDistribution(PARTONS::GPDType::Ht).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()+
                            gpdResult_min_VGG99.getPartonDistribution(PARTONS::GPDType::Ht).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
                double Htd_VGG99 = gpdResult_VGG99.getPartonDistribution(PARTONS::GPDType::Ht).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution()+
                            gpdResult_min_VGG99.getPartonDistribution(PARTONS::GPDType::Ht).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
                double Etu_VGG99 = gpdResult_VGG99.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()+
                            gpdResult_min_VGG99.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
                double Etd_VGG99 = gpdResult_VGG99.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution()+
                            gpdResult_min_VGG99.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();

                double Hu_MMS13 = gpdResult_MMS13.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()+
                            gpdResult_min_MMS13.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
                double Hd_MMS13 = gpdResult_MMS13.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution()+
                            gpdResult_min_MMS13.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
                double Eu_MMS13 = gpdResult_MMS13.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()+
                            gpdResult_min_MMS13.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
                double Ed_MMS13 = gpdResult_MMS13.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution()+
                            gpdResult_min_MMS13.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();

                double HT_iv = 2.*(gpdT.getHT_vector()+gpdT_min.getHT_vector());
                double HT_is = 2.*(gpdT.getHT_singlet()+gpdT_min.getHT_singlet());
                double ET_iv_0 = 2.*(gpdT.getET_vector(0)+gpdT_min.getET_vector(0));
                double ET_is_0 = 2.*(gpdT.getET_singlet(0)+gpdT_min.getET_singlet(0));
                double ET_iv_1 = 2.*(gpdT.getET_vector(1)+gpdT_min.getET_vector(1));
                double ET_is_1 = 2.*(gpdT.getET_singlet(1)+gpdT_min.getET_singlet(1));
                double ET_iv_2 = 2.*(gpdT.getET_vector(2)+gpdT_min.getET_vector(2));
                double ET_is_2 = 2.*(gpdT.getET_singlet(2)+gpdT_min.getET_singlet(2));

                cout << xi << " " << x << " " << Hu_GK16 << " " << Hd_GK16 << " " << Hu_VGG99 << " " << Hd_VGG99 << " " << Hu_MMS13 << " " << Hd_MMS13 << " " 
                                            << Eu_GK16 << " " << Ed_GK16 << " " << Eu_VGG99 << " " << Ed_VGG99 << " " << Eu_MMS13 << " " << Ed_MMS13 << " " 
                                            << Htu_GK16 << " " << Htd_GK16 << " " << Htu_VGG99 << " " << Htd_VGG99 << " " 
                                            << Etu_GK16 << " " << Etd_GK16 << " " << Etu_VGG99 << " " << Etd_VGG99 << " " 
                                            << HT_iv << " " << HT_is << " " << ET_iv_0 << " " << ET_is_0 << " " << ET_iv_1 << " " 
                                            << ET_is_1 << " " << ET_iv_2 << " " << ET_is_2 << " " << 
                                            Deut_Conv_GPD_T::getGPD_odd_nucl(-1,1,xi,mandelstam_t*1.06,mandelstam_t*1.06,0.,0,1,1,gpdT+gpdT_min).real()*2. << " "<<
                                            Deut_Conv_GPD_T::getGPD_odd_nucl(-1,1,xi,mandelstam_t*1.06,mandelstam_t*1.06,0.,1,1,1,gpdT+gpdT_min).real()*2. << " " <<
                                            Deut_Conv_GPD_T::getGPD_odd_nucl(-1,1,xi,mandelstam_t*1.06,mandelstam_t*1.06,0.,2,1,1,gpdT+gpdT_min).real()*2. << endl;
            }
            cout << endl << endl;
        }
        // Remove pointer references
        // Module pointers are managed by PARTONS
        PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
                pGK16, 0);
        pGK16 = 0;
        PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
                pMMS13, 0);
        pMMS13 = 0;
        PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(
                pVGG99, 0);
        pVGG99 = 0;
        

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
