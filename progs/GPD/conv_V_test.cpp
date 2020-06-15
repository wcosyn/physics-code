#include <iostream>
#include <cstdlib>
#include <vector>

using namespace std;


#include <iostream>
#include <cmath>
#include "constants.hpp"

#include <Deut_Conv_GPD_V.hpp>
#include <GPD_V_Nucl_grid.hpp>

#include <NucleonEMOperator.hpp>

#include <ElementaryUtils/logger/CustomException.h>
#include <ElementaryUtils/logger/LoggerManager.h>
#include <partons/Partons.h>
#include <partons/services/automation/AutomationService.h>
#include <partons/modules/running_alpha_strong/RunningAlphaStrongStandard.h>
#include <partons/ServiceObjectRegistry.h>
#include <QtCore/qcoreapplication.h>
#include <string>


#include <partons/ModuleObjectFactory.h>
#include <partons/services/GPDService.h>
#include <partons/services/ObservableService.h>
#include <partons/ServiceObjectRegistry.h>

#include <partons/modules/gpd/GPDMMS13.h>
#include <partons/modules/gpd/GPDVGG99.h>
#include <partons/modules/gpd/GPDGK16.h>
#include <partons/modules/gpd/GPDGK16Numerical.h>


int main(int argc, char *argv[]){
    std::string wf=argv[1];
    double xi=atof(argv[3]);
    double t=-atof(argv[2]); //[GeV^2] positive input!
    //int ERBL = atoi(argv[4]);
    double scale = sqrt(2.);//atof(argv[5]); //[GeV]
    bool setS = atoi(argv[4]);
    bool setD = atoi(argv[5]);
    bool setH = atoi(argv[6]);
    bool setE = atoi(argv[7]);
    string GPD_model = argv[8];  // "GPDMMS13" or "GPDVGG99" or "GPDGK16" or "GPDGK16Numerical"

    // cout << "xi " << xi << " t " << t << " model " << model << endl;
    // cout << -4.*MASSD*MASSD*xi*xi/(1-xi*xi)-t << endl;


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
        PARTONS::GPDModule* pGPDModel = NULL;
        unsigned int ID;
        if(!GPD_model.compare("GPDMMS13")){
                pGPDModel = PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(PARTONS::GPDMMS13::classId);
                ID = PARTONS::GPDMMS13::classId;
                }
        else if(!GPD_model.compare("GPDVGG99")){
                pGPDModel = PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(PARTONS::GPDVGG99::classId);
                ID = PARTONS::GPDVGG99::classId;
        }
        else if(!GPD_model.compare("GPDGK16Numerical")){
                pGPDModel = PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(PARTONS::GPDGK16Numerical::classId);
                ID = PARTONS::GPDGK16Numerical::classId;
        }
        else if(!GPD_model.compare("GPDGK16")){
                pGPDModel = PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(PARTONS::GPDGK16::classId);
                ID = PARTONS::GPDGK16::classId;
        }
        else {cout << "Invalid GPD model, exiting" << endl; exit(1);}
        PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule = PARTONS::Partons::getInstance()->getModuleObjectFactory()->newRunningAlphaStrongModule(
                    PARTONS::RunningAlphaStrongStandard::classId);

        //nucleon GPD test
        // for(int i=0;i<=200;i++){
        //     for(int j=0;j<=100;j++){
        //         PARTONS::GPDKinematic gpdKinematic(0.01*(i-100)+(i==100? 1.E-04:0),0.01*(j),t, 1., 1.);
        //         //PARTONS::GPDKinematic gpdKinematic(0.1, 0.2, -0.1, 1., 1.);
        //         PARTONS::GPDResult gpdResult = pGPDService->computeGPDModel(gpdKinematic,
        //                 pGPDModel);
        //         double H=0.5*(gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
        //             +gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution());
        //         double E=0.5*(gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
        //             +gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution());
        //         cout << 0.01*(i-100) << " " << 0.01*j << " " << gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution() <<
        //         " " << gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution()
        //          << " " << gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution() << " " 
        //          << gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution() << endl;

        //     }
        // }


        //cout << "alpha_s " << alpha_s << endl;

        Deut_Conv_GPD_V test=Deut_Conv_GPD_V(pGPDService,pGPDModel,wf,ID);
        test.getWf()->setD(setD);
        test.getWf()->setS(setS);
        test.setH(setH);
        test.setE(setE);
        //-99 to 99

        //test.getDeut_GPD_V_set(0.,xi,t,scale, ERBL).getAmp_00();
    
        //non rel test
        // for (int i=0;i<=200;i++){
        //     double Q = i*2./200.;
        //     double t = -Q*Q;
        //     double eta = -t/4./MASSD/MASSD*1.E06;
        //     vector<double> FF = test.calc_NR_ffs(t);
        //     double G1 = FF[0]-2./3.*eta*FF[2];
        //     double G3 = (FF[2]-G1+FF[1])/(1.+eta);
        //     cout << Q << " " << FF[0] << " " << FF[1] << " " << FF[2] << " " << G1 << " " << G3 << endl;
        // }
        // exit(1);

        // for (int i=0;i<=40;i++){
        //     double Q = i*2./40.;
        //     t = -Q*Q;
        //     cout << Q << " ";
        //     //NucleonEMOperator proton(-t*1.E06,1,0), neutron(-t*1.E06,0,0);
        //     //cout << proton.getGE()+neutron.getGE() << " " << proton.getGM()+neutron.getGM() << " " << proton.getF1() + neutron.getF1() << " " << proton.getF2() + neutron.getF2() << endl;
            
            
            
        //     vector< complex<double> > out = test.FF_conv(xi,t);
        //     for(int j=0;j<5;j++) cout << out[j].real() << " ";
            
        //      vector< complex<double> > hels = {out[0],out[2],out[4]};
        //      vector< complex<double> > FFs = Deut_Conv_GPD_V::helamps_to_FFs_V(xi,t,hels);
            
        //      double FM = FFs[1].real();
        //      double FQ = (FFs[0].real()-FFs[1].real()+(1.-t/4./MASSD/MASSD*1.E06)*FFs[2].real());
        //      double FC = (FFs[0].real()-t/4./MASSD/MASSD*1.E06*FQ*2./3.);
            
        //     for(int j=0;j<3;j++) cout << FFs[j].real() << " ";
        //     cout << FC << " "<< FM << " " << FQ << " ";
            
            
            
        //     vector< complex<double> > gpd = Deut_Conv_GPD_V::helamps_to_gpds_V(xi,t,out);
        //     // vector< complex<double> > hel = Deut_Conv_GPD_V::gpds_to_helamps_V(xi,t,gpd);

        //     double GM = gpd[1].real();
        //     double GQ = (gpd[0].real()-gpd[1].real()+(1.-t/4./MASSD/MASSD*1.E06)*gpd[2].real());
        //     double GC = (gpd[0].real()-t/4./MASSD/MASSD*1.E06*GQ);
            
        //     //for(int j=0;j<5;j++) cout << out[j].real() << " ";
        //     for(int j=0;j<5;j++) cout << gpd[j].real() << " ";
        //     cout << GC << " " << GM << " " << GQ << endl;
        // }
        double t0=-4.*MASSD_G*MASSD_G*xi*xi/(1-xi*xi); // GeV^2
        double D=(t0-t)/(4.*MASSD_G*MASSD_G); // GeV^2
        cout << "ximax " << sqrt(-t)/(sqrt(4.*MASSD_G*MASSD_G-t)) << " D " << D << endl;
        NucleonEMOperator proton(-t*1.E06,1,0), neutron(-t*1.E06,0,0);
        //cout << proton.getF1() << " " << proton.getF2() << " " << neutron.getF1() << " " << neutron.getF2() << endl;
        cout << (proton.getF1()+neutron.getF1())/2. << " " <<  (proton.getF2()+neutron.getF2())/2. << endl << endl << endl;
        for(int i=-99;i<=99;i++){
            double x=i*0.01;

            if(t>t0){
                cout << x << " ";
                for(int j=0;j<5;j++) cout << 0. << " ";
                for(int j=0;j<5;j++) cout << 0. << " ";
                // for(int j=0;j<5;j++) cout << hel[j].real() << " ";
                cout << endl;
            }
            else{
                vector< complex<double> > out = test.gpd_conv(xi,x,t,scale);
                vector< complex<double> > gpd = Deut_Conv_GPD_V::helamps_to_gpds_V(xi,t,out);
                // vector< complex<double> > hel = Deut_Conv_GPD_V::gpds_to_helamps_V(xi,t,gpd);
                
                cout << x << " ";
                for(int j=0;j<5;j++) cout << out[j].real() << " ";
                for(int j=0;j<5;j++) cout << gpd[j].real() << " ";
                // for(int j=0;j<5;j++) cout << hel[j].real() << " ";
                cout << endl;
            }
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

