#include <iostream>
#include <cstdlib>
#include <vector>

using namespace std;


#include <iostream>
#include <cmath>
#include "constants.hpp"

#include <Deut_Conv_GPD_V.hpp>
#include "TransGPD_set.hpp"

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


int main(int argc, char *argv[]){
    std::string wf=argv[1];
    double xi=atof(argv[3]);
    double t=atof(argv[2])*-1.E06;
    // cout << "xi " << xi << " t " << t << " model " << model << endl;
    // cout << -4.*MASSD*MASSD*xi*xi/(1-xi*xi)-t << endl;


    QCoreApplication a(argc, argv);
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
        Deut_Conv_GPD_V test=Deut_Conv_GPD_V(pGPDService,pGPDModel,wf);

        for(int i=-99;i<=99;i++){
            double x=i*0.01;
            vector< complex<double> > out = test.gpd_conv(xi,x,t);
            vector< complex<double> > gpd = Deut_Conv_GPD_V::helamps_to_gpds_V(xi,t,out);
            vector< complex<double> > hel = Deut_Conv_GPD_V::gpds_to_helamps_V(xi,t,gpd);
            
            cout << x << " ";
            for(int j=0;j<5;j++) cout << out[j] << " ";
            for(int j=0;j<5;j++) cout << gpd[j] << " ";
            cout << endl;
        
        }
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

