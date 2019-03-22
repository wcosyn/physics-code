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
#include <partons/Partons.h>
#include <partons/services/ConvolCoeffFunctionService.h>
#include <partons/services/GPDService.h>
#include <partons/services/ObservableService.h>
#include <partons/ServiceObjectRegistry.h>

#include <partons/modules/gpd/GPDMMS13.h>

//#include "../include/examples.h"
//headers for integration 
#include <numint/numint.hpp>
#include <constants.hpp>

using namespace std;

//needed for integrandum of amplitude
struct Ftor_2vector {

  static void exec(const numint::array<double,2> &x, void *param, numint::vector_d &ret) {
    Ftor_2vector &p = * (Ftor_2vector *) param;
    p.f(ret,x[0],x[1],*p.pGPDService, *p.pGPDModel, p.xi, p.t, p.psq, p.Qsq);
  }
  PARTONS::GPDModule* pGPDModel;
  PARTONS::GPDService* pGPDService;
  double xi;
  double t;
  double psq;
  double Qsq;
  void (*f)(numint::vector_d &, double u, double z, PARTONS::GPDService &gpdService, PARTONS::GPDModule &gpdModel, double xi, double t, double psq, double Qsq);

};

// integrandum for longitudinal polarized photon on nucleon, two times rho_L in final state
void integrandum_L(numint::vector_d &, double u, double z, PARTONS::GPDService &gpdService, PARTONS::GPDModule &gpdModel, double xi, double t, double psq, double Qsq);

// integrandum for transversily polarized photon on nucleon, two times rho_L in final state
void integrandum_T(numint::vector_d &, double u, double z, PARTONS::GPDService &gpdService, PARTONS::GPDModule &gpdModel, double xi, double t, double psq, double Qsq);

/*
 * Main function.
 */
int main(int argc, char** argv) {

    // Init Qt4
    QCoreApplication a(argc, argv);
    PARTONS::Partons* pPartons = 0;

    double xi=atof(argv[1]);
    double Qsq = atof(argv[2]); // virtual photon scale squared [GeV^2]

    double t= -4.*MASSP*MASSP*xi*xi/(1-xi*xi)*1.E-06;//tmin [GeV^2]

    


    double frhoplus = 0.198; //rho+ decay constant [GeV]
    double frho0 = 0.216; //rho0 decay constant [GeV]



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

        double alpha_s = pRunningAlphaStrongModule->compute(1.); //scale is 1 GeV^2
        //cout << "alpha_s " << alpha_s << endl;
        int Nc=3; //number of colors
        double CF=(Nc*Nc-1)/2./Nc;

        Ftor_2vector F;
        F.pGPDService = pGPDService;
        F.pGPDModel = pGPDModel;
        F.xi=xi;
        F.t=t;
        F.Qsq=Qsq;

        numint::mdfunction<numint::vector_d,2> mdf;
        mdf.func = &Ftor_2vector::exec;
        mdf.param = &F;

        numint::array<double,2> lower = {{0.,0.}};
        numint::array<double,2> upper = {{1.,1.}};
        for(int i=0;i<=16;i++){
        
            double psq = 2.+i*0.5; //pomeron scale squared [GeV^2]
            F.psq=psq;

            std::vector<double> integral(1,0.); //integrandum result array

            //gamma_L calculation
            F.f=integrandum_L;
            unsigned count=0;
            // integration over u and z
            if(Qsq>1.E-02) numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,2E05,integral,count,0);
            integral[0]*=36.; //prefactor DA (twice, squared)
            
            //cout << count << " " << integral[0] << endl;
            
            integral[0] *= 16.*PI*PI*alpha_s*frhoplus*xi*sqrt((1-xi)*(1+xi))*CF/Nc/psq/psq; //prefactors Eq (14) Enberg et al.
            integral[0] *= 2.*PI*alpha_s*sqrt(Qsq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI); //prefactors Eq (15) Enberg et al.

            double result_L=pow(integral[0],2.)/256./pow(PI,3.)/xi/(1.-xi); // Enberg et al Eq (12)
            //cout << Qsq << " " << xi << " " << psq << " " << result_L*0.389379E06 << endl;  // factor to go from GeV-6 to nb GeV-4
            
            //gamma_T calculation
            F.f=integrandum_T;
            count=0;
            integral[0]=0.;
            numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-02,1E02,2E03,integral,count,0);
            integral[0]*=36.; //prefactor DA (twice, squared)
            //cout << count << " " << integral[0] << endl;
 
            integral[0] *= 16.*PI*PI*alpha_s*frhoplus*xi*sqrt((1-xi)*(1+xi))*CF/Nc/psq/psq; //prefactors Eq (14) Enberg et al.
            integral[0] *= PI*alpha_s*sqrt(psq)/Nc/sqrt(2.)*frho0*sqrt(ALPHA*4.*PI); //prefactors Eq (17) Enberg et al.

            double result_T=pow(integral[0],2.)/256./pow(PI,3.)/xi/(1.-xi); // Enberg et al Eq (12)

 
            cout << Qsq << " " << xi << " " << psq << " " << result_T*0.389379E06 << " " << result_L*0.389379E06 << endl;
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

// integrandum for gamma_L amplitude
void integrandum_L(numint::vector_d &result, double u, double z, PARTONS::GPDService &gpdService, 
                    PARTONS::GPDModule &gpdModel, double xi, double t, double psq, double Qsq){
    result = vector<double>(1,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) result[0]=0.;
    else{ 
        
        PARTONS::GPDKinematic gpdKinematic(xi*(2.*u-1),xi,t, 1., 1.);

        PARTONS::GPDResult gpdResult = gpdService.computeGPDModel(gpdKinematic,
            &gpdModel);

        // (Hq^u-Hq^d) * z^2 * barz^2 / u / baru * P [z*barz,u*baru from 2 DA taken into account, factors 6 afterwards]
        result[0]=(gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                    -gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution())
                    *z*z*(1.-z)*(1.-z)/u/(1.-u)
                        *(1./(z*z*psq+Qsq*z*(1.-z)) + 1./((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) - 1./((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - 1./((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) );

    }
    //cout << u << " " << z << " " << result[0] << endl;
    return;

}

// integrandum for gamma_T amplitude
void integrandum_T(numint::vector_d &result, double u, double z, PARTONS::GPDService &gpdService, 
                    PARTONS::GPDModule &gpdModel, double xi, double t, double psq, double Qsq){
    result = vector<double>(1,0.);
    //limits are well behaved
    if(u==0||u==1||z==0||z==1) result[0]=0.;
    else{ 
        
        PARTONS::GPDKinematic gpdKinematic(xi*(2.*u-1),xi,t, 1., 1.);

        PARTONS::GPDResult gpdResult = gpdService.computeGPDModel(gpdKinematic,
            &gpdModel);

        // (Hq^u-Hq^d) * (2z-1) * z * barz / u / baru * Q    [z*barz,u*baru from 2 DA taken into account, factors 6 afterwards] [factor p_perp from nominator taken out]
        result[0]=(gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                    -gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution())
                    *(2.*z-1)*z*(1.-z)/u/(1.-u)
                        *(z/(z*z*psq+Qsq*z*(1.-z)) - (1.-z)/((1.-z)*(1.-z)*psq+Qsq*z*(1.-z)) + (u-z)/((u-z)*(u-z)*psq+Qsq*z*(1.-z)) - (u-1+z)/((u-1.+z)*(u-1.+z)*psq+Qsq*z*(1.-z)) );

    }
     return;

}