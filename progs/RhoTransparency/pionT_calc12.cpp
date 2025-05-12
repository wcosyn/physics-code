#include <iostream>
#include <cstdlib>

using namespace std;

#include <PionTCross.hpp>
#include <Utilfunctions.hpp>
#include <string>



int main(int argc, char *argv[])
{

    string arg_names[5]={"exec name", "nucleus name","beam energy [GeV]","scattered energy [GeV]","scattered angle [deg]"};
    if(argc!=5){
        cout << "Usage: " << argv[0] << " <nucleus name> <beam energy [GeV]> <scattered energy [GeV]> <scattered angle [deg]>" << endl;
        return 1;
    }
    std::cout << "Called from file: " << __FILE__ << std::endl;
    Bookkeep(argc,argv,arg_names);  



    string nucleus_name=argv[1];  //nucleus name: options are He, C, O, Fe, Pb, Al, Cu, Au
    if(nucleus_name!="He" && nucleus_name!="C" && nucleus_name!="O" && nucleus_name!="Fe" && nucleus_name!="Pb" && nucleus_name!="Al" && nucleus_name!="Cu" && nucleus_name!="Au"){
        cout << "Nucleus name not recognized. Options are He, C, O, Fe, Pb, Al, Cu, Au." << endl;
        return 1;
    }
    double Ebeam  = atof(argv[2])*1.E03;  //Beam energy in GeV, converted to MeV
    double E_out = atof(argv[3])*1.E03;  // Scattered electron energy in GeV, converted to MeV
    double theta_e = atof(argv[4])*DEGRTORAD; // Scattered electron angle in degrees input converted to radians
    
    double Q2 = 2.*Ebeam*E_out*(1.-cos(theta_e)); //Q2 in MeV^2

    PionTCross pionT_object = PionTCross(MeanFieldNucleus::TypeNames.at(nucleus_name),
                                  400.,HOMEDIR,0,0.,1.E-05,2,20000);
    double *results = new double[pionT_object.getNrofcross()];
    pionT_object.getCross(results, Ebeam,E_out,theta_e);

    cout << Q2*1.E-06 << " ";
    for(int i=0;i<pionT_object.getNrofcross();i++){
        cout << results[i] << " ";
    }
    cout << endl;
    delete [] results;
    return 0;
}