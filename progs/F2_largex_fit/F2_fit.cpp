//program used to fit the rescattering parameters in the FSI to the deeps data

#include <TVirtualFitter.h>
#include <TFitter.h>
#include <TMatrixT.h>
#include <TMinuit.h>
// #include <TMath.h>

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <complex>
#include <cstdarg>
#include <time.h>

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <signal.h>
#include <iomanip>
#include <vector>
using namespace std;

#include <constants.hpp>
#include "fomin_data.hpp"

vector<double> xi_data; //nachtmann
vector<double> x_data; //bjorken x
vector<double> F2_data;
vector<double> F2_err;
vector<double> Q2_data;
int dof=0;

void plotfit(double *par){
  double xvalues[]={0.75,0.85,0.95,1.05,1.15,1.25};
  double adam[]={0.52577737,  3.12260242, -8.5547532, 7.89292279,  2.22163175, -0.76036478};
  
  for(int x_it=0;x_it<6;x_it++){
    ofstream myfile;
    string filename = "x"+to_string(xvalues[x_it]);
    myfile.open(filename.c_str());
    for(unsigned i=0; i < xi_data.size(); i++){
      if(abs(x_data[i]-xvalues[x_it])<.01){
        double nu=Q2_data[i]/(2.*11.17793/12*xvalues[x_it]);
        double xi=2.*xvalues[x_it]/(1+sqrt(1.+Q2_data[i]/nu/nu));
        myfile << x_data[i] << " " << Q2_data[i] << " " << F2_data[i] << " " << F2_err[i] << " " << exp(par[0]+par[1]*xi+par[2]*pow(xi,2.))*
              (1.+par[3]*pow(xi,par[4])*(1.+xi*par[5])/Q2_data[i])
              << " " << exp(adam[0]+adam[1]*xi+adam[2]*pow(xi,2.))*
              (1.+adam[3]*pow(xi,adam[4])*(1.+xi*adam[5])/Q2_data[i]) << " " << 
              exp(par[0]+par[1]*xi_data[i]+par[2]*pow(xi_data[i],2.))*
              (1.+par[3]*pow(xi_data[i],par[4])*(1.+xi_data[i]*par[5])/Q2_data[i])
              << " " << exp(adam[0]+adam[1]*xi_data[i]+adam[2]*pow(xi_data[i],2.))*
              (1.+adam[3]*pow(xi_data[i],adam[4])*(1.+xi_data[i]*adam[5])/Q2_data[i]) << endl;
      }
    }
    // for(int i=1;i<=20;i++){
    //   double Q2=10.*30./20.*i;
    //   double nu=Q2/(2.*11177.93/12*xvalues[x_it]);
    //   double xi=2.*xvalues[x_it]/(1+sqrt(1.+Q2/nu/nu));
    //   myfile << xvalues[x_it] << " " << Q2 << " " << "nan" << " " << "nan" << " " << exp(par[0]+par[1]*xi+par[2]*pow(xi,2.))*
    //         (1.+par[3]*pow(xi,par[4])*(1.+xi*par[5])/Q2)
    //         << " " << exp(adam[0]+adam[1]*xi+adam[2]*pow(xi,2.))*
    //         (1.+adam[3]*pow(xi,adam[4])*(1.+xi*adam[5])/Q2) << endl;

    // }
    myfile.close();
  }



}

void Process_Data(){
  xi_data.clear(); F2_data.clear(); F2_err.clear(); Q2_data.clear(); x_data.clear();
  
  for(int i=0;i<fomin_data::datapoints;i++){
    if(fomin_data::A[i]==12){
      double nu = (fomin_data::Ein[i]-fomin_data::Eout[i])*1.E03;
      double Q2 = 4.*fomin_data::Ein[i]*fomin_data::Eout[i]*1.E06*pow(sin(fomin_data::theta[i]*DEGRTORAD*0.5),2.);
      double x= Q2/(2.*11177.93/12*nu);
      double nachtmann = 2.*x/(1.+sqrt(1+Q2/nu/nu));
      //cout << nu << " " << Q2*1.E-06 << " " << x  << " " << nachtmann << endl;
      if(nachtmann>0.5){
        xi_data.push_back(nachtmann);
        Q2_data.push_back(Q2*1.E-06);
        F2_data.push_back(fomin_data::F2[i]*1.E06);
        F2_err.push_back(fomin_data::F2_stat[i]*1.E06);
        x_data.push_back(x);
      }
    }

  }


}



void Fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  f=0.;
  dof=0;
  cout << npar << endl;
  //DeuteronCross DeepsCross("paris",proton,"SLAC",par[0],par[1],epsilon,betaoff,lambdain,offshellset,looplimit);
  cout << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4] << " " << par[5] << " " << endl;

  for(unsigned i=0; i < xi_data.size(); i++){
    f+=pow((F2_data[i]-exp(par[0]+par[1]*xi_data[i]+par[2]*pow(xi_data[i],2.))*
      (1.+par[3]*pow(xi_data[i],par[4])*(1.+xi_data[i]*par[5])/Q2_data[i]))/F2_err[i],2.);
    dof++;
  }
  cout << "intermediate chi2 " << f/(dof-npar) << " " << dof << " " << npar << " " << f << endl;
}


int main(int argc, char *argv[])
{
  
  Process_Data();
  // for (unsigned i=0; i < xi_data.size(); i++){
  //   cout << i << " " << xi_data[i] << " " << x_data[i] << " " << Q2_data[i]*1.E-06 << " " << F2_data[i] << endl;
  //  }
  // exit(1);

  int testing = 0;
  int fNDim = 6; // number of dimensions
  double fLo[] = {-100.,-100.,-100.,-100.,-100.,-1000.}; // lower limits of params
  double fHi[] = {100.,100.,100.,100.,100.,1000.}; // upper limits of params
  char* fName[] = {"LT_norm","LT_xi","LT_xi2", "HT_xiff", "HT_xipow", "HT_linxi"};
  
  int fBound[] = {1,1,1,1,1,1};

  TVirtualFitter *gMinuit = TVirtualFitter::Fitter ( 0, fNDim );

  // Start values of parameters
  // If you have a starting individual, you can simply use 
  double minuitIndividual[] = {0.,0.,0.,0.,0.,0.}; // FIXME insert your starting individual ( double array) here
  
  std::cout << "done" << endl;

  //------------------------------------------------------------------
  // The Minimisation section...

  std::cout << "Start minimizing...\n";

 // Initialize TMinuit via generic fitter (Minuit)
  //TVirtualFitter::SetDefaultFitter("Minuit2"); // if you want to use minuit2
  gMinuit->SetFCN (  (void (*)(Int_t &, Double_t *, Double_t &, Double_t *, Int_t))   &Fcn );  // which function to minimise

  // ignore until ---------
  Double_t arglist[fNDim]; // For interface to Minuit routines
  Int_t ierflg = 0;  // ditto

  // Let Minuit produce minimal output...
  arglist[0] = 0;
  gMinuit->ExecuteCommand ( "SET PRINT", arglist, 1 );

  // Set Minuit to the medium strategy
  arglist[0] = 1;
  gMinuit->ExecuteCommand ( "SET STR", arglist, 1 );
  
  // If seed is given, set random seed.
  // This can be any integer between 10000 and 900000000.
//     arglist[0] = (Double_t) fRandomseed; 
//     while (arglist[0] >= 900000000)  arglist[0] -= 900000000; 
//     while (arglist[0] <= 10000)  arglist[0] += 10000;
//     gMinuit->ExecuteCommand ( "SET RAN", arglist, 1 );
//   
  // Set the limits for each parameter...
  double step; // step size
  for ( int i = 0 ; i < fNDim ; i++ )
{
    step = 0.01 * minuitIndividual[i];  // step size is 1% of value - so don't set value to 0 ;-)
    if ( step < 0 ) step *= -1; // make it positive
    else if ( step == 0.0 ) step = 1.0e-3;
    gMinuit->SetParameter ( i, fName[i], minuitIndividual[i], step, fLo[i], fHi[i] );
}

  // For testing, fix all but three parameters...
  if (testing)
{
    std::cout << "Running reduced capability for testing...\n";
    std::cout << "Fixing all but three parameters...\n";
    for ( int i = 3 ; i < gMinuit->GetNumberTotalParameters() ; ++i )
    {
	  gMinuit->FixParameter ( i );
    }
}

  std::cout << "Fixing some parameters...\n";
//   if(offshellset==2){
//     for ( int i = fNDim-fixparam-1 ; i < gMinuit->GetNumberTotalParameters()-1 ; ++i ) gMinuit->FixParameter ( i );    
//   }
//   else{
//     for ( int i = fNDim-fixparam ; i < gMinuit->GetNumberTotalParameters() ; ++i ) gMinuit->FixParameter ( i );
//   }
  
  // Print out to verify...
  std::cout << "No. of parameters:\t\t"
  << gMinuit->GetNumberTotalParameters() << endl;
  std::cout << "No. of free parameters:\t\t"
  << gMinuit->GetNumberFreeParameters() << endl;
  std::cout << "Max. No. of iterations:\t\t"
  << gMinuit->GetMaxIterations() << endl;

  /* **************** *
  * START MINIMIZING *
  * **************** */

  Double_t value, verr; // value of parameter and its error
  Int_t nFreePars = gMinuit->GetNumberFreeParameters();// # of free parameters
  Int_t someParsFixed = 0; // # of fixed parameters

  bool outerMinimizationLoop = true; /* Stay in outside loop untill this
				     * flag is turned to false */
  bool innerMinimizationLoop = true; /* Stay in inside loop untill this
				     * flag is turned to false */

  /* The outer minimization loop:
  * - Using command: MIGRAD
  * - SET STRATEGY 1
  * - max. 7000 function calls
  * untill convergence with all
  * parameters away from limits */
  while ( outerMinimizationLoop )
{

    /* The inner minimization loop:
    * - Using command: MIGRAD (and SIMPLEX)
    * - SET STRATEGY 1
    * - max. 7000 function calls
    * untill convergence */
    while ( innerMinimizationLoop )
{

      // Run minimization...
      // ...untill convergence is reached (max 4 times)
      arglist[0] = 7000;
      gMinuit->ExecuteCommand ( "MIGRAD", arglist, 1 );
      for ( int i = 0; i < 3 &&
            static_cast<TFitter*> ( gMinuit->GetFitter() )->
            GetMinuit()->fCstatu.CompareTo ( "CONVERGED " );
            i++ )
        gMinuit->ExecuteCommand ( "MIGRAD", arglist, 1 );

      // When MIGRAD converged...
      if ( !static_cast<TFitter*> ( gMinuit->GetFitter() )->
           GetMinuit()->fCstatu.CompareTo ( "CONVERGED " ) )
{
        // ...try to IMPROVE the chi-squared
        arglist[0] = 7000;
        gMinuit->ExecuteCommand ( "IMPROVE", arglist, 1 );

        // we exit the inner loop
        innerMinimizationLoop = false;
}
      // Else...
      else
{
        // ...try to find a better spot in parameter space
        arglist[0] = 7000;
        arglist[1] = 0.01;
        gMinuit->ExecuteCommand ( "SIMPLEX", arglist, 2 );

        // if SIMPLEX did not converge, we give up!
        if ( static_cast<TFitter*> ( gMinuit->GetFitter() )->
             GetMinuit()->fCstatu.CompareTo ( "CONVERGED " ) )
{
          std::cout << endl
          << "**********************************************" << endl
          << "MIGRAD did not manage to converge: we give up!" << endl
          << "**********************************************" << endl;
          exit ( 1 );
}
}
} // end inner loop

    // Check to see if any parameters are at their limits. This is
    // done by checking whether they are within half an error unit of
    // either limit. If they are, then fix them and run
    // the first minimization loop again.
    // Obviously, we do not need to perform this check for
    // parameters without limits
    // When MINUIT issues a warning about the covariance matrix, it
    // is questionable to rely on the errors. Therefor we only trust
    // the errors when the status of the error matrix is set to
    // accurate or forced pos-def.
    if ( !static_cast<TFitter*> ( gMinuit->GetFitter() )->GetMinuit()->
         fCovmes[static_cast<TFitter*> ( gMinuit->GetFitter() )->
                 GetMinuit()->fISW[1]].CompareTo ( "ERROR MATRIX ACCURATE " ) ||
         !static_cast<TFitter*> ( gMinuit->GetFitter() )->GetMinuit()->
         fCovmes[static_cast<TFitter*> ( gMinuit->GetFitter() )->
                 GetMinuit()->fISW[1]].CompareTo ( "ERR MATRIX NOT POS-DEF" ) )
{
      for ( int i = 0 ; i < fNDim ; ++i )    // loop over parameters
{
	if (fBound[i])
        //if ( GetTCalc()->GetLimit ( i ).bound )   // check if bound
{
          value = gMinuit->GetParameter ( i );
          verr  = gMinuit->GetParError ( i );

          // If the parameter is within half an error bar of one of
          // the limits, it should be fixed, as it may be that the
          // true minimum lies beyond the defined physical limit
          if ( TMath::Abs ( value - fLo[i] ) < verr*0.5 || // half a sigma
               TMath::Abs ( value - fHi[i] ) < verr*0.5 )
{
            gMinuit->FixParameter ( i );
            someParsFixed++;
            std::cout << "Parameter no." << i + 1
            << " has been fixed: Too close to limits!\n";
}
} // if parameter bound
} // end loop over parameters
} // status check error matrix

    // If some parameters had to be fixed, re-run the loop...
    if ( someParsFixed )
{
      someParsFixed = 0; // reset
      std::cout << "Some parameters fixed, rerun Migrad()...\n";
}
    else
{
      outerMinimizationLoop = false; // exit outer loop
}
} // end outer loop

  /* ****************************** *
  * Determine precise error matrix *
  * ****************************** */

  /*    // Remove the limits on all parameters (if there are any)
  // That way the error matrix can be determined correctly.
  if (!static_cast<TFitter*>(gMinuit->GetFitter())->GetMinuit()->fLnolim) {
  std::cout << "Removing limits...\n";
  gMinuit->ExecuteCommand("SET LIM", arglist, 0);

    // Minimize one last time...
    // ...untill convergence is reached (max 5 times)
  std::cout << "Refit with unbound parameters...\n";
  arglist[0] = 7000;
  gMinuit->ExecuteCommand("MIGRAD", arglist, 1);
  for(int i=0; i<4 &&
  static_cast<TFitter*>(gMinuit->GetFitter())->
  GetMinuit()->fCstatu.CompareTo("CONVERGED ");
  i++)
  gMinuit->ExecuteCommand("MIGRAD", arglist, 1);

    // When MIGRAD does not converge... give up!
  if (static_cast<TFitter*>(gMinuit->GetFitter())->
  GetMinuit()->fCstatu.CompareTo("CONVERGED "))
{
  std::cout << endl
  << "**********************************************" << endl
  << "MIGRAD did not manage to converge: we give up!" << endl
  << "**********************************************" << endl;
  exit(1);
}
} // end remove limits
  */
  // Calculate the covariance matrix again...
  gMinuit->ExecuteCommand ( "MINOS", arglist, 0 );
  gMinuit->ExecuteCommand ( "SHO COV", arglist, 0 );

  // Extract the covariance matrix (method GetCovarianceMatrix only
  // returns a double*)
  nFreePars = gMinuit->GetNumberFreeParameters();
  TMatrixT<double> matrix ( nFreePars, nFreePars );

  for ( int i = 0 ; i < nFreePars ; ++i )
{
    for ( int j = 0 ; j < nFreePars ; ++j )
{
      matrix ( i, j ) = gMinuit->GetCovarianceMatrixElement ( i, j );
}
}
  cout << nFreePars << endl;
  double eplus[nFreePars],eminus[nFreePars],eparab[nFreePars],globcc[nFreePars], params[nFreePars];
  for(unsigned int i=0;i<nFreePars;++i){
    gMinuit->GetErrors(i,eplus[i],eminus[i],eparab[i],globcc[i]);
    params[i]=gMinuit->GetParameter(i);
  }
  //double paramin[3]={params[0],nFreePars>1?params[1]:8.,-0.5};
  cout << "bla" << endl;
  int n=nFreePars;
  double f;
  Fcn(n, &f, f, params, n);
  //cout << Warray[Windex] << " ";
  cout << "parameter values + parabolic errors" << endl;
  for(unsigned int i=0;i<nFreePars;++i) cout << params[i] << " " << eparab[i] << endl;
  cout << endl << endl << "chi2 per dof, chi2, dof, #freeparams" << endl;
  cout << f/(dof-nFreePars) << " " << f << " " << dof << " " << nFreePars << endl;

  double adam[]={0.52577737,  3.12260242, -8.5547532, 7.89292279,  2.22163175, -0.76036478};
  Fcn(n, &f, f, adam, n);

  cout << "Adam values" << endl << "chi2 per dof, chi2, dof, #freeparams" << endl;
  cout << f/(dof-nFreePars) << " " << f << " " << dof << " " << nFreePars << endl;
  
  delete gMinuit;
//   cout << "\nCleaning up..." << endl;
//   DeuteronCross::maint_deepsarray(deepsdata);

  plotfit(params);

  return 0;
}
