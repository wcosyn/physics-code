//program used to fit the EMC vs a2 data

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
#include <cassert>
using namespace std;

#include <constants.hpp>

vector<double> EMC_data={-0.207,
-0.326,
-0.285,
-0.222,
-0.283,
-0.322,
-0.391,
-0.391,
-0.511,
-0.34,
-0.347,
-0.472,
-0.539
};
vector<double> EMC_error={0.025,
0.026,
0.026,
0.045,
0.028,
0.033,
0.025,
0.025,
0.03,
0.022,
0.022,
0.023,
0.02}; 
vector<double> m1factor={3.05433370846,
3.14661535135,
3.48440239145,
3.05433370846,
3.14661535135,
3.48440239145,
3.96941579834,
3.89401669771,
4.12119329745,
3.48440239145,
3.75785387162,
3.96941579834,
4.11068824918};

vector<double> m2factor={0,
-0.262521632186,
0,
0,
-0.262521632186,
0,
-0.209628370758,
-0.285416235795,
-0.838865756409,
0,
-0.104697839081,
-0.209628370758,
-0.662374991039};
int dof=0;
bool SLAC;


void Fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag); //fit function



int main(int argc, char *argv[])
{
  
  
  

  int testing = 0;
  int fNDim = 2; // number of dimensions
  double fLo[] = {-100.,-100.}; // lower limits of params
  double fHi[] = {100.,100.}; // upper limits of params
  char* fName[] = {"m1","m2"};
  
  int fBound[] = {1,1};

  TVirtualFitter *gMinuit = TVirtualFitter::Fitter ( 0, fNDim );

  // Start values of parameters
  // If you have a starting individual, you can simply use 
  double minuitIndividual[] = {0.,0.}; // FIXME insert your starting individual ( double array) here
  
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

  gMinuit->FixParameter ( 1 );
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
      arglist[0] = 10000;
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
        arglist[0] = 10000;
        gMinuit->ExecuteCommand ( "IMPROVE", arglist, 1 );

        // we exit the inner loop
        innerMinimizationLoop = false;
}
      // Else...
      else
{
        // ...try to find a better spot in parameter space
        arglist[0] = 10000;
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
          // exit ( 1 );
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

  cout << "Covariance matrix" << endl;
  for ( int i = 0 ; i < nFreePars ; ++i )
{
    for ( int j = 0 ; j < nFreePars ; ++j )
{
      matrix ( i, j ) = gMinuit->GetCovarianceMatrixElement ( i, j );
      cout << matrix(i,j) << " ";
}
cout << endl;
}
  cout << nFreePars << endl;
  double eplus[nFreePars],eminus[nFreePars],eparab[nFreePars],globcc[nFreePars], params[nFreePars];
  for(unsigned int i=0;i<nFreePars;++i){
    gMinuit->GetErrors(i,eplus[i],eminus[i],eparab[i],globcc[i]);
    params[i]=gMinuit->GetParameter(i);
  }
  //double paramin[3]={params[0],nFreePars>1?params[1]:8.,-0.5};

  //Summary of my fit
  int n=nFreePars;
  double f;
  Fcn(n, &f, f, params, n);
  //cout << Warray[Windex] << " ";
  cout << "parameter values + parabolic errors" << endl;
  for(unsigned int i=0;i<nFreePars;++i) cout << params[i] << " " << eparab[i] << endl;
  cout << endl << endl << "chi2 per dof, chi2, dof, #freeparams" << endl;
  cout << f/(dof-nFreePars) << " " << f << " " << dof << " " << nFreePars << endl;
  cout << endl << endl;

  cout << "data | fit" << endl;
for(unsigned i=0; i < EMC_data.size(); i++){
  cout << -EMC_data[i] << " " << params[0]*m1factor[i] + params[1]*m2factor[i] << endl;
}
  
  delete gMinuit;
//   cout << "\nCleaning up..." << endl;

  
  return 0;
}



//fit function
void Fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  f=0.;
  dof=0;
  //cout << npar << endl;
  //DeuteronCross DeepsCross("paris",proton,"SLAC",par[0],par[1],epsilon,betaoff,lambdain,offshellset,looplimit);
  // cout << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4] << " " << par[5] << " " << endl;



  for(unsigned i=0; i < EMC_data.size(); i++){
    f+=pow((-EMC_data[i]-par[0]*m1factor[i]-par[1]*m2factor[i])/EMC_error[i],2.);
    dof++;
  }
  // cout << "intermediate chi2 " << f/(dof-npar) << " " << dof << " " << npar << " " << f << endl;
}
