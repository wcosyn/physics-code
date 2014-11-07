
#include <TVirtualFitter.h>
#include <TFitter.h>
#include <TMatrixT.h>
#include <TMinuit.h>
#include <TMath.h>

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
using namespace std;

#include <constants.hpp>
#include <DeuteronCross.hpp>
#include "bonusdata.h"
#include "bonusfits.h"

//double sigmainput=40.;
// int phiavg=1;
// int which_dx=0;
// int F_param=0;
int Qindex=-1;
int Windex=-1;
int startset=0;
int stopset=4;
int offshellset=0;
string dir;
int looplimit=-1;
bool lc=0;

double get_normfit_bonus(int beamindex, int Qindex, int Windex, int psindex, int offshell, bool lc);

void Fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  f=0.;
  int dof=0.;
  bool proton=0;
  double epsilon=-0.5;
  double betaoff=8.;
  double lambdain=1.2;
//   cout << npar << endl;
  if(offshellset==1) lambdain=par[2];
  if(offshellset==2) betaoff=par[2];
  else epsilon=par[2];  //we keep epsilon fixed, didn't improve fit!
  DeuteronCross DeepsCross("paris",proton,"CB",par[0],par[1],epsilon,betaoff,lambdain,offshellset,1E03);
  cout << "bla " << par[0] << " " << par[1] << " " << par[2]  << endl;
  for(int i=startset;i<stopset;i++){
    for(int j=0;j<10;j++){
      double error,result,costheta;
      if(!((Qindex==2)&&(Windex==4))){
	error=data::bonusdata4[Qindex][Windex][i][j][2];
	result=data::bonusdata4[Qindex][Windex][i][j][1];
	costheta=data::bonusdata4[Qindex][Windex][i][j][0];
	if(result!=0.){
	  double bonusMC=0., pw=0.,fsi=0.;
	  DeepsCross.getBonusMCresult(bonusMC, pw, fsi, 0.5*(data::Q2[Qindex]+data::Q2[Qindex+1]),0.5*(data::W[Windex]+data::W[Windex+1]), 
				      data::Ebeam[0], 0.5*(data::ps[i]+data::ps[i+1]), costheta, proton, 0, lc);
	  if(!isnan(fsi)&&!(bonusMC==0.)){ f+=pow((get_normfit_bonus(0,Qindex,Windex,i,offshellset,lc)
						  *fsi/bonusMC-result)/error,2.); dof++;}
	  cout << costheta << " " << bonusMC << " " << pw << " " << fsi << " " << result << " " << get_normfit_bonus(0,Qindex,Windex,i,offshellset,lc)
						  *fsi/bonusMC << " " << f << endl;
	}
      }
      if(Qindex>0){
	error=data::bonusdata5[Qindex-1][Windex][i][j][2];
	result=data::bonusdata5[Qindex-1][Windex][i][j][1];
	costheta=data::bonusdata5[Qindex-1][Windex][i][j][0];
	if(result!=0.){
	  double bonusMC=0., pw=0.,fsi=0.;
	  DeepsCross.getBonusMCresult(bonusMC, pw, fsi, 0.5*(data::Q2[Qindex]+data::Q2[Qindex+1]),0.5*(data::W[Windex]+data::W[Windex+1]), 
				      data::Ebeam[1], 0.5*(data::ps[i]+data::ps[i+1]), costheta, proton, 0, lc);
	  if(!isnan(fsi)&&!(bonusMC==0.)){ f+=pow((get_normfit_bonus(1,Qindex,Windex,i,offshellset,lc)
						  *fsi/bonusMC-result)/error,2.); dof++;}
	  cout << costheta << " " << bonusMC << " " << pw << " " << fsi << " " << result << " " << get_normfit_bonus(1,Qindex,Windex,i,offshellset,lc)
						  *fsi/bonusMC << " " << f << endl;
	}
      }
    }
  }
  // Function to minimize (chi^2)
  //  f = GetChiSquaredOfVertex(par) // your fitness function goes here: typically ~ sum_i {(model(par,i)-data(i))^2 / error(i)^2} 
  f/=(dof-npar);
  cout << f << " " << dof-npar << endl;
}


int main(int argc, char *argv[])
{
  
  Qindex = atoi(argv[2]); // parse from argv or something
  Windex = atoi(argv[3]); // parse from argv or something
  int fixparam = atoi(argv[1]); // parse from argv or something
//   startset=atoi(argv[4]);
//   stopset=atoi(argv[5]);
  offshellset=atoi(argv[4]);
//   looplimit = atoi(argv[7]);
  lc=atoi(argv[5]); //lc or vna density
  
  int testing = 0;
  int fNDim = 3; // number of dimensions
  double fLo[] = {0.,1.,-1.}; // lower limits of params
  double fHi[] = {100.,20.,1.}; // upper limits of params
  char* fName[] = {"sigma_tot","beta","epsilon"};
  
  if(offshellset==1){
    fLo[2]=0.5;
    fHi[2]=2.;
    fName[2]="Lambda";
  }
  if(offshellset==2){
    fLo[2]=1.;
    fHi[2]=20.;
    fName[2]="beta_off";
  }
  
  int fBound[] = {1,1,1};

  TVirtualFitter *gMinuit = TVirtualFitter::Fitter ( 0, fNDim );

  // Start values of parameters
  // If you have a starting individual, you can simply use 
  double minuitIndividual[] = {40.,8.,-0.5}; // FIXME insert your starting individual ( double array) here
  if(offshellset==1) minuitIndividual[2]=1.2;
  if(offshellset==2) minuitIndividual[2]=8.;
  
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
  if(offshellset==2){
    for ( int i = fNDim-fixparam-1 ; i < gMinuit->GetNumberTotalParameters()-1 ; ++i ) gMinuit->FixParameter ( i );    
  }
  else{
    for ( int i = fNDim-fixparam ; i < gMinuit->GetNumberTotalParameters() ; ++i ) gMinuit->FixParameter ( i );
  }
  
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
  delete gMinuit;
  return 0;
}



double get_normfit_bonus(int beamindex, int Qindex, int Windex, int psindex, int offshell, bool lc){
  if(lc){
    if(offshell==3){
      if(beamindex==0) return bonusfits::normfits_fix3_off3_lc1_beam4[Qindex][Windex][psindex];
      else return bonusfits::normfits_fix3_off3_lc1_beam5[Qindex-1][Windex][psindex];
    }
    if(offshell==4){
      if(beamindex==0) return bonusfits::normfits_fix3_off4_lc1_beam4[Qindex][Windex][psindex];
      else return bonusfits::normfits_fix3_off4_lc1_beam5[Qindex-1][Windex][psindex];
    }
  }
  else{
    if(offshell==3){
      if(beamindex==0) return bonusfits::normfits_fix3_off3_lc0_beam4[Qindex][Windex][psindex];
      else return bonusfits::normfits_fix3_off3_lc0_beam5[Qindex-1][Windex][psindex];
    }
    if(offshell==4){
      if(beamindex==0) return bonusfits::normfits_fix3_off4_lc0_beam4[Qindex][Windex][psindex];
      else return bonusfits::normfits_fix3_off4_lc0_beam5[Qindex-1][Windex][psindex];
    }
  }
  return 0./0.;
}
