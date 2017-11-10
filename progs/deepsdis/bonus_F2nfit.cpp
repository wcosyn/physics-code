//program to do on-shell extrapolation of F2n using the bonus data
//we use renormalization for the data derived from our fits of the FSI model at certain kinematics
//different strategies have been tried

//sumfile etc. is for an independent fit of the results over all angles, similar to using weighted averages

#include <TVirtualFitter.h>
#include <TFitter.h>
#include <TMatrixT.h>
#include <TMinuit.h>
#include <TMath.h>


#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <cmath>
#include <cassert>

#include "bonusdata.h"
#include "bonusfits2.h"

#include <DeuteronCross.hpp>
#include <NuclStructure.hpp>

using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::abs;
using std::ofstream;

vector<double> tarray;
vector<double> data_array;
vector<double> error_array;

vector<double> thetavals;
vector<double> thetasigma;

double get_normfit_array(string wf, string struc, int offshell, bool lc, bool q2dep, bool pw, int beamindex,
			 int Qindex, int Windex, int psindex);
double sigmaparam(double W, double Q2, bool Q2dep);
void get_data(double &data, double& error, int beamindex, int Qindex, int Windex, int psindex, int cosindex);
void get_norm_final(bool q2dep, int beamindex,int psindex,string struc, double &norm, double &errornorm);


void Fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // Function to minimize (chi^2)
  f= 0.; // easily minimizable function as an example...
  //  f = GetChiSquaredOfVertex(par) // your fitness function goes here: typically ~ sum_i {(model(par,i)-data(i))^2 / error(i)^2} 
  for(unsigned int i=0;i<tarray.size();i++) f+=pow((par[0]*pow(tarray[i],2.)-par[1]*tarray[i]
					  +par[2]-data_array[i])/error_array[i],2.);
//   f/=tarray.size()==npar?1:(tarray.size()-npar);
}

void Fcn2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // Function to minimize (chi^2)
  f= 0.; // easily minimizable function as an example...
  //  f = GetChiSquaredOfVertex(par) // your fitness function goes here: typically ~ sum_i {(model(par,i)-data(i))^2 / error(i)^2} 
  for(unsigned int i=0;i<thetavals.size();i++) f+=pow((par[0]-thetavals[i])/thetasigma[i],2.);
//   f/=tarray.size()==npar?1:(tarray.size()-npar);
}


int main(int argc, char *argv[])
{

  int offshellset = atoi(argv[1]);
  bool q2dep = atoi(argv[2]);

  bool lc=atoi(argv[4]);
  double betain=8.;  
  bool proton=0;
  double phi=0.;
  string wf=argv[5];
  string strucname=argv[6];
  bool pw=atoi(argv[3]);
  
  
  string outdir="/home/wim/Calculations/Bonus/results/parabfitsApr16/";
  ofstream generalfile;
  ofstream sumfile;
  generalfile.open(outdir+wf+"."+strucname+".pw"+std::to_string(pw)+
	  ".lc"+std::to_string(lc)+".q2dep"+std::to_string(q2dep)+".off"+std::to_string(offshellset)+".fitsrig.dat");
  sumfile.open(outdir+wf+"."+strucname+".pw"+std::to_string(pw)+
	  ".lc"+std::to_string(lc)+".q2dep"+std::to_string(q2dep)+".off"+std::to_string(offshellset)+".fitsrig_indepfit.dat");
  for(int Qindex=0;Qindex<3;Qindex++){
    double Q2 = 0.5*(data::Q2[Qindex]+data::Q2[Qindex+1]);
    for(int Windex=0;Windex<5;Windex++){
      double Wprime = 0.5*(data::W[Windex]+data::W[Windex+1]);
      DeuteronCross test(wf,proton,strucname,sigmaparam(Wprime,Q2,q2dep),betain,-0.5,8.,1.2,offshellset,1E03);
      double xref=Q2/(Wprime*Wprime-MASSN*MASSN+Q2);
      NuclStructure strfunction(proton,Q2,xref,0,strucname);
      NuclStructure strfunction2(!proton,Q2,xref,0,strucname);
      double F2refn, F2refp;
      F2refn=strfunction.getF2(); 
      F2refp=strfunction2.getF2();
      thetavals.clear();
      thetasigma.clear();
      for(int cosindex=0;cosindex<10;cosindex++){
	double costhetar=-0.9+cosindex*0.2;
	tarray.clear();
	data_array.clear();
	error_array.clear();
	for(int psindex=0;psindex<4;psindex++){
	  double ps = 0.5*(data::ps[psindex]+data::ps[psindex+1]);
	  
	  double Er=sqrt(ps*ps+MASSP*MASSP);
	  double Einoff=MASSD-Er;
	  double massoff=sqrt(Einoff*Einoff-ps*ps);
	  double tprime = MASSN*MASSN-massoff*massoff;
	  
	  int beamindex = 0;
	  double Ebeam = data::Ebeam[beamindex]; //[0,1]
	 
	  double result, error, norm, errornorm;
// 	  double norm=get_normfit_array(wf,strucname, 
// 					offshellset,lc,q2dep,pw,beamindex,Qindex,Windex,psindex);
	  get_norm_final(q2dep, beamindex,psindex,strucname,norm,errornorm);
// 	  cout << norm << endl;
	  get_data(result, error, beamindex,Qindex,Windex,psindex,cosindex);
	  double ratio=test.getBonus_extrapratio(Q2,Wprime,Ebeam,ps,costhetar,proton,lc);
	  if(abs(ratio)>1.E-09){ //unphysical points yield ratio 0
// 	    cout << result << " " << error << " " << error/result << " " << norm << " " << errornorm << " " << errornorm/norm << " ";
// 	    error= error/result+errornorm/norm;
	    result*=ratio/norm;
	    error*=ratio/norm;
// 	    cout << result <<  " " << error << " " << error/result << endl;
	    tarray.push_back(tprime*1.E-06);
	    data_array.push_back(result);
	    error_array.push_back(error);
	  }
// 	  else cout << ratio << endl;
	  if(Qindex>0){
	    beamindex=1;
	    Ebeam=data::Ebeam[beamindex];
// 	    norm=get_normfit_array(wf,strucname, 
// 					offshellset,lc,q2dep,pw,beamindex,Qindex,Windex,psindex);
	    get_norm_final(q2dep, beamindex,psindex,strucname,norm,errornorm); 
	    get_data(result, error, beamindex,Qindex,Windex,psindex,cosindex);
	    ratio=test.getBonus_extrapratio(Q2,Wprime,Ebeam,ps,costhetar,proton,lc);
	    if(abs(ratio)>1.E-09){ //unphysical points yield ratio 0
// 	      cout << result << " " << error << " " << error/result << " " << norm << " " << errornorm << " " << errornorm/norm << " ";
// 	      error= error/result+errornorm/norm;
	      result*=ratio/norm;
	      error*=ratio/norm;
// 	      cout << result <<  " " << error << " " << error/result << endl;
	      tarray.push_back(tprime*1.E-06);
	      data_array.push_back(result);
	      error_array.push_back(error);
	    }
	  }
	}
	if(tarray.size()>3){ //should have enough points for a fit!
	  int testing = 0;
	  int fNDim = 3; // number of dimensions
	  int fRandomseed = 12345; // parse from argv or something
	  double fLo[] = {0.,0.,0.}; // lower limits of params
	  double fHi[] = {100.,10000.,3.}; // upper limits of params
	  char* fName[] = {"a","b","c"};
	  
	  int fBound[] = {1,1,1};
	
	    
	  TVirtualFitter *gMinuit = TVirtualFitter::Fitter ( 0, fNDim );
	
	  // Start values of parameters
	  // If you have a starting individual, you can simply use 
	  double minuitIndividual[] = {10.,1.,0.1}; // FIXME insert your starting individual ( double array) here
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
	  arglist[0] = (Double_t) fRandomseed; 
	  while (arglist[0] >= 900000000)  arglist[0] -= 900000000; 
	  while (arglist[0] <= 10000)  arglist[0] += 10000;
	  gMinuit->ExecuteCommand ( "SET RAN", arglist, 1 );
	  
	  // Set the limits for each parameter...
	  double step; // step size
	  for ( int i = 0 ; i < fNDim ; i++ ){
	    step = 0.01 * minuitIndividual[i];  // step size is 1% of value - so don't set value to 0 ;-)
	    if ( step < 0 ) step *= -1; // make it positive
	    else if ( step == 0.0 ) step = 1.0e-3;
	    gMinuit->SetParameter ( i, fName[i], minuitIndividual[i], step, fLo[i], fHi[i] );
	  }
	  
	  //   gMinuit->FixParameter ( 1 );
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
	  while ( outerMinimizationLoop ){
	
	    /* The inner minimization loop:
	    * - Using command: MIGRAD (and SIMPLEX)
	    * - SET STRATEGY 1
	    * - max. 7000 function calls
	    * untill convergence */
	    while ( innerMinimizationLoop ){
	
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
	      if(!static_cast<TFitter*> ( gMinuit->GetFitter() )->
		    GetMinuit()->fCstatu.CompareTo ( "CONVERGED " ) ){
		// ...try to IMPROVE the chi-squared
		arglist[0] = 7000;
		gMinuit->ExecuteCommand ( "IMPROVE", arglist, 1 );
	
		// we exit the inner loop
		innerMinimizationLoop = false;
	      }
	      // Else...
	      else{
		// ...try to find a better spot in parameter space
		arglist[0] = 7000;
		arglist[1] = 0.01;
		gMinuit->ExecuteCommand ( "SIMPLEX", arglist, 2 );
	
		// if SIMPLEX did not converge, we give up!
		if ( static_cast<TFitter*> ( gMinuit->GetFitter() )->
		      GetMinuit()->fCstatu.CompareTo ( "CONVERGED " ) ){
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
			  GetMinuit()->fISW[1]].CompareTo ( "ERR MATRIX NOT POS-DEF" ) ){
	      for ( int i = 0 ; i < fNDim ; ++i ){    // loop over parameters
		if (fBound[i]){
		//if ( GetTCalc()->GetLimit ( i ).bound )   // check if bound
		  value = gMinuit->GetParameter ( i );
		  verr  = gMinuit->GetParError ( i );
	
		  // If the parameter is within half an error bar of one of
		  // the limits, it should be fixed, as it may be that the
		  // true minimum lies beyond the defined physical limit
		  if ( TMath::Abs ( value - fLo[i] ) < verr*0.5 || // half a sigma
			TMath::Abs ( value - fHi[i] ) < verr*0.5 ){
		    gMinuit->FixParameter ( i );
		    someParsFixed++;
		    std::cout << "Parameter no." << i + 1
		    << " has been fixed: Too close to limits!\n";
		  }
		} // if parameter bound
	      } // end loop over parameters
	    } // status check error matrix
	
	    // If some parameters had to be fixed, re-run the loop...
	    if ( someParsFixed ){
	      someParsFixed = 0; // reset
	      std::cout << "Some parameters fixed, rerun Migrad()...\n";
	    }
	    else outerMinimizationLoop = false; // exit outer loop
	  } // end outer loop
      
	  // Calculate the covariance matrix again...
	  gMinuit->ExecuteCommand ( "MINOS", arglist, 0 );
	  gMinuit->ExecuteCommand ( "SHO COV", arglist, 0 );
	
	  // Extract the covariance matrix (method GetCovarianceMatrix only
	  // returns a double*)
	  nFreePars = gMinuit->GetNumberFreeParameters();
	  TMatrixT<double> matrix ( nFreePars, nFreePars );
	
	  for ( int i = 0 ; i < nFreePars ; ++i ){
	    for ( int j = 0 ; j < nFreePars ; ++j ){
	      matrix ( i, j ) = gMinuit->GetCovarianceMatrixElement ( i, j );
	    }
	  }
	  cout << "\nCleaning up..." << endl;
// 	  cout << Qindex << " " << Windex << " " << cosindex << " " << gMinuit->GetParameter(0) << " "<< gMinuit->GetParameter(1) << " " << gMinuit->GetParameter(2) << endl;
	  double eplus[nFreePars],eminus[nFreePars],eparab[nFreePars],globcc[nFreePars], params[nFreePars];
	  for(unsigned int i=0;i<nFreePars;++i){
	    gMinuit->GetErrors(i,eplus[i],eminus[i],eparab[i],globcc[i]);
	    params[i]=gMinuit->GetParameter(i);
// 	    cout << gMinuit->GetParameter(i) << " " << eplus[i] << " " << eminus[i] << 
// 	    " " << eparab[i] << " " << globcc[i] << endl;
	  }
	  ofstream myfile;
	  myfile.open(outdir+wf+"."+strucname+".pw"+std::to_string(pw)+
	  ".lc"+std::to_string(lc)+".q2dep"+std::to_string(q2dep)+".off"+std::to_string(offshellset)+
	  ".Q"+std::to_string(Qindex)
	    +".W"+std::to_string(Windex)+".th"+std::to_string(cosindex)+".rigfit.dat");
// 	  for(unsigned int i=0;i<nFreePars;++i) myfile << gMinuit->GetParameter(i) << " " << eparab[i] << " ";
// 	  myfile << F2refn << " " << F2refp << endl << endl << endl;
          myfile << "x\ty\tdy\t#points for pole extrapolation t[GeV-2] F2(ps)[]" << endl;
	  for(unsigned int i=0;i<tarray.size();i++) myfile << tarray[i] << " " << data_array[i] << " " << error_array[i] << endl;	  
	  myfile.close();
	  int n=3;
	  double f;
	  Fcn(n, &f, f, params, n);
	  generalfile << xref << " " << Q2*1.E-06 << " " << Wprime << " " << costhetar << " " 
	    << gMinuit->GetParameter(0) << " " << eparab[0] << " "
	    << gMinuit->GetParameter(1) << " " << eparab[1] << " "
	    << gMinuit->GetParameter(2) << " " << eparab[2] << " " << F2refn << " " << F2refp << " "
	    << f/(tarray.size()-3) << endl;
	  thetavals.push_back(gMinuit->GetParameter(2));
	  thetasigma.push_back(eparab[2]);
	  delete gMinuit;
	}
      }
      
      
      int testing = 0;
      int fNDim = 1; // number of dimensions
      int fRandomseed = 3454; // parse from argv or something
      double fLo[] = {0.}; // lower limits of params
      double fHi[] = {1.}; // upper limits of params
      char* fName[] = {"F2n"};
      
      int fBound[] = {1,1,1};
    
	
      TVirtualFitter *gMinuit2 = TVirtualFitter::Fitter ( 0, fNDim );
    
      // Start values of parameters
      // If you have a starting individual, you can simply use 
      double minuitIndividual[] = {0.2}; // FIXME insert your starting individual ( double array) here
      std::cout << "done" << endl;
    
      //------------------------------------------------------------------
      // The Minimisation section...
    
      std::cout << "Start minimizing...\n";
    
      // Initialize TMinuit via generic fitter (Minuit)
      //TVirtualFitter::SetDefaultFitter("Minuit2"); // if you want to use minuit2
      gMinuit2->SetFCN (  (void (*)(Int_t &, Double_t *, Double_t &, Double_t *, Int_t))   &Fcn2 );  // which function to minimise
    
      // ignore until ---------
      Double_t arglist[fNDim]; // For interface to Minuit routines
      Int_t ierflg = 0;  // ditto
    
      // Let Minuit produce minimal output...
      arglist[0] = 0;
      gMinuit2->ExecuteCommand ( "SET PRINT", arglist, 1 );
    
      // Set Minuit to the medium strategy
      arglist[0] = 1;
      gMinuit2->ExecuteCommand ( "SET STR", arglist, 1 );
      
      // If seed is given, set random seed.
      // This can be any integer between 10000 and 900000000.
      arglist[0] = (Double_t) fRandomseed; 
      while (arglist[0] >= 900000000)  arglist[0] -= 900000000; 
      while (arglist[0] <= 10000)  arglist[0] += 10000;
      gMinuit2->ExecuteCommand ( "SET RAN", arglist, 1 );
      
      // Set the limits for each parameter...
      double step; // step size
      for ( int i = 0 ; i < fNDim ; i++ ){
	step = 0.01 * minuitIndividual[i];  // step size is 1% of value - so don't set value to 0 ;-)
	if ( step < 0 ) step *= -1; // make it positive
	else if ( step == 0.0 ) step = 1.0e-3;
	gMinuit2->SetParameter ( i, fName[i], minuitIndividual[i], step, fLo[i], fHi[i] );
      }
      
      //   gMinuit2->FixParameter ( 1 );
      // Print out to verify...
      std::cout << "No. of parameters:\t\t"
      << gMinuit2->GetNumberTotalParameters() << endl;
      std::cout << "No. of free parameters:\t\t"
      << gMinuit2->GetNumberFreeParameters() << endl;
      std::cout << "Max. No. of iterations:\t\t"
      << gMinuit2->GetMaxIterations() << endl;
    
      /* **************** *
      * START MINIMIZING *
      * **************** */
    
      Double_t value, verr; // value of parameter and its error
      Int_t nFreePars = gMinuit2->GetNumberFreeParameters();// # of free parameters
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
      while ( outerMinimizationLoop ){
    
	/* The inner minimization loop:
	* - Using command: MIGRAD (and SIMPLEX)
	* - SET STRATEGY 1
	* - max. 7000 function calls
	* untill convergence */
	while ( innerMinimizationLoop ){
    
	  // Run minimization...
	  // ...untill convergence is reached (max 4 times)
	  arglist[0] = 7000;
	  gMinuit2->ExecuteCommand ( "MIGRAD", arglist, 1 );
	  for ( int i = 0; i < 3 &&
		static_cast<TFitter*> ( gMinuit2->GetFitter() )->
		GetMinuit()->fCstatu.CompareTo ( "CONVERGED " );
		i++ )
	    gMinuit2->ExecuteCommand ( "MIGRAD", arglist, 1 );
    
	  // When MIGRAD converged...
	  if(!static_cast<TFitter*> ( gMinuit2->GetFitter() )->
		GetMinuit()->fCstatu.CompareTo ( "CONVERGED " ) ){
	    // ...try to IMPROVE the chi-squared
	    arglist[0] = 7000;
	    gMinuit2->ExecuteCommand ( "IMPROVE", arglist, 1 );
    
	    // we exit the inner loop
	    innerMinimizationLoop = false;
	  }
	  // Else...
	  else{
	    // ...try to find a better spot in parameter space
	    arglist[0] = 7000;
	    arglist[1] = 0.01;
	    gMinuit2->ExecuteCommand ( "SIMPLEX", arglist, 2 );
    
	    // if SIMPLEX did not converge, we give up!
	    if ( static_cast<TFitter*> ( gMinuit2->GetFitter() )->
		  GetMinuit()->fCstatu.CompareTo ( "CONVERGED " ) ){
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
	if ( !static_cast<TFitter*> ( gMinuit2->GetFitter() )->GetMinuit()->
	      fCovmes[static_cast<TFitter*> ( gMinuit2->GetFitter() )->
		      GetMinuit()->fISW[1]].CompareTo ( "ERROR MATRIX ACCURATE " ) ||
	      !static_cast<TFitter*> ( gMinuit2->GetFitter() )->GetMinuit()->
	      fCovmes[static_cast<TFitter*> ( gMinuit2->GetFitter() )->
		      GetMinuit()->fISW[1]].CompareTo ( "ERR MATRIX NOT POS-DEF" ) ){
	  for ( int i = 0 ; i < fNDim ; ++i ){    // loop over parameters
	    if (fBound[i]){
	    //if ( GetTCalc()->GetLimit ( i ).bound )   // check if bound
	      value = gMinuit2->GetParameter ( i );
	      verr  = gMinuit2->GetParError ( i );
    
	      // If the parameter is within half an error bar of one of
	      // the limits, it should be fixed, as it may be that the
	      // true minimum lies beyond the defined physical limit
	      if ( TMath::Abs ( value - fLo[i] ) < verr*0.5 || // half a sigma
		    TMath::Abs ( value - fHi[i] ) < verr*0.5 ){
		gMinuit2->FixParameter ( i );
		someParsFixed++;
		std::cout << "Parameter no." << i + 1
		<< " has been fixed: Too close to limits!\n";
	      }
	    } // if parameter bound
	  } // end loop over parameters
	} // status check error matrix
    
	// If some parameters had to be fixed, re-run the loop...
	if ( someParsFixed ){
	  someParsFixed = 0; // reset
	  std::cout << "Some parameters fixed, rerun Migrad()...\n";
	}
	else outerMinimizationLoop = false; // exit outer loop
      } // end outer loop
  
      // Calculate the covariance matrix again...
      gMinuit2->ExecuteCommand ( "MINOS", arglist, 0 );
      gMinuit2->ExecuteCommand ( "SHO COV", arglist, 0 );
    
      // Extract the covariance matrix (method GetCovarianceMatrix only
      // returns a double*)
      nFreePars = gMinuit2->GetNumberFreeParameters();
      TMatrixT<double> matrix ( nFreePars, nFreePars );
    
      for ( int i = 0 ; i < nFreePars ; ++i ){
	for ( int j = 0 ; j < nFreePars ; ++j ){
	  matrix ( i, j ) = gMinuit2->GetCovarianceMatrixElement ( i, j );
	}
      }
      cout << "\nCleaning up..." << endl;
// 	  cout << Qindex << " " << Windex << " " << cosindex << " " << gMinuit2->GetParameter(0) << " "<< gMinuit2->GetParameter(1) << " " << gMinuit2->GetParameter(2) << endl;
      double eplus[nFreePars],eminus[nFreePars],eparab[nFreePars],globcc[nFreePars], params[nFreePars];
      for(unsigned int i=0;i<nFreePars;++i){
	gMinuit2->GetErrors(i,eplus[i],eminus[i],eparab[i],globcc[i]);
	params[i]=gMinuit2->GetParameter(i);
// 	    cout << gMinuit2->GetParameter(i) << " " << eplus[i] << " " << eminus[i] << 
// 	    " " << eparab[i] << " " << globcc[i] << endl;
      }
      int n=1;
      double f;
      Fcn2(n, &f, f, params, n);
      double avgform=0.,sigmaform=0.;
      for(size_t i=0;i<thetavals.size();i++){
	avgform+=thetavals[i]/thetasigma[i]/thetasigma[i];
	sigmaform+=1./thetasigma[i]/thetasigma[i];
      }
      avgform/=sigmaform;
      sigmaform=sqrt(1./sigmaform);
      sumfile << xref << " " << Q2*1.E-06 << " " << Wprime << " "
	<< gMinuit2->GetParameter(0) << " " << eparab[0] << " " << F2refn << " " << F2refp << " "
	<< f/(thetavals.size()-1.) << " " << avgform << " " << sigmaform << endl;
      delete gMinuit2;
      
      
      
      generalfile << endl << endl;
    }
    generalfile << endl << endl;
  }
  generalfile.close();
}
  
  

//   return 0;
// }

void get_norm_final(bool q2dep, int beamindex,int psindex, string struc, double & norm, double &errornorm){
  norm=1.;
  errornorm=0.;
  return;
  if(!struc.compare("CB")){
    if(q2dep){
      if(beamindex==0) {
	norm=bonusfits2::normfitsrig_final_q2dep1_beam4_setCB[psindex];
	errornorm=bonusfits2::errornormfitsrig_final_q2dep1_beam4_setCB[psindex];
      }
      else{
	norm=bonusfits2::normfitsrig_final_q2dep1_beam5_setCB[psindex];
	errornorm=bonusfits2::errornormfitsrig_final_q2dep1_beam5_setCB[psindex];
      }
    }
    else{
      if(beamindex==0) {
	norm=bonusfits2::normfits_final_q2dep0_beam4_setCB[psindex];
	errornorm=bonusfits2::errornormfits_final_q2dep0_beam4_setCB[psindex];
      }
      else{
	norm=bonusfits2::normfits_final_q2dep0_beam5_setCB[psindex];
	errornorm=bonusfits2::errornormfits_final_q2dep0_beam5_setCB[psindex];
      }	      
    } 
  }
  if(!struc.compare("SLAC")){
    if(q2dep){
      if(beamindex==0) {
	norm=bonusfits2::normfits_final_q2dep1_beam4_setSLAC[psindex];
	errornorm=bonusfits2::errornormfits_final_q2dep1_beam4_setSLAC[psindex];
      }
      else{
	norm=bonusfits2::normfits_final_q2dep1_beam5_setSLAC[psindex];
	errornorm=bonusfits2::errornormfits_final_q2dep1_beam5_setSLAC[psindex];
      }
    }
    else{
      if(beamindex==0) {
	norm=bonusfits2::normfits_final_q2dep0_beam4_setSLAC[psindex];
	errornorm=bonusfits2::errornormfits_final_q2dep0_beam4_setSLAC[psindex];
      }
      else{
	norm=bonusfits2::normfits_final_q2dep0_beam5_setSLAC[psindex];
	errornorm=bonusfits2::errornormfits_final_q2dep0_beam5_setSLAC[psindex];
      }	      
    } 
  }
}

//sigma in mbarn
double sigmaparam(double W, double Q2, bool Q2dep){
  if(!Q2dep||Q2<1.8E06){
      if(abs(W*W-1.232*1.232E06)<2.5E5) return 65;
  return (25.3*1.E-06*2.3+53*(sqrt(W*W>5.76E06?5.76E06:W*W)-MASSP)*1.E-03)
	/(1.8);
  }
  else{ 
    if(abs(W*W-1.232*1.232E06)<2.5E5) return 65.*1.8E06/Q2;
    return (25.3*1.E-06*2.3+53*(sqrt(W*W>5.76E06?5.76E06:W*W)-MASSP)*1.E-03)
	  /(1.E-06*Q2);
  }
}

void  get_data(double &data, double& error, int beamindex, int Qindex, int Windex, int psindex, int costhetaindex){
  if(beamindex==0){ 
    error=data::bonusdata4[Qindex][Windex][psindex][costhetaindex][2];
    data=data::bonusdata4[Qindex][Windex][psindex][costhetaindex][1];
  }
  else{
    error=data::bonusdata5[Qindex-1][Windex][psindex][costhetaindex][2];
    data=data::bonusdata5[Qindex-1][Windex][psindex][costhetaindex][1];
  }
}
