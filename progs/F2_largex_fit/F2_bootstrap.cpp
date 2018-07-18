//program used to fit the high x F2 JLab data
// takes two parameters: 1) prefix of filenames that have output of fits and data 2) include sys errors or not (1/0)
//run for instance  ~/physics-code/bin/F2_fit speedyevofitstat 0

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
#include <random>
#include <chrono>
using namespace std;

#include <constants.hpp>
#include "fomin_data.hpp"
#include "slac_data.hpp"


vector<double> xi_data; //nachtmann
vector<double> x_data; //bjorken x
vector<double> F2_data;
vector<double> F2_err;
vector<double> F2_sys;
vector<double> Q2_data;
int dof=0;
bool SLAC;


//from fortran evolution code!
extern "C"{
    void evolve_ns_adam(int *nx, double xx[], double HNS[], double *Q2i, double *Q2f, double *xi, double *A);

}

void Process_Data(bool sys, bool SLAC, bool hixcut, default_random_engine & generator); //selects the data to fit to, only C12 and x>0.5
void Fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag); //fit function

/**
 * @brief Computes the evolution of the F2 structure function, using the standard parametrization, which is assumed to have an initial scale of Q2i
 * 
 * @param x momentum fraction, or Nachtmann variable if computing from DIS kinematics
 * @param Q2f [GeV^2]  scale at which you want to evaluate F2
 * @param Q2i [GeV^2] initial scale for which the LT parametrization is valid
 * @param params parameters of the LT and HT fit
 * @return double F2(x,Q2f)
 */
double F2_evo(double x, double Q2f, double Q2i, double params[]);


/**
 * @brief interpolates precomputed LT grid with evolution
 * 
 * @param F2_grid LT part with evolution, precomputed
 * @param x_evo grid with x values
 * @param x momentum fraction, or Nachtmann variable
 * @param Q2 [GeV^2] scale
 * @return double 
 */
double F2_interp(double **F2_grid, double x_evo[], double x, double Q2, bool SLAC);

/**
 * @brief fill a grid with the LT Q2 evolution to interpolate in the fit function (speed up)
 * 
 * @param F2_grid index in Q2 (from 2 to 10, per 1), and x
 * @param x_evo values of x we use, momentum fraction or Nachtmann variable
 * @param params fit function parameters
 */
void fill_F2evo_grid(double **F2_grid, double x_evo[], double params[], bool SLAC);

void bootstrap_evo(double *par, string fileprefix, double reduced_chi2);

int main(int argc, char *argv[])
{
  
  bool sys=0;
  SLAC=0;
  bool HT=0;
  bool hixcut=0;
  
  string fn_prefix=argv[1];
  
//   if(string(argv[1])=="stat") { sys=0; fn_prefix=string(argv[1])+'.';}
//   else if(string(argv[1])=="sys") { sys=1; fn_prefix=string(argv[1])+'.';}
//   else { cerr << "invalid first commandline argument (sys or stat): " << argv[1]<<  endl; assert(1==0); }

//   if(string(argv[2])=="noSLAC") { SLAC=0; fn_prefix+=string(argv[2])+'.';}
//   else if(string(argv[2])=="wSLAC") { SLAC=1; fn_prefix+=string(argv[2])+'.';}
//   else { cerr << "invalid second commandline argument (noSLAC or wSLAC): "  << argv[2]<< endl; assert(1==0); }

//   if(string(argv[3])=="noHT") { HT=0; fn_prefix+=string(argv[3])+'.';}
//   else if(string(argv[3])=="wHT") { HT=1; fn_prefix+=string(argv[3])+'.';}
//   else { cerr << "invalid third commandline argument (noHT or wHT): "  << argv[3] << endl; assert(1==0); }

//   if(string(argv[4])=="noxcut") { hixcut=0; fn_prefix+=string(argv[4])+'.';}
//   else if(string(argv[4])=="xcut") { hixcut=1; fn_prefix+=string(argv[4])+'.';}
//   else { cerr << "invalid third commandline argument (noxcut or xcut): " << argv[4]<< endl; assert(1==0); }

//   if(SLAC && sys){
//     cerr << "When including SLAC data you should only run with stat errors as I don't have the sys ones" << endl;
//     assert(1==0);
//   }
//   cout << "prefix is " << fn_prefix << endl;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
  // for (unsigned i=0; i < xi_data.size(); i++){
  //   cout << i << " " << xi_data[i] << " " << x_data[i] << " " << Q2_data[i] << " " << F2_data[i] << " " << F2_err[i] << " " << F2_sys[i] << endl;
  //  }
  // exit(1);

  // int nn;
  // double ff;
  // double adamm[]={0.22969463,  4.44739662, -9.1538221, 0., 0., 0.};
  // Fcn(nn, &ff, ff, adamm, nn);
  for(int i=0;i<10000;i++){
    Process_Data(sys,SLAC,hixcut, generator);
    int testing = 0;
    int fNDim = 6; // number of dimensions
    double fLo[] = {-100.,-100.,-100.,-10.,-5.,-10.}; // lower limits of params
    double fHi[] = {100.,100.,100.,10.,5.,10.}; // upper limits of params
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

    // When HT is disabled, only fit LT parameters!
    if (!HT)
    {
        std::cout << "Running reduced capability for testing...\n";
        std::cout << "Fixing all but three parameters...\n";
        for ( int i = 3 ; i < gMinuit->GetNumberTotalParameters() ; ++i )
        {
        gMinuit->FixParameter ( i );
        }
    }

    // std::cout << "Fixing some parameters...\n";
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

    //Summary of my fit
    int n=nFreePars;
    double f;
    Fcn(n, &f, f, params, n);
    //cout << Warray[Windex] << " ";
    cout << "parameter values + parabolic errors" << endl;
    for(unsigned int i=0;i<nFreePars;++i) cout << params[i] << " " << eparab[i] << endl;
    cout << endl << endl << "chi2 per dof, chi2, dof, #freeparams" << endl;
    cout << f/(dof-nFreePars) << " " << f << " " << dof << " " << nFreePars << endl;

  
    delete gMinuit;
    bootstrap_evo(params,fn_prefix, f/(dof-nFreePars));
  }
  
  return 0;
}


// output bootstrap evoluted values for error estimates
void bootstrap_evo(double *par, string fileprefix, double reduced_chi2){
  double xvalues[8]={0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25};
  ofstream paramfile;
  paramfile.open(fileprefix+"param",std::ofstream::out | std::ofstream::app);
  paramfile << par[0] << " " << par[1] << " " << par[2] << " " << reduced_chi2 << endl;
  for(int x_it=0;x_it<8;x_it++){
    ofstream myfile;
    string filename = fileprefix+"x"+to_string(xvalues[x_it]);
    myfile.open(filename.c_str(),std::ofstream::out | std::ofstream::app);
    for(int i=0;i<=20;i++){
      double Q2=2.*pow(150.,i/20.);
      double nu=Q2/(2.*11.174862401 /12*xvalues[x_it]);
      double xi=2.*xvalues[x_it]/(1+sqrt(1.+Q2/nu/nu));
      //cout << "extra " << xvalues[x_it] << " " << xi << " " << Q2 << endl;
      myfile << F2_evo(xi,Q2,sqrt(18.),par) << " ";
    }
    myfile << endl;
    myfile.close();
  }



}

//selects the data to fit to, only C12 and x>0.5
void Process_Data(bool sys, bool SLAC, bool hixcut, default_random_engine & generator){
  xi_data.clear(); F2_data.clear(); F2_err.clear(); Q2_data.clear(); x_data.clear();
  
  for(int i=0;i<fomin_data::datapoints;i++){
    if(fomin_data::A[i]==12){
      double nu = (fomin_data::Ein[i]-fomin_data::Eout[i]);
      double Q2 = 4.*fomin_data::Ein[i]*fomin_data::Eout[i]*pow(sin(fomin_data::theta[i]*DEGRTORAD*0.5),2.);
      double x= Q2/(2.*11.174862401/12*nu);
      double nachtmann = 2.*x/(1.+sqrt(1+Q2/nu/nu));
      //cout << nu << " " << Q2*1.E-06 << " " << x  << " " << nachtmann << endl;
      if(x>0.5 && ((!hixcut)||x<1.7)){
        xi_data.push_back(nachtmann);
        Q2_data.push_back(Q2);
        normal_distribution<double> distribution(fomin_data::F2[i]*1.E06,fomin_data::F2_stat[i]*1.E06);
        F2_data.push_back(distribution(generator));
        F2_err.push_back(fomin_data::F2_stat[i]*1.E06);
        if(sys) F2_sys.push_back(fomin_data::F2_sys[i]*1.E06);  //change to zero if you do not want sys errors in the fit
        else F2_sys.push_back(0.);  //change to zero if you do not want sys errors in the fit
        x_data.push_back(x);
      }
    }

  }

  if(SLAC){
    for(int i=0;i<SLAC_data::datapoints;i++){
      double nu = SLAC_data::Q2[i]/(2.*11.174862401/12*SLAC_data::x[i]);
      double nachtmann = 2.*SLAC_data::x[i]/(1.+sqrt(1+SLAC_data::Q2[i]/nu/nu));
      if(SLAC_data::x[i]>0.5 && SLAC_data::Q2[i]>2. && ((!hixcut)||SLAC_data::x[i]<1.7)){
        xi_data.push_back(nachtmann);
        Q2_data.push_back(SLAC_data::Q2[i]);
        F2_data.push_back(SLAC_data::F2[i]*1.E06);
        F2_err.push_back(SLAC_data::F2_stat[i]*1.E06);
        F2_sys.push_back(0.);  //change to zero if you do not want sys errors in the fit
        x_data.push_back(SLAC_data::x[i]);
      }
    }
  }


}


//fit function
void Fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  f=0.;
  dof=0;
  //cout << npar << endl;
  //DeuteronCross DeepsCross("paris",proton,"SLAC",par[0],par[1],epsilon,betaoff,lambdain,offshellset,looplimit);
  // cout << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4] << " " << par[5] << " " << endl;

  int n_x=512;
  int nQ=SLAC?26*2:9*2;
  double **F2_grid= new double *[nQ];
  for(int i=0;i<nQ;i++) F2_grid[i]=new double[n_x];

  double x_evo[n_x];
  fill_F2evo_grid(F2_grid,x_evo,par,SLAC);


  for(unsigned i=0; i < xi_data.size(); i++){
    
    // f+=pow((F2_data[i]-exp(par[0]+par[1]*xi_data[i]+par[2]*pow(xi_data[i],2.))*
    //   (1.+par[3]*pow(xi_data[i],par[4])*(1.+xi_data[i]*par[5])/Q2_data[i]))/sqrt(pow(F2_err[i],2.)+pow(F2_sys[i],2.)),2.);
    
    // f+=pow((F2_data[i]-F2_evo(xi_data[i],Q2_data[i],sqrt(18.),par))/sqrt(pow(F2_err[i],2.)+pow(F2_sys[i],2.)),2.);
    // cout << xi_data[i] << " " << Q2_data[i] << " " << F2_interp(F2_grid,x_evo,xi_data[i],Q2_data[i],SLAC) << " " << F2_evo(xi_data[i],Q2_data[i],sqrt(18.),par) << endl;
    f+=pow((F2_data[i]-F2_interp(F2_grid,x_evo,xi_data[i],Q2_data[i],SLAC)*
            (1.+par[3]*pow(xi_data[i],par[4])*(1.+xi_data[i]*par[5])/Q2_data[i]))/sqrt(pow(F2_err[i],2.)+pow(F2_sys[i],2.)),2.);
    dof++;
  }
  for(int i=0;i<nQ;i++) delete [] F2_grid[i];
  delete [] F2_grid;
  // cout << "intermediate chi2 " << f/(dof-npar) << " " << dof << " " << npar << " " << f << endl;
}

//computes F2(x,Q2f) starting from fit at LT Q2i with params
double F2_evo(double x, double Q2f, double Q2i, double params[]){
  int n_x=512;
  double x_evo[n_x];
  double pdf[n_x];
  //double pdf_ini[n_x];
  double xi=0.;
  double A=12;
  double xi_low=0.3;
  for(int i=0;i<n_x;i++){
    x_evo[i]=xi_low*pow(12./xi_low,double(i)/(n_x-1.));
    pdf[i]=exp(params[0]+params[1]*x_evo[i]+params[2]*pow(x_evo[i],2.))/x_evo[i]; //pdf \propto F2/x, only LT part!!    
  }
  // for(int i=0;i<100;i++) cout << x_evo[i] << " " << pdf[i] << endl;
  evolve_ns_adam(&n_x, x_evo, pdf, &Q2i, &Q2f, &xi, &A);
  // cout << "after " << endl;

  // for(int i=0;i<100;i++) cout << x_evo[i] << " " << pdf[i] << endl;
  int index=log(x/xi_low)/log(12./xi_low)*(n_x-1.);
  // cout << x << " " << x_evo[index] << " " << x_evo[index+1] << " " << pdf[index]*x_evo[index] << " " << pdf[index+1]*x_evo[index+1] << endl;
  double result = ((x-x_evo[index])*pdf[index+1]+(x_evo[index+1]-x)*pdf[index])/(x_evo[index+1]-x_evo[index])*x
          *(1.+params[3]*pow(x,params[4])*(1.+x*params[5])/Q2f);  //multiply by x + HT part to get F2!

  if(std::isnan(result)) result=0.;
  if(std::isinf(result)) result=0.;
  return result;
}


double F2_interp(double **F2_grid, double x_evo[], double x, double Q2, bool SLAC){
  int n_x=512;
  double xi_low=0.3;
  int x_index=log(x/xi_low)/log(12./xi_low)*(n_x-1.);
  int Q_index=(Q2-2)*2;
  if(x_index==(n_x-1)) x_index--;
  if(Q_index==((SLAC?26*2:2*9)-1)) Q_index--;
  // cout << "index " << Q_index << " " << Q2 << " " << x_index << " " << x_evo[x_index] << " " << x_evo[x_index+1] << " " << x << endl;

  //2d interpolation
  return ((Q2-Q_index*0.5-2.)*((x-x_evo[x_index])*F2_grid[Q_index+1][x_index+1]+(x_evo[x_index+1]-x)*F2_grid[Q_index+1][x_index])+
          (Q_index*0.5+2.5-Q2)*((x-x_evo[x_index])*F2_grid[Q_index][x_index+1]+(x_evo[x_index+1]-x)*F2_grid[Q_index][x_index]))/(x_evo[x_index+1]-x_evo[x_index])/0.5;

}

void fill_F2evo_grid(double **F2_grid, double x_evo[], double params[], bool SLAC){
  int n_x=512;
  double xi=0.;
  double A=12;
  double xi_low=0.3;
  double base_F2[n_x];
  for(int i=0;i<n_x;i++){
    x_evo[i]=xi_low*pow(12./xi_low,double(i)/(n_x-1.));
    base_F2[i]=exp(params[0]+params[1]*x_evo[i]+params[2]*pow(x_evo[i],2.))/x_evo[i]; //pdf \propto F2/x, only LT part!!    
  }
  //here the F2_grid contains pdf values, and we do evolution per 1 GeV^2
  int nQ=9*2;
  if(SLAC) nQ=26*2;
  // double Q2i=sqrt(18);
  // double Q2f=2.;
  // evolve_ns_adam(&n_x, x_evo, F2_grid[0], &Q2i, &Q2f, &xi, &A);

  for(int i=0;i<nQ;i++){
    for(int j=0;j<n_x;j++) F2_grid[i][j]=base_F2[j];
    double Q2i=sqrt(18.);
    double Q2f=2.+i*0.5;
    evolve_ns_adam(&n_x, x_evo, F2_grid[i], &Q2i, &Q2f, &xi, &A);
  }
  //normalize so we have F2 again
  for(int i=0;i<nQ;i++){
    for(int j=0;j<n_x;j++) {
      F2_grid[i][j]*=x_evo[j];
      if(std::isnan(F2_grid[i][j])) F2_grid[i][j]=0.;
      if(F2_grid[i][j]>1.E10) F2_grid[i][j]=0.;
      // cout << i+2 << " " << x_evo[j] << " " << F2_grid[i][j] << endl;
    }
  }

  return;
}