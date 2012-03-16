#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include "FsiCorrelator.hpp"
#include <Utilfunctions.hpp>


#define TOTALGRIDPNTS (getCorr_rgrid()+1)*(getCorr_cthgrid()+1)

//constructor
FsiCorrelator::FsiCorrelator(MeanFieldNucleusThick *inputnucleus, const int rgrid, const int cthgrid, const string dir):
corr_rgrid(rgrid),corr_cthgrid(cthgrid),corr_phigrid(6),pnucleus(inputnucleus)
{
  cout << "Initializing Object for correlated glauber" << endl;
  //initialization
  invrstep=corr_rgrid/pnucleus->getRange();
  invcthstep=corr_cthgrid/2.;
  setCorrFilename(dir);
  corrgridfull=new double*[getCorr_rgrid()+1];
  corrgridproton=new double*[getCorr_rgrid()+1];
  corrgridneutron=new double*[getCorr_rgrid()+1];
  for(int i=0;i<=getCorr_rgrid();i++){
    corrgridfull[i]=new double[getCorr_cthgrid()+1];
    corrgridproton[i]=new double[getCorr_cthgrid()+1];
    corrgridneutron[i]=new double[getCorr_cthgrid()+1];
  }
  functionmatrix=NULL;
  x=NULL;
  ifstream infile(corrfilename.c_str(),ios::in|ios::binary);
  //check if object has been created sometime earlier and read it in
  if(infile.is_open()){
    cout << "Reading in correlation grids from memory: " << corrfilename << endl;
    for(int i=0;i<=corr_rgrid;i++){
      for(int j=0;j<=corr_cthgrid;j++){
	infile.read(reinterpret_cast<char *>(&corrgridfull[i][j]),sizeof(double));
	infile.read(reinterpret_cast<char *>(&corrgridproton[i][j]),sizeof(double));
	infile.read(reinterpret_cast<char *>(&corrgridneutron[i][j]),sizeof(double));
      }
    }
    infile.close();
    //printCorrGridAll();
  }
  //else construct it now!
  else{
    cout << "Constructing grids" << endl;
    constructCorrgrid(pnucleus->getTotalDensity(),corrgridfull);
    constructCorrgrid(pnucleus->getProtonDensity(),corrgridproton);
    constructCorrgrid(pnucleus->getNeutronDensity(),corrgridneutron);
    //printCorrGridAll();
    ofstream outfile(corrfilename.c_str(),ios::out|ios::binary);
    if(outfile.is_open()){
      cout << "Writing grids to file :" << corrfilename << endl;
      for(int i=0;i<=corr_rgrid;i++){
	for(int j=0;j<=corr_cthgrid;j++){
	  outfile.write(reinterpret_cast<char *>(&corrgridfull[i][j]),sizeof(double));
	  outfile.write(reinterpret_cast<char *>(&corrgridproton[i][j]),sizeof(double));
	  outfile.write(reinterpret_cast<char *>(&corrgridneutron[i][j]),sizeof(double));
	}
      }
      outfile.close();
    }
    else{
      cerr << "could not open file for writing corrgrid output: " << corrfilename << endl;
    }
  }

}

//mem cleanup
FsiCorrelator::~FsiCorrelator(){
  cout << "Cleaning up correlator class object" << endl;
  for(int i=0;i<=corr_rgrid;i++){
    delete []corrgridfull[i];
    delete []corrgridproton[i];
    delete []corrgridneutron[i];
  }
  delete [] corrgridfull;
  delete [] corrgridproton;
  delete [] corrgridneutron;    
}

//get corrfilename
const string FsiCorrelator::getCorrFilename() const{
  return corrfilename;
}

//get corr_rgrid
int FsiCorrelator::getCorr_rgrid() const{
  return corr_rgrid;
}

//get corr_thgrid
int FsiCorrelator::getCorr_cthgrid() const{
  return corr_cthgrid;
}

//get rstep
double FsiCorrelator::getInvRstep() const{
  return invrstep;
}

//get thstep
double FsiCorrelator::getInvCthstep() const{
  return invcthstep;
}

//get s_interp
double FsiCorrelator::getS_interp() const{
  return s_interp;
}

//get t_interp
double FsiCorrelator::getT_interp() const{
  return t_interp;
}

//get comp_s_interp
double FsiCorrelator::getComp_s_interp() const{
  return comp_s_interp;
}

//get comp_t_interp
double FsiCorrelator::getComp_t_interp() const{
  return comp_t_interp;
}

//get rindex
int FsiCorrelator::getRindex() const{
  return rindex;
}

//get thindex
int FsiCorrelator::getCthindex() const{
  return cthindex;
}

//set up inteprolation help variables for r
void FsiCorrelator::setRinterp(const double r){
  if(r>pnucleus->getRange()){
    cerr << "r out of range in getWave_Rpart: " << r << endl;
    exit(1);
  }
  rindex = int(floor(r*getInvRstep()));
  if(rindex==getCorr_rgrid()) rindex-=1;
  //if(rindex==(getCorr_rgrid()-1)) rindex-=1;
  s_interp = (r*getInvRstep() - (rindex));
  comp_s_interp=1.-s_interp;
}

//set up interpolation help variables for theta
void FsiCorrelator::setCthinterp(double costheta){
  if(costheta<0.) costheta=-costheta;
  if(abs(costheta)>1.) costheta=1.;
  cthindex = int(floor((1.-costheta)*getInvCthstep()));
  if(cthindex==getCorr_cthgrid()) cthindex-=1;
  t_interp = ((1.-costheta)*getInvCthstep() - (cthindex));  
  comp_t_interp=1.-t_interp;
}

//set r, theta and then get interpolation value of the grid
double FsiCorrelator::getCorrGridFull_interp(const double r, const double costheta){
  setRinterp(r); 
  setCthinterp(costheta);
  return getCorrGridFull_interp();
}

//set r, theta and then get interpolation value of the grid
double FsiCorrelator::getCorrGridProton_interp(const double r, const double costheta){
  setRinterp(r); 
  setCthinterp(costheta);
  return getCorrGridProton_interp();
}


//set r, theta and then get interpolation value of the grid
double FsiCorrelator::getCorrGridNeutron_interp(const double r, const double costheta){
  setRinterp(r); 
  setCthinterp(costheta);
  return getCorrGridNeutron_interp();
}

//set theta and then get interpolation value of the grid
double FsiCorrelator::getCorrGridFull_interp(const double costheta){
  setCthinterp(costheta);
  return getCorrGridFull_interp();
}

//set theta and then get interpolation value of the grid
double FsiCorrelator::getCorrGridProton_interp(const double costheta){
  setCthinterp(costheta);
  return getCorrGridProton_interp();
}


//set theta and then get interpolation value of the grid
double FsiCorrelator::getCorrGridNeutron_interp(const double costheta){
  setCthinterp(costheta);
  return getCorrGridNeutron_interp();
}

double FsiCorrelator::getCorrGrid_interp(const double r, const double costheta, bool proton){
 return proton? getCorrGridProton_interp(r,costheta) : getCorrGridNeutron_interp(r,costheta);
}

double FsiCorrelator::getCorrGrid_interp(const double costheta, bool proton){
 return proton? getCorrGridProton_interp(costheta) : getCorrGridNeutron_interp(costheta);
}

// get interpolation value of the grid
double FsiCorrelator::getCorrGridFull_interp() const{
  return comp_t_interp*(comp_s_interp*corrgridfull[rindex][cthindex]
	    + s_interp*corrgridfull[rindex+1][cthindex])
	    + t_interp*(comp_s_interp*corrgridfull[rindex][cthindex+1]
	    + s_interp*corrgridfull[rindex+1][cthindex+1]);  
}


// get interpolation value of the grid
double FsiCorrelator::getCorrGridProton_interp() const{
  return comp_t_interp*(comp_s_interp*corrgridproton[rindex][cthindex]
	    + s_interp*corrgridproton[rindex+1][cthindex])
	    + t_interp*(comp_s_interp*corrgridproton[rindex][cthindex+1]
	    + s_interp*corrgridproton[rindex+1][cthindex+1]);  
}


// get interpolation value of the grid
double FsiCorrelator::getCorrGridNeutron_interp() const{
//  cout << rindex << " " << cthindex << " " << s_interp << " " << t_interp << endl;
  return comp_t_interp*(comp_s_interp*corrgridneutron[rindex][cthindex]
	    + s_interp*corrgridneutron[rindex+1][cthindex])
	    + t_interp*(comp_s_interp*corrgridneutron[rindex][cthindex+1]
	    + s_interp*corrgridneutron[rindex+1][cthindex+1]);  
}
  
double FsiCorrelator::getCorrGrid_interp(bool proton) const{
 return proton? getCorrGridProton_interp() : getCorrGridNeutron_interp();
}
  
  
  
  
//set filename where the class object is stored or written to
void FsiCorrelator::setCorrFilename(const string dir){
  corrfilename=dir+"/input/CorrGrid."+pnucleus->getNucleusName()+".tot."+to_string(corr_rgrid)+".c"+to_string(corr_cthgrid)+".dat";  
}

void FsiCorrelator::printCorrGridFull() const{
  for(int i=0;i<=corr_rgrid;i++){
    for(int j=0;j<=corr_cthgrid;j++){
      cout << pnucleus->getRange()*i/corr_rgrid << " " << 1.-j/invcthstep << " " << corrgridfull[i][j] << endl;
    }
  }
  cout << endl << endl;
}
  
void FsiCorrelator::printCorrGridProton() const{
  for(int i=0;i<=corr_rgrid;i++){
    for(int j=0;j<=corr_cthgrid;j++){
      cout << pnucleus->getRange()*i/corr_rgrid << " " << 1.-j/invcthstep << " " << corrgridproton[i][j] << endl;
    }
  }
  cout << endl << endl;
}
  
void FsiCorrelator::printCorrGridNeutron() const{
  for(int i=0;i<=corr_rgrid;i++){
    for(int j=0;j<=corr_cthgrid;j++){
      cout << pnucleus->getRange()*i/corr_rgrid << " " << 1.-j/invcthstep << " " << corrgridneutron[i][j] << endl;
    }
  }
  cout << endl << endl;
}
  
void FsiCorrelator::printCorrGridAll() const{
  for(int i=0;i<=corr_rgrid;i++){
    for(int j=0;j<=corr_cthgrid;j++){
      cout << pnucleus->getRange()*i/corr_rgrid << " " << 1.-j/invcthstep << " " << corrgridfull[i][j]
      << " " << corrgridproton[i][j] << " " << corrgridneutron[i][j] << endl;
    }
  }
  cout << endl << endl;
}
  




void FsiCorrelator::constructCorrgrid(const double *density, double **corrgrid){

  cout << "Constructing correlation grid..." << endl;
  functionmatrix = new double *[TOTALGRIDPNTS];//TOTALGRIDPNTS= (getCorr_rgrid()+1(*(getCorr_thgrid()+1)
  for(int i=0;i<TOTALGRIDPNTS;i++) functionmatrix[i] = new double[TOTALGRIDPNTS];
  cout << "Construct Matrix of the functions..." << endl;
  constructMatrixSet(density);
  x = new double[TOTALGRIDPNTS];
  for(int i=0;i<TOTALGRIDPNTS;i++) x[i]=1.;
  int *check= new int;
  newt(check);
  cout << ReturnFunction() << endl;
  for(int i=0;i<TOTALGRIDPNTS;i++) delete [] functionmatrix[i];
  delete [] functionmatrix;
  if(*check>0) {
    delete [] x; 
    delete check;
    cerr << "no solution found for integration equation" << endl;
    exit(1);
  }
  else{
    for(int i=0;i<TOTALGRIDPNTS;i++){
      int j = i%(corr_cthgrid+1);
      int ii = i/(corr_cthgrid+1);
//       cout << i << " " << ii << " " << j << " " << x[i] << endl;
      corrgrid[ii][j]=x[i];
    }
    cout << endl << endl;
    delete [] x;
    delete check;
    return;
  }  
}

//////////////////////
//construct the matrix needed to solve the non-linear set of equations
////////////////////////

void FsiCorrelator::constructMatrixSet(const double *density){

  double VOLEL=pnucleus->getRange()/corr_rgrid*2./corr_cthgrid*PI/corr_phigrid; //phasespace
  //create grid, phi dependence is ignored as it is "very small"
  for(int rindex1=0; rindex1<=corr_rgrid; rindex1++){
    double r1 = pnucleus->getRange()*rindex1/corr_rgrid;
    for(int cthetaindex1=0; cthetaindex1<=corr_cthgrid; cthetaindex1++){
      double costheta1 = 1.-2.*cthetaindex1/corr_cthgrid;
/*      double costheta1 = cos(theta1);
      double sintheta1 = sin(theta1);*/
      //double costheta1,sintheta1;
      //sincos(theta1,&sintheta1,&costheta1);
      double sintheta1 = sqrt(1.-costheta1*costheta1);
      int index1 = (rindex1)*(corr_cthgrid+1)+cthetaindex1;
      for(int rindex2=0; rindex2<=corr_rgrid; rindex2++){
	double r2 = pnucleus->getRange()*rindex2/corr_rgrid;
	for(int cthetaindex2=0; cthetaindex2<=corr_cthgrid; cthetaindex2++){
	  double costheta2 = 1.-2.*cthetaindex2/corr_cthgrid;
// 	  cout << costheta1 << " " << costheta2 << endl;
	  //double costheta2,sintheta2;
	  //sincos(theta2,&sintheta2,&costheta2);
	  double sintheta2 = sqrt(1.-costheta2*costheta2);
	  int index2 = (rindex2)*(corr_cthgrid+1)+cthetaindex2;	  
	  functionmatrix[index2][index1] = sintheta1*pnucleus->getDensity(r1,density)*VOLEL;
	  double gsum = 0.;
	  for(int phiindex1=0; phiindex1<2*corr_phigrid; phiindex1++){
	    double phi1 = PI*phiindex1/corr_phigrid;
	    double cosphi1, sinphi1;
	    sincos(phi1,&sinphi1,&cosphi1);
	    for(int phiindex2=0; phiindex2<2*corr_phigrid; phiindex2++){
	      double phi2 = PI*phiindex2/corr_phigrid;
	      double cosphi2, sinphi2;
	      sincos(phi2,&sinphi2,&cosphi2);
	      gsum +=correlation(normr(r1, costheta1, sintheta1, cosphi1, sinphi1, r2, costheta2, sintheta2, cosphi2, sinphi2));
// 	      gsum +=correlation(normr(r1, -costheta1, sintheta1, cosphi1, sinphi1, r2, costheta2, sintheta2, cosphi2, sinphi2));
// 	      gsum +=correlation(normr(r1, costheta1, sintheta1, cosphi1, sinphi1, r2, -costheta2, sintheta2, cosphi2, sinphi2));
// 	      gsum +=correlation(normr(r1, -costheta1, sintheta1, cosphi1, sinphi1, r2, -costheta2, sintheta2, cosphi2, sinphi2));
		//begin and endpoint trapezium integration coefficients
		/*if(phiindex1==0) fillmatrix[index1][index2]/=2.;		
		if(phiindex1==2.*PHIGRID) fillmatrix[index1][index2]/=2.;
		if(thetaindex1==0) fillmatrix[index1][index2]/=2.;
		if(thetaindex1==THETAGRID) fillmatrix[index1][index2]/=2.;
		if(rindex1==0) fillmatrix[index1][index2]/=2.;
		if(rindex1==RGRID)  fillmatrix[index1][index2]/=2.;*/
	      }
	  }
	  functionmatrix[index2][index1] *=gsum;
	}
      }
    }
  }
}



/////////////////
//G-D correlation function
///////////////////////
  
double FsiCorrelator::correlation(double r) const{

  double a[11]={0.99994,0.99600E-01,-2.9385,28.513,-132.53,238.33,-151.54,-69.845,159.94,-88.043,16.930};
  if(r>1.24) return 1+0.1204071*exp(-2.5002*(r-1.24)*(r-1.24));
  else{
    double result=0.;
    for(int i=0;i<11;i++) result += a[i]*pow(r,i);
    return 1-result;
    }
  //return 1-exp(-r*r*1.18147448);
}

////
//some more num recipes shit i copied
///////


void FsiCorrelator::newt(int *check){

  int MAXITS=200;
  double TOLF = 1.0e-7, TOLMIN = 1.0e-06, TOLX = 1.0e-07, STPMX = 100.0;

  double *functionvalues = new double[TOTALGRIDPNTS];
  double function = ReturnFunction(functionvalues);
  double test=0.;
  //test if initial guess is zero point
  for(int i=0;i<TOTALGRIDPNTS;i++) if(abs(functionvalues[i])>test) test=abs(functionvalues[i]);
  if(test<TOLF) {
    delete [] functionvalues;
    cout << "initial point is zero " << test << endl;
    return;
  }
 
  double sum = 0.;
  for(int i=0;i<TOTALGRIDPNTS;i++) sum+=x[i]*x[i];
  double stpmax = STPMX*(sqrt(sum)>float(TOTALGRIDPNTS)? sqrt(sum):float(TOTALGRIDPNTS));
  for(int its=0;its<MAXITS;its++){
    cout << "Entering iteration " << its << endl;
    double **Jacobian = new double*[TOTALGRIDPNTS];
    for(int i=0;i<TOTALGRIDPNTS;i++) Jacobian[i] = new double[TOTALGRIDPNTS];
    double *grad = new double[TOTALGRIDPNTS];
    double *xold = new double[TOTALGRIDPNTS];
    double functionold;
    double *step = new double[TOTALGRIDPNTS];
    cout << "Constructing Jacobian...." << endl;
    ConstructJacobian(Jacobian);
    
    for(int i=0; i<TOTALGRIDPNTS; i++){
      sum=0.;
      for(int j=0;j<TOTALGRIDPNTS;j++) sum+=functionvalues[j]*Jacobian[j][i];
      grad[i] = sum;
    }
    for (int i=0;i<TOTALGRIDPNTS;i++) {xold[i] = x[i]; step[i] = -functionvalues[i];}
    functionold = function;
    cout << "Solving for step... " << TOTALGRIDPNTS << endl;
    LUsolve(Jacobian,step,TOTALGRIDPNTS);
    // for(int i=0;i<TOTALGRIDPNTS;i++) cout << step[i] << " ";
    cout << "Doing linesearch..." << endl;
    lnsrch(xold,functionold,grad,step,&function,stpmax,check,functionvalues);
    test=0.;
    for(int i=0;i<TOTALGRIDPNTS;i++) if(abs(functionvalues[i])>test) test = abs(functionvalues[i]);
    if(test<TOLF) {
      cout << "Convergence reached on function values " << endl;
      corrmaintenance(TOTALGRIDPNTS, functionvalues, Jacobian, grad, xold, step);
      *check=0;
      return;
    }
    
    if(*check){  //check for grad f equal to zero
      test=0.;
      double den=(function>0.5*TOTALGRIDPNTS? function : 0.5*TOTALGRIDPNTS);
      for(int i=0;i<TOTALGRIDPNTS;i++){
	double temp=abs(grad[i]) * (abs(x[i])>1.? abs(x[i]) : 1.)/den;
	if(temp>test) test=temp;
      }
      *check=(test < TOLMIN ? 1:0);
      cout << "grad f was equal to zero, spurious thingy" << endl;
      corrmaintenance(TOTALGRIDPNTS, functionvalues, Jacobian, grad, xold, step);
      return;
    }
    test=0.;
    for(int i=0;i<TOTALGRIDPNTS;i++){
      double temp=(abs(x[i]-xold[i]))/(abs(x[i])>1.? abs(x[i]) : 1.);
      if(temp>test) test=temp;
    }
    if(test < TOLX) {
      corrmaintenance(TOTALGRIDPNTS, functionvalues, Jacobian, grad, xold, step); 
      cout << "Convergence reached on dx" << endl;
      return;
    }
    delete [] grad;
    delete [] xold;
    delete [] step;
    for(int i=0;i<TOTALGRIDPNTS;i++) delete [] Jacobian[i];
    delete [] Jacobian;

  }
  
  cout << "Max iterations reached" << endl;
  *check = 2;
  delete [] functionvalues;  
  return;

}

/////////////////
//Jacobian is contructed, what did i need it for again? :)
/////////////////

void FsiCorrelator::ConstructJacobian(double **Jacobian){

  for(int i=0;i<TOTALGRIDPNTS;i++){
    for(int j=0;j<TOTALGRIDPNTS;j++){
      Jacobian[i][j]=functionmatrix[i][j]*x[i];
      if(i==j) {
	for(int k=0;k<TOTALGRIDPNTS;k++) Jacobian[i][i]+=functionmatrix[i][k]*x[k];
      }
    }
  }
}


///////////////
//i think this is another matrix needed to solve the non-lin set of equations :)
//aha, i think i remember, this is the function that needs minimizing!
/////////////////

double FsiCorrelator::ReturnFunction(double *functionvalues){
  for(int i=0;i<TOTALGRIDPNTS;i++){
    functionvalues[i]=-/*2.**/2*corr_phigrid;
    for(int j=0;j<TOTALGRIDPNTS;j++){
      double coeff = 1.;
      if(j%(corr_cthgrid+1)==0||j%(corr_cthgrid+1)==corr_cthgrid) coeff*=0.5;
      if(j/(corr_cthgrid+1)==0||j/(corr_cthgrid+1)==corr_rgrid) coeff*=0.5;
      functionvalues[i]+=functionmatrix[i][j]*x[i]*x[j]*coeff;
    }
  }
  double f=0.0;
  for(int i=0;i<TOTALGRIDPNTS;i++) f+=functionvalues[i]*functionvalues[i];
  f/=2.;
  return f;
}

double FsiCorrelator::ReturnFunction(){
  double *functionvalues = new double[TOTALGRIDPNTS];
  for(int i=0;i<TOTALGRIDPNTS;i++){
    functionvalues[i]=-1.*2*corr_phigrid;
    for(int j=0;j<TOTALGRIDPNTS;j++){
      functionvalues[i]+=functionmatrix[i][j]*x[i]*x[j];
    }
  }
  double f=0.0;
  for(int i=0;i<TOTALGRIDPNTS;i++) f+=functionvalues[i]*functionvalues[i];
  delete [] functionvalues;
  f/=2.;
  return f;
}

///////////////////
//memory maintenance
/////////////////////

void FsiCorrelator::corrmaintenance(int n, double *functionvalues, double **Jacobian, double *grad, double *xold, double *step){

  delete [] functionvalues;
  delete [] grad;
  delete [] xold;
  delete [] step;
  for(int i=0;i<n;i++) delete [] Jacobian[i];
  delete [] Jacobian;
  return;
}

///////////////
//taken from numerical recipes, needed in the procedure to solve non-lin set, by minimalizing a testfunction
///////////////////////

void FsiCorrelator::lnsrch(double *xold, double functionold, double *grad, double *step, double *functionnew,
	    double stpmax, int *check, double *functionvalues){

  double TOLX = 1.0e-07;
  double ALF = 1.0e-04; 

  *check=0;
  float sum=0.0;
  for(int i=0;i<TOTALGRIDPNTS;i++) sum += step[i]*step[i];
  sum=sqrt(sum);
  if(sum > stpmax) for(int i=0;i<TOTALGRIDPNTS;i++) step[i] *= stpmax/sum;  //scaling
  
  double slope=0.0;
  for(int i=0;i<TOTALGRIDPNTS;i++) slope += grad[i]*step[i];

  double test=0.;
  for(int i=0;i<TOTALGRIDPNTS;i++){
    double temp=abs(step[i])/(fabs(xold[i])>1.0? xold[i]:1.0);
    if (temp>test) test=temp;
  }
  double alamin=TOLX/test;  //compute minimum value of lambda
  double alam=1.0;  //start with full step

  for(;;){

    double alam2=0., functionnew2=0., functionold2=0., tmplam=0.;

    for(int i=0;i<TOTALGRIDPNTS;i++) x[i]=xold[i]+alam*step[i];
    *functionnew = ReturnFunction(functionvalues);
    if(alam <alamin) {
      for(int i=0;i<TOTALGRIDPNTS;i++) x[i] = xold[i];
      *check=1;
      return;
    }
    else if (*functionnew <= functionold + ALF*alam*slope) return;
    else{

      if(alam == 1.) tmplam = -slope/(2.0*(*functionnew-functionold-slope));
      else{
	double rhs1 = *functionnew-functionold-alam*slope;
	double rhs2 = functionnew2-functionold2-alam2*slope;
	double a = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	double b = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	if (a==0.) tmplam = -slope/(2.0*b);
	else{
	  double disc = b*b-3.0*a*slope;
	  if(disc<0.) cout << "roundoff problem in lnsrch" << endl;
	  else tmplam = (-b+sqrt(disc))/(3.0*a);
	}
	if(tmplam>0.5*alam) tmplam = 0.5*alam;
      }
    }
    alam2=alam;
    functionnew2=*functionnew;
    functionold2=functionold;
    alam=(tmplam>0.1*alam? tmplam : 0.1*alam);
  }
}


