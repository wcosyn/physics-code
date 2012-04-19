#include<cmath>
#include<iostream>
#include<cstdlib>

using namespace std;

#include "Utilfunctions.hpp"

////////////////////////////////////////
//calculates points between whom      //  
//the interpolation is performed.     //  
//lowest order interpolation          //
////////////////////////////////////////

double interpolate(const double *array, double r, double rstep,int lines, int offset){
  //interpolation points is double of order
  //rstep is set distance between known points
  //offset indicates times rstep from where the interpolation array starts
  //lines is number of elements in array

  if(array==NULL) {
    cerr <<"Array has not been initialized!!!!" << endl;
    exit(1);
  }
  int index=int(floor(r/rstep))-offset;

  if(index==lines-1) index--; // if index is at the end of the array.
  if(index<0) index=0;

  return ((r-rstep*(index+1+offset))*array[index] + (rstep*(index+offset)-r)*array[index+1]) / -rstep;
}

///////////////////////////////////
//some functions to calculate powers
//these are faster than regular pow function...
///////////////////////////////////////

double power(double x, int y){
  switch(y){
  case 0:
    return 1.;
    break;
  case 1:
    return x;
    break;
  case 2:
    return x*x;
    break;
  case 3:
    return x*x*x;
    break;
  default:
    if(y%2==0) {
      double temp = power(x, y/2);
      return temp*temp;
    }
    else{
      double temp = power(x, (y-1)/2);
      return x*temp*temp;
    }
  }
}


int power(int x, int y){
  switch(y){
  case 0:
    return 1;
    break;
  case 1:
    return x;
    break;
  case 2:
    return x*x;
    break;
  case 3:
    return x*x*x;
    break;
  default:
    if(y%2==0) {
      int temp = power(x, y/2);
      return temp*temp;
    }
    else{
      int temp = power(x, (y-1)/2);
      return x*temp*temp;
    }
  }
}

////////////////
//distance between 2 points
//////////////////

double normr(double r1, double costheta1, double sintheta1, double cosphi1, double sinphi1,
	     double r2, double costheta2, double sintheta2, double cosphi2, double sinphi2){
  return sqrt((r1*costheta1-r2*costheta2)*(r1*costheta1-r2*costheta2)+(r1*sintheta1*cosphi1-r2*sintheta2*cosphi2)*(r1*sintheta1*cosphi1-r2*sintheta2*cosphi2)+
    (r1*sintheta1*sinphi1-r2*sintheta2*sinphi2)*(r1*sintheta1*sinphi1-r2*sintheta2*sinphi2));
}


//////////////////////////////////////
//routine to solve matrix equation  //
//order is dim of square matrix     //
//numrec p 46                       //
//////////////////////////////////////


void LUsolve(double **a, double *b, int order){
  double sum,big,temp,dum;
  int imax=0;
  int *index=new int[order];
  //determine scaling for each row
  double *scaling = new double[order];
  for(int i=0;i<order;i++){
    big=0.0;
    for(int j=0;j<order;j++){
      temp=fabs(a[i][j]);
      if(temp > big) big=temp;
    }
    if (big==0.0) cout << "Singular Matrix!" << endl;
    scaling[i]=1.0/big;
  }
  //Crout's Method
  for (int j=0; j<order;j++){
    for(int i=0; i<j; i++){
      sum=a[i][j];
      for(int k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for(int i=j;i<order;i++){
      sum=a[i][j];
      for(int k=0;k<j;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      dum=scaling[i]*abs(sum);
      if(dum>=big){
	big=dum;
	imax=i;
      }
    }
    //  cout << a[25][65] << endl;
   
    //interchanges needed?
    if (j != imax){
      for(int k=0;k<order;k++){
	double dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      scaling[imax]=scaling[j];
    }
    index[j]=imax;
    if(a[j][j]==0)    a[j][j] = 1.0e-20;

    if(j!=order){
      for(int i=j+1;i<order;i++) a[i][j] /= a[j][j];
    }      
  }

  delete []scaling;

     
  //back and forward substitutions
  int ii=-1,ip;
  for (int i=0;i<order;i++) {
    ip=index[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii!=-1)
      for (int j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }

   
  for (int i=order-1;i>=0;i--) {
    sum=b[i];
    for (int j=i+1;j<order;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }

  delete [] index;
}


// void sincos(const double x, double * sin, double * cos){
//   if(fabs(x) < 1.E-03) {
//     *sin = x;
//     *cos = 1.0 - x*x/2.0;
//     return;
//   }
//   double s, c;
//   sincos(x/2., &s, &c);
//   *sin = 2.*s*c;
//   *cos = c*c - s*s;
// }


/////////////////////////////////
//algorythm to fit data.       //
//according to numrecipes p671 //
//base functions can be chosen,//
//here polynoms are taken.     //
//chi squared is also returned //
//via reference                //
/////////////////////////////////


double* fitdata(double *xarray, double *yarray, double *sigma, int lower, int upper, int order, double & chi2){
  int N = upper-lower+1; //# of data points to fit
  if (N<=order) cout << "Number of data points should be larger than total base functions" << endl;
  //initialise alpha and beta matrix
  double **alpha = matrix(order,order);
  for(int i=0;i<order;i++){
    for(int j=0;j<order;j++){
      alpha[i][j]=0.;
    }
  }
  double *beta = new double[order];
  for(int i=0;i<order;i++) beta[i]=0.;
  //fill them
  for(int i=0;i<order;i++){
    for(int j=0;j<order;j++){
      for(int k=lower;k<=upper;k++){
	alpha[i][j] += fitfunction(xarray[k],i)*fitfunction(xarray[k],j) / sigma[k]/sigma[k];
      }
    }
  }
  for(int i=0;i<order;i++){
    for(int j=lower;j<=upper;j++){
      beta[i] += yarray[j]*fitfunction(xarray[j],i) / sigma[j]/sigma[j];
    }
  }
  //solve the set
  LUsolve(alpha,beta,order);

  //calculation of the chi squared of the fit
  chi2=0.;
  for(int i=lower;i<=upper;i++){
    double temp=yarray[i];
    for(int j=0;j<order;j++) //temp-=beta[j]*cos(xarray[i]*j/range);
      temp-=beta[j]*fitfunction(xarray[i],j);
    temp/=sigma[i];
    chi2+=pow(temp,2);
  }

  
  //maintenance
  freematrix(alpha,order);
  return beta;
}

/////////////////////////////////
//set of base functions used in fit
///////////////////////////////////

double fitfunction(double x, int i){
  //if(i==order) return exp(-a*x);
   return pow(x,i);
}

//////////////////////
//create n*m matrix //
//////////////////////

double **matrix(int n, int m){
  double **matrix=new double*[n];
  for(int i=0;i<n;i++) matrix[i]=new double[m];
  return matrix;
}

////////////////
//free matrix //
////////////////

void freematrix(double **matrix, int n){
  for(int i=0;i<n;i++) delete [] matrix[i];
  delete[] matrix;
}



void intTest(double r, complex<double> *results, va_list ap){
}


complex<double> Interp3d(complex<double> ***grid, double s, double t, double u, 
			  double comps, double compt, double compu, int sindex, int tindex, int uindex){
  return complex<double>(compu*(compt*(comps*real(grid[sindex][tindex][uindex])
	+ s*real(grid[sindex+1][tindex][uindex])) 
	+ t*(comps*real(grid[sindex][tindex+1][uindex] )
	+ s*real(grid[sindex+1][tindex+1][uindex])))
	+  u*(compt*(comps*real(grid[sindex][tindex][uindex+1])
	+ s*real(grid[sindex+1][tindex][uindex+1])) 
	+ t*(comps*real(grid[sindex][tindex+1][uindex+1])
	+ s*real(grid[sindex+1][tindex+1][uindex+1]))),
	compu*(compt*(comps*imag(grid[sindex][tindex][uindex])
	+ s*imag(grid[sindex+1][tindex][uindex])) 
	+ t*(comps*imag(grid[sindex][tindex+1][uindex] )
	+ s*imag(grid[sindex+1][tindex+1][uindex])))
	+  u*(compt*(comps*imag(grid[sindex][tindex][uindex+1])
	+ s*imag(grid[sindex+1][tindex][uindex+1])) 
	+ t*(comps*imag(grid[sindex][tindex+1][uindex+1])
	+ s*imag(grid[sindex+1][tindex+1][uindex+1]))));
}
double Interp3d(double ***grid, double s, double t, double u, 
			  double comps, double compt, double compu, int sindex, int tindex, int uindex){
  return compu*(compt*(comps*grid[sindex][tindex][uindex]
	+ s*grid[sindex+1][tindex][uindex]) 
	+ t*(comps*grid[sindex][tindex+1][uindex] 
	+ s*grid[sindex+1][tindex+1][uindex]))
	+  u*(compt*(comps*grid[sindex][tindex][uindex+1]
	+ s*grid[sindex+1][tindex][uindex+1]) 
	+ t*(comps*grid[sindex][tindex+1][uindex+1]
	+ s*grid[sindex+1][tindex+1][uindex+1]));
}
