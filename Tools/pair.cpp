// See header for more information about class and functions
#include "pair.h"

Pair::Pair(string wfname, string path, int A, int n1, int l1, int two_ms1, int two_j1, int two_mj1, int two_t1, int n2, int l2, int two_ms2, int two_j2, int two_mj2, int two_t2 )
	: n1( n1), l1( l1), two_ms1( two_ms1), two_j1( two_j1), two_mj1( two_mj1), two_t1( two_t1),
	  n2( n2), l2( l2), two_ms2( two_ms2), two_j2( two_j2), two_mj2( two_mj2), two_t2( two_t2),
	  mosh(n1,l1,n2,l2,path), A( A), coeflistmade( false)
{
	double hbaromega =45.*pow(A, -1./3.) - 25 * pow( A, -2./3.); //MeV
	nu = 938.*hbaromega/2./197./197.;
	deuteronwf = TDeuteron::Wavefunction::CreateWavefunction(wfname);
	if( coeflistmade == false) makecoeflist();
}

Pair::~Pair()
{
	for( int i = 0; i < coeflist.size(); i++ )
		delete coeflist[i];
	delete deuteronwf;
}


//Warning!!! makes no sense anymore!!!

// double Pair::getRelPair( int l ) const
// {
// 	
// 	double result = 0;
// 	for( int i = 0; i < coeflist.size(); i++ )
// 	{
// 		Newcoef* coef = coeflist[i];
// 		if( coef->getl() != l && l >  -1 ) continue;
// 		double value = coef->getCoef();
// 		result += value*value;
// 	}
// 	return result;
// }


double Pair::getcoef( int n, int l, int S, int j, int mj, int T, int MT, int N, int L, int ML )
{
	double result = 0;
	Newcoef* coef = new Newcoef( n1, l1, two_j1, two_mj1, two_t1, n2, l2, two_j2, two_mj2, two_t2, &mosh, N, L, ML, n, l, S, j, mj, T, MT, two_ms1, two_ms2 ); 
	result = coef->getCoef();
	delete coef;
	return result;
}

complex<double> Pair::getwf(bool corr, double r1, double costh1, double phi1, double r2, double costh2, double phi2) const
{

  if(!corr){
//     cout << costh1 << " " << phi1 << " " << costh2 << " " << phi2 << endl;
//     cout << l1 << " " << two_j1 << " " << two_mj1 << " " << two_ms1 << " " << l2 << " " <<two_j2 << " " <<two_mj2 << " " << two_ms2 << endl;
    return uncorrelatedradialwf(n1, l1, r1)*uncorrelatedradialwf(n2, l2, r2)
	*angularwfsplit(l1,1,two_ms1, two_j1, two_mj1, costh1, phi1)*angularwfsplit(l2,1,two_ms2, two_j2, two_mj2, costh2, phi2);
  }
  else{
    double R, costhR, phiR;
    getcomcoord( r1, costh1, phi1, r2, costh2, phi2, &R, &costhR, &phiR );
    double r, costhr, phir;
    getrelcoord( r1, costh1, phi1, r2, costh2, phi2, &r, &costhr, &phir );
    
    // hold sum over real and imag part
    complex<double>  result = 0;
    // start for over coefs in coeflist
    for( int i = 0; i < coeflist.size(); i ++ )
    {
	Newcoef* coef = coeflist[ i ];
	if(coef->getCoef()!=0){
	  // get com wf: radial * angular
	  complex<double> comwf = uncorrelatedradialwf(coef->getN(), coef->getL(), R)
				  * angularwf( coef->getL(), coef->getML(), costhR, phiR);
	  
	  //get rel wf: radial, real ang and imag ang part
	  int jj = coef->getj();
	  for( int la = (jj==0?0:jj-1); la <= jj+1; la++ )
	  {			
		  // sum = coeff * com wave function * rel radial * rel angular
		  //cout << la << " " << coef->getS() << " " << jj << " " << coef->getmj() << endl;
		  result += coef->getCoef() * comwf 
			    * radialwf(coef->getn(), coef->getl(), la, coef->getS(), jj, coef->getT(), r)
			    * angularwf( la, 2*coef->getS(), 2*jj, 2*coef->getmj(),  costhr, phir);			
	  }
			
	}
    }
    // return results
    return result;
  }
}



complex<double> Pair::getcomwf( double r1, double costh1, double phi1, double r2, double costh2, double phi2, int n, int l, int S, int j, int mj, int T) const
{
	double R, costhR, phiR;
	getcomcoord( r1, costh1, phi1, r2, costh2, phi2, &R, &costhR, &phiR );
	//cout << R << " " << costhR << " " << phiR << endl;
	
	// hold sum over real and imag part
	complex<double> result=0.;
	// start for over coefs in coeflist
	for( int i = 0; i < coeflist.size(); i ++ )
	{
		Newcoef* coef = coeflist[ i ];
		if( n != coef->getn() ) continue;
		if( l != coef->getl() ) continue;
		if( S != coef->getS() ) continue;
		if( j != coef->getj() ) continue;
		if( mj != coef->getmj() ) continue;
		if( T != coef->getT() ) continue;

		
		// sum = coeff * radialcom*angularcom
		result += coef->getCoef() * uncorrelatedradialwf(coef->getN(), coef->getL(), R) 
			  * angularwf( coef->getL(), coef->getML(), costhR, phiR);
	// end for loop
	}
	// return results
	return result;
  
}

complex<double> Pair::getrelwf( double r1, double costh1, double phi1, double r2, double costh2, double phi2, int n, int l, int la, int S, int j, int mj, int T) const
{
	double r, costh, phi;
	getrelcoord( r1, costh1, phi1, r2, costh2, phi2, &r, &costh, &phi );

	
	// ! the coef part is in the com wf !

	return radialwf(n, l, la, S, j, T, r) * angularwf( la, 2*S, 2*j, 2*mj, costh, phi);
}

// complex<double> Pair::getrelwf(double r, double costh, double phi, int n, int l, int S, int j, int mj, int T) const
// {
// 	
// 	if(n==0&&l==0&&S==1&&T==0) return deuteronwf->
// 	else return uncorrelatedradialwf(n, l, r) * angularwf( l, 2*S, 2*j, 2*mj, costh, phi)*correlation(r);
// }

//la is l after tensor operator (so l or l+2)
//j is coupling (l,S).
double Pair::radialwf( int n, int l, int la, int S, int j, int T, double r) const
{
	if(n==0&&l==0&&S==1&&T==0) return (la==l? deuteronwf->GetUr(r)/r: deuteronwf->GetWr(r)/r);
	else return (la==l? uncorrelatedradialwf(n,l,r)*correlation(r):0.);

}

complex<double> Pair::angularwf( int l, int ml, double costh, double phi) const
{
	double legendre;
	if( ml < 0 )
	{
		legendre = pow(-1.,-ml) * gsl_sf_legendre_sphPlm(l,-ml,costh);
	}
	else	
	{
		legendre = gsl_sf_legendre_sphPlm(l, ml, costh);
	}
	gsl_complex wf = gsl_complex_polar( legendre, ml*phi );
	return complex<double>(GSL_REAL( wf ), GSL_IMAG(wf));
}

complex<double> Pair::angularwf( int l, int two_S, int two_j, int two_mj, double costh, double phi) const
{
	complex<double> resultsum=0.;
	for( int two_mS = -two_S; two_mS <= two_S; two_mS+=2 )
	{
	    // clebsch * angularwf
	    double clebsch = gsl_sf_coupling_3j( 2*l, two_S, two_j, two_mj-two_mS, two_mS, -two_mj );
	    if(abs(clebsch)>1.E-03){
	      resultsum+= pow(-1., (2*l-two_S+two_mj)/2) * sqrt(two_j+1) * clebsch
			  * angularwf( l, (two_mj-two_mS)/2, costh, phi);			
			   }
	}

	return resultsum;
}

complex<double> Pair::angularwfsplit( int l, int two_S, int two_ms, int two_j, int two_mj, double costh, double phi) const
{
  if(abs((two_mj-two_ms)/2)>l) { /*cout << "blaaa " << (two_mj-two_ms)/2 << " " << two_mj << " " << two_ms << endl ;*/ return 0.;}
  else return  pow(-1., (2*l-two_S+two_mj)/2) * sqrt(two_j+1) 
		*gsl_sf_coupling_3j( 2*l, two_S, two_j, two_mj-two_ms, two_ms, -two_mj )
		    * angularwf( l, (two_mj-two_ms)/2, costh, phi);			
}


void Pair::getcomcoord( double r1, double costh1, double phi1, double r2, double costh2, double phi2, double* R, double* costhR, double* phiR ) const
{
	double sinth1 = sqrt( 1- costh1*costh1);
	double sinth2 = sqrt( 1- costh2*costh2);

	gsl_sf_result sinphi1;
	gsl_sf_result cosphi1;
	gsl_sf_result sinphi2;
	gsl_sf_result cosphi2;
	int status = gsl_sf_sin_e( phi1, &sinphi1);
	if( status ) cerr << "gsl_sf_sin failed, gsl_errno = " << status;
	status = gsl_sf_cos_e( phi1, &cosphi1);
	if( status ) cerr << "gsl_sf_cos failed, gsl_errno = " << status;
	status = gsl_sf_sin_e( phi2, &sinphi2);
	if( status ) cerr << "gsl_sf_sin failed, gsl_errno = " << status;
	status = gsl_sf_cos_e( phi2, &cosphi2);
	if( status ) cerr << "gsl_sf_cos failed, gsl_errno = " << status;

	double cosphi1mphi2 = cosphi1.val* cosphi2.val+ sinphi1.val* sinphi2.val;
	double Rsq = 0.5*r1*r1+ 0.5*r2*r2 + r1*r2* (sinth1*sinth2* cosphi1mphi2 + costh1* costh2 );
	*R = sqrt(Rsq);

	*costhR = (r1*costh1+ r2*costh2)/ sqrt(2.)/ *R;

	double tany= r1* sinth1* sinphi1.val+ r2* sinth2* sinphi2.val;
	double tanx= r1* sinth1* cosphi1.val+ r2* sinth2* cosphi2.val;

	*phiR = atan2( tany, tanx);
}

void Pair::getrelcoord( double r1, double costh1, double phi1, double r2, double costh2, double phi2, double* r, double* costhr, double* phir ) const
{
	double sinth1 = sqrt( 1- costh1*costh1);
	double sinth2 = sqrt( 1- costh2*costh2);

	gsl_sf_result sinphi1;
	gsl_sf_result cosphi1;
	gsl_sf_result sinphi2;
	gsl_sf_result cosphi2;
	int status = gsl_sf_sin_e( phi1, &sinphi1);
	if( status ) cerr << "gsl_sf_sin failed, gsl_errno = " << status;
	status = gsl_sf_cos_e( phi1, &cosphi1);
	if( status ) cerr << "gsl_sf_cos failed, gsl_errno = " << status;
	status = gsl_sf_sin_e( phi2, &sinphi2);
	if( status ) cerr << "gsl_sf_sin failed, gsl_errno = " << status;
	status = gsl_sf_cos_e( phi2, &cosphi2);
	if( status ) cerr << "gsl_sf_cos failed, gsl_errno = " << status;

	double cosphi1mphi2 = cosphi1.val* cosphi2.val+ sinphi1.val* sinphi2.val;
	double rsq = 0.5*r1*r1+ 0.5*r2*r2 - r1*r2* (sinth1*sinth2* cosphi1mphi2 + costh1* costh2 );
	*r = sqrt(rsq);

	*costhr = (r1*costh1- r2*costh2)/ sqrt(2.)/ *r;

	double tany= r1* sinth1* sinphi1.val+ r2* sinth2* sinphi2.val;
	double tanx= r1* sinth1* cosphi1.val+ r2* sinth2* cosphi2.val;

	*phir = atan2( tany, tanx);
}

void Pair::makecoeflist()
{
	if( coeflistmade == true ) return;
	// Next line can be added, but the following should eliminate this result by itself
	//if( n1 == n2 && l1 == l2 && two_j1 == two_j2 && two_mj1 == two_mj2 ) return;
	
	int Smin = 0;
	int Smax = 1;
	int Tmin = 0;
	int Tmax = 1;
	
	int totalEnergy = 2*n1+l1+2*n2+l2;
	double sum = 0;
	// Summation over all quantum numbers
	for( int S = Smin; S <= Smax; S++)
	{
		for( int T = Tmin; T <= Tmax; T++)
		{
			for (int MT = -T; MT <= T; MT++ )
			{
				for( int n = 0; 2*n <= totalEnergy; n++)
				{
					int lmax = totalEnergy - 2*n;
					int lmin = 0;
					for( int l = lmin; l <= lmax ; l++ )
					{
						for( int N = 0; 2*N <= totalEnergy-2*n-l; N++)
						{
							int L;
							L = totalEnergy-2*n-l-2*N;
							for( int ML = -L; ML <= L; ML++)
							{
								for( int j = fabs( l-S ); j <= l+S; j++ )
								{
									for( int mj = -j; mj <= j; mj++ )
									{
										// Create coeff
										Newcoef* coeff = new Newcoef( n1, l1, two_j1, two_mj1, two_t1, n2, l2, two_j2, two_mj2, two_t2, &mosh, N, L, ML, n, l, S, j, mj, T, MT, two_ms1, two_ms2);
										double val = coeff->getCoef();
										if( fabs(val) > 1e-6 )
										{
											coeflist.push_back( coeff );
											sum += val*val;
										}
										else
										{
											delete coeff;
										}
										
										
									}
								}
							}
						}
					}
				}
			}
		}
	}
	//cout << "sum is " << sum << endl;
	coeflistmade = true;
}


double Pair::uncorrelatedradialwf(int n, int l, double r) const
{
	
	gsl_sf_result exp;
	int status = gsl_sf_exp_e(-1.*nu*r*r, &exp);
	if (status)
	{
		if(status == GSL_EUNDRFLW)
		{
			return 0;
		}
		else cerr << "failed, gsl_errno = " << status << endl;
	}

	double radial = pow(r, l) * exp.val * gsl_sf_laguerre_n(n, l+0.5, 2.*nu*r*r);
	double norm = sqrt( 2*gamma( 2*(n+1)) / gamma( 2*n+2*l+3) );
	norm *= pow( 2*nu, l/2.+3./4.);
	
	return norm*radial;
}



double Pair::exponent(double x) const
{
	gsl_sf_result exp;
	int status = gsl_sf_exp_e( x, &exp);
	if (status)
	{
		if(status == GSL_EUNDRFLW)
		{
			return 0;
		}
		else cerr << "gsl_sf_exp failed, gsl_errno = " << status;
		cerr << ": " << x << endl;
	}
	return exp.val;
}

double Pair::log(double x) const
{
	gsl_sf_result log;
	int status = gsl_sf_log_e( x, &log);
	if (status)
	{
		cerr << "gsl_sf_exp failed, gsl_errno = " << status;
		cerr << ": " << x << endl;
	}
	return log.val;
}

double Pair::gamma( int two_n) const
{
	if( two_n == 2 ) return 1;
	if( two_n == 1 ) return sqrt(M_PI);
	if( two_n == 0 ) cerr << "two_n == 0 !" << endl;
	return (two_n-2)/2. * gamma( two_n-2);
}

/////////////////
//G-D correlation function
///////////////////////
  
double Pair::correlation(double r) const{

  double a[11]={0.99994,0.99600E-01,-2.9385,28.513,-132.53,238.33,-151.54,-69.845,159.94,-88.043,16.930};
  if(r>1.24) return 1+0.1204071*exp(-2.5002*(r-1.24)*(r-1.24));
  else{
    double result=0.;
    for(int i=0;i<11;i++) result += a[i]*pow(r,i);
    return 1-result;
    }
  //return 1-exp(-r*r*1.18147448);
}
