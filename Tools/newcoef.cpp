// Documentation in header file.
#include "newcoef.h"

  
Newcoef::Newcoef( int n1, int l1, int two_j1, int two_mj1, int two_t1, int n2, int l2, int two_j2, int two_mj2, int two_t2, RecMosh* mosh, 
			int N, int L, int ML, int n, int l, int S, int j, int mj, int T, int MT, int two_ms1, int two_ms2)
	: N(N),
	  L(L),
	  ML(ML),
	  n(n),
	  l(l),
	  S( S ),
	  j( j),
	  mj( mj ),
	  T( T ),
	  MT( MT )
{
	
	double value = 0;
	// Term due to anti-symmetrisation (1 - (-1)^(S+l+T) )
	if( (S+l+T)%2 == 0 )
	{
		coeff = 0;
		return;
	}
	else
	{
	  int two_ml1= two_mj1- two_ms1;
	  int two_ml2= two_mj2- two_ms2;
	  double c1= gsl_sf_coupling_3j( 2*l1, 1, two_j1, two_ml1, two_ms1, -two_mj1);
	  double c2= gsl_sf_coupling_3j( 2*l2, 1, two_j2, two_ml2, two_ms2, -two_mj2);
	  double factor= (two_j1+1)* (two_j2+1)* c1*c1* c2*c2;
	  if(abs(factor)<1.E-04) {
	    coeff=0.;
	    return;
	  }
	  else{
	    for( int two_J = fabs(two_j1-two_j2); two_J <= two_j1+two_j2; two_J +=2)
	    {
	      int J = two_J/2;
	      for( int MJ = -J; MJ <= J; MJ++)
	      {
		      double threej1 = pow( -1., ( two_j1- two_j2+ 2*MJ)/2 ) * sqrt( two_J+1) * gsl_sf_coupling_3j( two_j1, two_j2, two_J, two_mj1, two_mj2, 2*MJ ); 
		      double threej2 = pow( -1., j-L+MJ) * sqrt( 2*J+1) * gsl_sf_coupling_3j( 2*j, 2*L, 2*J, 2*mj, 2*ML, -2*MJ);
		      double isospin = sqrt( 2*T + 1.) * pow(-1., MT) * gsl_sf_coupling_3j( 1, 1, 2*T, two_t1, two_t2, -2*MT);
		      for( int Lambda = fabs(l1-l2); Lambda <= l1+l2; Lambda++)
		      {
			      double moshme = mosh->getCoefficient(n, l, N, L, Lambda);
			      if( moshme == 0 )
			      {
				      continue;
			      }
			      double ninej = sqrt(two_j1+1) * sqrt(two_j2+1) * sqrt(2*S+1) * sqrt(2*Lambda+1) 
			      * gsl_sf_coupling_9j(2*l1, 1, two_j1, 2*l2, 1, two_j2, 2*Lambda, 2*S, 2*J);
			      double sixj = sqrt(2*Lambda+1) * sqrt(2*j+1) * pow(-1., j+Lambda+S+L) * gsl_sf_coupling_6j(2*j, 2*L, 2*J, 2*Lambda, 2*S, 2*l);
			      value += factor* ninej* sixj*threej1* threej2* isospin* 2*moshme;
		      }
	      }
	    }
	  }
	  // Normalisation: See notes about two and three body transformation, antisymmetrisation and normalisation.
	  double norm = sqrt(2.);
  //  	if( n1 == n2 && l1 == l2 && two_j1 == two_j2 && two_t1 == two_t2 )
  //  		norm = 2.;
	  if( value*value < 0.0001)
		  value = 0;
	  
	  coeff=value/norm;		
	}
		
}
