#include <cassert>

using namespace std;

#include "constants.hpp"
#include "OldDeuteron.hpp"

#define NORMFACTOR 11.4081948 //1/HBARC^{3/2}

OldDeuteron::OldDeuteron(const std::string &wf_name){
    if(!wf_name.compare("Paris")){
        mass=vector<double>(13,0.);
        S_pre=vector<double>(13,0.);
        D_pre=vector<double>(13,0.);
        S_pre[0]=0.88688076;                                                   
        S_pre[1]=-0.34717093;                                                  
        S_pre[2]=-3.050238;                                                   
        S_pre[3]=56.207766;                                                    
        S_pre[4]=-749.57334;                                                   
        S_pre[5]=5336.5279;                                                    
        S_pre[6]=-22706.863;                                                   
        S_pre[7]=60434.4690;                                                  
        S_pre[8]=-102920.58 ; 
        S_pre[9]=112233.57;                                                
        S_pre[10]=-75925.226;                                       
        S_pre[11]=29059.715; 
        S_pre[12]=0.;
        for(int i=0;i<12;i++) S_pre[12]-=S_pre[i];

        D_pre[0]=0.023135193;
        D_pre[1]=-0.85604572;
        D_pre[2]=5.6068193  ;
        D_pre[3]=-69.462922 ;
        D_pre[4]=416.31118  ;
        D_pre[5]=-1254.6621 ;
        D_pre[6]=1238.783   ;
        D_pre[7]=3373.9172  ;
        D_pre[8]=-13041.151 ;
        D_pre[9]=19512.524 ;
        
        for(int i=0;i<13;i++) mass[i]=0.23162461+i;
        
        double a=0.,b=0.,cc=0.;
        for(int i=0;i<10;i++){
            a+=D_pre[i]/(mass[i]*mass[i]);
            b+=D_pre[i];
            cc+=D_pre[i]*(mass[i]*mass[i]);
        }
        D_pre[10]=mass[10]*mass[10]/(mass[12]*mass[12]-mass[10]*mass[10])/(mass[11]*mass[11]-mass[10]*mass[10])
            *(-mass[11]*mass[11]*mass[12]*mass[12]*a+(mass[11]*mass[11]+mass[12]*mass[12])*b-cc);  
        D_pre[11]=mass[11]*mass[11]/(mass[10]*mass[10]-mass[11]*mass[11])/(mass[12]*mass[12]-mass[11]*mass[11])
            *(-mass[12]*mass[12]*mass[10]*mass[10]*a+(mass[12]*mass[12]+mass[10]*mass[10])*b-cc);  
        D_pre[12]=mass[12]*mass[12]/(mass[11]*mass[11]-mass[12]*mass[12])/(mass[10]*mass[10]-mass[12]*mass[12])
            *(-mass[10]*mass[10]*mass[11]*mass[11]*a+(mass[10]*mass[10]+mass[11]*mass[11])*b-cc);  
        norm=0.79788456/sqrt(4.*PI)*pow(INVHBARC,3./2.);
    }
    else if(!wf_name.compare("AV18")){
        mass=vector<double>(12,0.);
        S_pre=vector<double>(12,0.);
        D_pre=vector<double>(12,0.);
        S_pre[0]  =  0.706699E+00;                                                
        S_pre[1]  = -0.169743E+00;                                             
        S_pre[2]  =  0.112368E+01;                                             
        S_pre[3]  = -0.852995E+01;                                             
        S_pre[4]  =  0.195033E+02;                                             
        S_pre[5]  = -0.757831E+02;                                             
        S_pre[6]  =  0.283739E+03;                                             
        S_pre[7]  = -0.694734E+03;                                             
        S_pre[8]  =  0.885257E+03;                                             
        S_pre[9] = -0.720739E+03 ;                                            
        S_pre[10] =  0.412969E+03;
        S_pre[11] = -0.103336E+03;                                             

        D_pre[0]  =  0.176655E-01;                                               
        D_pre[1]  = -0.124551E+00;                                            
        D_pre[2]  = -0.108815E+01;                                            
        D_pre[3]  =  0.384848E+01;                                            
        D_pre[4]  = -0.852442E+01;                                            
        D_pre[5]  =  0.209435E+02;                                            
        D_pre[6]  = -0.490728E+02;                                            
        D_pre[7]  =  0.577382E+02;                                            
        D_pre[8]  = -0.127114E+01;                                            
        D_pre[9] = -0.628361E+02;
        D_pre[10] =  0.581016E+02;
        D_pre[11] = -0.177062E+02;                                            

        mass[0]  = 0.2316; 
        mass[1]  = 1.0;
        mass[2]  = 1.5;
        mass[3]  = 2.0;
        mass[4]  = 2.5;
        mass[5]  = 3.5;
        mass[6]  = 4.5;
        mass[7]  = 5.5;
        mass[8]  = 6.5;
        mass[9] = 8.0;
        mass[10] = 9.5;
        mass[11] = 11.0;
        norm=1./sqrt(4.*PI)*pow(INVHBARC,3./2.);
        
    }
    else if(!wf_name.compare("CDBonn")){
        mass=vector<double>(11,0.);
        S_pre=vector<double>(11,0.);
        D_pre=vector<double>(11,0.);
        S_pre[0] =   0.88472985;
        S_pre[1] = - 0.26408759;
        S_pre[2] = - 0.44114404e-01;
        S_pre[3] = - 0.14397512e+02;
        S_pre[4] =   0.85591256e+02;
        S_pre[5] = - 0.31876761e+03;
        S_pre[6] =   0.70336701e+03;
        S_pre[7] = - 0.90049586e+03;
        S_pre[8] =   0.66145441e+03;
        S_pre[9]= - 0.25958894e+03;
        S_pre[10]=0.;
        for(int i=0;i<10;i++) S_pre[10]-=S_pre[i];
        
        D_pre[0] =   0.22623762e-01;
        D_pre[1] = - 0.50471056e+00;
        D_pre[2] =   0.56278897e+00;
        D_pre[3] = - 0.16079764e+02;
        D_pre[4] =   0.11126803e+03;
        D_pre[5] = - 0.44667490e+03;
        D_pre[6] =   0.10985907e+04;
        D_pre[7] = - 0.16114995e+04;

        for(int i=0;i<11;i++) mass[i]=0.2315380 + 0.9*i;
        
        double tm0=0.,tm1=0.,tm2=0.;
        for(int i=0;i<8;i++){
        tm0+=D_pre[i];
        tm1+=D_pre[i]/(mass[i]*mass[i]);
        tm2+=D_pre[i]*(mass[i]*mass[i]);
        }
        D_pre[8] = mass[8]*mass[8]/(mass[10]*mass[10]-mass[8]*mass[8])/(mass[9]*mass[9]-mass[8]*mass[8]) *
                ( -mass[9]*mass[9]*mass[10]*mass[10]*tm1 + (mass[9]*mass[9]+mass[10]*mass[10])*tm0-tm2);
        D_pre[9] = mass[9]*mass[9]/(mass[10]*mass[10]-mass[9]*mass[9])/(mass[8]*mass[8]-mass[9]*mass[9]) *
                ( -mass[10]*mass[10]*mass[8]*mass[8]*tm1 + (mass[10]*mass[10]+mass[8]*mass[8])*tm0-tm2);
        D_pre[10] = mass[10]*mass[10]/(mass[9]*mass[9]-mass[10]*mass[10])/(mass[8]*mass[8]-mass[10]*mass[10]) *
                ( -mass[8]*mass[8]*mass[9]*mass[9]*tm1 + (mass[8]*mass[8]+mass[9]*mass[9])*tm0-tm2);
        
        norm=1./(sqrt(2.)*PI)*pow(INVHBARC,3./2.);
    }

    else if(!wf_name.compare("AV18b")){
        mass=vector<double>(12,0.);
        S_pre=vector<double>(12,0.);
        D_pre=vector<double>(12,0.);
        mass[0] = 0.232500e+00;
        mass[1] = 0.500000e+00;
        mass[2] = 0.800000e+00;
        mass[3] = 0.120000e+01;
        mass[4] = 0.160000e+01;
        mass[5] = 0.200000e+01;
        mass[6] = 0.400000e+01;
        mass[7] = 0.600000e+01;
        mass[8] = 0.100000e+02;
        mass[9] = 0.140000e+02;
        mass[10] = 0.180000e+02;
        mass[11] = 0.220000e+02;

        S_pre[0] =  0.105252223e+02;
        S_pre[1] =  0.124352529e+02;
        S_pre[2] = -0.687541641e+02;
        S_pre[3] =  0.239111042e+03;
        S_pre[4] = -0.441014422e+03;
        S_pre[5] =  0.300140328e+03;
        S_pre[6] = -0.230639939e+03;
        S_pre[7] =  0.409671540e+03;
        S_pre[8] = -0.733453611e+03;
        S_pre[9]=  0.123506081e+04;
        S_pre[10]= -0.120520606e+04;
        S_pre[11]=0.;
        for(int i=0;i<11;i++) S_pre[11]-=S_pre[i];

        D_pre[0] =  0.280995496e+00;
        D_pre[1] =  0.334117629e-01;
        D_pre[2] = -0.727192237e+00;
        D_pre[3] = -0.302809607e+01;
        D_pre[4] = -0.903824982e+01;
        D_pre[5] =  0.496045967e+01;
        D_pre[6] = -0.271985613e+02;
        D_pre[7] =  0.125334598e+03;
        D_pre[8] = -0.346742235e+03;
        
        double sp2=0.,sm=0.,sm2=0.;
        for(int i=0;i<9;i++){
        sp2+=D_pre[i]/(mass[i]*mass[i]);
        sm+=D_pre[i];
        sm2+=D_pre[i]*(mass[i]*mass[i]);
        }
        double a,b,cc;
        a = mass[11]*mass[11];
        b = mass[10]*mass[10];
        cc = mass[9]*mass[9];
        D_pre[9] = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2);
    
        a = mass[9]*mass[9];
        b = mass[11]*mass[11];
        cc = mass[10]*mass[10];
        D_pre[10] = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2);
    
        a = mass[10]*mass[10];
        b = mass[9]*mass[9];
        cc = mass[11]*mass[11];
        D_pre[11] = (cc/((a-cc)*(b-cc)))*(-b*a*sp2+(b+a)*sm-sm2);
        
        double fact=4.*PI*sqrt(PI/2.);
        for(int i=0;i<12;i++){
        S_pre[i]/=fact;
        D_pre[i]/=fact;
        }
        norm=1./sqrt(4.*PI)*pow(INVHBARC,3./2.);
    }
    else{
        cerr << " invalid deuteron wf name " << endl;
        assert(1==0);
    }
    

}



double OldDeuteron::Ufront(int M, int spinp, int spinn){
  double front=0.;
  if(M==2&&spinp==1&&spinn==1) front=1.;
  if(M==-2&&spinp==-1&&spinn==-1) front=1.;
  if(M==0&&spinp==1&&spinn==-1) front=1./sqrt(2.);
  if(M==0&&spinp==-1&&spinn==1) front=1./sqrt(2.);
  
  return front;
}

std::complex<double> OldDeuteron::get_stensor(int dspin, int spinp, int spinn, double theta, double phi){

   complex<double> stensor;
   double costheta,sintheta;
   sincos(theta,&sintheta,&costheta);
   double cosphi,sinphi;
   sincos(phi,&sinphi,&cosphi);
   double cos2phi,sin2phi;
   sincos(2.*phi,&sin2phi,&cos2phi);
   
   switch(dspin){
    case -2:
      switch(spinp){
	  case -1:
	  switch(spinn){
	    case -1:
	      stensor=complex<double>(3.0*costheta*costheta-1.0,0.);
	      break;
	    case 1:
	      stensor=complex<double>(-3.0*costheta*sintheta*cosphi,3.0*costheta*sintheta*sinphi);
	      break;
	    default:
	      cout << "invalid neutron spin" << endl;
	      exit(1);
	  }
	  break;
	case 1:
	  switch(spinn){
	    case -1:
	      stensor=complex<double>(-3.0*costheta*sintheta*cosphi,3.0*costheta*sintheta*sinphi);
	      break;
	    case 1:
	      stensor=complex<double>(3.0*sintheta*sintheta*cos2phi,-3.0*sintheta*sintheta*sin2phi);
	      break;
	    default:
	      cout << "invalid neutron spin" << endl;
	      exit(1);
	  }
	  break;
	default:
	  cout << "invalid proton spin" << endl;
	  exit(1);
      }
      break;
    case 0:
      switch(spinp){
	  case -1:
	  switch(spinn){
	    case -1:
	      stensor=complex<double>(-3.0*costheta*sintheta*cosphi,-3.0*costheta*sintheta*sinphi);
	      break;
	    case 1:
	      stensor=complex<double>(-(3.0*costheta*costheta-1.0),0.);
	      break;
	    default:
	      cout << "invalid neutron spin" << endl;
	      exit(1);
	  }
	  break;
	  case 1:
	  switch(spinn){
	    case -1:
	      stensor=complex<double>(-(3.0*costheta*costheta-1.0),0.);
	      break;
	    case 1:
	      stensor=complex<double>(3.0*costheta*sintheta*cosphi,-3.0*costheta*sintheta*sinphi);
	      break;
	    default:
	      cout << "invalid neutron spin" << endl;
	      exit(1);
	  }
	  break;
	default:
	  cout << "invalid proton spin" << endl;
	  exit(1);
      }
      stensor*=sqrt(2.);
      break;
    case 2:
      switch(spinp){
	  case -1:
	  switch(spinn){
	    case -1:
	      stensor=complex<double>(3.0*sintheta*sintheta*cos2phi,3.0*sintheta*sintheta*sin2phi);
	      break;
	    case 1:
	      stensor=complex<double>(3.0*costheta*sintheta*cosphi,3.0*costheta*sintheta*sinphi);
	      break;
	    default:
	      cout << "invalid neutron spin" << endl;
	      exit(1);
	  }
	  break;
	  case 1:
	  switch(spinn){
	    case -1:
	      stensor=complex<double>(3.0*costheta*sintheta*cosphi,3.0*costheta*sintheta*sinphi);
	      break;
	    case 1:
	      stensor=complex<double>(3.0*costheta*costheta-1.0,0.);
	      break;
	    default:
	      cout << "invalid neutron spin" << endl;
	      exit(1);
	  }
	  break;
	default:
	  cout << "invalid proton spin" << endl;
	  exit(1);
      }
      break;
    default:
      cout << "invalid deuteron spin" << endl;
      exit(1);
  } 
  return stensor;
 

}


double OldDeuteron::U(double p){
    double result=0.;
    double x=p*INVHBARC;
    for(int i=0;i<S_pre.size();i++) result+=S_pre[i]/(x*x+mass[i]*mass[i]);
    
    return result*norm;


}


double OldDeuteron::W(double p){
    double result=0.;
    double x=p*INVHBARC;
    for(int i=0;i<D_pre.size();i++) result+=D_pre[i]/(x*x+mass[i]*mass[i]);
    
    return result*norm;


}

std::complex< double> OldDeuteron::deuteronwf(int dspin, int proton, int spinp, int spinn, double p, double theta, double phi){
    complex<double> result = Ufront(dspin,spinp,spinn)*U(p) + W(p)*get_stensor(dspin,spinp,spinn,theta,phi)/sqrt(8.);
    return (proton? result:-result);
}