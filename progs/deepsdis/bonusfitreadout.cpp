//write out some of the fit parameters



#include <cstdlib>
#include <cmath>
#include <complex>
#include <cstdarg>
#include <time.h>

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <signal.h>
#include <iomanip>
// using namespace std;

using std::cout;
using std::endl;

#include <constants.hpp>
#include <DeuteronCross.hpp>
#include "bonusfits.h"
#include <NuclStructure.hpp>
#include "bonusdata.h"
#include "bonusfits2.h"

int main(int argc, char *argv[])
{
  
//   for(int i=0;i<3;i++){
//     double Q2=(data::Q2[i]+data::Q2[i+1])/2.;
//     for(int j=1;j<60;j++){
//       double Wsq=2.8E03*2.8E03-j*11.E04;
//       NuclStructure s1_n(0,Q2,Wsq,1,"CB");
//       NuclStructure s2_n(0,Q2,Wsq,1,"SLAC");
//       NuclStructure s1_p(1,Q2,Wsq,1,"CB");
//       NuclStructure s2_p(1,Q2,Wsq,1,"SLAC");
//       cout << Q2 << " " << sqrt(Wsq) << " " << s1_n.getF2() << " " << s2_n.getF2() 
//       << " " << s1_p.getF2() << " " << s2_p.getF2() << endl;
//       
//     }
//     cout << endl << endl;
//   }
//   exit(1);
  
  //norms
  cout <<"#rows: pm setting" << endl 
    << "#column1: Q2=3.6,W=2,Ein=4" << endl
    << "#column2: Q2=3.6,W=2,Ein=5" << endl
    << "#column3: Q2=3.6,W=2.4,Ein=5" << endl
    << "#column4: Q2=1.6,W=2,Ein=4" << endl
    << "#column5: Q2=1.6,W=2.4,Ein=4" << endl
    << "#column6: Q2=1.6,W=2,Ein=5" << endl
    << "#column7: Q2=1.6,W=2.4,Ein=5" << endl<<endl;
    
  cout << "#paris CB" << endl;
  for(int i=0;i<4;i++){ //spectator momentum index
    cout << i << " ";
    cout << bonusfits2::normfits_paris_CB_off3_lc0_fsi_q2dep1_beam4[0][3][i] << " "; //Q2=3.6,W=2,Ein=4
    cout << bonusfits2::normfits_paris_CB_off3_lc0_fsi_q2dep1_beam5[0][3][i] << " "; //Q2=3.6,W=2,Ein=5
    cout << bonusfits2::normfits_paris_CB_off3_lc0_fsi_q2dep1_beam5[0][4][i] << " "; //Q2=3.6,W=2.4,Ein=5
    cout << bonusfits2::normfits_paris_CB_off3_lc0_fsi_q2dep0_beam4[1][3][i] << " ";//Q2=1.6,W=2,Ein=4
    cout << bonusfits2::normfits_paris_CB_off3_lc0_fsi_q2dep0_beam4[1][4][i] << " ";//Q2=1.6,W=2.4,Ein=4
    cout << bonusfits2::normfits_paris_CB_off3_lc0_fsi_q2dep0_beam5[0][3][i] << " ";//Q2=1.6,W=2,Ein=5
    cout << bonusfits2::normfits_paris_CB_off3_lc0_fsi_q2dep0_beam5[0][4][i] << endl;//Q2=1.6,W=2,Ein=5	
  }
  cout << endl << endl << endl;

  cout << "#AV18 CB" << endl;
  for(int i=0;i<4;i++){ //spectator momentum index
    cout << i << " ";
    cout << bonusfits2::normfits_AV18_CB_off3_lc0_fsi_q2dep1_beam4[0][3][i] << " ";
    cout << bonusfits2::normfits_AV18_CB_off3_lc0_fsi_q2dep1_beam5[0][3][i] << " ";
    cout << bonusfits2::normfits_AV18_CB_off3_lc0_fsi_q2dep1_beam5[0][4][i] << " ";
    cout << bonusfits2::normfits_AV18_CB_off3_lc0_fsi_q2dep0_beam4[1][3][i] << " ";
    cout << bonusfits2::normfits_AV18_CB_off3_lc0_fsi_q2dep0_beam4[1][4][i] << " ";
    cout << bonusfits2::normfits_AV18_CB_off3_lc0_fsi_q2dep0_beam5[0][3][i] << " ";
    cout << bonusfits2::normfits_AV18_CB_off3_lc0_fsi_q2dep0_beam5[0][4][i] << endl;	
  }
  cout << endl << endl << endl;
  
  cout << "#CDBonn CB" << endl;
  for(int i=0;i<4;i++){ //spectator momentum index
    cout << i << " ";
    cout << bonusfits2::normfits_CDBonn_CB_off3_lc0_fsi_q2dep1_beam4[0][3][i] << " ";
    cout << bonusfits2::normfits_CDBonn_CB_off3_lc0_fsi_q2dep1_beam5[0][3][i] << " ";
    cout << bonusfits2::normfits_CDBonn_CB_off3_lc0_fsi_q2dep1_beam5[0][4][i] << " ";
    cout << bonusfits2::normfits_CDBonn_CB_off3_lc0_fsi_q2dep0_beam4[1][3][i] << " ";
    cout << bonusfits2::normfits_CDBonn_CB_off3_lc0_fsi_q2dep0_beam4[1][4][i] << " ";
    cout << bonusfits2::normfits_CDBonn_CB_off3_lc0_fsi_q2dep0_beam5[0][3][i] << " ";
    cout << bonusfits2::normfits_CDBonn_CB_off3_lc0_fsi_q2dep0_beam5[0][4][i] << endl;	
  }
  cout << endl << endl << endl;

  cout << "#GrossWJC1 CB" << endl;    
  for(int i=0;i<4;i++){ //spectator momentum index
    cout << i << " ";
    cout << bonusfits2::normfits_GrossWJC1_CB_off3_lc0_fsi_q2dep1_beam4[0][3][i] << " ";
    cout << bonusfits2::normfits_GrossWJC1_CB_off3_lc0_fsi_q2dep1_beam5[0][3][i] << " ";
    cout << bonusfits2::normfits_GrossWJC1_CB_off3_lc0_fsi_q2dep1_beam5[0][4][i] << " ";
    cout << bonusfits2::normfits_GrossWJC1_CB_off3_lc0_fsi_q2dep0_beam4[1][3][i] << " ";
    cout << bonusfits2::normfits_GrossWJC1_CB_off3_lc0_fsi_q2dep0_beam4[1][4][i] << " ";
    cout << bonusfits2::normfits_GrossWJC1_CB_off3_lc0_fsi_q2dep0_beam5[0][3][i] << " ";
    cout << bonusfits2::normfits_GrossWJC1_CB_off3_lc0_fsi_q2dep0_beam5[0][4][i] << endl;	
  }
  cout << endl << endl << endl;

  cout << "#paris SLAC" << endl;
  for(int i=0;i<4;i++){ //spectator momentum index
    cout << i << " ";
    cout << bonusfits2::normfits_paris_SLAC_off3_lc0_fsi_q2dep1_beam4[0][3][i] << " ";
    cout << bonusfits2::normfits_paris_SLAC_off3_lc0_fsi_q2dep1_beam5[0][3][i] << " ";
    cout << bonusfits2::normfits_paris_SLAC_off3_lc0_fsi_q2dep1_beam5[0][4][i] << " ";
    cout << bonusfits2::normfits_paris_SLAC_off3_lc0_fsi_q2dep0_beam4[1][3][i] << " ";
    cout << bonusfits2::normfits_paris_SLAC_off3_lc0_fsi_q2dep0_beam4[1][4][i] << " ";
    cout << bonusfits2::normfits_paris_SLAC_off3_lc0_fsi_q2dep0_beam5[0][3][i] << " ";
    cout << bonusfits2::normfits_paris_SLAC_off3_lc0_fsi_q2dep0_beam5[0][4][i] << endl;	
  }
  cout << endl << endl << endl;

  //errors
  cout <<"#rows: pm setting" << endl 
    << "#column1: Q2=3.6,W=2,Ein=4" << endl
    << "#column2: Q2=3.6,W=2,Ein=5" << endl
    << "#column3: Q2=3.6,W=2.4,Ein=5" << endl
    << "#column4: Q2=1.6,W=2,Ein=4" << endl
    << "#column5: Q2=1.6,W=2.4,Ein=4" << endl
    << "#column6: Q2=1.6,W=2,Ein=5" << endl
    << "#column7: Q2=1.6,W=2.4,Ein=5" << endl<<endl;
    
  cout << "#paris CB" << endl;
  for(int i=0;i<4;i++){ //spectator momentum index
    cout << i << " ";
    cout << bonusfits2::normfiterrors_paris_CB_off3_lc0_fsi_q2dep1_beam4[0][3][i] << " "; //Q2=3.6,W=2,Ein=4
    cout << bonusfits2::normfiterrors_paris_CB_off3_lc0_fsi_q2dep1_beam5[0][3][i] << " "; //Q2=3.6,W=2,Ein=5
    cout << bonusfits2::normfiterrors_paris_CB_off3_lc0_fsi_q2dep1_beam5[0][4][i] << " "; //Q2=3.6,W=2.4,Ein=5
    cout << bonusfits2::normfiterrors_paris_CB_off3_lc0_fsi_q2dep0_beam4[1][3][i] << " ";//Q2=1.6,W=2,Ein=4
    cout << bonusfits2::normfiterrors_paris_CB_off3_lc0_fsi_q2dep0_beam4[1][4][i] << " ";//Q2=1.6,W=2.4,Ein=4
    cout << bonusfits2::normfiterrors_paris_CB_off3_lc0_fsi_q2dep0_beam5[0][3][i] << " ";//Q2=1.6,W=2,Ein=5
    cout << bonusfits2::normfiterrors_paris_CB_off3_lc0_fsi_q2dep0_beam5[0][4][i] << endl;//Q2=1.6,W=2,Ein=5	
  }
  cout << endl << endl << endl;

  cout << "#AV18 CB" << endl;
  for(int i=0;i<4;i++){ //spectator momentum index
    cout << i << " ";
    cout << bonusfits2::normfiterrors_AV18_CB_off3_lc0_fsi_q2dep1_beam4[0][3][i] << " ";
    cout << bonusfits2::normfiterrors_AV18_CB_off3_lc0_fsi_q2dep1_beam5[0][3][i] << " ";
    cout << bonusfits2::normfiterrors_AV18_CB_off3_lc0_fsi_q2dep1_beam5[0][4][i] << " ";
    cout << bonusfits2::normfiterrors_AV18_CB_off3_lc0_fsi_q2dep0_beam4[1][3][i] << " ";
    cout << bonusfits2::normfiterrors_AV18_CB_off3_lc0_fsi_q2dep0_beam4[1][4][i] << " ";
    cout << bonusfits2::normfiterrors_AV18_CB_off3_lc0_fsi_q2dep0_beam5[0][3][i] << " ";
    cout << bonusfits2::normfiterrors_AV18_CB_off3_lc0_fsi_q2dep0_beam5[0][4][i] << endl;	
  }
  cout << endl << endl << endl;
  
  cout << "#CDBonn CB" << endl;
  for(int i=0;i<4;i++){ //spectator momentum index
    cout << i << " ";
    cout << bonusfits2::normfiterrors_CDBonn_CB_off3_lc0_fsi_q2dep1_beam4[0][3][i] << " ";
    cout << bonusfits2::normfiterrors_CDBonn_CB_off3_lc0_fsi_q2dep1_beam5[0][3][i] << " ";
    cout << bonusfits2::normfiterrors_CDBonn_CB_off3_lc0_fsi_q2dep1_beam5[0][4][i] << " ";
    cout << bonusfits2::normfiterrors_CDBonn_CB_off3_lc0_fsi_q2dep0_beam4[1][3][i] << " ";
    cout << bonusfits2::normfiterrors_CDBonn_CB_off3_lc0_fsi_q2dep0_beam4[1][4][i] << " ";
    cout << bonusfits2::normfiterrors_CDBonn_CB_off3_lc0_fsi_q2dep0_beam5[0][3][i] << " ";
    cout << bonusfits2::normfiterrors_CDBonn_CB_off3_lc0_fsi_q2dep0_beam5[0][4][i] << endl;	
  }
  cout << endl << endl << endl;

  cout << "#GrossWJC1 CB" << endl;    
  for(int i=0;i<4;i++){ //spectator momentum index
    cout << i << " ";
    cout << bonusfits2::normfiterrors_GrossWJC1_CB_off3_lc0_fsi_q2dep1_beam4[0][3][i] << " ";
    cout << bonusfits2::normfiterrors_GrossWJC1_CB_off3_lc0_fsi_q2dep1_beam5[0][3][i] << " ";
    cout << bonusfits2::normfiterrors_GrossWJC1_CB_off3_lc0_fsi_q2dep1_beam5[0][4][i] << " ";
    cout << bonusfits2::normfiterrors_GrossWJC1_CB_off3_lc0_fsi_q2dep0_beam4[1][3][i] << " ";
    cout << bonusfits2::normfiterrors_GrossWJC1_CB_off3_lc0_fsi_q2dep0_beam4[1][4][i] << " ";
    cout << bonusfits2::normfiterrors_GrossWJC1_CB_off3_lc0_fsi_q2dep0_beam5[0][3][i] << " ";
    cout << bonusfits2::normfiterrors_GrossWJC1_CB_off3_lc0_fsi_q2dep0_beam5[0][4][i] << endl;	
  }
  cout << endl << endl << endl;

  cout << "#paris SLAC" << endl;
  for(int i=0;i<4;i++){ //spectator momentum index
    cout << i << " ";
    cout << bonusfits2::normfiterrors_paris_SLAC_off3_lc0_fsi_q2dep1_beam4[0][3][i] << " ";
    cout << bonusfits2::normfiterrors_paris_SLAC_off3_lc0_fsi_q2dep1_beam5[0][3][i] << " ";
    cout << bonusfits2::normfiterrors_paris_SLAC_off3_lc0_fsi_q2dep1_beam5[0][4][i] << " ";
    cout << bonusfits2::normfiterrors_paris_SLAC_off3_lc0_fsi_q2dep0_beam4[1][3][i] << " ";
    cout << bonusfits2::normfiterrors_paris_SLAC_off3_lc0_fsi_q2dep0_beam4[1][4][i] << " ";
    cout << bonusfits2::normfiterrors_paris_SLAC_off3_lc0_fsi_q2dep0_beam5[0][3][i] << " ";
    cout << bonusfits2::normfiterrors_paris_SLAC_off3_lc0_fsi_q2dep0_beam5[0][4][i] << endl;	
  }
  cout << endl << endl << endl;


  
}
