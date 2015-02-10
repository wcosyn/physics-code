//prog to fit the normalization of the bonus data
//sigma parameters used from deeps fit
//take scattering parameters fixed, fit norms (other strategies do not yield significant better results_)



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

int main(int argc, char *argv[])
{
  for(int i=0;i<4;i++){
    cout << i << " ";
    for(int j=0;j<3;j++){
      for(int k=0;k<5;k++){
	cout << bonusfits::normfits_off3_paris_CB_off3_lc1_fsi_q2dep0_beam4[j][k][i] << " ";
	if(j>0) cout << bonusfits::normfits_off3_paris_CB_off3_lc1_fsi_q2dep0_beam5[j-1][k][i] << " ";
      }
    }
    cout << endl;
  }
  cout << endl << endl << endl;

  for(int i=0;i<4;i++){
    cout << i << " ";
    for(int j=0;j<3;j++){
      for(int k=0;k<5;k++){
	cout << bonusfits::normfits_off3_paris_CB_off3_lc1_fsi_q2dep1_beam4[j][k][i] << " ";
	if(j>0) cout << bonusfits::normfits_off3_paris_CB_off3_lc1_fsi_q2dep1_beam5[j-1][k][i] << " ";
      }
    }
    cout << endl;
  }
  cout << endl << endl << endl;

  for(int i=0;i<4;i++){
    cout << i << " ";
    for(int j=0;j<3;j++){
      for(int k=0;k<5;k++){
	cout << bonusfits::normfits_off3_paris_CB_off3_lc0_pw_beam4[j][k][i] << " ";
	if(j>0) cout << bonusfits::normfits_off3_paris_CB_off3_lc0_pw_beam5[j-1][k][i] << " ";
      }
    }
    cout << endl;
  }
  cout << endl << endl << endl;
  
  for(int i=0;i<4;i++){
    cout << i << " ";
    for(int j=0;j<3;j++){
      for(int k=0;k<5;k++){
	cout << bonusfits::normfits_off3_AV18_CB_off3_lc0_fsi_q2dep1_beam4[j][k][i] << " ";
	if(j>0) cout << bonusfits::normfits_off3_AV18_CB_off3_lc0_fsi_q2dep1_beam5[j-1][k][i] << " ";
      }
    }
    cout << endl;
  }
  cout << endl << endl << endl;


  for(int i=0;i<4;i++){
    cout << i << " ";
    for(int j=0;j<3;j++){
      for(int k=0;k<5;k++){
	cout << bonusfits::normfits_off3_paris_SLAC_off3_lc0_fsi_q2dep1_beam4[j][k][i] << " ";
	if(j>0) cout << bonusfits::normfits_off3_paris_SLAC_off3_lc0_fsi_q2dep1_beam5[j-1][k][i] << " ";
      }
    }
    cout << endl;
  }
  cout << endl << endl << endl;

  for(int i=0;i<4;i++){
    cout << i << " ";
    for(int j=0;j<3;j++){
      for(int k=0;k<5;k++){
	cout << bonusfits::normfits_off3_paris_SLAC_off3_lc0_pw_beam4[j][k][i] << " ";
	if(j>0) cout << bonusfits::normfits_off3_paris_SLAC_off3_lc0_pw_beam5[j-1][k][i] << " ";
      }
    }
    cout << endl;
  }
  cout << endl << endl << endl;
  
}
