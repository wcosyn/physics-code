#ifndef BONUSFITS_H
#define BONUSFITS_H

struct bonusfits{

const static double normfits_fix3_off3_lc0_beam4[3][5][4];
const static double normfits_fix3_off3_lc0_beam5[2][5][4];
const static double normfits_fix3_off3_lc1_beam4[3][5][4];
const static double normfits_fix3_off3_lc1_beam5[2][5][4];
const static double normfits_fix3_off4_lc0_beam4[3][5][4];
const static double normfits_fix3_off4_lc0_beam5[2][5][4];
const static double normfits_fix3_off4_lc1_beam4[3][5][4];
const static double normfits_fix3_off4_lc1_beam5[2][5][4];

const static double normfits_sigmadeeps_fix3_off3_lc0_q2dep0_beam4[3][5][4];
const static double normfits_sigmadeeps_fix3_off3_lc0_q2dep0_beam5[2][5][4];
const static double normfits_sigmadeeps_fix3_off3_lc0_q2dep1_beam4[3][5][4];
const static double normfits_sigmadeeps_fix3_off3_lc0_q2dep1_beam5[2][5][4];
const static double normfits_sigmadeeps_fix3_off4_lc0_q2dep0_beam4[3][5][4];
const static double normfits_sigmadeeps_fix3_off4_lc0_q2dep0_beam5[2][5][4];
const static double normfits_sigmadeeps_fix3_off4_lc0_q2dep1_beam4[3][5][4];
const static double normfits_sigmadeeps_fix3_off4_lc0_q2dep1_beam5[2][5][4];


const static double normfits_off3_paris_CB_off3_lc0_pw_beam4[3][5][4];
const static double normfits_off3_paris_CB_off3_lc0_pw_beam5[2][5][4];
const static double normfits_off3_paris_CB_off3_lc1_pw_beam4[3][5][4];
const static double normfits_off3_paris_CB_off3_lc1_pw_beam5[2][5][4];
const static double normfits_off3_AV18_CB_off3_lc0_pw_beam4[3][5][4];
const static double normfits_off3_AV18_CB_off3_lc0_pw_beam5[2][5][4];
const static double normfits_off3_paris_SLAC_off3_lc0_pw_beam4[3][5][4];
const static double normfits_off3_paris_SLAC_off3_lc0_pw_beam5[2][5][4];
const static double normfits_off3_paris_CB_off3_lc0_fsi_q2dep0_beam4[3][5][4];
const static double normfits_off3_paris_CB_off3_lc0_fsi_q2dep0_beam5[2][5][4];
const static double normfits_off3_paris_CB_off3_lc0_fsi_q2dep1_beam4[3][5][4];
const static double normfits_off3_paris_CB_off3_lc0_fsi_q2dep1_beam5[2][5][4];
const static double normfits_off3_paris_CB_off3_lc1_fsi_q2dep0_beam4[3][5][4];
const static double normfits_off3_paris_CB_off3_lc1_fsi_q2dep0_beam5[2][5][4];
const static double normfits_off3_paris_CB_off3_lc1_fsi_q2dep1_beam4[3][5][4];
const static double normfits_off3_paris_CB_off3_lc1_fsi_q2dep1_beam5[2][5][4];
const static double normfits_off3_AV18_CB_off3_lc0_fsi_q2dep1_beam4[3][5][4];
const static double normfits_off3_AV18_CB_off3_lc0_fsi_q2dep1_beam5[2][5][4];
const static double normfits_off3_paris_SLAC_off3_lc0_fsi_q2dep1_beam4[3][5][4];
const static double normfits_off3_paris_SLAC_off3_lc0_fsi_q2dep1_beam5[2][5][4];
const static double normfits_off3_paris_Alekhin_off3_lc0_fsi_q2dep1_beam4[3][5][4];
const static double normfits_off3_paris_Alekhin_off3_lc0_fsi_q2dep1_beam5[2][5][4];





};
//fit parameters for data normalization with all other parameters fixed, no offshell and VNA
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_fix3_off3_lc0_beam4[3][5][4]={ 
{ 
{ 2.07487 , 2.19835 , 2.66858 , 3.04783 }
,
{ 1.42899 , 1.44783 , 1.52886 , 1.63854 }
,
{ 1.25514 , 1.22437 , 1.31444 , 1.40801 }
,
{ 1.27716 , 1.27665 , 1.35659 , 1.50608 }
,
{ 1.19976 , 1.16977 , 1.23044 , 1.36489 }

},
{ 
{ 2.2721 , 2.42396 , 2.8341 , 3.41139 }
,
{ 1.566 , 1.52146 , 1.68302 , 1.8164 }
,
{ 1.25177 , 1.2504 , 1.30427 , 1.32778 }
,
{ 1.35496 , 1.32783 , 1.42987 , 1.55267 }
,
{ 1.2879 , 1.23995 , 1.35205 , 1.51351 }

},
{ 
{ 2.50028 , 2.51602 , 2.93264 , 2.80703 }
,
{ 1.7345 , 1.82157 , 1.81671 , 1.75733 }
,
{ 1.23362 , 1.14081 , 1.23171 , 1.3066 }
,
{ 1.4927 , 1.40698 , 1.58191 , 1.73406 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_fix3_off3_lc0_beam5[2][5][4]={ 
{ 
{ 2.08457 , 2.39333 , 2.67035 , 3.35753 }
,
{ 1.47213 , 1.54952 , 1.63679 , 1.89244 }
,
{ 1.21796 , 1.23697 , 1.27735 , 1.42833 }
,
{ 1.26332 , 1.29684 , 1.35119 , 1.52857 }
,
{ 1.12609 , 1.14418 , 1.19374 , 1.3631 }

},
{ 
{ 2.18341 , 2.56946 , 2.59072 , 3.10616 }
,
{ 1.69735 , 1.72793 , 1.87347 , 2.11095 }
,
{ 1.15925 , 1.17993 , 1.19815 , 1.29663 }
,
{ 1.31273 , 1.38583 , 1.40873 , 1.68692 }
,
{ 1.23162 , 1.27959 , 1.28961 , 1.56226 }

},
};


//fit parameters for data normalization with all other parameters fixed, full offshell and VNA
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_fix3_off4_lc0_beam4[3][5][4]={ 
{ 
{ 2.05684 , 2.16853 , 2.61006 , 2.91154 }
,
{ 1.41709 , 1.4262 , 1.49197 , 1.56164 }
,
{ 1.24507 , 1.20681 , 1.28234 , 1.34085 }
,
{ 1.26743 , 1.25822 , 1.32373 , 1.43603 }
,
{ 1.19211 , 1.15624 , 1.20764 , 1.31858 }

},
{ 
{ 2.25283 , 2.39123 , 2.77278 , 3.28113 }
,
{ 1.55406 , 1.50129 , 1.64544 , 1.73807 }
,
{ 1.24177 , 1.23296 , 1.27329 , 1.26853 }
,
{ 1.34522 , 1.31006 , 1.39759 , 1.48709 }
,
{ 1.28143 , 1.23279 , 1.3386 , 1.49413 }

},
{ 
{ 2.48254 , 2.48878 , 2.89705 , 2.73917 }
,
{ 1.7245 , 1.8021 , 1.7799 , 1.71036 }
,
{ 1.22709 , 1.12812 , 1.20715 , 1.25955 }
,
{ 1.48712 , 1.40051 , 1.57165 , 1.71518 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_fix3_off4_lc0_beam5[2][5][4]={ 
{ 
{ 2.0666 , 2.35821 , 2.60857 , 3.22484 }
,
{ 1.46009 , 1.52868 , 1.59931 , 1.80977 }
,
{ 1.20799 , 1.21929 , 1.24705 , 1.36356 }
,
{ 1.25372 , 1.27862 , 1.31977 , 1.45844 }
,
{ 1.11839 , 1.12988 , 1.16962 , 1.31118 }

},
{ 
{ 2.1686 , 2.53928 , 2.54128 , 3.01317 }
,
{ 1.68667 , 1.70792 , 1.83336 , 2.03578 }
,
{ 1.15122 , 1.16412 , 1.16878 , 1.24146 }
,
{ 1.30405 , 1.36913 , 1.38048 , 1.62046 }
,
{ 1.22687 , 1.27309 , 1.28107 , 1.54279 }

},
};


//fit parameters for data normalization with all other parameters fixed, no offshell and VNA
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_fix3_off3_lc1_beam4[3][5][4]={ 
{ 
{ 2.05407 , 2.17418 , 2.61588 , 2.91775 }
,
{ 1.41516 , 1.42197 , 1.48866 , 1.557 }
,
{ 1.24289 , 1.20408 , 1.2772 , 1.3312 }
,
{ 1.26542 , 1.25389 , 1.31686 , 1.42619 }
,
{ 1.1942 , 1.16136 , 1.21622 , 1.33293 }

},
{ 
{ 2.24925 , 2.39217 , 2.78171 , 3.32232 }
,
{ 1.55553 , 1.50187 , 1.64681 , 1.74195 }
,
{ 1.24072 , 1.23115 , 1.26862 , 1.26196 }
,
{ 1.34503 , 1.31002 , 1.39437 , 1.48341 }
,
{ 1.28873 , 1.25791 , 1.37644 , 1.57079 }

},
{ 
{ 2.48532 , 2.50873 , 2.96196 , 2.84325 }
,
{ 1.73444 , 1.81395 , 1.783 , 1.74166 }
,
{ 1.23515 , 1.13559 , 1.2155 , 1.27328 }
,
{ 1.50366 , 1.43811 , 1.62792 , 1.81455 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_fix3_off3_lc1_beam5[2][5][4]={ 
{ 
{ 2.06246 , 2.35402 , 2.60552 , 3.25257 }
,
{ 1.45884 , 1.52944 , 1.59862 , 1.80033 }
,
{ 1.20627 , 1.21713 , 1.24286 , 1.35285 }
,
{ 1.25308 , 1.27597 , 1.31545 , 1.44547 }
,
{ 1.11932 , 1.13152 , 1.17242 , 1.31632 }

},
{ 
{ 2.17395 , 2.54886 , 2.56193 , 3.06017 }
,
{ 1.69293 , 1.71832 , 1.83677 , 2.04985 }
,
{ 1.15112 , 1.16342 , 1.16268 , 1.23783 }
,
{ 1.30718 , 1.3729 , 1.38548 , 1.62568 }
,
{ 1.24019 , 1.3044 , 1.32793 , 1.62899 }

},
};


//fit parameters for data normalization with all other parameters fixed, full offshell and VNA
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_fix3_off4_lc1_beam4[3][5][4]={ 
{ 
{ 2.03408 , 2.13931 , 2.54505 , 2.75151 }
,
{ 1.40215 , 1.39727 , 1.44518 , 1.46582 }
,
{ 1.23181 , 1.18409 , 1.23984 , 1.25284 }
,
{ 1.25477 , 1.23308 , 1.27876 , 1.34434 }
,
{ 1.18582 , 1.14562 , 1.18846 , 1.27436 }

},
{ 
{ 2.22789 , 2.35424 , 2.70731 , 3.15467 }
,
{ 1.54224 , 1.47851 , 1.60165 , 1.64663 }
,
{ 1.22985 , 1.21111 , 1.23225 , 1.1917 }
,
{ 1.33411 , 1.28961 , 1.35625 , 1.40479 }
,
{ 1.28155 , 1.24936 , 1.35956 , 1.54459 }

},
{ 
{ 2.46547 , 2.47588 , 2.91224 , 2.74225 }
,
{ 1.72297 , 1.79097 , 1.73809 , 1.67517 }
,
{ 1.22789 , 1.12048 , 1.18545 , 1.21346 }
,
{ 1.49723 , 1.42999 , 1.61511 , 1.7886 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_fix3_off4_lc1_beam5[2][5][4]={ 
{ 
{ 2.04263 , 2.31362 , 2.531 , 3.08495 }
,
{ 1.44552 , 1.50528 , 1.5538 , 1.70026 }
,
{ 1.19551 , 1.19703 , 1.20718 , 1.27637 }
,
{ 1.2426 , 1.2551 , 1.27856 , 1.36363 }
,
{ 1.11095 , 1.115 , 1.14349 , 1.25259 }

},
{ 
{ 2.15718 , 2.51287 , 2.50038 , 2.93294 }
,
{ 1.68064 , 1.69408 , 1.78847 , 1.95163 }
,
{ 1.14205 , 1.14522 , 1.12857 , 1.17104 }
,
{ 1.29749 , 1.35296 , 1.35122 , 1.54398 }
,
{ 1.23486 , 1.29649 , 1.31721 , 1.6022 }

},
};

//fit parameters for data normalization with all other parameters fixed,
//no offshell and VNA, sigma taken from deeps parametrization, no Q2 dependence
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_sigmadeeps_fix3_off3_lc0_q2dep0_beam4[3][5][4]={ 
{ 
{ 2.2037 , 2.38073 , 2.97073 , 3.53992 }
,
{ 1.34979 , 1.34461 , 1.38529 , 1.42484 }
,
{ 1.2063 , 1.1627 , 1.22731 , 1.27878 }
,
{ 1.25311 , 1.24573 , 1.31311 , 1.43833 }
,
{ 1.20858 , 1.18105 , 1.24653 , 1.39108 }

},
{ 
{ 2.41221 , 2.62368 , 3.15212 , 3.97714 }
,
{ 1.47924 , 1.41239 , 1.52429 , 1.57699 }
,
{ 1.20316 , 1.18762 , 1.21828 , 1.20467 }
,
{ 1.32955 , 1.29555 , 1.38379 , 1.48116 }
,
{ 1.29748 , 1.2525 , 1.37089 , 1.54589 }

},
{ 
{ 2.65697 , 2.72873 , 3.28565 , 3.31683 }
,
{ 1.637 , 1.68795 , 1.64426 , 1.50843 }
,
{ 1.18475 , 1.08228 , 1.14802 , 1.17925 }
,
{ 1.46409 , 1.37057 , 1.52498 , 1.64047 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_sigmadeeps_fix3_off3_lc0_q2dep0_beam5[2][5][4]={ 
{ 
{ 2.21317 , 2.58837 , 2.96341 , 3.91359 }
,
{ 1.3911 , 1.43861 , 1.48308 , 1.64547 }
,
{ 1.17066 , 1.17512 , 1.19305 , 1.2968 }
,
{ 1.23963 , 1.2655 , 1.3078 , 1.45988 }
,
{ 1.13432 , 1.15508 , 1.20908 , 1.38851 }

},
{ 
{ 2.32065 , 2.78377 , 2.88701 , 3.64775 }
,
{ 1.60224 , 1.60267 , 1.69635 , 1.82625 }
,
{ 1.114 , 1.12069 , 1.11996 , 1.17597 }
,
{ 1.28803 , 1.35202 , 1.36256 , 1.60842 }
,
{ 1.24097 , 1.29254 , 1.30812 , 1.59548 }

},
};


//fit parameters for data normalization with all other parameters fixed, full offshell and VNA,
//sigma taken from deeps parametrization, no Q2 dependence
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_sigmadeeps_fix3_off4_lc0_q2dep0_beam4[3][5][4]={ 
{ 
{ 2.16981 , 2.32127 , 2.84396 , 3.21137 }
,
{ 1.34573 , 1.33753 , 1.37414 , 1.40418 }
,
{ 1.20098 , 1.15383 , 1.21182 , 1.24877 }
,
{ 1.24563 , 1.23196 , 1.28902 , 1.38876 }
,
{ 1.2002 , 1.16612 , 1.22112 , 1.33859 }

},
{ 
{ 2.37603 , 2.55839 , 3.01996 , 3.65561 }
,
{ 1.47513 , 1.40584 , 1.51296 , 1.55619 }
,
{ 1.19794 , 1.17874 , 1.20332 , 1.17847 }
,
{ 1.32196 , 1.28218 , 1.36011 , 1.43517 }
,
{ 1.2904 , 1.24461 , 1.3559 , 1.52389 }

},
{ 
{ 2.62344 , 2.67369 , 3.20393 , 3.13536 }
,
{ 1.63353 , 1.68168 , 1.63327 , 1.4969 }
,
{ 1.18136 , 1.07586 , 1.13625 , 1.1587 }
,
{ 1.45979 , 1.3657 , 1.5174 , 1.62747 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_sigmadeeps_fix3_off4_lc0_q2dep0_beam5[2][5][4]={ 
{ 
{ 2.17911 , 2.51838 , 2.83045 , 3.58695 }
,
{ 1.38698 , 1.43185 , 1.47177 , 1.62351 }
,
{ 1.16553 , 1.16625 , 1.1784 , 1.26804 }
,
{ 1.23232 , 1.25173 , 1.28476 , 1.4104 }
,
{ 1.12593 , 1.13925 , 1.18224 , 1.32989 }

},
{ 
{ 2.29257 , 2.72292 , 2.77898 , 3.40681 }
,
{ 1.59863 , 1.59607 , 1.68451 , 1.80703 }
,
{ 1.10972 , 1.1127 , 1.10582 , 1.15164 }
,
{ 1.28136 , 1.33922 , 1.34198 , 1.56199 }
,
{ 1.23563 , 1.28539 , 1.29835 , 1.57362 }

},
};


//fit parameters for data normalization with all other parameters fixed, no offshell and VNA,
//sigma taken from deeps parametrization, Q2 dependence
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_sigmadeeps_fix3_off3_lc0_q2dep1_beam4[3][5][4]={ 
{ 
{ 2.2037 , 2.38073 , 2.97073 , 3.53992 }
,
{ 1.34979 , 1.34461 , 1.38529 , 1.42484 }
,
{ 1.2063 , 1.1627 , 1.22731 , 1.27878 }
,
{ 1.25311 , 1.24573 , 1.31311 , 1.43833 }
,
{ 1.20858 , 1.18105 , 1.24653 , 1.39108 }

},
{ 
{ 2.41221 , 2.62368 , 3.15212 , 3.97714 }
,
{ 1.47924 , 1.41239 , 1.52429 , 1.57699 }
,
{ 1.20316 , 1.18762 , 1.21828 , 1.20467 }
,
{ 1.32955 , 1.29555 , 1.38379 , 1.48116 }
,
{ 1.29748 , 1.2525 , 1.37089 , 1.54589 }

},
{ 
{ 2.46814 , 2.47299 , 2.86276 , 2.70894 }
,
{ 1.60865 , 1.6497 , 1.59556 , 1.44084 }
,
{ 1.15465 , 1.04667 , 1.09794 , 1.10494 }
,
{ 1.41188 , 1.30536 , 1.42531 , 1.48277 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_sigmadeeps_fix3_off3_lc0_q2dep1_beam5[2][5][4]={ 
{ 
{ 2.21317 , 2.58837 , 2.96341 , 3.91359 }
,
{ 1.3911 , 1.43861 , 1.48308 , 1.64547 }
,
{ 1.17066 , 1.17512 , 1.19305 , 1.2968 }
,
{ 1.23963 , 1.2655 , 1.3078 , 1.45988 }
,
{ 1.13432 , 1.15508 , 1.20908 , 1.38851 }

},
{ 
{ 2.15522 , 2.52618 , 2.53191 , 3.00028 }
,
{ 1.57476 , 1.56672 , 1.64654 , 1.74797 }
,
{ 1.0861 , 1.0846 , 1.07295 , 1.10496 }
,
{ 1.2434 , 1.29117 , 1.28127 , 1.47257 }
,
{ 1.18156 , 1.20998 , 1.19308 , 1.38978 }

},
};


//fit parameters for data normalization with all other parameters fixed, full offshell and VNA,
//sigma taken from deeps parametrization, Q2 dependence
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_sigmadeeps_fix3_off4_lc0_q2dep1_beam4[3][5][4]={ 
{ 
{ 2.16981 , 2.32127 , 2.84396 , 3.21137 }
,
{ 1.34573 , 1.33753 , 1.37414 , 1.40418 }
,
{ 1.20098 , 1.15383 , 1.21182 , 1.24877 }
,
{ 1.24563 , 1.23196 , 1.28902 , 1.38876 }
,
{ 1.2002 , 1.16612 , 1.22112 , 1.33859 }

},
{ 
{ 2.37603 , 2.55839 , 3.01996 , 3.65561 }
,
{ 1.47513 , 1.40584 , 1.51296 , 1.55619 }
,
{ 1.19794 , 1.17874 , 1.20332 , 1.17847 }
,
{ 1.32196 , 1.28218 , 1.36011 , 1.43517 }
,
{ 1.2904 , 1.24461 , 1.3559 , 1.52389 }

},
{ 
{ 2.45325 , 2.45049 , 2.83417 , 2.65632 }
,
{ 1.60689 , 1.64657 , 1.59015 , 1.43554 }
,
{ 1.15296 , 1.04354 , 1.09237 , 1.09586 }
,
{ 1.4098 , 1.30305 , 1.42184 , 1.47742 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_sigmadeeps_fix3_off4_lc0_q2dep1_beam5[2][5][4]={ 
{ 
{ 2.17911 , 2.51838 , 2.83045 , 3.58695 }
,
{ 1.38698 , 1.43185 , 1.47177 , 1.62351 }
,
{ 1.16553 , 1.16625 , 1.1784 , 1.26804 }
,
{ 1.23232 , 1.25173 , 1.28476 , 1.4104 }
,
{ 1.12593 , 1.13925 , 1.18224 , 1.32989 }

},
{ 
{ 2.14271 , 2.50107 , 2.49162 , 2.92859 }
,
{ 1.57286 , 1.56343 , 1.64071 , 1.73901 }
,
{ 1.08401 , 1.08077 , 1.06638 , 1.09395 }
,
{ 1.24025 , 1.28535 , 1.27204 , 1.45304 }
,
{ 1.17909 , 1.20674 , 1.18901 , 1.38145 }

},
};

//fit parameters for data normalization, no offshell and VNA
//paris wf, Christy Bosted F2, PLANE WAVE fit
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_paris_CB_off3_lc0_pw_beam4[3][5][4]={ 
{ 
{ 1.86672 , 1.90922 , 2.25352 , 2.37885 }
,
{ 1.30798 , 1.29062 , 1.31424 , 1.30922 }
,
{ 1.1477 , 1.09482 , 1.11997 , 1.13039 }
,
{ 1.16273 , 1.1314 , 1.15543 , 1.19954 }
,
{ 1.09137 , 1.03383 , 1.04126 , 1.06884 }

},
{ 
{ 2.05105 , 2.11068 , 2.40222 , 2.63315 }
,
{ 1.43381 , 1.35977 , 1.43667 , 1.44682 }
,
{ 1.14836 , 1.11442 , 1.11659 , 1.06002 }
,
{ 1.23362 , 1.17639 , 1.21672 , 1.23057 }
,
{ 1.17044 , 1.08933 , 1.133 , 1.15812 }

},
{ 
{ 2.23411 , 2.18937 , 2.44007 , 2.1517 }
,
{ 1.58807 , 1.61328 , 1.5576 , 1.37929 }
,
{ 1.12818 , 1.01793 , 1.05047 , 1.03286 }
,
{ 1.35587 , 1.23617 , 1.32175 , 1.3265 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//fit parameters for data normalization, no offshell and VNA
//paris wf, Christy Bosted F2, PLANE WAVE fit
//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_paris_CB_off3_lc0_pw_beam5[2][5][4]={ 
{ 
{ 1.8817 , 2.09013 , 2.27026 , 2.59111 }
,
{ 1.34868 , 1.38508 , 1.39831 , 1.51091 }
,
{ 1.11759 , 1.10325 , 1.09369 , 1.14147 }
,
{ 1.15073 , 1.14958 , 1.15046 , 1.21766 }
,
{ 1.02492 , 1.01252 , 1.01275 , 1.07375 }

},
{ 
{ 1.95143 , 2.24338 , 2.16798 , 2.39726 }
,
{ 1.5561 , 1.53324 , 1.60754 , 1.67477 }
,
{ 1.0614 , 1.05473 , 1.02822 , 1.03579 }
,
{ 1.19518 , 1.22652 , 1.19581 , 1.33407 }
,
{ 1.11855 , 1.12423 , 1.07803 , 1.19516 }

},
};


//fit parameters for data normalization, no offshell and LC density
//paris wf, Christy Bosted F2, PLANE WAVE fit
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_paris_CB_off3_lc1_pw_beam4[3][5][4]={ 
{ 
{ 1.86543 , 1.90653 , 2.24942 , 2.37045 }
,
{ 1.30747 , 1.28872 , 1.31078 , 1.30213 }
,
{ 1.14741 , 1.09331 , 1.11677 , 1.12286 }
,
{ 1.16258 , 1.12988 , 1.15193 , 1.19177 }
,
{ 1.09071 , 1.0327 , 1.04049 , 1.06856 }

},
{ 
{ 2.05004 , 2.10787 , 2.39846 , 2.63297 }
,
{ 1.43268 , 1.35793 , 1.43382 , 1.44168 }
,
{ 1.14786 , 1.11292 , 1.11346 , 1.05436 }
,
{ 1.23316 , 1.1749 , 1.21379 , 1.22508 }
,
{ 1.16899 , 1.08867 , 1.13588 , 1.17062 }

},
{ 
{ 2.23217 , 2.18672 , 2.44528 , 2.16948 }
,
{ 1.58566 , 1.61159 , 1.55529 , 1.38491 }
,
{ 1.1264 , 1.01678 , 1.04935 , 1.03259 }
,
{ 1.35298 , 1.23532 , 1.3266 , 1.34303 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//fit parameters for data normalization, no offshell and LC density
//paris wf, Christy Bosted F2, PLANE WAVE fit
//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_paris_CB_off3_lc1_pw_beam5[2][5][4]={ 
{ 
{ 1.88082 , 2.08702 , 2.2651 , 2.58791 }
,
{ 1.34793 , 1.38315 , 1.39522 , 1.50319 }
,
{ 1.11721 , 1.10169 , 1.0907 , 1.13455 }
,
{ 1.15029 , 1.14806 , 1.14746 , 1.20941 }
,
{ 1.02442 , 1.01134 , 1.01123 , 1.07112 }

},
{ 
{ 1.94915 , 2.24065 , 2.16671 , 2.40288 }
,
{ 1.55428 , 1.53134 , 1.60498 , 1.67507 }
,
{ 1.06093 , 1.05334 , 1.02502 , 1.03138 }
,
{ 1.19421 , 1.22505 , 1.19426 , 1.33075 }
,
{ 1.11627 , 1.1235 , 1.0821 , 1.20921 }

},
};


//fit parameters for data normalization, no offshell and VNA
//paris wf, SLAC F2, PLANE WAVE fit
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_paris_SLAC_off3_lc0_pw_beam4[3][5][4]={ 
{ 
{ 3.51926 , 3.60019 , 4.25044 , 4.48915 }
,
{ 1.48727 , 1.46748 , 1.49446 , 1.4891 }
,
{ 1.18412 , 1.12943 , 1.1551 , 1.16545 }
,
{ 1.24042 , 1.20644 , 1.23144 , 1.27747 }
,
{ 1.201 , 1.13736 , 1.14498 , 1.17417 }

},
{ 
{ 3.95262 , 4.06864 , 4.63239 , 5.0 }
,
{ 1.61801 , 1.53461 , 1.62169 , 1.63382 }
,
{ 1.19669 , 1.16125 , 1.16344 , 1.10467 }
,
{ 1.25696 , 1.19848 , 1.23919 , 1.25293 }
,
{ 1.27948 , 1.19139 , 1.23908 , 1.26684 }

},
{ 
{ 5.0 , 5.0 , 5.0 , 5.0 }
,
{ 1.78495 , 1.81356 , 1.75238 , 1.55349 }
,
{ 1.2421 , 1.12172 , 1.15833 , 1.14024 }
,
{ 1.37193 , 1.25024 , 1.33702 , 1.34198 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//fit parameters for data normalization, no offshell and VNA
//paris wf, SLAC F2, PLANE WAVE fit
//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_paris_SLAC_off3_lc0_pw_beam5[2][5][4]={ 
{ 
{ 3.61755 , 4.019 , 4.36664 , 4.98495 }
,
{ 1.52682 , 1.5681 , 1.58343 , 1.71123 }
,
{ 1.1736 , 1.15855 , 1.14854 , 1.19868 }
,
{ 1.18645 , 1.18509 , 1.18585 , 1.25462 }
,
{ 1.12924 , 1.11555 , 1.11563 , 1.18239 }

},
{ 
{ 5.0 , 5.0 , 5.0 , 5.0 }
,
{ 1.75699 , 1.73171 , 1.81618 , 1.8936 }
,
{ 1.15235 , 1.14555 , 1.11771 , 1.12653 }
,
{ 1.1962 , 1.22797 , 1.19763 , 1.33699 }
,
{ 1.22729 , 1.23243 , 1.18122 , 1.30903 }

},
};


//fit parameters for data normalization, no offshell and VNA
//AV18 wf, Christy Bosted F2, PLANE WAVE fit
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_AV18_CB_off3_lc0_pw_beam4[3][5][4]={ 
{ 
{ 1.86672 , 1.90922 , 2.25352 , 2.37885 }
,
{ 1.30798 , 1.29062 , 1.31424 , 1.30922 }
,
{ 1.1477 , 1.09482 , 1.11997 , 1.13039 }
,
{ 1.16273 , 1.1314 , 1.15543 , 1.19954 }
,
{ 1.09137 , 1.03383 , 1.04126 , 1.06884 }

},
{ 
{ 2.05105 , 2.11068 , 2.40222 , 2.63315 }
,
{ 1.43381 , 1.35977 , 1.43667 , 1.44682 }
,
{ 1.14836 , 1.11442 , 1.11659 , 1.06002 }
,
{ 1.23362 , 1.17639 , 1.21672 , 1.23057 }
,
{ 1.17044 , 1.08933 , 1.133 , 1.15812 }

},
{ 
{ 2.23411 , 2.18937 , 2.44007 , 2.1517 }
,
{ 1.58807 , 1.61328 , 1.5576 , 1.37929 }
,
{ 1.12818 , 1.01793 , 1.05047 , 1.03286 }
,
{ 1.35587 , 1.23617 , 1.32175 , 1.3265 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//fit parameters for data normalization, no offshell and VNA
//AV18 wf, Christy Bosted F2, PLANE WAVE fit
//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_AV18_CB_off3_lc0_pw_beam5[2][5][4]={ 
{ 
{ 1.8817 , 2.09013 , 2.27026 , 2.59111 }
,
{ 1.34868 , 1.38508 , 1.39831 , 1.51091 }
,
{ 1.11759 , 1.10325 , 1.09369 , 1.14147 }
,
{ 1.15073 , 1.14958 , 1.15046 , 1.21766 }
,
{ 1.02492 , 1.01252 , 1.01275 , 1.07375 }

},
{ 
{ 1.95143 , 2.24338 , 2.16798 , 2.39726 }
,
{ 1.5561 , 1.53324 , 1.60754 , 1.67477 }
,
{ 1.0614 , 1.05473 , 1.02822 , 1.03579 }
,
{ 1.19518 , 1.22652 , 1.19581 , 1.33407 }
,
{ 1.11855 , 1.12423 , 1.07803 , 1.19516 }

},
};


//fit parameters for data normalization, no offshell and VNA
//paris wf, Christy Bosted F2, FSI fit, deeps parameters, no Q2 dependence
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_paris_CB_off3_lc0_fsi_q2dep0_beam4[3][5][4]={ 
{ 
{ 2.17895 , 2.33995 , 2.96154 , 3.49905 }
,
{ 1.35692 , 1.35336 , 1.40003 , 1.4324 }
,
{ 1.21117 , 1.17391 , 1.22879 , 1.28857 }
,
{ 1.25321 , 1.24584 , 1.31316 , 1.43855 }
,
{ 1.20865 , 1.1811 , 1.24654 , 1.39114 }

},
{ 
{ 2.39114 , 2.58473 , 3.15214 , 3.89988 }
,
{ 1.4874 , 1.42622 , 1.53087 , 1.58413 }
,
{ 1.21164 , 1.1948 , 1.22443 , 1.20978 }
,
{ 1.32954 , 1.29596 , 1.38382 , 1.48132 }
,
{ 1.29736 , 1.25252 , 1.37082 , 1.54601 }

},
{ 
{ 2.61096 , 2.69332 , 3.25984 , 3.28565 }
,
{ 1.64832 , 1.69416 , 1.66078 , 1.52034 }
,
{ 1.19167 , 1.09318 , 1.15478 , 1.18599 }
,
{ 1.46427 , 1.3708 , 1.52498 , 1.64097 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//fit parameters for data normalization, no offshell and VNA
//paris wf, Christy Bosted F2, FSI fit, deeps parameters, no Q2 dependence
//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_paris_CB_off3_lc0_fsi_q2dep0_beam5[2][5][4]={ 
{ 
{ 2.19348 , 2.55392 , 2.96343 , 3.83474 }
,
{ 1.39875 , 1.45254 , 1.48965 , 1.65306 }
,
{ 1.17913 , 1.18242 , 1.19924 , 1.30194 }
,
{ 1.23972 , 1.26592 , 1.30786 , 1.45994 }
,
{ 1.13435 , 1.15511 , 1.2091 , 1.38853 }

},
{ 
{ 2.28098 , 2.75281 , 2.85854 , 3.60978 }
,
{ 1.61478 , 1.60921 , 1.71328 , 1.83806 }
,
{ 1.12023 , 1.13084 , 1.12644 , 1.18241 }
,
{ 1.28816 , 1.35229 , 1.3626 , 1.6089 }
,
{ 1.24086 , 1.29258 , 1.30789 , 1.59573 }

},
};


//fit parameters for data normalization, no offshell and VNA
//paris wf, Christy Bosted F2, FSI fit, deeps parameters, Q2 dependence
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_paris_CB_off3_lc0_fsi_q2dep1_beam4[3][5][4]={ 
{ 
{ 2.17895 , 2.33995 , 2.96154 , 3.49905 }
,
{ 1.35692 , 1.35336 , 1.40003 , 1.4324 }
,
{ 1.21117 , 1.17391 , 1.22879 , 1.28857 }
,
{ 1.25321 , 1.24584 , 1.31316 , 1.43855 }
,
{ 1.20865 , 1.1811 , 1.24654 , 1.39114 }

},
{ 
{ 2.39114 , 2.58473 , 3.15214 , 3.89988 }
,
{ 1.4874 , 1.42622 , 1.53087 , 1.58413 }
,
{ 1.21164 , 1.1948 , 1.22443 , 1.20978 }
,
{ 1.32954 , 1.29596 , 1.38382 , 1.48132 }
,
{ 1.29736 , 1.25252 , 1.37082 , 1.54601 }

},
{ 
{ 2.42519 , 2.44108 , 2.8399 , 2.68707 }
,
{ 1.61977 , 1.65572 , 1.61152 , 1.45233 }
,
{ 1.16141 , 1.05718 , 1.10444 , 1.11129 }
,
{ 1.41205 , 1.30558 , 1.4253 , 1.48322 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//fit parameters for data normalization, no offshell and VNA
//paris wf, Christy Bosted F2, FSI fit, deeps parameters, Q2 dependence
//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_paris_CB_off3_lc0_fsi_q2dep1_beam5[2][5][4]={ 
{ 
{ 2.19348 , 2.55392 , 2.96343 , 3.83474 }
,
{ 1.39875 , 1.45254 , 1.48965 , 1.65306 }
,
{ 1.17913 , 1.18242 , 1.19924 , 1.30194 }
,
{ 1.23972 , 1.26592 , 1.30786 , 1.45994 }
,
{ 1.13435 , 1.15511 , 1.2091 , 1.38853 }

},
{ 
{ 2.11836 , 2.49819 , 2.50659 , 2.97396 }
,
{ 1.587 , 1.57311 , 1.66286 , 1.75939 }
,
{ 1.0922 , 1.0944 , 1.07915 , 1.11107 }
,
{ 1.24353 , 1.29143 , 1.28131 , 1.47301 }
,
{ 1.18158 , 1.20999 , 1.19309 , 1.38982 }

},
};


//fit parameters for data normalization, no offshell and LC density
//paris wf, Christy Bosted F2, FSI fit, deeps parameters, no Q2 dependence
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_paris_CB_off3_lc1_fsi_q2dep0_beam4[3][5][4]={ 
{ 
{ 2.17778 , 2.33629 , 2.95167 , 3.46008 }
,
{ 1.3564 , 1.35132 , 1.39589 , 1.42238 }
,
{ 1.21092 , 1.17223 , 1.22469 , 1.27697 }
,
{ 1.25311 , 1.24409 , 1.30831 , 1.42457 }
,
{ 1.20799 , 1.17973 , 1.24471 , 1.38632 }

},
{ 
{ 2.39035 , 2.58089 , 3.14239 , 3.87093 }
,
{ 1.48628 , 1.42425 , 1.52732 , 1.57604 }
,
{ 1.21117 , 1.19313 , 1.22038 , 1.20052 }
,
{ 1.32911 , 1.29423 , 1.37954 , 1.47013 }
,
{ 1.29582 , 1.25176 , 1.37403 , 1.56199 }

},
{ 
{ 2.60906 , 2.68969 , 3.26257 , 3.29267 }
,
{ 1.64586 , 1.69234 , 1.65776 , 1.52451 }
,
{ 1.18983 , 1.0919 , 1.15307 , 1.18346 }
,
{ 1.46119 , 1.36986 , 1.53053 , 1.66099 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//fit parameters for data normalization, no offshell and LC density
//paris wf, Christy Bosted F2, FSI fit, deeps parameters, no Q2 dependence
//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_paris_CB_off3_lc1_fsi_q2dep0_beam5[2][5][4]={ 
{ 
{ 2.19284 , 2.54973 , 2.95152 , 3.8016 }
,
{ 1.39803 , 1.45048 , 1.48583 , 1.64182 }
,
{ 1.17878 , 1.1807 , 1.19536 , 1.29094 }
,
{ 1.23932 , 1.26415 , 1.30354 , 1.44529 }
,
{ 1.13386 , 1.15366 , 1.20627 , 1.37999 }

},
{ 
{ 2.27895 , 2.74897 , 2.85268 , 3.59263 }
,
{ 1.61294 , 1.60716 , 1.70994 , 1.83531 }
,
{ 1.11977 , 1.1293 , 1.12235 , 1.17459 }
,
{ 1.28717 , 1.35058 , 1.35995 , 1.60014 }
,
{ 1.23837 , 1.29171 , 1.31275 , 1.61392 }

},
};


//fit parameters for data normalization, no offshell and LC density
//paris wf, Christy Bosted F2, FSI fit, deeps parameters, Q2 dependence
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_paris_CB_off3_lc1_fsi_q2dep1_beam4[3][5][4]={ 
{ 
{ 2.17778 , 2.33629 , 2.95167 , 3.46008 }
,
{ 1.3564 , 1.35132 , 1.39589 , 1.42238 }
,
{ 1.21092 , 1.17223 , 1.22469 , 1.27697 }
,
{ 1.25311 , 1.24409 , 1.30831 , 1.42457 }
,
{ 1.20799 , 1.17973 , 1.24471 , 1.38632 }

},
{ 
{ 2.39035 , 2.58089 , 3.14239 , 3.87093 }
,
{ 1.48628 , 1.42425 , 1.52732 , 1.57604 }
,
{ 1.21117 , 1.19313 , 1.22038 , 1.20052 }
,
{ 1.32911 , 1.29423 , 1.37954 , 1.47013 }
,
{ 1.29582 , 1.25176 , 1.37403 , 1.56199 }

},
{ 
{ 2.42326 , 2.43796 , 2.84416 , 2.70145 }
,
{ 1.61734 , 1.65396 , 1.60885 , 1.45725 }
,
{ 1.1596 , 1.05596 , 1.10302 , 1.10993 }
,
{ 1.40923 , 1.30471 , 1.43036 , 1.5015 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//fit parameters for data normalization, no offshell and LC density
//paris wf, Christy Bosted F2, FSI fit, deeps parameters, Q2 dependence
//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_paris_CB_off3_lc1_fsi_q2dep1_beam5[2][5][4]={ 
{ 
{ 2.19284 , 2.54973 , 2.95152 , 3.8016 }
,
{ 1.39803 , 1.45048 , 1.48583 , 1.64182 }
,
{ 1.17878 , 1.1807 , 1.19536 , 1.29094 }
,
{ 1.23932 , 1.26415 , 1.30354 , 1.44529 }
,
{ 1.13386 , 1.15366 , 1.20627 , 1.37999 }

},
{ 
{ 2.11633 , 2.49489 , 2.50332 , 2.97084 }
,
{ 1.58518 , 1.57114 , 1.6599 , 1.75818 }
,
{ 1.09173 , 1.09293 , 1.0755 , 1.105 }
,
{ 1.24255 , 1.28985 , 1.27923 , 1.46712 }
,
{ 1.17918 , 1.20922 , 1.19756 , 1.40587 }

},
};


//fit parameters for data normalization, no offshell and VNA
//AV18 wf, Christy Bosted F2, FSI fit, deeps parameters, Q2 dependence
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_AV18_CB_off3_lc0_fsi_q2dep1_beam4[3][5][4]={ 
{ 
{ 2.18076 , 2.34264 , 2.96641 , 3.50907 }
,
{ 1.35717 , 1.35369 , 1.40052 , 1.43321 }
,
{ 1.21147 , 1.17431 , 1.22938 , 1.28961 }
,
{ 1.25363 , 1.24641 , 1.31402 , 1.44015 }
,
{ 1.20919 , 1.18182 , 1.24762 , 1.39324 }

},
{ 
{ 2.39295 , 2.58745 , 3.15696 , 3.91047 }
,
{ 1.48766 , 1.42656 , 1.53138 , 1.58501 }
,
{ 1.21193 , 1.1952 , 1.22501 , 1.21074 }
,
{ 1.32998 , 1.29654 , 1.3847 , 1.48293 }
,
{ 1.29792 , 1.25327 , 1.372 , 1.54839 }

},
{ 
{ 2.42603 , 2.44221 , 2.84176 , 2.69025 }
,
{ 1.61992 , 1.65592 , 1.61179 , 1.45273 }
,
{ 1.16156 , 1.05737 , 1.10471 , 1.11174 }
,
{ 1.41229 , 1.30589 , 1.42578 , 1.48407 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//fit parameters for data normalization, no offshell and VNA
//AV18 wf, Christy Bosted F2, FSI fit, deeps parameters, Q2 dependence
//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_AV18_CB_off3_lc0_fsi_q2dep1_beam5[2][5][4]={ 
{ 
{ 2.19515 , 2.55663 , 2.96797 , 3.84502 }
,
{ 1.399 , 1.45288 , 1.49015 , 1.65397 }
,
{ 1.17942 , 1.18282 , 1.19981 , 1.30297 }
,
{ 1.24013 , 1.26649 , 1.3087 , 1.46151 }
,
{ 1.13484 , 1.1558 , 1.21013 , 1.39059 }

},
{ 
{ 2.11905 , 2.49933 , 2.50827 , 2.9774 }
,
{ 1.58714 , 1.5733 , 1.66314 , 1.75989 }
,
{ 1.09234 , 1.09459 , 1.07941 , 1.11152 }
,
{ 1.24375 , 1.29173 , 1.28174 , 1.47383 }
,
{ 1.18184 , 1.21036 , 1.19361 , 1.39044 }

},
};



//fit parameters for data normalization, no offshell and VNA
//paris wf, SLAC F2, FSI fit, deeps parameters, Q2 dependence
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_paris_SLAC_off3_lc0_fsi_q2dep1_beam4[3][5][4]={ 
{ 
{ 2.17895 , 2.33995 , 2.96154 , 3.49905 }
,
{ 1.35692 , 1.35336 , 1.40003 , 1.4324 }
,
{ 1.21117 , 1.17391 , 1.22879 , 1.28857 }
,
{ 1.25321 , 1.24584 , 1.31316 , 1.43855 }
,
{ 1.20865 , 1.1811 , 1.24654 , 1.39114 }

},
{ 
{ 2.39114 , 2.58473 , 3.15214 , 3.89988 }
,
{ 1.4874 , 1.42622 , 1.53087 , 1.58413 }
,
{ 1.21164 , 1.1948 , 1.22443 , 1.20978 }
,
{ 1.32954 , 1.29596 , 1.38382 , 1.48132 }
,
{ 1.29736 , 1.25252 , 1.37082 , 1.54601 }

},
{ 
{ 2.42519 , 2.44108 , 2.8399 , 2.68707 }
,
{ 1.61977 , 1.65572 , 1.61152 , 1.45233 }
,
{ 1.16141 , 1.05718 , 1.10444 , 1.11129 }
,
{ 1.41205 , 1.30558 , 1.4253 , 1.48322 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//fit parameters for data normalization, no offshell and VNA
//paris wf, SLAC F2, FSI fit, deeps parameters, Q2 dependence
//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_paris_SLAC_off3_lc0_fsi_q2dep1_beam5[2][5][4]={ 
{ 
{ 2.19348 , 2.55392 , 2.96343 , 3.83474 }
,
{ 1.39875 , 1.45254 , 1.48965 , 1.65306 }
,
{ 1.17913 , 1.18242 , 1.19924 , 1.30194 }
,
{ 1.23972 , 1.26592 , 1.30786 , 1.45994 }
,
{ 1.13435 , 1.15511 , 1.2091 , 1.38853 }

},
{ 
{ 2.11836 , 2.49819 , 2.50659 , 2.97396 }
,
{ 1.587 , 1.57311 , 1.66286 , 1.75939 }
,
{ 1.0922 , 1.0944 , 1.07915 , 1.11107 }
,
{ 1.24353 , 1.29143 , 1.28131 , 1.47301 }
,
{ 1.18158 , 1.20999 , 1.19309 , 1.38982 }

},
};


//fit parameters for data normalization, no offshell and VNA
//paris wf, Alekhin F2, FSI fit, deeps parameters, Q2 dependence
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_paris_Alekhin_off3_lc0_fsi_q2dep1_beam4[3][5][4]={ 
{ 
{ 4.8608 , 5.0 , 5.0 , 5.0 }
,
{ 1.36961 , 1.3665 , 1.41394 , 1.44753 }
,
{ 1.26834 , 1.22955 , 1.28735 , 1.35061 }
,
{ 1.18719 , 1.18054 , 1.24461 , 1.36401 }
,
{ 1.22497 , 1.19717 , 1.26394 , 1.4119 }

},
{ 
{ 5.0 , 5.0 , 5.0 , 5.0 }
,
{ 2.0353 , 1.95269 , 2.09666 , 2.17145 }
,
{ 1.49241 , 1.47224 , 1.50963 , 1.49282 }
,
{ 1.34649 , 1.31279 , 1.40239 , 1.5021 }
,
{ 1.37933 , 1.3306 , 1.45618 , 1.64117 }

},
{ 
{ 5.0 , 5.0 , 5.0 , 5.0 }
,
{ 4.49434 , 4.59645 , 4.48381 , 4.02556 }
,
{ 2.22585 , 2.02863 , 2.12293 , 2.13709 }
,
{ 1.66925 , 1.53766 , 1.67737 , 1.74094 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//fit parameters for data normalization, no offshell and VNA
//paris wf, Christy Bosted F2, FSI fit, deeps parameters, Q2 dependence
//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits::normfits_off3_paris_Alekhin_off3_lc0_fsi_q2dep1_beam5[2][5][4]={ 
{ 
{ 5.0 , 5.0 , 5.0 , 5.0 }
,
{ 1.88054 , 1.95312 , 2.00363 , 2.22504 }
,
{ 1.43261 , 1.43698 , 1.45779 , 1.58355 }
,
{ 1.24401 , 1.27055 , 1.3129 , 1.46631 }
,
{ 1.19336 , 1.21546 , 1.27244 , 1.46197 }

},
{ 
{ 5.0 , 5.0 , 5.0 , 5.0 }
,
{ 4.22004 , 4.18315 , 4.42724 , 4.68276 }
,
{ 1.97692 , 1.98231 , 1.95842 , 2.01696 }
,
{ 1.40454 , 1.45935 , 1.44849 , 1.66744 }
,
{ 1.35701 , 1.38728 , 1.36683 , 1.59073 }

},
};




#endif