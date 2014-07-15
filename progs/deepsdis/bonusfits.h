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


//fit parameters for data normalization with all other parameters fixed, no offshell and VNA
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


//fit parameters for data normalization with all other parameters fixed, no offshell and VNA
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



#endif