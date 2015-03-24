#ifndef BONUSFITS2_H
#define BONUSFITS2_H

struct bonusfits2{

const static double normfits_paris_CB_off3_lc0_fsi_q2dep0_beam4[3][5][4];
const static double normfits_paris_CB_off3_lc0_fsi_q2dep0_beam5[2][5][4];
const static double normfits_paris_CB_off3_lc0_fsi_q2dep1_beam4[1][5][4];
const static double normfits_paris_CB_off3_lc0_fsi_q2dep1_beam5[1][5][4];
//
const static double normfits_AV18_CB_off3_lc0_fsi_q2dep0_beam4[3][5][4];
const static double normfits_AV18_CB_off3_lc0_fsi_q2dep0_beam5[2][5][4];
const static double normfits_AV18_CB_off3_lc0_fsi_q2dep1_beam4[1][5][4];
const static double normfits_AV18_CB_off3_lc0_fsi_q2dep1_beam5[1][5][4];
//
const static double normfits_CDBonn_CB_off3_lc0_fsi_q2dep0_beam4[3][5][4];
const static double normfits_CDBonn_CB_off3_lc0_fsi_q2dep0_beam5[2][5][4];
const static double normfits_CDBonn_CB_off3_lc0_fsi_q2dep1_beam4[1][5][4];
const static double normfits_CDBonn_CB_off3_lc0_fsi_q2dep1_beam5[1][5][4];
//
const static double normfits_GrossWJC1_CB_off3_lc0_fsi_q2dep0_beam4[3][5][4];
const static double normfits_GrossWJC1_CB_off3_lc0_fsi_q2dep0_beam5[2][5][4];
const static double normfits_GrossWJC1_CB_off3_lc0_fsi_q2dep1_beam4[1][5][4];
const static double normfits_GrossWJC1_CB_off3_lc0_fsi_q2dep1_beam5[1][5][4];
//
const static double normfits_paris_SLAC_off3_lc0_fsi_q2dep0_beam4[3][5][4];
const static double normfits_paris_SLAC_off3_lc0_fsi_q2dep0_beam5[2][5][4];
const static double normfits_paris_SLAC_off3_lc0_fsi_q2dep1_beam4[1][5][4];
const static double normfits_paris_SLAC_off3_lc0_fsi_q2dep1_beam5[1][5][4];

};
//fit parameters for data normalization with all other parameters fixed, no offshell and VNA


//fit parameters for data normalization, no offshell and VNA, no q2 dependence in fits
//paris wf, Christy Bosted F2, FSI fit
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_paris_CB_off3_lc0_fsi_q2dep0_beam4[3][5][4]={ 
{ 
{ 2.17845 , 2.33953 , 2.96072 , 3.49764 }
,
{ 1.35644 , 1.35278 , 1.39942 , 1.43166 }
,
{ 1.21039 , 1.17304 , 1.22768 , 1.28714 }
,
{ 1.25199 , 1.24442 , 1.31129 , 1.43587 }
,
{ 1.20535 , 1.1769 , 1.24072 , 1.3812 }

},
{ 
{ 2.39092 , 2.58448 , 3.15181 , 3.89969 }
,
{ 1.48725 , 1.42616 , 1.53078 , 1.5839 }
,
{ 1.21103 , 1.19404 , 1.22349 , 1.20876 }
,
{ 1.3279 , 1.29394 , 1.38097 , 1.47706 }
,
{ 1.2883 , 1.25014 , 1.36711 , 1.52716 }

},
{ 
{ 2.61151 , 2.69423 , 3.26238 , 3.28726 }
,
{ 1.65139 , 1.69766 , 1.66309 , 1.52242 }
,
{ 1.1928 , 1.09385 , 1.1551 , 1.18637 }
,
{ 1.47094 , 1.37421 , 1.54539 , 1.6548 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//fit parameters for data normalization, no offshell and VNA, no q2 dependence in fits
//paris wf, Christy Bosted F2, FSI fit
//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_paris_CB_off3_lc0_fsi_q2dep0_beam5[2][5][4]={ 
{ 
{ 2.193 , 2.55313 , 2.96253 , 3.83351 }
,
{ 1.39839 , 1.45218 , 1.48921 , 1.65232 }
,
{ 1.17854 , 1.18176 , 1.19845 , 1.3009 }
,
{ 1.23878 , 1.26466 , 1.30634 , 1.45765 }
,
{ 1.13277 , 1.15316 , 1.20651 , 1.38448 }

},
{ 
{ 2.28064 , 2.75214 , 2.8576 , 3.60744 }
,
{ 1.61588 , 1.61046 , 1.71424 , 1.83826 }
,
{ 1.12035 , 1.13104 , 1.12651 , 1.18279 }
,
{ 1.28691 , 1.3503 , 1.36017 , 1.60437 }
,
{ 1.23266 , 1.29629 , 1.31709 , 1.59312 }

},
};


//fit parameters for data normalization, no offshell and VNA, q2 dependence in fits
//paris wf, Christy Bosted F2, FSI fit
//Bonus data beam 4GeV
//first index is Q^2={3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_paris_CB_off3_lc0_fsi_q2dep1_beam4[1][5][4]={ 
{ 
{ 2.42564 , 2.44212 , 2.84258 , 2.68928 }
,
{ 1.62284 , 1.65913 , 1.61387 , 1.45452 }
,
{ 1.16256 , 1.05794 , 1.10489 , 1.11203 }
,
{ 1.4187 , 1.3101 , 1.44694 , 1.49761 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

},
};


//fit parameters for data normalization, no offshell and VNA, q2 dependence in fits
//paris wf, Christy Bosted F2, FSI fit
//Bonus data beam 5GeV
//first index is Q^2={3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_paris_CB_off3_lc0_fsi_q2dep1_beam5[1][5][4]={ 
{ 
{ 2.11797 , 2.49768 , 2.50604 , 2.97271 }
,
{ 1.58806 , 1.57432 , 1.66374 , 1.75988 }
,
{ 1.09235 , 1.0946 , 1.07926 , 1.11155 }
,
{ 1.24235 , 1.28975 , 1.27909 , 1.46928 }
,
{ 1.17439 , 1.21542 , 1.20469 , 1.39134 }

},
};


//fit parameters for data normalization, no offshell and VNA, no q2 dependence in fits
//AV18 wf, Christy Bosted F2, FSI fit
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_AV18_CB_off3_lc0_fsi_q2dep0_beam4[3][5][4]={ 
{ 
{ 2.18026 , 2.34223 , 2.9656 , 3.50764 }
,
{ 1.35668 , 1.35311 , 1.39991 , 1.43248 }
,
{ 1.21069 , 1.17343 , 1.22827 , 1.28818 }
,
{ 1.25241 , 1.24499 , 1.31214 , 1.43747 }
,
{ 1.20588 , 1.17762 , 1.24179 , 1.38328 }

},
{ 
{ 2.39274 , 2.58705 , 3.15662 , 3.91028 }
,
{ 1.48751 , 1.4265 , 1.53129 , 1.58478 }
,
{ 1.21133 , 1.19443 , 1.22407 , 1.20972 }
,
{ 1.32834 , 1.29452 , 1.38185 , 1.47866 }
,
{ 1.28886 , 1.25089 , 1.36829 , 1.52949 }

},
{ 
{ 2.61326 , 2.69668 , 3.26672 , 3.29553 }
,
{ 1.65166 , 1.69805 , 1.66362 , 1.52323 }
,
{ 1.19309 , 1.09421 , 1.15565 , 1.18731 }
,
{ 1.47142 , 1.37483 , 1.54637 , 1.65668 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//fit parameters for data normalization, no offshell and VNA, no q2 dependence in fits
//AV18 wf, Christy Bosted F2, FSI fit
//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_AV18_CB_off3_lc0_fsi_q2dep0_beam5[2][5][4]={ 
{ 
{ 2.19467 , 2.55584 , 2.96708 , 3.84379 }
,
{ 1.39864 , 1.45253 , 1.48971 , 1.65323 }
,
{ 1.17883 , 1.18215 , 1.19902 , 1.30193 }
,
{ 1.23919 , 1.26523 , 1.30718 , 1.45922 }
,
{ 1.13327 , 1.15385 , 1.20754 , 1.38654 }

},
{ 
{ 2.28218 , 2.75465 , 2.86145 , 3.6162 }
,
{ 1.61615 , 1.61082 , 1.71479 , 1.83925 }
,
{ 1.12063 , 1.13142 , 1.12704 , 1.18373 }
,
{ 1.28733 , 1.35091 , 1.36103 , 1.60613 }
,
{ 1.23319 , 1.29706 , 1.31821 , 1.59461 }

},
};


//fit parameters for data normalization, no offshell and VNA, q2 dependence in fits
//AV18 wf, Christy Bosted F2, FSI fit
//Bonus data beam 4GeV
//first index is Q^2={3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_AV18_CB_off3_lc0_fsi_q2dep1_beam4[1][5][4]={ 
{ 
{ 2.42648 , 2.44325 , 2.84445 , 2.69247 }
,
{ 1.62298 , 1.65932 , 1.61414 , 1.45492 }
,
{ 1.16271 , 1.05813 , 1.10517 , 1.11248 }
,
{ 1.41894 , 1.31041 , 1.44742 , 1.49847 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

},
};


//fit parameters for data normalization, no offshell and VNA, q2 dependence in fits
//AV18 wf, Christy Bosted F2, FSI fit
//Bonus data beam 5GeV
//first index is Q^2={3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_AV18_CB_off3_lc0_fsi_q2dep1_beam5[1][5][4]={ 
{ 
{ 2.11871 , 2.49884 , 2.50772 , 2.97615 }
,
{ 1.5882 , 1.57451 , 1.66402 , 1.76037 }
,
{ 1.09249 , 1.0948 , 1.07952 , 1.112 }
,
{ 1.24256 , 1.29005 , 1.27952 , 1.4701 }
,
{ 1.17465 , 1.21579 , 1.20521 , 1.39197 }

},
};


//fit parameters for data normalization, no offshell and VNA, no q2 dependence in fits
//CDBonn wf, Christy Bosted F2, FSI fit
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_CDBonn_CB_off3_lc0_fsi_q2dep0_beam4[3][5][4]={ 
{ 
{ 2.19614 , 2.36513 , 3.00493 , 3.57604 }
,
{ 1.35905 , 1.35624 , 1.40422 , 1.43875 }
,
{ 1.21377 , 1.17738 , 1.2338 , 1.29646 }
,
{ 1.25668 , 1.25063 , 1.32012 , 1.44997 }
,
{ 1.21126 , 1.18459 , 1.25164 , 1.39887 }

},
{ 
{ 2.40802 , 2.60917 , 3.19282 , 3.97501 }
,
{ 1.48987 , 1.4295 , 1.53557 , 1.59102 }
,
{ 1.21423 , 1.19822 , 1.22924 , 1.21697 }
,
{ 1.33272 , 1.30018 , 1.38995 , 1.49101 }
,
{ 1.29449 , 1.25813 , 1.37887 , 1.54605 }

},
{ 
{ 2.62708 , 2.71526 , 3.29552 , 3.33529 }
,
{ 1.65399 , 1.70117 , 1.66769 , 1.52822 }
,
{ 1.19571 , 1.09738 , 1.16008 , 1.19366 }
,
{ 1.47601 , 1.38049 , 1.55478 , 1.66908 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//fit parameters for data normalization, no offshell and VNA, no q2 dependence in fits
//CDBonn wf, Christy Bosted F2, FSI fit
//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_CDBonn_CB_off3_lc0_fsi_q2dep0_beam5[2][5][4]={ 
{ 
{ 2.2087 , 2.57759 , 3.00135 , 3.90762 }
,
{ 1.40087 , 1.45559 , 1.49388 , 1.65978 }
,
{ 1.18165 , 1.1859 , 1.20409 , 1.30975 }
,
{ 1.24328 , 1.27077 , 1.31484 , 1.47143 }
,
{ 1.13822 , 1.16056 , 1.21695 , 1.40192 }

},
{ 
{ 2.2942 , 2.77385 , 2.8878 , 3.66277 }
,
{ 1.61843 , 1.6138 , 1.71898 , 1.84544 }
,
{ 1.1231 , 1.13471 , 1.13142 , 1.19018 }
,
{ 1.29137 , 1.35651 , 1.36857 , 1.61872 }
,
{ 1.23844 , 1.30437 , 1.32807 , 1.61124 }

},
};


//fit parameters for data normalization, no offshell and VNA, q2 dependence in fits
//CDBonn wf, Christy Bosted F2, FSI fit
//Bonus data beam 4GeV
//first index is Q^2={3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_CDBonn_CB_off3_lc0_fsi_q2dep1_beam4[1][5][4]={ 
{ 
{ 2.4331 , 2.45185 , 2.85703 , 2.70819 }
,
{ 1.62419 , 1.66093 , 1.61621 , 1.4574 }
,
{ 1.16406 , 1.05973 , 1.10738 , 1.11556 }
,
{ 1.42127 , 1.31323 , 1.4515 , 1.50422 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

},
};


//fit parameters for data normalization, no offshell and VNA, q2 dependence in fits
//CDBonn wf, Christy Bosted F2, FSI fit
//Bonus data beam 5GeV
//first index is Q^2={3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_CDBonn_CB_off3_lc0_fsi_q2dep1_beam5[1][5][4]={ 
{ 
{ 2.12447 , 2.50773 , 2.51932 , 2.9947 }
,
{ 1.58938 , 1.57604 , 1.66616 , 1.76346 }
,
{ 1.09376 , 1.09647 , 1.08172 , 1.11513 }
,
{ 1.24461 , 1.29284 , 1.28319 , 1.47601 }
,
{ 1.17726 , 1.21934 , 1.20984 , 1.39932 }

},
};


//fit parameters for data normalization, no offshell and VNA, no q2 dependence in fits
//GrossWJC1 wf, Christy Bosted F2, FSI fit
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_GrossWJC1_CB_off3_lc0_fsi_q2dep0_beam4[3][5][4]={ 
{ 
{ 2.18192 , 2.34323 , 2.96395 , 3.49563 }
,
{ 1.35679 , 1.35308 , 1.3995 , 1.431 }
,
{ 1.2108 , 1.17335 , 1.2277 , 1.28619 }
,
{ 1.2526 , 1.24493 , 1.31138 , 1.4345 }
,
{ 1.20618 , 1.17755 , 1.24076 , 1.37901 }

},
{ 
{ 2.39557 , 2.58967 , 3.15749 , 3.90103 }
,
{ 1.48772 , 1.42658 , 1.53103 , 1.58341 }
,
{ 1.2115 , 1.19444 , 1.22363 , 1.208 }
,
{ 1.32858 , 1.2945 , 1.38111 , 1.47565 }
,
{ 1.2892 , 1.2508 , 1.36703 , 1.52401 }

},
{ 
{ 2.61802 , 2.70181 , 3.27145 , 3.29216 }
,
{ 1.65209 , 1.6984 , 1.66368 , 1.52223 }
,
{ 1.19339 , 1.09437 , 1.15543 , 1.18582 }
,
{ 1.47176 , 1.37484 , 1.54546 , 1.65241 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//fit parameters for data normalization, no offshell and VNA, no q2 dependence in fits
//GrossWJC1 wf, Christy Bosted F2, FSI fit
//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_GrossWJC1_CB_off3_lc0_fsi_q2dep0_beam5[2][5][4]={ 
{ 
{ 2.19726 , 2.55843 , 2.96801 , 3.83497 }
,
{ 1.39884 , 1.45262 , 1.48946 , 1.65184 }
,
{ 1.179 , 1.18216 , 1.19859 , 1.30011 }
,
{ 1.23942 , 1.26522 , 1.30649 , 1.45642 }
,
{ 1.13358 , 1.15386 , 1.20666 , 1.38258 }

},
{ 
{ 2.28635 , 2.75989 , 2.86564 , 3.61349 }
,
{ 1.61656 , 1.61116 , 1.71485 , 1.83816 }
,
{ 1.12091 , 1.13159 , 1.12688 , 1.18238 }
,
{ 1.28766 , 1.351 , 1.36044 , 1.60305 }
,
{ 1.23355 , 1.29701 , 1.31701 , 1.58876 }

},
};


//fit parameters for data normalization, no offshell and VNA, q2 dependence in fits
//GrossWJC1 wf, Christy Bosted F2, FSI fit
//Bonus data beam 4GeV
//first index is Q^2={3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_GrossWJC1_CB_off3_lc0_fsi_q2dep1_beam4[1][5][4]={ 
{ 
{ 2.42879 , 2.44569 , 2.84665 , 2.69119 }
,
{ 1.6232 , 1.65951 , 1.61417 , 1.45442 }
,
{ 1.16287 , 1.05821 , 1.10506 , 1.11176 }
,
{ 1.41912 , 1.31042 , 1.44699 , 1.49651 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

},
};


//fit parameters for data normalization, no offshell and VNA, q2 dependence in fits
//GrossWJC1 wf, Christy Bosted F2, FSI fit
//Bonus data beam 5GeV
//first index is Q^2={3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_GrossWJC1_CB_off3_lc0_fsi_q2dep1_beam5[1][5][4]={ 
{ 
{ 2.12074 , 2.50134 , 2.50968 , 2.97507 }
,
{ 1.58841 , 1.57469 , 1.66405 , 1.75983 }
,
{ 1.09264 , 1.09489 , 1.07945 , 1.11134 }
,
{ 1.24273 , 1.2901 , 1.27923 , 1.46863 }
,
{ 1.17484 , 1.21578 , 1.20467 , 1.38942 }

},
};


//fit parameters for data normalization, no offshell and VNA, no q2 dependence in fits
//paris wf, SLAC F2, FSI fit
//Bonus data beam 4GeV
//first index is Q^2={0.93,1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_paris_SLAC_off3_lc0_fsi_q2dep0_beam4[3][5][4]={ 
{ 
{ 4.10435 , 4.40867 , 5.0 , 5.0 }
,
{ 1.54301 , 1.53892 , 1.59211 , 1.62898 }
,
{ 1.25137 , 1.21258 , 1.26879 , 1.32965 }
,
{ 1.34061 , 1.33186 , 1.40263 , 1.53422 }
,
{ 1.33468 , 1.30259 , 1.37225 , 1.52519 }

},
{ 
{ 4.60444 , 4.97855 , 6.07388 , 7.52073 }
,
{ 1.67991 , 1.61101 , 1.72953 , 1.79006 }
,
{ 1.26485 , 1.24701 , 1.27763 , 1.26214 }
,
{ 1.35777 , 1.32277 , 1.41116 , 1.50849 }
,
{ 1.41161 , 1.36999 , 1.49791 , 1.67316 }

},
{ 
{ 5.0 , 5.0 , 5.0 , 5.0 }
,
{ 1.85853 , 1.91128 , 1.8731 , 1.71673 }
,
{ 1.30773 , 1.20016 , 1.26863 , 1.30446 }
,
{ 1.48355 , 1.38606 , 1.55919 , 1.66929 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

}
};


//fit parameters for data normalization, no offshell and VNA, no q2 dependence in fits
//paris wf, SLAC F2, FSI fit
//Bonus data beam 5GeV
//first index is Q^2={1.66,3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_paris_SLAC_off3_lc0_fsi_q2dep0_beam5[2][5][4]={ 
{ 
{ 4.21435 , 4.90756 , 5.6958 , 7.37344 }
,
{ 1.58394 , 1.64503 , 1.68713 , 1.8722 }
,
{ 1.23915 , 1.2425 , 1.26001 , 1.36758 }
,
{ 1.27962 , 1.30612 , 1.34893 , 1.50436 }
,
{ 1.24978 , 1.27209 , 1.33066 , 1.52623 }

},
{ 
{ 5.0 , 5.0 , 5.0 , 5.0 }
,
{ 1.8259 , 1.82026 , 1.93802 , 2.07953 }
,
{ 1.21343 , 1.22555 , 1.22163 , 1.2836 }
,
{ 1.2855 , 1.34931 , 1.35972 , 1.60519 }
,
{ 1.347 , 1.41597 , 1.43842 , 1.73883 }

},
};


//fit parameters for data normalization, no offshell and VNA, q2 dependence in fits
//paris wf, SLAC F2, FSI fit
//Bonus data beam 4GeV
//first index is Q^2={3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_paris_SLAC_off3_lc0_fsi_q2dep1_beam4[1][5][4]={ 
{ 
{ 5.0 , 5.0 , 5.0 , 5.0 }
,
{ 1.8264 , 1.86791 , 1.81769 , 1.6402 }
,
{ 1.27457 , 1.16071 , 1.21336 , 1.22243 }
,
{ 1.43085 , 1.32136 , 1.45979 , 1.51064 }
,
{ 1.0 , 1.0 , 1.0 , 1.0 }

},
};


//fit parameters for data normalization, no offshell and VNA, q2 dependence in fits
//paris wf, SLAC F2, FSI fit
//Bonus data beam 5GeV
//first index is Q^2={3.38}
//second index is W={1.17,1.48,1.73,2.03,2.44}
//third index is p_s[MeV]={78,93,110,135}
const double bonusfits2::normfits_paris_SLAC_off3_lc0_fsi_q2dep1_beam5[1][5][4]={ 
{ 
{ 5.0 , 5.0 , 5.0 , 5.0 }
,
{ 1.79447 , 1.77943 , 1.88095 , 1.9909 }
,
{ 1.18308 , 1.18603 , 1.17031 , 1.20607 }
,
{ 1.24098 , 1.28875 , 1.27857 , 1.46975 }
,
{ 1.2833 , 1.32754 , 1.31555 , 1.51841 }

},
};

#endif