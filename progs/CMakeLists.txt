project(progs)

add_subdirectory(Obstest) #prog used to compute A(e,e'p) observables
add_subdirectory(twonucl) #old code that was work in progress on a 2n knockout model w correlated pairs
# add_subdirectory(threeres) #old code that was never properly used, resonance model for incl DIS
add_subdirectory(T_for_Or)  #program used to do some transparency calculations for a jlab proposal
add_subdirectory(SCX_test) # test program to check charge exchange densities etc. Camille wrote it, doesn't seem to be running fine on my system (roundoff errors in gsl)