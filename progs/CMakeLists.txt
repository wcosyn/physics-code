project(progs)

add_subdirectory(Obstest) #prog used to compute A(e,e'p) observables
add_subdirectory(twonucl) #old code that was work in progress on a 2n knockout model w correlated pairs
# add_subdirectory(threeres) #old code that was never properly used, resonance model for incl DIS
add_subdirectory(onenucl)  #programs used to compute A(e,e'p) cross sections and transparencies
add_subdirectory(SCX_test) # test program to check charge exchange densities etc. Camille wrote it, doesn't seem to be running fine on my system (roundoff errors in gsl)
add_subdirectory(RhoTransparency)  #several programs used to compute rho transparencies
add_subdirectory(poldeut) #program for testing lightfront deuteron spectator tagging code
add_subdirectory(pionprog) #program that wrote out pion glauber parameters
add_subdirectory(pionestimate) #program that estimated uncertainties caused by rho decay to 2 pions in rho transparency calculations