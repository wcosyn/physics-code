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
add_subdirectory(neutrinoQE) #set of programs to compute neutrino QE double differential and single differential cross sections, relativistic mean field, relativistic fermi gas, nucleon and nucleus
add_subdirectory(momdistr) #small program that computed momentum distributions with relativistic MF nuclear wf
add_subdirectory(LDRD_test) #small test program for the LDRD project, tagging with He3, never properly used
add_subdirectory(integrator) #small test program for integration routines, calculates volume of part of a sphere
add_subdirectory(incldis) #contains programs that were used to compute fsi in inclusive deuteron dis [unpolarized and Azz]
add_subdirectory(gammatest) #test program for gamma matrix identities and to check some trace calculations explicitly
add_subdirectory(DQEinc) #inclusive deuteron QE FSI calculations, work in progress

add_subdirectory(cuba_test) #test prog for the cuba library
add_subdirectory(camille) #some test programs that camille used, don't know the details...