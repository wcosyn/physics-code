project(spectagdis)

#program that calculates deeps-like reduced structure functions
add_executable(deepstest deepstest.cpp )  
target_link_libraries(deepstest DIS MePhys)

#program used to fit the FSI parameters to the deeps data
add_executable(deepsfit deepsfit.cpp )  
target_link_libraries(deepsfit DIS MePhys)

#program used output bonus-like ratios (but output is nominator and denom separately).
# Used for the plots in the bonus paper.  Uses renormalization of the data with the bonus MC results.
# this is the one to use [not bonusetrapolate(2)]
add_executable(bonustest bonustest.cpp )  
target_link_libraries(bonustest DIS MePhys)

#program used to refit rescattering parameters but no sensible results (too low spectator momenta) 
add_executable(bonusscatterfit bonusscatterfit.cpp )  
target_link_libraries(bonusscatterfit DIS MePhys)

#small program that wrote out some of the data renormalization fit values
add_executable(bonusfitreadout bonusfitreadout.cpp )  
target_link_libraries(bonusfitreadout DIS MePhys)

#small program that wrote out some of the data renormalization fit values
add_executable(bonusfit_rig bonusfit_rig.cpp )  
target_link_libraries(bonusfit_rig DIS MePhys)

#calculates bonus like F2n ratios, used in the python script fitscript.py for the data points
add_executable(bonusextrapolate bonusextrapolate.cpp )  
target_link_libraries(bonusextrapolate DIS MePhys)

#calculates bonus like F2n ratios, used in the python script fitscript.py for the extrapolation procedure
add_executable(bonusextrapolate2 bonusextrapolate2.cpp )  
target_link_libraries(bonusextrapolate2 DIS MePhys)

#program that dit the pole extrapolation for the bonus data (renormalized and al), also independent fit of all angles which is the same as below 
add_executable(bonus_F2nfit bonus_F2nfit.cpp )  
target_link_libraries(bonus_F2nfit DIS MePhys)

#program that dit the pole extrapolation for the bonus data (renormalized and al), this was used to compute the systematic errors
add_executable(bonus_F2nfit_sys bonus_F2nfit_sys.cpp )  
target_link_libraries(bonus_F2nfit_sys DIS MePhys)

#program that dit the pole extrapolation for the bonus data (renormalized and al), by fitting all angles at the same time, NOT USED no good results
add_executable(bonus_F2nfit_allangles bonus_F2nfit_allangles.cpp )  
target_link_libraries(bonus_F2nfit_allangles DIS MePhys)

#small program that checked hermes b1 kinematics and possible problems with mixing of structure functions
add_executable(hermes_comp hermes_comp.cpp )  
target_link_libraries(hermes_comp DIS MePhys)

#program used to compute the convolution b1 VNA, also has the comparison with Shunzo Kumano's model
add_executable(b1_calc b1_calc.cpp )  
target_link_libraries(b1_calc DIS MePhys)

#program used to compute the convolution b1 VNA, also has the comparison with Shunzo Kumano's model
add_executable(b1_calc_comp b1_calc_comp.cpp )  
target_link_libraries(b1_calc_comp DIS MePhys)

add_executable(band_calc band_calc.cpp )  
target_link_libraries(band_calc DIS MePhys)

if(NOT BUILD_PROGS)
set_target_properties(deepstest PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(deepsfit PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(bonustest PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(bonusscatterfit PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(bonusfitreadout PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(bonusfit_rig PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(bonusextrapolate PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(bonusextrapolate2 PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(bonus_F2nfit PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(bonus_F2nfit_sys PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(bonus_F2nfit_allangles PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(hermes_comp PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(b1_calc PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(b1_calc_comp PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(band_calc PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS band_calc deepstest deepsfit bonustest bonusscatterfit bonusfitreadout bonusfit_rig hermes_comp 
            bonusextrapolate bonusextrapolate2 bonus_F2nfit bonus_F2nfit_sys bonus_F2nfit_allangles b1_calc b1_calc_comp  DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()