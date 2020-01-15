#ifndef GPD_V_NUCL_GRID
#define GPD_V_NUCL_GRID

#include <vector>
#include <cmath>
#include <iostream>
#include <cstdarg>
#include <cstdlib>
#include <complex>
#include <numint/numint.hpp>

#include <partons/ModuleObjectFactory.h>
#include <partons/services/GPDService.h>
#include <partons/services/ObservableService.h>
#include <partons/ServiceObjectRegistry.h>




/**
 * @brief contains a grid with isoscalar chiral even nucleon GPDs obtained from partons
 * 
 */
class GPD_V_Nucl_grid{

public:
/**
 * @brief Construct a new gpd_V nucl grid object
 * 
 * @param 
 */
GPD_V_Nucl_grid(PARTONS::GPDService *pGPDService, PARTONS::GPDModule *pGPDModel);

~GPD_V_Nucl_grid(){;}



/**
 * @brief Obtains the isoscalar chiral even nucleon gpds for a certain kinematics (interpolated from a grid)
 * 
 * @param x average lf momentum fraction quark
 * @param xi skewness nucleon
 * @param t [MeV^2] momentum transfer sq.
 * @param scale [GeV] factorization = renormalization scale
 * @return [0] H, [1] E
 */
std::vector<double> getVectorGPDSet(const double x, const double xi, const double t, const double scale);


private:

double t_grid; ///< [MeV^2] momentum transfer sq value the grid has
bool grid_set; ///< is the grid with transv gpds set or not
bool ERBL_set; ///< is the grid only ERBL or not
double grid[201][101][2];

PARTONS::GPDService* pGPDService;
PARTONS::GPDModule* pGPDModel;

};


#endif