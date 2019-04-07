#include "GPD_V_Nucl_grid.hpp"

#include <constants.hpp>
#include <cassert>
#include <numint/numint.hpp>


using namespace std;

GPD_V_Nucl_grid::GPD_V_Nucl_grid(PARTONS::GPDService *pGPDService, PARTONS::GPDModule *pGPDModel):
pGPDService(pGPDService),pGPDModel(pGPDModel){
    grid_set=false;
    t_grid=0.;
    ERBL_set=false;

}


vector<double> GPD_V_Nucl_grid::getVectorGPDSet(const double x, const double xi, const double t, const bool ERBL){
    if(xi<0) return getVectorGPDSet(x,-xi,t,ERBL);  //both GPDs are even
    //make a grid in x,xi since the integrals to compute the chiral odd gpds take some time, t is normally constant for a computation
    if(t!=t_grid||grid_set==false||ERBL!=ERBL_set){
        cout << "constructing chiral even gpd grid" << endl;
        for(int i=0;i<=200;i++){
            for(int j=0;j<=100;j++){
                PARTONS::GPDKinematic gpdKinematic(0.01*(i-100)*(ERBL? abs(xi): 1.),0.01*(j),t, 1., 1.);
                PARTONS::GPDResult gpdResult = pGPDService->computeGPDModel(gpdKinematic,
                    pGPDModel);
                double H=0.5*(gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                    +gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution());
                double E=0.5*(gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                    +gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution());

                grid[i][j][0]=H;
                grid[i][j][1]=E;
            }
        }
        cout << "construction done" << endl;
        grid_set=true;
        t_grid=t;
        ERBL_set=ERBL;
   }
   //interpolation
    double index_i=0.,index_j=0.;
    double frac_i=0.,frac_j=0.;
    frac_i=modf(x*100/(ERBL? abs(xi): 1.)+100,&index_i);
    frac_j=modf(xi*100,&index_j);

    vector<double> result=vector<double>(2,0.);
    result[0]=grid[int(index_i)][int(index_j)][0]*(1.-frac_i)*(1.-frac_j)
                                        +grid[int(index_i)+1][int(index_j)][0]*(frac_i)*(1.-frac_j)
                                        +grid[int(index_i)][int(index_j)+1][0]*(1.-frac_i)*(frac_j)
                                        +grid[int(index_i)+1][int(index_j)+1][0]*(frac_i)*(frac_j);
    result[1]=grid[int(index_i)][int(index_j)][0]*(1.-frac_i)*(1.-frac_j)
                                        +grid[int(index_i)+1][int(index_j)][1]*(frac_i)*(1.-frac_j)
                                        +grid[int(index_i)][int(index_j)+1][1]*(1.-frac_i)*(frac_j)
                                        +grid[int(index_i)+1][int(index_j)+1][1]*(frac_i)*(frac_j);
                                    
                        


    // TransGPD_set gpd_nucl=getGK_param(x,xi,t);

    // cout << "gpd " << gpd_nucl_grid.getHTu() << " " << gpd_nucl.getHTu() << " " << gpd_nucl_grid.getHTd() << " " << gpd_nucl.getHTd() << " "
    //     << gpd_nucl_grid.getEbarTd() << " " << gpd_nucl.getEbarTd() << " " << gpd_nucl_grid.getEbarTu() << " " << gpd_nucl.getEbarTu() << endl;
    return result;
}

