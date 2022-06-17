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


vector<double> GPD_V_Nucl_grid::getVectorGPDSet(const double x, const double xi, const double t, const double scale, const bool gpdvector){
    if(xi<0) return getVectorGPDSet(x,-xi,t,scale, gpdvector);  //both GPDs are even
    //make a grid in x,xi since the integrals to compute the chiral odd gpds take some time, t is normally constant for a computation
    if(t!=t_grid||grid_set==false||gpdvector != grid_vector){
        //cout << "constructing chiral even gpd grid " << t_grid << " " << t << " " << grid_set << " " << grid_vector << endl;
        for(int i=0;i<=200;i++){
            for(int j=0;j<=100;j++){
                double x = 0.01*(i-100)+(i==100?1.E-04:0.);
                PARTONS::GPDKinematic gpdKinematic(x,0.01*(j)+(j==0?1.E-04:0.),t*1.E-06, scale*scale, scale*scale);
                PARTONS::GPDResult gpdResult = pGPDService->computeSingleKinematic(gpdKinematic,
                    pGPDModel);
                double H=0.,E=0.;
                if(gpdvector){
                    H=0.5*(gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                        +gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution());
                    E=0.5*(gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                        +gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution());
                }
                else{
                    H=0.5*(gpdResult.getPartonDistribution(PARTONS::GPDType::Ht).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                        +gpdResult.getPartonDistribution(PARTONS::GPDType::Ht).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution());
                    E=0.5*(gpdResult.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                        +gpdResult.getPartonDistribution(PARTONS::GPDType::Et).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution());
                }
                //if(isnan(H)||isnan(E)) cout << x << " " << 0.01*j << endl;
                grid[i][j][0]=H;
                grid[i][j][1]=E;
                if(isnan(H)||isnan(E)) grid[i][j][0]=grid[i][j][1]=0.;
                // grid[i][j][2]=0.5*(gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                //     -gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution());
                // grid[i][j][3]=0.5*(gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution()
                //     -gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution());
                // grid[i][j][0]=gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
                // grid[i][j][1]=gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::UP).getQuarkDistribution();
                // grid[i][j][2]=gpdResult.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
                // grid[i][j][3]=gpdResult.getPartonDistribution(PARTONS::GPDType::E).getQuarkDistribution(PARTONS::QuarkFlavor::DOWN).getQuarkDistribution();
                
                //cout << x << " " << 0.01*j << " " << grid[i][j][0] << " " << grid[i][j][1] << endl;
            }
        }
        //test nucleon

        // for(int j=0;j<=20;j++){
        //     double sumEs=0,sumHs=0;
        //     // sumEv=0., sumHv=0.;
        //     for(int i=1;i<199;i++) {sumHs+=grid[i][j][0]; sumEs+=grid[i][j][1];/* sumHv+=grid[i][j][2]; sumEv+=grid[i][j][3];*/}
        //     cout << j*0.01 << " " /* << (sumHs+3.*sumHv)/3.*0.01 << " " << (sumEs+3.*sumEv)/3.*0.01 << " " << (sumHs-3.*sumHv)/3.*0.01 << " " << (sumEs-3.*sumEv)/3.*0.01 << 
        //     " "*/ << sumHs/3.*0.01 << " " << sumEs/3.*0.01 /* << " " << sumHv*0.01 << " " << sumEv*0.01*/ <<  endl;
        //     // double sumEu=0,sumHu=0, sumEd=0., sumHd=0.;
        //     // for(int i=0;i<200;i++) {sumHu+=grid[i][j][0]; sumEu+=grid[i][j][1];sumHd+=grid[i][j][2]; sumEd+=grid[i][j][3];}
        //     // cout << j << " " << (2.*sumHu-sumHd)/3.*0.01 << " " << (2.*sumEu-sumEd)/3.*0.01 << " " << (2.*sumHd-sumHu)/3.*0.01 << " " << (2.*sumEd-sumEu)/3.*0.01 << endl;
        //  }
        //  cout << endl;
        // // cout << "construction done" << endl;
        grid_set=true;
        t_grid=t;
        grid_vector=gpdvector;
   }
   //interpolation
    double index_i=0.,index_j=0.;
    double frac_i=0.,frac_j=0.;
    frac_i=modf(x*100+100,&index_i);
    frac_j=modf(xi*100,&index_j);
    if(int(index_j)==0.) {index_j=1.; frac_j-=1.;} //xi=0 gives nan values from partons usually


    vector<double> result=vector<double>(2,0.);
    result[0]=grid[int(index_i)][int(index_j)][0]*(1.-frac_i)*(1.-frac_j)
                                        +grid[int(index_i)+1][int(index_j)][0]*(frac_i)*(1.-frac_j)
                                        +grid[int(index_i)][int(index_j)+1][0]*(1.-frac_i)*(frac_j)
                                        +grid[int(index_i)+1][int(index_j)+1][0]*(frac_i)*(frac_j);
    result[1]=grid[int(index_i)][int(index_j)][1]*(1.-frac_i)*(1.-frac_j)
                                        +grid[int(index_i)+1][int(index_j)][1]*(frac_i)*(1.-frac_j)
                                        +grid[int(index_i)][int(index_j)+1][1]*(1.-frac_i)*(frac_j)
                                        +grid[int(index_i)+1][int(index_j)+1][1]*(frac_i)*(frac_j);
                                    
                        

    
    // TransGPD_set gpd_nucl=getGK_param(x,xi,t);

    // cout << "gpd " << gpd_nucl_grid.getHTu() << " " << gpd_nucl.getHTu() << " " << gpd_nucl_grid.getHTd() << " " << gpd_nucl.getHTd() << " "
    //     << gpd_nucl_grid.getEbarTd() << " " << gpd_nucl.getEbarTd() << " " << gpd_nucl_grid.getEbarTu() << " " << gpd_nucl.getEbarTu() << endl;
    return result;
}

