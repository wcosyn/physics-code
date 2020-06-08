#include "GPD_T_Nucl_grid.hpp"
#include "TransGPD_set.hpp"

#include <constants.hpp>
#include <mstwpdf.h>
#include <cassert>
#include <numint/numint.hpp>


using namespace std;

GPD_T_Nucl_grid::GPD_T_Nucl_grid(const std::string &pdf_name){
   if(!pdf_name.compare("MSTW")){
        string file= string(HOMEDIR)+"/mstw_grids/mstw2008lo.00.dat";
        mstw = new c_mstwpdf(file,false,true);
    }
    else {
        mstw=NULL;
        cerr << "you have not chosen a valid pdf set in Deut_Conv_GPD_T::GPD()" << endl;
        assert(1==0);
        exit(1);
    }
    grid_set=false;
    t_grid=0.;
    ERBL_set=false;

}

GPD_T_Nucl_grid::~GPD_T_Nucl_grid(){
    if(mstw) delete mstw;
}

TransGPD_set GPD_T_Nucl_grid::getTransGPDSet(const double x, const double xi, const double t, const double scale){
    if(xi<0) return getTransGPDSet(x,-xi,t,scale);  //Etilde is always zero in the models so we have only even GPDs in xi
    //make a grid in x,xi since the integrals to compute the chiral odd gpds take some time, t is normally constant for a computation
    if(t!=t_grid||grid_set==false){
        //cout << "constructing chiral odd gpd grid " << t_grid << " " << t << " " << xi << " " << grid_set << endl;
        for(int i=0;i<=200;i++){
            for(int j=int(xi*100);j<=int(xi*100)+1;j++){
                grid[i][j]=getGK_param(0.01*(i-100),0.01*(j),t, scale);
        //         cout << 0.01*(i-100) << " " << 0.01*(j) << " " << grid[i][j].getHTu() << " " << grid[i][j].getHTd() << " "
        // << grid[i][j].getEbarTu() <<  " " << grid[i][j].getEbarTd() << endl;
            }
        }
        // cout << "construction done" << endl;
        grid_set=true;
        t_grid=t;
   }
   //interpolation
    double index_i=0.,index_j=0.;
    double frac_i=0.,frac_j=0.;
    frac_i=modf(x*100+100,&index_i);
    frac_j=modf(xi*100,&index_j);

    TransGPD_set gpd_nucl_grid=grid[int(index_i)][int(index_j)]*(1.-frac_i)*(1.-frac_j)+grid[int(index_i)+1][int(index_j)]*(frac_i)*(1.-frac_j)
                                +grid[int(index_i)][int(index_j)+1]*(1.-frac_i)*(frac_j)+grid[int(index_i)+1][int(index_j)+1]*(frac_i)*(frac_j);


    // TransGPD_set gpd_nucl=getGK_param(x,xi,t);

    // cout << "gpd " << gpd_nucl_grid.getHTu() << " " << gpd_nucl.getHTu() << " " << gpd_nucl_grid.getHTd() << " " << gpd_nucl.getHTd() << " "
    //     << gpd_nucl_grid.getEbarTd() << " " << gpd_nucl.getEbarTd() << " " << gpd_nucl_grid.getEbarTu() << " " << gpd_nucl.getEbarTu() << endl;
    return gpd_nucl_grid;
}


TransGPD_set GPD_T_Nucl_grid::getGK_param(const double x, const double xi, const double t, const double scale){

    if(fabs(x)>1.) return TransGPD_set(0.,0.,0.,0.);
    GPD_T_Nucl_grid::Ftor_doubledistr F;
    F.x=x;
    F.xi=xi;
    F.t=t;
    F.scale=scale;
    F.gpdobject=this;
    numint::mdfunction<numint::vector_d,1> mdf;
    mdf.func = &Ftor_doubledistr::exec;
    mdf.param = &F;
    numint::vector_d ret(4,0.);
    F.f=GPD_T_Nucl_grid::DD_int_rho;
    double low=0.,hi=1.;
    if(xi>=0.){
        low=max(0.,(x-xi)/(1.-xi));
        hi=min(1.,(x+xi)/(1.+xi));
    }
    else{
        low=max(0.,(x+xi)/(1.+xi));
        hi=min(1.,(x-xi)/(1.-xi));
    }

    if(low>hi) return TransGPD_set(0.,0.,0.,0.);

    //  cout << "bounds " << x << " " << low << " " << hi << endl;
    if(low<1.E-09) low=1.E-09;
    if((1.-hi)<1.E-05) hi=1.-1.E-05;
    numint::array<double,1> lower = {{low}};
    numint::array<double,1> upper = {{hi}};

    unsigned count;
    numint::cube_adaptive(mdf,lower,upper,1.E-08,1.E-03,1E02,6E04,ret,count,0);
    // cout << "c " << ret[0] << " " << ret[1] << " " << ret[2] << " " << ret[3] << " " << count << endl;
    return TransGPD_set(ret[0],ret[1],ret[2],ret[3]);


}


void GPD_T_Nucl_grid::DD_int_rho(numint::vector_d & res, double rho, double x, double xi, double t, double scale, GPD_T_Nucl_grid &gpdobject){

    res=numint::vector_d(4,0.);
    double eta=(x-rho)/xi;
    double temp=3./4.*((1.-rho)*(1.-rho)-eta*eta)/pow(1.-rho,3.)/abs(xi);
    double HTdfront,HTufront, EbarTdfront, EbarTufront;
    gpdobject.getHTfront(rho,HTdfront,HTufront, scale);
    gpdobject.getEbarTfront(rho,EbarTdfront, EbarTufront);

    double exp1=exp(-0.45*log(rho)*t*1.E-06);
    double exp2=exp1*exp(0.5*t*1.E-06);
    // cout << rho << " " << x<< " "<< xi << " " << t << " " << exp1 << " " << exp2 << " " << temp << " " << HTdfront << " " << HTufront << " " << EbarTdfront << " " << EbarTufront << endl;
    res[0]=exp1*temp*HTdfront;
    res[1]=exp1*temp*HTufront;
    res[2]=exp2*temp*EbarTdfront;
    res[3]=exp2*temp*EbarTufront;
    // cout << rho << " " << x << " " << xi << " " << t << " " << res[0] << " " << res[1] << " " << res[2] << " " << res[3] << endl;

}

void GPD_T_Nucl_grid::getHTfront(const double x, double &HTdfront, double &HTufront, const double scale) const{
    mstw->update(x,2.);


    double d=mstw->parton(1,x,scale)/x;
    double dbar=mstw->parton(-1,x,scale)/x;
    double u=mstw->parton(2,x,scale)/x;
    double ubar=mstw->parton(-2,x,scale)/x;

    double Dd=-0.7*sqrt(x)*d;
    double Ddbar=-0.3*pow(x,0.4)*dbar;
    double Du=Dd/d/-0.7*u;
    double Dubar=Ddbar/dbar*ubar;


    HTdfront=-1.01*sqrt(x)*(1.-x)*(d-dbar+Dd-Ddbar);
    HTufront=0.78*sqrt(x)*(1.-x)*(u-ubar+Du-Dubar);
    return;
}

void GPD_T_Nucl_grid::getEbarTfront(const double x, double &EbarTdfront, double &EbarTufront) const{
    EbarTdfront=5.05*pow(x,-0.3)*pow(1.-x,5.);
    EbarTufront=EbarTdfront/(1.-x)/5.05*6.83;
    return;
}
