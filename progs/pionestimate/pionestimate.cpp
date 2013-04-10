#include <iostream>
#include <cstdlib>

using namespace std;

#include <MeanFieldNucleusThick.hpp>
#include <FsiCorrelator.hpp>
#include <FastParticle.hpp>
#include <constants.hpp>
#include <GlauberGrid.hpp>
#include <OneGlauberGrid.hpp>
#include <GlauberGridThick.hpp>
#include <TKinematics2to2.h>
#include <TElectronKinematics.h>
#include <Cross.hpp>
#include <TRotation.h>
#include <DistMomDistrGrid.hpp>

//run ./pionest [nucleus] [one doublepion=0, two pions=1] 
//estimates effect of taking collinear pions in rho transparency calculations


void intpm(numint::vector_d & res, double pm, double costheta, double phi, DistMomDistrGrid **pgrid, int totshells){
  res=numint::vector_d(1,0.);
  for(int i=0;i<totshells;i++) res[0]+=pgrid[i]->getRhoGridFull_interp3(0, pm, costheta, phi);
}

void intpm_pw(numint::vector_d & res, double pm, double costheta, double phi, DistMomDistrGrid **pgrid, int totshells){
  res=numint::vector_d(1,0.);
  for(int i=0;i<totshells;i++) res[0]+=pgrid[i]->getRhopwGridFull_interp(pm);
}


int main(int argc, char *argv[])
{
  int nucleus=atoi(argv[1]);
  double prec=(nucleus==1? 1.E-05: 1.E-06);//atof(argv[7]);
  int integrator=2;//atoi(argv[8]);
  int thick=1;//atoi(argv[9]);
  int maxEval=20000;//atoi(argv[10]);
  bool ppion=atoi(argv[2]);
  
  string homedir=argv[3];

  MeanFieldNucleusThick nucl(nucleus,homedir);
  DistMomDistrGrid *pgrid[nucl.getTotalLevels()];
  double prho=3E03;
  double ppi=sqrt((MASSRHO*MASSRHO+prho*prho)/4.-MASSPI*MASSPI);
  double ppiz=prho/2.;
  double thetapi=2.*acos(ppiz/ppi);
  cout << thetapi*RADTODEGR << endl;
  
  
  if(ppion){
    GlauberGridThick grid = GlauberGridThick(60,18,5,&nucl,prec,integrator,homedir);
    FastParticle pion1(2, 0, ppi,0.,0.,3.,0.,homedir);
    FastParticle pion2(3, 0, ppi,thetapi,0.,3.,0.,homedir);
    grid.clearParticles();
    if(ppion==1||ppion==2) grid.addParticle(pion1);
    if(ppion==1||ppion==3) grid.addParticle(pion2);
    grid.updateGrids();
    grid.clearKnockout();
    for(int i=0;i<nucl.getTotalLevels();i++){
      pgrid[i] = new DistMomDistrGrid(i, 300., 30,20,5,&grid,1.E-03,2,5E04,0.,homedir);
      TRotation rot;
      rot.Rotate(0.,TVector3(0.,1.,0.));
      pgrid[i]->updateGrids(&grid,i,rot);
      pgrid[i]->printRho_grid(0);
    }
  }
  else{
    OneGlauberGrid grid = OneGlauberGrid(60,18,&nucl,prec,integrator,homedir);
    FastParticle rhopion(7, 0, ppi,0.,0.,3.,0.,homedir);
    grid.clearParticles();
    grid.addParticle(rhopion);
    grid.updateGrids();
    grid.clearKnockout();
    for(int i=0;i<nucl.getTotalLevels();i++){
      pgrid[i] = new DistMomDistrGrid(i, 300., 30,20,5,&grid,1.E-03,2,5E04,0.,homedir);
      TRotation rot;
      rot.Rotate(0.,TVector3(0.,1.,0.));
      pgrid[i]->updateGrids(&grid,i,rot);
      pgrid[i]->printRho_grid(0);
    }
  }
  
    struct Ftor_grid {

      /*! integrandum function */
      static void exec(const numint::array<double,3> &x, void *param, numint::vector_d &ret) {
	Ftor_grid &p = * (Ftor_grid *) param;
	p.f(ret,x[0],x[1],x[2],p.pgrid,p.shells);
      }
      DistMomDistrGrid **pgrid;/*!< pointer to the grid where the integration is performed */
      int shells;
      /*! integrandum 
      * \param res results
      * \param pm first integration variable
      * \param costheta second integration variable
      * \param phi third integration variable
      * \param pgrid the DistMomDistrGrid instance
      */
      void (*f)(numint::vector_d & res, double pm, double costheta, double phi, DistMomDistrGrid **pgrid, int totshells);
    };
    int res=90;
    unsigned count=0;
    numint::vector_d ret(1,0.);
    numint::array<double,3> lower = {{0.,-1.,0.}};
    numint::array<double,3> upper = {{300.,1.,2.*PI}};
    
    Ftor_grid F;
    F.pgrid = pgrid;
    F.shells=nucl.getTotalLevels();
    numint::mdfunction<numint::vector_d,3> mdf;
    mdf.func = &Ftor_grid::exec;
    mdf.param = &F;
    F.f=intpm_pw;
    res = numint::cube_adaptive(mdf,lower,upper,1.E-12,1.E-03,20000,ret,count,0);
    cout << res << " " << count << endl;
    numint::vector_d ret2(1,0.);
    F.f=intpm;
    res = numint::cube_adaptive(mdf,lower,upper,1.E-12,1.E-03,20000,ret2,count,0);
    cout << res << " " << count << endl << endl;
    
    cout << ret2[0]/ret[0] << endl;
  

  for(int i=0;i<nucl.getTotalLevels();i++) delete pgrid[i];
  
}