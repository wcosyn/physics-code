#include <iostream>
using std::cout;
using std::endl;

#include <TMFSpinor.hpp>
#include <TSpinor.h>
#include <FourVector.h>
#include <constants.hpp>
#include <MeanFieldNucleus.hpp>
#include <string>
using std::string;
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <fstream>
using std::ofstream;
#include <numint/numint.hpp>
#include <complex>
using std::complex;
#include <iomanip>
using std::setw;
#include <vector>
using std::vector;
using std::pair;
#include <fstream>
using std::ofstream;
#include <cmath>
#include <sstream>
#include <cassert>

#define FORCE_POS_ENERGY_PROJ 1 // to check WIM's previous results set to True, otherwise set false

void printArray(MeanFieldNucleus&, const int*);
void printNucleusInfo(MeanFieldNucleus&);
void integrandum(complex<double> &, const numint::array<double,3>&, TVector3&, TSpinor*, TSpinor*, MeanFieldNucleus*, int shellindex1,int m1, int shellindex2, int m2);


void genPcmList(TVector3 [] , int nsteps); // pass a TVector pointer, will contain list, second argument is number of steps in each perpendicular direction, so length of list is 3*nsteps

double dist2bodymom(MeanFieldNucleus&, int,int,TVector3&,TVector3&, TVector3& );
//double dist2bodymom(double,int,int,TVector3& P);

struct F {
	static void exec(const numint::array<double,3> &x, void *param, complex<double> &ret){
		F &p = * (F *) param; // cast void pointer to an struct F pointer
		p.f(ret,x,p.P,p.u_pm,p.u_ps,p.nuc,p.shellindex1,p.m1,p.shellindex2,p.m2);
	}
	
	MeanFieldNucleus* nuc;
	TSpinor *u_pm,*u_ps;
	int shellindex1,m1,shellindex2,m2;
	TVector3 P;

	void (*f)(complex<double>& , const numint::array<double,3>&, TVector3& P,TSpinor*, TSpinor*, MeanFieldNucleus*, int, int,int,int);
};


void note_kinematics(MeanFieldNucleus&); // kinematics as in my notes
void symmetric_kinematics(MeanFieldNucleus&); // symmetric kinematics like WIM pg 65


// some util function for my 2bodymomdistributionmatrix
vector< pair<TVector3,double> >** initRho(unsigned n){ // init function -> pass rho as dynamically sized 2D array
       	vector< pair<TVector3,double> > **rho;
	rho = new vector< pair<TVector3,double> >*[n];
	for (unsigned i=0;i<n;i++){ rho[i]=new vector< pair<TVector3,double> >[n]; for (unsigned j=0; j<n; j++) rho[i][j] = vector< pair<TVector3,double> >();}
	return rho;
}	
void delRho(  vector<pair<TVector3,double> >** rho, unsigned n){ // del funtion -> pass rho as dynamically sized 2D array	
	for (unsigned i=0;i<n;i++){ delete[] rho[i];}
	delete[] rho;
}
void writeRho(string filename, MeanFieldNucleus& nuc, vector< pair<TVector3,double> >** rho); // write my vector pair matrix shizzle to a file

// util print function

string vecString(TVector3& v)
{
	std::stringstream ss;
	ss << "[" << v.X() << ", " << v.Y() << ", " << v.Z() << "]";
	return ss.str();
}


int main()
{
	string homedir = "/home/camille/Code/share";
	// -- MeanFieldNucleus
	MeanFieldNucleus nuc = MeanFieldNucleus(3,homedir);
	if (FORCE_POS_ENERGY_PROJ) { cout << "Forcing positive energy contributions in spinors!!!" << endl; }
	symmetric_kinematics(nuc);	
	return 0;
}

void symmetric_kinematics(MeanFieldNucleus& nuc)
{
	unsigned nsteps = 400;
	double PcmRange = 600; // [-...,+...] MeV
	TVector3 Pcm[nsteps];
	for (unsigned i=0; i<nsteps; i++)
	{
		Pcm[i] = TVector3(0.,0.,2.*PcmRange/(nsteps-1)*i-PcmRange); // zrange from -Range to Range
	}
	TVector3 q = TVector3(0.,0.,3000.); // q along z
	vector< pair<TVector3, double> > **rho = initRho(nuc.getPLevels()); //[nuc.getPLevels()][nuc.getPLevels()];
	for (int shellindex1 = 0; shellindex1 < nuc.getPLevels() ; shellindex1++)
	{
		cout << shellindex1+1 << " to " << nuc.getPLevels() << endl;
		for (int shellindex2 = shellindex1; shellindex2 < nuc.getPLevels() ; shellindex2++)
		{
			cout <<">" << shellindex2+1 << " to " << nuc.getPLevels() << endl;
			for (unsigned i=0; i<nsteps;i++)
			{
				double EP = sqrt(4.*MASSP*MASSP+Pcm[i]*Pcm[i]); // energy in COM
				double kz = (Pcm[i].Z()+q.Z())/2.; // because symmetric kinematics: kz is P/2+q, -- Z compos, everything is z --
				double K = sqrt(pow((q.Mag()+EP)/2.,2)-MASSP*MASSP); // kinetic energy of ej nucleons  = (q+P)/2 - rest mass (total com energy at coll = P+q)
				double kperp = sqrt(K*K-kz*kz); // perp component of k of ej nucleons
				double theta = atan2(kperp,kz); // theta of ej nucleons
				TVector3 pf(0.,0.,K), ps(0.,0.,K); // create vectors of ej nucleons
				pf.RotateY(theta); ps.RotateY(-theta); //now rotate them into xz plane
				TVector3 pm = pf-q; // calculate pm
				//TVector3 pm = Pcm[i]-ps;
				cout << i << " to " << nsteps << endl;
				//cout << "P " << vecString(Pcm[i]) << " pm " << vecString(pm) << " ps " << vecString(ps) << endl;	
				double res = dist2bodymom(nuc,shellindex1,shellindex2,Pcm[i],pm,ps);
				rho[shellindex1][shellindex2].push_back(pair<TVector3,double>(Pcm[i],res));
			}
		}
		cout << endl;
	}
	if (FORCE_POS_ENERGY_PROJ) writeRho("rho_symm_posE.dat",nuc,rho);
	else writeRho("rho_symm.dat",nuc,rho);
	delRho(rho,nuc.getPLevels());
}	

void note_kinematics(MeanFieldNucleus& nuc)
{

	// create Pcm list
	//TVector3 *Pcm = NULL; // explicit NULL initialization, this is good common practise for creating empty pointers!
	int nsteps = 100;
	TVector3 Pcm[3*nsteps];
	genPcmList(Pcm,nsteps);
	// -- Vectors
	TVector3 q = TVector3(0.,0.,1400.); // see distribution of Maarten
	TVector3 pf = 0.8*q;
	pf.RotateY(DEGRTORAD*5); // pf now in xz plane
	TVector3 pm = pf-q; // pm also in xz plane
	vector< pair<TVector3, double> >**  rho = NULL; // two shellindices of the nucleons, pair contains P vector and result
	rho = initRho(nuc.getPLevels());
	vector<TVector3> ps_list;	
	
	// make function that constructs whole rho matrix, or list for different quantum number n1,k1,n2,k2
	
	for (int shellindex1 = 0; shellindex1 < nuc.getPLevels() ; shellindex1++)
	{
		cout << shellindex1+1 << " to " << nuc.getPLevels() << endl;
		for (int shellindex2 = shellindex1; shellindex2 < nuc.getPLevels() ; shellindex2++)
		{
			cout <<">" << shellindex2+1 << " to " << nuc.getPLevels() << endl;
	
			for (int i=0; i<3*nsteps; i++)
			{
				cout << ">>" << i+1 << " to " << 3*nsteps << endl;
				TVector3 ps = Pcm[i]-pm;
				double res = dist2bodymom(nuc,shellindex1,shellindex2,Pcm[i],pm,ps);
				rho[shellindex1][shellindex2].push_back(pair<TVector3,double>(Pcm[i],res));
				ps_list.push_back(ps);
			}
		}
	}

	writeRho("rho_notes.dat",nuc,rho);
	delRho(rho,nuc.getPLevels());
}

void writeRho(string filename, MeanFieldNucleus& nuc, vector< pair<TVector3,double> > **rho)
{
	// write to file
	ofstream data(filename.c_str());
	for (int shellindex1=0; shellindex1 < nuc.getPLevels(); shellindex1++)
	{
		for (int shellindex2=shellindex1; shellindex2 < nuc.getPLevels();shellindex2++)
		{
			data << "**n1=" << nuc.getN_array()[shellindex1] << "|k1="<< nuc.getKappas()[shellindex1] << "|l1=" << nuc.getL_array()[shellindex1] << "|j1="<< nuc.getJ_array()[shellindex1];
			data <<  "|n2=" << nuc.getN_array()[shellindex2] << "|k2="<< nuc.getKappas()[shellindex2] << "|l2=" << nuc.getL_array()[shellindex2] << "|j2="<< nuc.getJ_array()[shellindex2] << endl;
			for (vector< pair<TVector3,double> >::iterator it=rho[shellindex1][shellindex2].begin(); it != rho[shellindex1][shellindex2].end(); it++) {
				data << it->first.X() << "\t" << it->first.Y() << "\t" << it->first.Z() << "\t" << it->second << endl;
			}
		}
	}
	data.close();
}
double dist2bodymom(MeanFieldNucleus& nuc, int shellindex1, int shellindex2, TVector3& P,TVector3& pm, TVector3& ps)
{
	double intgrl_res = 0.; // result of the integral we will return at the end
	TVector3 diff = P-(pm+ps);
	if (diff.Mag() > 1e-12)
	{
		cerr << " com momentum mismatch; difference between P and pm+ps should be zero " << endl;
		cerr << " P  = " << vecString(P) << endl;
		cerr << " pm = " << vecString(pm) << endl;
		cerr << " ps = " << vecString(ps) << endl;
		exit(1);
	}
	// init fourvectors to construct spinors
	FourVector<double> f_pm = FourVector<double>(sqrt(MASSP*MASSP + pm.Dot(pm)) ,pm.X(),pm.Y(),pm.Z());
	FourVector<double> f_ps = FourVector<double>(sqrt(MASSP*MASSP + ps.Dot(ps)), ps.X(), ps.Y(), ps.Z());
	
	if (FORCE_POS_ENERGY_PROJ) // to check Wim's previous results we take a look at pos energy terms (neglect lower part of spinor, so set momentum->zero)
	{
		f_pm = FourVector<double>(sqrt(MASSP*MASSP+pm*pm), 0, 0, 0);
		f_ps = FourVector<double>(sqrt(MASSP*MASSP+ps*ps), 0, 0, 0);
	}

	// prepare integration
	numint::array<double,3> lowerb = {{0.,1.,0.}}; // r, cos(theta), phi
	numint::array<double,3> upperb = {{nuc.getRange(),-1,2.*PI}}; // r, cos(theta), phi
	numint::mdfunction<complex<double>,3> mdf;

	// make struct and set struct parameters
	F f;
	f.nuc = &nuc;
	f.f = integrandum;
	mdf.func = &F::exec;
	mdf.param = &f;
	f.shellindex1 = shellindex1;
	f.shellindex2 = shellindex2;
	f.P = P;
	
	TSpinor::Polarization::State pol[2] = {TSpinor::Polarization::kDown,TSpinor::Polarization::kUp}; // make a list to use my beloved basic forloops
	for (int s1=0; s1<2; s1++) // sum s_1
	{
		TSpinor u_pm = TSpinor(f_pm,MASSP,TSpinor::Polarization(0.,0.,pol[s1]),TSpinor::kUnity);
		f.u_pm = &u_pm;
		for (int s2=0; s2<2; s2++) // sum s_2
		{
			TSpinor u_ps = TSpinor(f_ps,MASSP,TSpinor::Polarization(0.,0.,pol[s2]),TSpinor::kUnity);
			f.u_ps = &u_ps;
			for (int m1 = -nuc.getJ_array()[shellindex1]; m1 <= nuc.getJ_array()[shellindex1]; m1+=2) // sum m1
			{
				// if shellindices are equal, m2 > m1 else we will have double counting!!!
				for (int m2 = (f.shellindex1 == f.shellindex2)? m1+2 : -nuc.getJ_array()[shellindex2] ; m2 <= nuc.getJ_array()[shellindex2]; m2+=2) // Pauli principle -- never forget --
				{
					if ((f.shellindex1==f.shellindex2) && (m1==m2)) m2+=2; // pauli!!!
					if (m2 > nuc.getJ_array()[shellindex2]) break; // we tried to change m2 because pauli but m2 is now 2 large
					//cout << "shellind 1 : " << setw(2) << f.shellindex1 << " shellind2 " << setw(2) <<f.shellindex2 <<" m1 : " << setw(2) << m1 << " , m2 : " << setw(2) << m2 << " pol1 : " << setw(2) << pol[s1] << " pol2: " << setw(2) << pol[s2] << endl;
					f.m1 = m1;
					f.m2 = m2;
					double prec=1e-08;
					complex<double> ret;
					int maxEval = 10000;
					unsigned count;
					int succes = numint::cube_adaptive(mdf,lowerb,upperb,1e-06,prec,maxEval,ret,count,0);	
					intgrl_res += norm(ret);
					assert(succes==0);
					//cout << succes << "  --  " << count << " -- " << ret << "    " << intgrl_res << endl;
				} // m2 loop
			} // m1 loop
		} // s2 loop
	} // s1 loop
	return (1./pow(2.*PI,3.))*intgrl_res;
	//return intgrl_res;
}


void integrandum(complex<double>& res, const numint::array<double,3>& x, TVector3& P, TSpinor* u_pm, TSpinor* u_ps, MeanFieldNucleus* nuc, int shellindex1, int m1, int shellindex2,int m2)
{
	double sinth = sqrt(1.-x[1]*x[1]);
	double PR = P*TVector3(x[0]*sinth*cos(x[2]),x[0]*sinth*sin(x[2]),x[0]*x[1]); // dot product of vector P with vector R
	//TVector3 R(); // convert spherical x[] array coordinates to a cartesian vector
	//R.SetMagThetaPhi(x[0],acos(x[1]),x[2]); // acos because x[1] is cos(theta)
	//double PR = P*R;
	res = exp(-std::I_UNIT*INVHBARC*PR)*TSpinor::Bar(*u_pm)*TMFSpinor(*nuc,shellindex1,m1,x[0],x[1],x[2])*TSpinor::Bar(*u_ps)*TMFSpinor(*nuc,shellindex2,m2,x[0],x[1],x[2]); // Wavefunctions already defined with factor r, so no r^2 that would appear from jacobian, no sin(theta) because we integrate over d(cos(theta))

}



void genPcmList(TVector3 Pcm[], int nsteps) // vary Pcm from -400 to 400 in n steps
{
	// construct a list for Pcm in three perp directions, take x y z, we can rotate later
	// vary Pcm from -400 to 400 [Mev]
	//Pcm = new TVector3[nsteps*3]; // times three because three perp directions
	double normPmin = -800.;
	double normPmax = 800.;
	double Pstep = (normPmax-normPmin)/(nsteps-1); // step value for Pcm
	for (int i=0; i<nsteps;i++) {Pcm[i]=TVector3(normPmin+i*Pstep,0.,0.);} // x-dir
	for (int i=nsteps; i<2*nsteps;i++) Pcm[i]=TVector3(0.,normPmin+(i-nsteps)*Pstep,0); //y-dir
	for (int i=2*nsteps; i<3*nsteps;i++) { Pcm[i]=TVector3(0.,0.,normPmin+(i-2*nsteps)*Pstep);} //z-dir
}










// --------------------- pure info print functions below ---------------------------- //
void printNucleusInfo(MeanFieldNucleus& nuc)
{
	cout << "nuc nucleus, with " << nuc.getZ() << " protons and " << nuc.getN() << " neutrons " << endl;
	cout << "Total number of shells " << nuc.getTotalLevels() << " of which " << nuc.getPLevels() << " are proton shells and " << nuc.getNLevels() << " are neutron shells " << endl;
	cout << "n array " << endl; printArray(nuc,nuc.getN_array());
	cout << "k array " << endl; printArray(nuc,nuc.getKappas());
	cout << "l array " << endl; printArray(nuc,nuc.getL_array());
	cout << "j array " << endl; printArray(nuc,nuc.getJ_array()); // returns j_max times two per level!
}

void printArray(MeanFieldNucleus &nuc,const int* arr)
{
	cout << " proton shells: " << endl;
	for (int i=0; i<nuc.getPLevels(); i++)
	{
		cout << "   "<< i << ": " << arr[i] << endl;
	}
	cout << " neutron shells: " << endl;
	for (int i=0; i<nuc.getNLevels(); i++)
	{
		cout << "   " << i << ": " << arr[i] << endl;
	}
}
