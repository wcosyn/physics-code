# include <pair.h>

#include <string>
using namespace std;


void test1();
void test2();
/*
 * Examples of usage Pair class.
 */
int main( int argc, char *argv[])
{
	test1();
	test2();

}

void test2()
{
	Pair* pn = new Pair("paris", "..", 2, 0, 0, 1, -1, 1, 0, 0, 1, -1, -1);
	double rmax= 5;
	double deltar = 0.05;
	int n= 0;
	int l= 0;
	int S= 1;
	int j= 1;
	int mj= 1;
	int T= 0;
	ofstream file( "test.dat", ofstream::out );

	for( double r1= 0; r1< rmax; r1+=deltar)
	{
		for( double r2= 0; r2< rmax; r2+=deltar)
		{
			complex<double> wf = pn->getcomwf( r1, 0, M_PI, r2, 0, M_PI, n, l, S, j, mj, T);
			file << r1 << "\t" << r2;
			file << "\t" << real(wf);
			file << "\t" << imag(wf);
			file << endl;
		}
		file << endl;
	}
	file.close();
	delete pn;
}

void test1()
{
	
	Pair* pp = new Pair("paris", "..", 2, 0, 0, 1, -1, 1, 0, 0, 1, -1, 1 );
	cout << pp->getRelPair( 0 ) << endl;	
	delete pp;
	pp = new Pair("paris", "..", 2, 0, 0, 1, -1, 1, 0, 0, 1, 1, 1 );
	cout << pp->getRelPair( 0 ) << endl;	
	delete pp;
	pp = new Pair("paris", "..", 2, 0, 0, 1, 1, 1, 0, 0, 1, -1, 1 );
	cout << pp->getRelPair( 0 ) << endl;	
	delete pp;
	pp = new Pair("paris", "..", 2, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1 );
	cout << pp->getRelPair( 0 ) << endl << endl;	
	delete pp;

	Pair* pn = new Pair("paris", "..", 2, 0, 0, 1, -1, 1, 0, 0, 1, -1, -1 );
	cout << pn->getRelPair( 0 ) << endl;	
	delete pn;
	pn = new Pair("paris", "..", 2, 0, 0, 1, -1, 1, 0, 0, 1, 1, -1 );
	cout << pn->getRelPair( 0 ) << endl;	
	delete pn;
	pn = new Pair("paris", "..", 2, 0, 0, 1, 1, 1, 0, 0, 1, -1, -1 );
	cout << pn->getRelPair( 0 ) << endl;	
	delete pn;
	pn = new Pair("paris", "..", 2, 0, 0, 1, 1, 1, 0, 0, 1, 1, -1 );
	cout << pn->getRelPair( 0 ) << endl << endl;	
	delete pn;

	Pair* nn = new Pair("paris", "..", 2, 0, 0, 1, -1, -1, 0, 0, 1, -1, -1 );
	cout << nn->getRelPair( 0 ) << endl;	
	delete nn;
	nn = new Pair("paris", "..", 2, 0, 0, 1, -1, -1, 0, 0, 1, 1, -1 );
	cout << nn->getRelPair( 0 ) << endl;	
	delete nn;
	nn = new Pair("paris", "..", 2, 0, 0, 1, 1, -1, 0, 0, 1, -1, -1 );
	cout << nn->getRelPair( 0 ) << endl;	
	delete nn;
	nn = new Pair("paris", "..", 2, 0, 0, 1, 1, -1, 0, 0, 1, 1, -1 );
	cout << nn->getRelPair( 0 ) << endl << endl;	
	delete nn;
}
