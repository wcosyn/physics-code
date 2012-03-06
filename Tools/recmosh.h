#ifndef RECMOSH_H
#define RECMOSH_H

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_coupling.h>
#include <cmath>
using std::fabs;
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <fstream>
using std::ofstream;
using std::ifstream;
#include <sstream>
using std::stringstream;
#include <map>
using std::map;

#include <string>
using namespace std;

class RecMosh
{
private: 	
	int n1;
	int l1;
	int n2;
	int l2;
	map< int, double > coefficients;
	string path;
	double Wc( int a, int b, int c, int d, int e, int f);
	double A( int l1, int l, int l2, int Lambda, int x );
	double calculate( int n, int l, int N, int Lambda, int n1, int l1, int n2, int l2, int L);
	double getMatrixElement( int n, int l, int N, int Lambda, int na, int la, int Na, int Lambdaa, int L, int f);
	void loadFile();
	void writeToFile();
public:
	RecMosh(int n1, int l1, int n2, int l2, string path);
	~RecMosh();
	double getCoefficient( int n, int l, int N, int Lambda, int L );
// 	RecMosh(int n1, int l1, int n2, int l2, int n, int l, int N, int Lambda, int L, bool reversed);
};

#endif // RECMOSH_H
