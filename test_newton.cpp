// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "fcn.hpp"

using namespace std;
// Function declaration
extern double newton(Fcn& f, Fcn& df, double x, int maxit, double tol, bool show_iterates);
// Root-finding residual function
class fcn : public Fcn {
public:
	double operator()(double x) {   // function evaluation
		return pow(x, 2.0) * (x - 3) * (x + 2);
	}
};

class dfcn : public Fcn {
public:
	double operator()(double x) {   // derivative function evaluation
		return 4 * pow(x, 3.0) - 3 * pow(x, 2.0) - 12 * x;
	}
};


// This routine tests the function newton.cpp on a nonlinear 
// function.  It uses a sevral initial guesses of x { -3, 1, 2 }, and solves for 
// the root to a different tolerance values { 10e-1, 10e-5, 10e-9 }
int main(int argc, char* argv[]) {

	// set the initial values
	double x0[] = { -3, 1, 2 };
	double tol[] = { 10e-1, 10e-5, 10e-9 };
	int maxit = 50;
	// create the functions objects
	fcn f;
	dfcn df;
	// Save initial x, tolerance, and result x values to files
	ofstream result_out("newton_result.txt");

	for( int i = 0; i < sizeof(x0) / sizeof(double); ++i ) {
		
		double initial_x = x0[i];
		
		for( int j = 0; j < sizeof(tol) / sizeof(double); ++j ) {
			double tolerance = tol[j];			
			cout << "Initial x=" << initial_x << ", tol=" << tolerance << endl;
			
			// call newton to compute the root for current initial guess and tolerance, and output result to screen
			double x = newton(f, df, initial_x, maxit, tolerance, true);
			cout << endl 
				<< " The approximate root is " << setprecision(16) 
				<< x << endl;
			cout << "========================================" << endl;
			result_out << initial_x << " " << tolerance << " " << x << endl;
		}
	}

	// Close result files
	result_out.close();

	return 0;
}