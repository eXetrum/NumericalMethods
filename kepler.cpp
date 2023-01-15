// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "fcn.hpp"
#include "matrix.hpp"
using namespace std;

// Object position 
struct Point { double x, y; };
// Root-finding residual function
class fcn : public Fcn {
	double e, M;
public:
	fcn(double e, double M) : e(e), M(M) {}

	double operator()(double E) {   // function evaluation
		return e * sin(E) - E + M; //return (M - E + e*sin(E));//w - e * sin(w) - t;
	}
};

class dfcn : public Fcn {
	double e;
public:
	dfcn(double e) : e(e) {}
	double operator()(double E) {   // derivative function evaluation
		return e * cos(E) - 1; //return 1 - e * cos(E);
	}
};

// Function declaration
extern double newton(Fcn& f, Fcn& df, double x, int maxit, double tol, bool show_iterates);
double getRadialPosition(double a, double b, double w);
Matrix generateTimeVec(double min_value, double max_value, double step);
Point getCartesian(double r, double w);

// This routine solve kepler problem
int main(int argc, char* argv[]) {
	// Newton method params
	const int maxiter = 6;
	const double tol = 10e-5;
	// Orbital params
	const double a = 2.0;
	const double b = 1.25;
	const double eccentricity = sqrt( 1.0 - (b * b) / (a * a) );
	// Generate time vector
	Matrix t = generateTimeVec(0.0, 10, 0.001);
	// Determine result vectors size
	int size = t.Size();
	// Create result vectors for x(t), y(t)
	Matrix x(size), y(size);
	// Calculate w, x(t), y(t) for each t[i]
	for(int i = 0; i < size; ++i) {
		// Get next time
		double t0 = t.data[0][i];
		// Create f object
		fcn f(eccentricity, t0);
		// Create df object
		dfcn df(eccentricity);
		// Set initial angle
		double w = t0;
		// Calculate w angle
		w = newton(f, df, w, maxiter, tol, false);
		// Get radial position
		double r = getRadialPosition(a, b, w);
		// Calculate cartesian coords
		Point p = getCartesian(r, w);
		// Save to the vectors
		x.data[0][i] = p.x;
		y.data[0][i] = p.y;
	}
	// Write results to the files
	t.Write("t.txt");
	x.Write("x.txt");
	y.Write("y.txt");

	return 0;
}


Matrix generateTimeVec(double min_value, double max_value, double step) {
	const double EPS = 10e-18;
	int size = (max_value - min_value) / step + 1 ;
	Matrix t( size );
	int i = 0;
	double value = min_value;
	while( abs(value - max_value) >= EPS && i < size ) {
		t.data[0][i++] = value;
		value += step;
	}

	return t;
}

double getRadialPosition(double a, double b, double w) {
	double bc = b * cos(w);
	double ac = a * sin(w);
	return a * b / sqrt( bc * bc + ac * ac );
}

Point getCartesian(double r, double w) {
	Point p;
	p.x = r * cos(w);
	p.y = r * sin(w);
	return p;
}