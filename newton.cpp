#include "fcn.hpp"

#include <iostream>
#include <cmath>
using namespace std;

double newton(Fcn& f, Fcn& df, double x, int maxit, double tol, bool show_iterates) {
	// Iteration counter
	int iter = 1;
	// Initial guess of 'x'
	double xk = x;
	// Absolute value update
	double h = 0;

	do {
		// If flag enabled -> print info for each step
		if(show_iterates) {
			cout << "Iter=" << iter 
				<< ", x=" << xk 
				<< ", h=" << h
				<< ", |f(x)|=" << abs(f(xk)) 
				<< endl;
		}
		// Calculate new value of 'x'
		double new_x = xk - f(xk) / df(xk);
		// Calculate diff between 'old_x' and 'new_x'
		h = abs(new_x - xk);
		// Set up new 'x' value
		xk = new_x;
		// Go to the next iteration
		iter++;

	}while( iter < maxit && h > tol ); // Stop if we reach maxiterations or if we reach precision
	// Return result
	return xk;
}