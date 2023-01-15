#include "matrix.hpp"
#include <iostream>
using std::cout;
using std::cin;
using std::endl;
#include <vector>
using std::vector;
#include <iomanip>

Matrix createVandermonde(Matrix &v);

int main(int argc, char *argv[]) {
	// All size values for Vandermond matrix	
	vector<double> nValues;
	nValues.push_back(5);
	nValues.push_back(9);
	nValues.push_back(17);
	nValues.push_back(33);
	nValues.push_back(65);
	int size = nValues.size();
	// Result vectors
	Matrix N(size);
	Matrix errors(size);
	Matrix residulas(size);

	for( int i = 0; i < size; ++i ) {
		int n = nValues[i];	
		// Generate vector 'v' with 'n' equaly-spaced entries between 0 and 1 (inclusive both bounds)
		Matrix v = Linspace(0, 1, n);
		// Create Vandermonde 'n' square matrix
		Matrix A = createVandermonde(v);
		// Generate random vector 'x'
		Matrix x = Random(n);
		// Calculate vector 'b' (right hand side of equation A * x = b)
		Matrix b = A * x;
		// Find the solution x^
		Matrix x_ = Solve(A, b);
		// Calculate difference betwen original 'x' and solution 'x^'
		Matrix error = x - x_;
		// Calculate resudula vector
		Matrix b_ = A * x_;
		Matrix residula = (b - b_);
		double error_norm = Norm(error);
		double residula_norm = Norm(residula);
		// Print error, residula 2-norm
		cout << "N = " << std::setw(4) << n 
			<< ", ||Error||: " << std::setw(16) << error_norm 
			<< ", ||Residuals||: " << residula_norm 
			<< endl;
		// Save results
		N.data[0][i]			= n;
		errors.data[0][i]		= error_norm;
		residulas.data[0][i]	= residula_norm;		
	}

	// Write results to the files
	N.Write("N.txt");
	errors.Write("errors.txt");
	residulas.Write("residulas.txt");
	
	return 0;
}

Matrix createVandermonde(Matrix &v) {
	int n = v.Rows();
	Matrix A(n, n);
	// Fill Vandemonde matrix
	for( int i = 0; i < n; ++i ) {
		// Get next value from 'v' vector
		double value = v.Column(0)[i];
		A.data[0][i] = 1;
		for( int j = 1; j < n; ++j )
			A.data[j][i] = value * A.data[j - 1][i];
	}
	// Return result
	return A;
}