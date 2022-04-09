#pragma once
#include "Geometry.h"

class PoissonSolver : public Geometry {
public:
	explicit PoissonSolver();
	~PoissonSolver();

	void ConfigureCoefficients();
private:
	// Poisson equation in Ax=b form:
	double* m_A;		// m_A[i] = the values of the nonzero elements
	int* m_Arow;		// m_Ar[i] = the row indices of the nonzero elements
	int* m_Acol;		// m_Acol[i] = the column indices of the nonzero elements
	double* m_b;		// b vector.
};