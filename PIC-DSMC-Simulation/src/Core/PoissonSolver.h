#pragma once
#include "Geometry.h"

class PoissonSolver : public Geometry {
public:
	explicit PoissonSolver();
	~PoissonSolver();

	void AllocateMemory();
	void ConfigureCoefficients();
	void ApplyBoundaryConditions();
	void SolvePoisson(int m_MaxIterations, double m_AbsoluteTolerance);
private:
	// Poisson equation in Ax=b form:
	int m_n;			// The order of the linear system.
	double* m_A;		// m_A[i] = the values of the nonzero elements

	int m_Alength;		// The number of nonzero elements.
	int* m_Arow;		// m_Ar[i] = the row indices of the nonzero elements
	int* m_Acol;		// m_Acol[i] = the column indices of the nonzero elements

	double* m_b;		// b vector
	double* m_x;		// The solution x, that satisfies Ax=b
};