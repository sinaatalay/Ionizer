#pragma once
#include <Eigen/Sparse>
#include <vector>

#include "Geometry.h"


namespace picdsmc {

	class PoissonSolver : public Geometry {
	public:
		PoissonSolver();
		explicit PoissonSolver(const Geometry& geometry);
		~PoissonSolver();

		void ConfigureMatrixA();
		void ApplyBoundaryConditions();
		void OutputMatrixA();

		void buildProblem(std::vector<T>& coefficients, Eigen::VectorXd& b, int n);
		void insertCoefficient(int id, int i, int j, int k, double w, std::vector<T>& coeffs, Eigen::VectorXd& b, int step_j, int step_i);

		void SolvePoisson(int m_MaxIterations, int KrylovSubspaceDimension, double m_AbsoluteTolerance);
	private:
		// Poisson equation in Ax=b form:
		int m_n;			// The order of the linear system
		int m_Alength;		// The number of nonzero elements

		double* m_A;		// m_A[i] = the values of the nonzero elements
		int* m_Arow;		// m_Ar[i] = the row indices of the nonzero elements
		int* m_Acol;		// m_Acol[i] = the column indices of the nonzero elements

		double* m_b;		// b vector
		double* m_x;		// The solution x, that satisfies Ax=b
	};

}