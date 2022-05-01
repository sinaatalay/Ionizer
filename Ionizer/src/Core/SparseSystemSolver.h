#pragma once

//C++ Libraries (Precompiled):
#include "Ionpch.h"

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

namespace Ionizer{

	class SparseSystemSolver{
	public:
		SparseSystemSolver();
		~SparseSystemSolver();

		// SetOrder() sets the order of the linear system to be solved:
		void SetOrder(int n);

		// InsertCoefficient() inserts a coefficient to matrix A or vector b into the specified index:
		void InsertCoefficient(uint8_t MatrixSelection, int row, int col, double value);

		// Solve() applies the biconjugate gradient stabilized method to the sparse linear system, Ax=b:
		void Solve();

		// OutputSolution() outputs the solution vector x to a *.txt file.
		void OutputSolution(const std::string& FileName);

	private:
		int m_n;												// The order of Ax=b.

		std::vector<Eigen::Triplet<double>> m_coefficients;		// The values of the nonzero elements of A of Ax=b.

		Eigen::VectorXd m_b;									// b vector of Ax=b.
		Eigen::SparseMatrix<double,Eigen::RowMajor> m_A;		// The sparse matrix A of Ax=b.
		Eigen::VectorXd m_x;									// The solution of Ax=b.
	};

}