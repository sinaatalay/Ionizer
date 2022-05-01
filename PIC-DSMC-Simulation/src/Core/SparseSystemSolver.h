#pragma once

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <vector>

namespace picdsmc{

	class SparseSystemSolver{
	public:
		SparseSystemSolver();
		~SparseSystemSolver();

		void SetOrder(int n);
		void SetbZero();

		void InsertCoefficient(unsigned char MatrixSelection, int row, int col, double value);
		void Solve();
		void OutputSolution();

	private:
		int m_n;											// The order of Ax=b.

		//A triplet is a simple object representing a non-zero entry as the triplet: row index, column index, value:
		std::vector<Eigen::Triplet<double>> m_coefficients;	// The values of the nonzero elements of A of Ax=b.

		Eigen::VectorXd m_b;									// b vector of Ax=b.
		Eigen::SparseMatrix<double,Eigen::RowMajor> m_A;						// The sparse matrix A of Ax=b.
		Eigen::VectorXd m_x;									// The solution of Ax=b.
	};

}