#include <iostream>
#include <fstream>

#include "SparseSystemSolver.h"
#include "Log.h"
#include "Timer.h"


namespace picdsmc {

	SparseSystemSolver::SparseSystemSolver() {}
	SparseSystemSolver::~SparseSystemSolver() {}

	void SparseSystemSolver::SetOrder(int n) {
		m_n = n;

		m_A.resize(m_n, m_n);
		m_A.setZero();

		m_b.resize(m_n);
		m_b.setZero();
	}

	void SparseSystemSolver::InsertCoefficient(unsigned char MatrixSelection, int row, int col, double value) {
		if (MatrixSelection == 0) { // Insert coefficient for matrix A.
			m_coefficients.push_back(Eigen::Triplet<double>(row, col, value));
		}
		else { // Inster coefficient for the vector b.
			m_b(row) = value;
		}

	}
	void SparseSystemSolver::Solve() {
		m_A.setFromTriplets(m_coefficients.begin(), m_coefficients.end());	 // Fill the sparse matrix.
		m_A.makeCompressed();

		LOG_CRITICAL("Threads {}", Eigen::nbThreads());

		LOG_INFO("Solving a {} x {} system.", m_n, m_n);
		Timer timer;
		Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> cg;

		cg.setMaxIterations(5000);
		cg.setTolerance(0.00003);
		cg.compute(m_A);
		m_x = cg.solve(m_b);
		LOG_INFO("Solved in {:.5} s.", timer.Elapsed());


		#if 1
		LOG_INFO("Number of iterations: {}", cg.iterations());
		if(cg.error()>1){
			LOG_CRITICAL("Estimated error: {:.5}", cg.error());
		}
		else {
			LOG_INFO("Estimated error: {:.5}", cg.error());
		}
		#endif
	}
	void SparseSystemSolver::OutputSolution(const std::string& FileName) {
		std::ofstream file(FileName + ".txt");

		for (int i = 0; i < m_n; i++) {
			file << m_x[i] << std::endl;
		}
		file.close();
	}

}