#//C++ Libraries (Precompiled):
#include "Ionpch.h"

#include "SparseSystemSolver.h"
#include "Log.h"
#include "Timer.h"


namespace Ionizer {

	SparseSystemSolver::SparseSystemSolver() {}
	SparseSystemSolver::~SparseSystemSolver() {}

	void SparseSystemSolver::SetOrder(int n) {
		m_n = n;

		m_A.resize(m_n, m_n);
		m_A.setZero();

		m_b.resize(m_n);
		m_b.setZero();
	}

	void SparseSystemSolver::InsertCoefficient(uint8_t MatrixSelection, int row, int col, double value) {
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

		LOG_INFO("Solving a {} x {} system.", m_n, m_n);
		LOG_INFO("The number of threads will be used: {}", Eigen::nbThreads());
		Timer timer;
		Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;

		solver.setMaxIterations(500);
		solver.setTolerance(0.00001);
		solver.compute(m_A);
		m_x = solver.solve(m_b);
		LOG_INFO("Solved in {:.5} s.", timer.Elapsed());

		LOG_INFO("Number of iterations: {}", solver.iterations());
		if (solver.error() > 1) {
			LOG_CRITICAL("Estimated error: {:.5}", solver.error());
		}
		else {
			LOG_INFO("Estimated error: {:.5}", solver.error());
		}
	}

	void SparseSystemSolver::OutputSolution(const std::string& FileName) {
		std::ofstream file;
		if (FileName.find(".txt") == std::string::npos) {
			file.open(FileName + ".txt");
		}
		else {
			file.open(FileName);
		}

		for (int i = 0; i < m_n; i++) {
			file << m_x[i] << std::endl;
		}
		file.close();
	}

	std::vector<double> SparseSystemSolver::GetSolution() {
		std::vector<double> result;
		result.resize(m_n);
		for (int i = 0; i < m_n; i++) {
			result[i] = m_x[i];
		}
		return result;
	}
}