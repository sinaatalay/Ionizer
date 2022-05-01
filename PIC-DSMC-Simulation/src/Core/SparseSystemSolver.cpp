#include <omp.h>

#include <iostream>
#include <fstream>
#include <chrono>

#include "SparseSystemSolver.h"
#include "Log.h"
#include <Eigen/src/SparseCore/SparseDenseProduct.h>

using namespace picdsmc;


SparseSystemSolver::SparseSystemSolver() {}
SparseSystemSolver::~SparseSystemSolver() {}
//===============================================
void SparseSystemSolver::SetOrder(int n) {
	m_n = n;

	m_A.resize(m_n, m_n);
	m_b.resize(m_n);
}
void SparseSystemSolver::SetbZero() {
	m_b.setZero();
}
//===============================================
void SparseSystemSolver::InsertCoefficient(unsigned char MatrixSelection, int row, int col, double value) {
	if (MatrixSelection == 0) { // Insert coefficient for matrix A.
		m_coefficients.push_back(Eigen::Triplet<double>(row, col, value));
	}
	else { // Inster coefficient for the vector b.
		m_b(row) = value;
	}

}
void SparseSystemSolver::Solve() {

	Eigen::setNbThreads(4);
	m_A.setFromTriplets(m_coefficients.begin(), m_coefficients.end());	 // Fill the sparse matrix.
	m_A.makeCompressed();

	LOG_INFO("Solving a {} x {} system.", m_n, m_n);
	auto start = std::chrono::steady_clock::now();
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::IncompleteLUT<double>> cg;
	cg.setMaxIterations(5000);
	cg.setTolerance(0.00001);
	cg.compute(m_A);
	m_x = cg.solve(m_b);
	auto end = std::chrono::steady_clock::now();
	std::chrono::milliseconds t = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	LOG_INFO("Solved in {:.5} s.", (float)t.count() / 1000);
	LOG_INFO("Number of iterations: {}", cg.iterations());
	LOG_INFO("Estimated error: {:.5}", cg.error());

	OutputSolution();
}
void SparseSystemSolver::OutputSolution() {
	std::ofstream file("PoissonMatrixA.txt");

	for (int i = 0; i < m_n; i++) {
		file << m_x[i] << std::endl;
	}
	file.close();
}
