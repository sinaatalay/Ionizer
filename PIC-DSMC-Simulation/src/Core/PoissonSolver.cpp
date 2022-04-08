#include "PoissonSolver.h"

PoissonSolver::PoissonSolver() {
	m_ThetaBegin = 0.0;
	m_AxialBegin = 0.0;
	m_RadialBegin = 0.0;
}
PoissonSolver::~PoissonSolver() {}

void PoissonSolver::ConfigureCoefficients() {
	m_A = new double[m_TotalNodeCount*6]; // Allocate memory for our m_A array for the first time.


	for (int i; m_RadialNodeCount; i++) {
		for (int j; m_ThetaNodeCount; j++) {
			for (int k; k < m_AxialNodeCount; k++) {
				m_A[k + m_AxialNodeCount*j + m_RadialNodeCount*m_AxialNodeCount*i]=
			}
		}
	}
}