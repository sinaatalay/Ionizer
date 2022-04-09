#include "PoissonSolver.h"

PoissonSolver::PoissonSolver() {
	m_ThetaBegin = 0.0;
	m_AxialBegin = 0.0;
	m_RadialBegin = 0.0;
}
PoissonSolver::~PoissonSolver() {}

void PoissonSolver::ConfigureCoefficients() {
	m_A = new double[m_TotalNodeCount * 7];		// Allocate dynamic memory for m_A array for the first time.
	m_Arow = new int[m_TotalNodeCount * 7];		// Allocate dynamic memory for m_Arow array for the first time.
	m_Acol = new int[m_TotalNodeCount * 7];		// Allocate dynamic memory for m_Acol array for the first time.
	m_b = new double[m_TotalNodeCount * 7];		// Allocate dynamic memory for m_Ab array for the first time.

	double r;
	int index;
	int step_i = m_RadialNodeCount * m_AxialNodeCount * 7;
	int step_j = m_AxialNodeCount * 7;
	int step_k = 7;

	// i: r nodes.
	// j: theta nodes.
	// k: z nodes.
	// Note that the nodes at the boundaries are untouched in the loop below.
	for (int i = 1; i < m_ThetaNodeCount - 1; i++) {
		for (int j = 1; j < m_RadialNodeCount - 1; j++) {
			for (int k = 1; k < m_AxialNodeCount - 1; k++) {
				r = m_RadialBegin + j * m_dr;
				index = i * step_i + j * step_j + k * step_k;

				m_A[index + 0] = 1.0 / (r * r * m_dtheta * m_dtheta);											// U_{i,	j-1,	k}
				m_A[index + 1] = 1.0 / (m_dr * m_dr) - 1.0 / (r * 2.0 * m_dr);									// U_{i-1,	j,		k}
				m_A[index + 2] = 1.0 / (m_dz * m_dz);															// U_{i,	j,		k-1}
				m_A[index + 3] = -2 / (m_dz * m_dz) - 2 / (m_dr * m_dr) - 2 / (r * r * m_dtheta * m_dtheta);	// U_{i,	j,		k}
				m_A[index + 4] = m_A[index + 2];																// U_{i,	j,		k+1}
				m_A[index + 5] = 1.0 / (m_dr * m_dr) + 1.0 / (r * 2.0 * m_dr);									// U_{i+1,	j,		k}
				m_A[index + 6] = m_A[index + 0];																// U_{i,	j+1,	k}

				m_Arow[index + 0] = index;
				m_Arow[index + 1] = index;
				m_Arow[index + 2] = index;
				m_Arow[index + 3] = index;
				m_Arow[index + 4] = index;
				m_Arow[index + 5] = index;
				m_Arow[index + 6] = index;

				m_Acol[index + 0] = index - step_j;
				m_Acol[index + 1] = index - step_i;
				m_Acol[index + 2] = index - step_k;
				m_Acol[index + 3] = index;
				m_Acol[index + 4] = index + step_k;
				m_Acol[index + 5] = index + step_i;
				m_Acol[index + 6] = index + step_j;

				m_b[index] = 0;
			}
		}
	}
}