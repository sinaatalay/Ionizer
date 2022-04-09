#include "PoissonSolver.h"
#include "Mgmres.h"
#include "Log.h"

PoissonSolver::PoissonSolver() {
	m_ThetaBegin = 0.0;
	m_AxialBegin = 0.0;
	m_RadialBegin = 0.0;
}
PoissonSolver::~PoissonSolver() {
	delete[] m_A;		// Release obtained memory area.
	delete[] m_Arow;	// Release obtained memory area.
	delete[] m_Acol;	// Release obtained memory area.
	delete[] m_b;		// Release obtained memory area.
	delete[] m_x;		// Release obtained memory area.
}

void PoissonSolver::AllocateMemory() {
	m_n = m_TotalNodeCount;
	m_Alength = m_n * 7;

	m_A = new double[m_Alength];		// Allocate dynamic memory for m_A array for the first time.
	m_Arow = new int[m_Alength];		// Allocate dynamic memory for m_Arow array for the first time.
	m_Acol = new int[m_Alength];		// Allocate dynamic memory for m_Acol array for the first time.
	m_b = new double[m_Alength];		// Allocate dynamic memory for m_Ab array for the first time.
	m_x = new double[m_Alength];		// Allocate dynamic memory for m_Ab array for the first time.
}

void PoissonSolver::ConfigureCoefficients() {
	double r;
	int index;
	int step_j = m_RadialNodeCount * m_AxialNodeCount * 7;
	int step_i = m_AxialNodeCount * 7;
	int step_k = 7;

	// i: r nodes.
	// j: theta nodes.
	// k: z nodes.
	// Note that the nodes at the boundaries are untouched in the loop below.
	for (int j = 1; j < m_ThetaNodeCount - 1; j++) {
		for (int i = 1; i < m_RadialNodeCount - 1; i++) {
			for (int k = 1; k < m_AxialNodeCount - 1; k++) {
				r = m_RadialBegin + i * m_dr;
				index = j * step_j + i * step_i + k * step_k;

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

void PoissonSolver::ApplyBoundaryConditions() {
	int index;
	int step_j = m_RadialNodeCount * m_AxialNodeCount * 7;
	int step_i = m_AxialNodeCount * 7;
	int step_k = 7;

	// i: r nodes.
	// j: theta nodes.
	// k: z nodes.

	// Screen Grid:
	for (int j = 0; j < m_ThetaNodeCount; j++) {
		for (int i = m_rScreenBeginNode; i < m_RadialNodeCount; i++) {
			for (int k = m_zScreenBeginNode; k <= m_zScreenEndNode; k++) {
				index = j * step_j + i * step_i + k * step_k;

				m_A[index + 0] = 0;			// U_{i,	j-1,	k}
				m_A[index + 1] = 0;			// U_{i-1,	j,		k}
				m_A[index + 2] = 0;			// U_{i,	j,		k-1}
				m_A[index + 3] = 1;			// U_{i,	j,		k}
				m_A[index + 4] = 0;			// U_{i,	j,		k+1}
				m_A[index + 5] = 0;			// U_{i+1,	j,		k}
				m_A[index + 6] = 0;			// U_{i,	j+1,	k}

				m_b[index] = m_VScreen;
			}
		}
	}

	// Acceleration Grid:
	for (int j = 0; j < m_ThetaNodeCount; j++) {
		for (int i = m_rScreenBeginNode; i < m_RadialNodeCount; i++) {
			for (int k = m_zScreenBeginNode; k <= m_zScreenEndNode; k++) {
				index = j * step_j + i * step_i + k * step_k;

				m_A[index + 0] = 0;			// U_{i,	j-1,	k}
				m_A[index + 1] = 0;			// U_{i-1,	j,		k}
				m_A[index + 2] = 0;			// U_{i,	j,		k-1}
				m_A[index + 3] = 1;			// U_{i,	j,		k}
				m_A[index + 4] = 0;			// U_{i,	j,		k+1}
				m_A[index + 5] = 0;			// U_{i+1,	j,		k}
				m_A[index + 6] = 0;			// U_{i,	j+1,	k}

				m_b[index] = m_VAccel;
			}
		}
	}

	// r = m_RadialLength surface:
	int i = m_RadialNodeCount - 1; // -1 because, indexing starts from 0.
	for (int j = 0; j < m_ThetaNodeCount; j++) {
		for (int k = 0; k < m_AxialNodeCount; k++) {
			index = j * step_j + i * step_i + k * step_k;

			m_A[index + 0] = 0;			// U_{i,	j-1,	k}
			m_A[index + 1] = -1;		// U_{i-1,	j,		k}
			m_A[index + 2] = 0;			// U_{i,	j,		k-1}
			m_A[index + 3] = 1;			// U_{i,	j,		k}
			m_A[index + 4] = 0;			// U_{i,	j,		k+1}
			m_A[index + 5] = 0;			// U_{i+1,	j,		k}
			m_A[index + 6] = 0;			// U_{i,	j+1,	k}

			m_b[index] = 0;
		}
	}

	// r = 0 surface:
	i = 0;
	for (int j = 0; j < m_ThetaNodeCount; j++) {
		for (int k = 0; k < m_AxialNodeCount; k++) {
			index = j * step_j + i * step_i + k * step_k;

			m_A[index + 0] = 0;			// U_{i,	j-1,	k}
			m_A[index + 1] = 0;			// U_{i-1,	j,		k}
			m_A[index + 2] = 0;			// U_{i,	j,		k-1}
			m_A[index + 3] = -1;		// U_{i,	j,		k}
			m_A[index + 4] = 0;			// U_{i,	j,		k+1}
			m_A[index + 5] = 1;			// U_{i+1,	j,		k}
			m_A[index + 6] = 0;			// U_{i,	j+1,	k}

			m_b[index] = 0;
		}
	}

	// z = 0 surface:
	int k = 0;
	for (int j = 0; j < m_ThetaNodeCount; j++) {
		for (int i = 0; i < m_RadialNodeCount; i++) {
			index = j * step_j + i * step_i + k * step_k;

			m_A[index + 0] = 0;			// U_{i,	j-1,	k}
			m_A[index + 1] = 0;			// U_{i-1,	j,		k}
			m_A[index + 2] = 0;			// U_{i,	j,		k-1}
			m_A[index + 3] = 1;			// U_{i,	j,		k}
			m_A[index + 4] = 0;			// U_{i,	j,		k+1}
			m_A[index + 5] = 0;			// U_{i+1,	j,		k}
			m_A[index + 6] = 0;			// U_{i,	j+1,	k}

			m_b[index] = m_VDischarge;
		}
	}

	// z = m_AxialLength surface:
	k = m_AxialNodeCount-1; // -1 because, indexing starts from 0.
	for (int j = 0; j < m_ThetaNodeCount; j++) {
		for (int i = 0; i < m_RadialNodeCount; i++) {
			index = j * step_j + i * step_i + k * step_k;

			m_A[index + 0] = 0;			// U_{i,	j-1,	k}
			m_A[index + 1] = 0;			// U_{i-1,	j,		k}
			m_A[index + 2] = 0;			// U_{i,	j,		k-1}
			m_A[index + 3] = 1;			// U_{i,	j,		k}
			m_A[index + 4] = 0;			// U_{i,	j,		k+1}
			m_A[index + 5] = 0;			// U_{i+1,	j,		k}
			m_A[index + 6] = 0;			// U_{i,	j+1,	k}

			m_b[index] = m_VPlume;
		}
	}

	// theta = 0 surface:
	int j = 0;
	for (int i = 0; i < m_RadialNodeCount; i++) {
		for (int k = 0; k < m_AxialNodeCount; k++) {
			index = j * step_j + i * step_i + k * step_k;

			m_A[index + 0] = 0;			// U_{i,	j-1,	k}
			m_A[index + 1] = 0;			// U_{i-1,	j,		k}
			m_A[index + 2] = 0;			// U_{i,	j,		k-1}
			m_A[index + 3] = -1;		// U_{i,	j,		k}
			m_A[index + 4] = 0;			// U_{i,	j,		k+1}
			m_A[index + 5] = 0;			// U_{i+1,	j,		k}
			m_A[index + 6] = 1;			// U_{i,	j+1,	k}

			m_b[index] = 0;
		}
	}

	// theta = m_ThetaLength surface:
	j = m_ThetaNodeCount-1;
	for (int i = 0; i < m_RadialNodeCount; i++) {
		for (int k = 0; k < m_AxialNodeCount; k++) {
			index = j * step_j + i * step_i + k * step_k;

			m_A[index + 0] = -1;		// U_{i,	j-1,	k}
			m_A[index + 1] = 0;			// U_{i-1,	j,		k}
			m_A[index + 2] = 0;			// U_{i,	j,		k-1}
			m_A[index + 3] = 1;			// U_{i,	j,		k}
			m_A[index + 4] = 0;			// U_{i,	j,		k+1}
			m_A[index + 5] = 0;			// U_{i+1,	j,		k}
			m_A[index + 6] = 0;			// U_{i,	j+1,	k}

			m_b[index] = 0;
		}
	}

}

void PoissonSolver::SolvePoisson(int MaxIterations, double AbsoluteTolerance) {
	//    Input, int N, the order of the linear system.
	//
	//    Input, int NZ_NUM, the number of nonzero matrix values.
	//
	//    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column indices
	//    of the matrix values.
	//
	//    Input, double A[NZ_NUM], the matrix values.
	//
	//    Input/output, double X[N]; on input, an approximation to
	//    the solution.  On output, an improved approximation.
	//
	//    Input, double RHS[N], the right hand side of the linear system.
	//
	//    Input, int ITR_MAX, the maximum number of (outer) iterations to take.
	//
	//    Input, int MR, the maximum number of (inner) iterations to take.
	//    MR must be less than N.
	//
	//    Input, double TOL_ABS, an absolute tolerance applied to the
	//    current residual.
	//
	//    Input, double TOL_REL, a relative tolerance comparing the
	//    current residual to the initial residual.

	for (int i = 0; i < m_Alength; i++) {
		m_x[i] = 10;
	}
	
	//mgmres_st(m_n, m_Alength, m_Arow, m_Acol, m_A, m_x, m_b, MaxIterations, m_n-500, AbsoluteTolerance, 0.0);

	int n = 3;
	int len = 5;
	double A[5] = { 3,2,1,4,6 };
	int row[5] = { 0,0,1,1,2 };
	int col[5] = { 1,2,0,2,0 };
	double x[3] = { 1,1,1 };
	double b[3] = { 0,0,4 };
	int nn = 3;

	mgmres_st(n, len, row, col, A, x, b, MaxIterations, nn, AbsoluteTolerance);
}