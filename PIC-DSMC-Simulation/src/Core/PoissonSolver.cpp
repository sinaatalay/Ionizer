#include <cstdlib>
#include <iostream>
#include <fstream>

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
	m_b = new double[m_n];				// Allocate dynamic memory for m_Ab array for the first time.
	m_x = new double[m_n];				// Allocate dynamic memory for m_Ab array for the first time.
}

void PoissonSolver::ConfigureMatrixA() {
	int step_j = m_RadialNodeCount * m_AxialNodeCount * 7;
	int step_i = m_AxialNodeCount * 7;
	int step_k = 7;

	// i: r nodes.
	// j: theta nodes.
	// k: z nodes.

	for (int j = 0; j < m_ThetaNodeCount; j++) {
		for (int i = 0; i < m_RadialNodeCount; i++) {
			for (int k = 0; k < m_AxialNodeCount; k++) {
				double r = m_RadialBegin + i * m_dr;
				int index = j * step_j + i * step_i + k * step_k;
				int realIndex = index / 7;

				m_A[index + 0] = 1.0 / (r * r * m_dtheta * m_dtheta);											// U_{i,	j-1,	k}
				m_A[index + 1] = 1.0 / (m_dr * m_dr) - 1.0 / (r * 2.0 * m_dr);									// U_{i-1,	j,		k}
				m_A[index + 2] = 1.0 / (m_dz * m_dz);															// U_{i,	j,		k-1}
				m_A[index + 3] = -2 / (m_dz * m_dz) - 2 / (m_dr * m_dr) - 2 / (r * r * m_dtheta * m_dtheta);	// U_{i,	j,		k}
				m_A[index + 4] = m_A[index + 2];																// U_{i,	j,		k+1}
				m_A[index + 5] = 1.0 / (m_dr * m_dr) + 1.0 / (r * 2.0 * m_dr);									// U_{i+1,	j,		k}
				m_A[index + 6] = m_A[index + 0];																// U_{i,	j+1,	k}

				m_Arow[index + 0] = realIndex;
				m_Arow[index + 1] = realIndex;
				m_Arow[index + 2] = realIndex;
				m_Arow[index + 3] = realIndex;
				m_Arow[index + 4] = realIndex;
				m_Arow[index + 5] = realIndex;
				m_Arow[index + 6] = realIndex;

				m_Acol[index + 0] = realIndex - step_j/7;
				m_Acol[index + 1] = realIndex - step_i/7;
				m_Acol[index + 2] = realIndex - step_k/7;
				m_Acol[index + 3] = realIndex;
				m_Acol[index + 4] = realIndex + step_k/7;
				m_Acol[index + 5] = realIndex + step_i/7;
				m_Acol[index + 6] = realIndex + step_j/7;

				m_b[realIndex] = 0;
			}
		}
	}
	ApplyBoundaryConditions();;
}

void PoissonSolver::ApplyBoundaryConditions() {
	int step_j = m_RadialNodeCount * m_AxialNodeCount * 7;
	int step_i = m_AxialNodeCount * 7;
	int step_k = 7;

	// i: r nodes.
	// j: theta nodes.
	// k: z nodes.

	//===============================================
	// z = 0 surface:
	int k = 0;
	for (int j = 0; j < m_ThetaNodeCount; j++) {
		for (int i = 0; i < m_RadialNodeCount; i++) {
			int index = j * step_j + i * step_i + k * step_k;
			int realIndex = index / 7;

			m_A[index + 0] = 0;			// U_{i,	j-1,	k}
			m_A[index + 1] = 0;			// U_{i-1,	j,		k}
			m_A[index + 2] = 0;			// U_{i,	j,		k-1}
			m_A[index + 3] = 1;			// U_{i,	j,		k}
			m_A[index + 4] = 0;			// U_{i,	j,		k+1}
			m_A[index + 5] = 0;			// U_{i+1,	j,		k}
			m_A[index + 6] = 0;			// U_{i,	j+1,	k}

			m_b[realIndex] = m_VDischarge;
		}
	}
	//===============================================
	// z = m_AxialLength surface:
	k = m_AxialNodeCount - 1; // -1 because, indexing starts from 0.
	for (int j = 0; j < m_ThetaNodeCount; j++) {
		for (int i = 0; i < m_RadialNodeCount; i++) {
			int index = j * step_j + i * step_i + k * step_k;
			int realIndex = index / 7;

			m_A[index + 0] = 0;			// U_{i,	j-1,	k}
			m_A[index + 1] = 0;			// U_{i-1,	j,		k}
			m_A[index + 2] = 0;			// U_{i,	j,		k-1}
			m_A[index + 3] = 1;			// U_{i,	j,		k}
			m_A[index + 4] = 0;			// U_{i,	j,		k+1}
			m_A[index + 5] = 0;			// U_{i+1,	j,		k}
			m_A[index + 6] = 0;			// U_{i,	j+1,	k}

			m_b[realIndex] = m_VPlume;
		}
	}
	//===============================================
	// theta = 0 surface:
	int j = 0;
	for (int i = 0; i < m_RadialNodeCount; i++) {
		for (int k = 0; k < m_AxialNodeCount; k++) {
			int index = j * step_j + i * step_i + k * step_k;
			int realIndex = index / 7;

			m_A[index + 0] = 0;			// U_{i,	j-1,	k}
			m_A[index + 1] = 0;			// U_{i-1,	j,		k}
			m_A[index + 2] = 0;			// U_{i,	j,		k-1}
			m_A[index + 3] = -1;		// U_{i,	j,		k}
			m_A[index + 4] = 0;			// U_{i,	j,		k+1}
			m_A[index + 5] = 0;			// U_{i+1,	j,		k}
			m_A[index + 6] = 1;			// U_{i,	j+1,	k}

			m_b[realIndex] = 0;
		}
	}
	//===============================================
	// theta = m_ThetaLength surface:
	j = m_ThetaNodeCount - 1;
	for (int i = 0; i < m_RadialNodeCount; i++) {
		for (int k = 0; k < m_AxialNodeCount; k++) {
			int index = j * step_j + i * step_i + k * step_k;
			int realIndex = index / 7;

			m_A[index + 0] = -1;		// U_{i,	j-1,	k}
			m_A[index + 1] = 0;			// U_{i-1,	j,		k}
			m_A[index + 2] = 0;			// U_{i,	j,		k-1}
			m_A[index + 3] = 1;			// U_{i,	j,		k}
			m_A[index + 4] = 0;			// U_{i,	j,		k+1}
			m_A[index + 5] = 0;			// U_{i+1,	j,		k}
			m_A[index + 6] = 0;			// U_{i,	j+1,	k}

			m_b[realIndex] = 0;
		}
	}
	//===============================================
	// r = 0 surface:
	int i = 0;
	for (int j = 0; j < m_ThetaNodeCount; j++) {
		for (int k = 0; k < m_AxialNodeCount; k++) {
			int index = j * step_j + i * step_i + k * step_k;
			int realIndex = index / 7;

			m_A[index + 0] = 0;			// U_{i,	j-1,	k}
			m_A[index + 1] = 0;			// U_{i-1,	j,		k}
			m_A[index + 2] = 0;			// U_{i,	j,		k-1}
			m_A[index + 3] = -1;		// U_{i,	j,		k}
			m_A[index + 4] = 0;			// U_{i,	j,		k+1}
			m_A[index + 5] = 1;			// U_{i+1,	j,		k}
			m_A[index + 6] = 0;			// U_{i,	j+1,	k}

			m_b[realIndex] = 0;
		}
	}
	//===============================================
	// r = m_RadialLength surface:
	i = m_RadialNodeCount - 1; // -1 because, indexing starts from 0.
	for (int j = 0; j < m_ThetaNodeCount; j++) {
		for (int k = 0; k < m_AxialNodeCount; k++) {
			int index = j * step_j + i * step_i + k * step_k;
			int realIndex = index / 7;

			m_A[index + 0] = 0;			// U_{i,	j-1,	k}
			m_A[index + 1] = -1;		// U_{i-1,	j,		k}
			m_A[index + 2] = 0;			// U_{i,	j,		k-1}
			m_A[index + 3] = 1;			// U_{i,	j,		k}
			m_A[index + 4] = 0;			// U_{i,	j,		k+1}
			m_A[index + 5] = 0;			// U_{i+1,	j,		k}
			m_A[index + 6] = 0;			// U_{i,	j+1,	k}

			m_b[realIndex] = 0;
		}
	}
	//===============================================
	// Screen Grid:
	for (int j = 1; j < m_ThetaNodeCount - 1; j++) {
		for (int i = m_rScreenBeginNode; i < m_RadialNodeCount - 1; i++) {
			for (int k = m_zScreenBeginNode; k <= m_zScreenEndNode; k++) {
				int index = j * step_j + i * step_i + k * step_k;
				int realIndex = index / 7;

				m_A[index + 0] = 0;			// U_{i,	j-1,	k}
				m_A[index + 1] = 0;			// U_{i-1,	j,		k}
				m_A[index + 2] = 0;			// U_{i,	j,		k-1}
				m_A[index + 3] = 1;			// U_{i,	j,		k}
				m_A[index + 4] = 0;			// U_{i,	j,		k+1}
				m_A[index + 5] = 0;			// U_{i+1,	j,		k}
				m_A[index + 6] = 0;			// U_{i,	j+1,	k}

				m_b[realIndex] = m_VScreen;
			}
		}
	}
	//===============================================
	// Acceleration Grid:
	for (int j = 1; j < m_ThetaNodeCount - 1; j++) {
		for (int i = m_rScreenBeginNode; i < m_RadialNodeCount - 1; i++) {
			for (int k = m_zScreenBeginNode; k <= m_zScreenEndNode; k++) {
				int index = j * step_j + i * step_i + k * step_k;
				int realIndex = index / 7;

				m_A[index + 0] = 0;			// U_{i,	j-1,	k}
				m_A[index + 1] = 0;			// U_{i-1,	j,		k}
				m_A[index + 2] = 0;			// U_{i,	j,		k-1}
				m_A[index + 3] = 1;			// U_{i,	j,		k}
				m_A[index + 4] = 0;			// U_{i,	j,		k+1}
				m_A[index + 5] = 0;			// U_{i+1,	j,		k}
				m_A[index + 6] = 0;			// U_{i,	j+1,	k}

				m_b[realIndex] = m_VAccel;
			}
		}
	}
}

void PoissonSolver::OutputMatrixA() {
	std::ofstream file("PoissonMatrixA.txt");

	for (int i = 0; i < m_n; i++) {
		file << m_x[i] << std::endl;
	}
	file.close();
}

void PoissonSolver::SolvePoisson(int MaxIterations, int KSD, double AbsoluteTolerance) {

}