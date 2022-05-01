#include <cmath>

#include "PoissonSolver.h"
#include "Log.h"
#include "Timer.h"

namespace picdsmc {

	PoissonSolver::PoissonSolver() {}
	PoissonSolver::PoissonSolver(const IonThrusterGeometry& geometry) {
		m_dtheta = geometry.Getdtheta();
		m_dr = geometry.Getdr();
		m_dz = geometry.Getdz();

		m_ThetaLength = geometry.GetThetaLength();
		m_RadialLength = geometry.GetRadialLength();
		m_AxialLength = geometry.GetAxialLength();

		m_wScreen = geometry.GetwScreen();
		m_rScreen = geometry.GetrScreen();
		m_VScreen = geometry.GetVScreen();

		m_wAccel = geometry.GetwAccel();
		m_rAccel = geometry.GetrAccel();
		m_VAccel = geometry.GetVAccel();

		m_zDischarge = geometry.GetzDischarge();
		m_zDistance = geometry.GetzDistance();

		m_VDischarge = geometry.GetVDischarge();
		m_VPlume = geometry.GetVPlume();

		GeometryCalculation();

		m_n = m_TotalNodeCount;
		m_Solver.SetOrder(m_n);
	}
	PoissonSolver::~PoissonSolver() {}

	double PoissonSolver::Normalizer(int i){
		return -1;
	}
	void PoissonSolver::ConfigureFDSystem() {
		enum MatrixSelection : unsigned char {
			A = 0, b = 1
		};

		// i: r nodes.
		// j: theta nodes.
		// k: z nodes.

		int step_j = m_RadialNodeCount * m_AxialNodeCount;
		int step_i = m_AxialNodeCount;

		for (int j = 0; j < m_ThetaNodeCount; j++) {
			for (int i = 0; i < m_RadialNodeCount; i++) {
				for (int k = 0; k < m_AxialNodeCount; k++) {
					int row = j * step_j + i * step_i + k;

					if (k == 0) {
						// z = 0 surface:
						m_Solver.InsertCoefficient(A, row, row, 1);
						m_Solver.InsertCoefficient(b, row, row, m_VDischarge);
					}
					else if (k == m_AxialNodeCount - 1) {
						// z = m_AxialLength surface:
						m_Solver.InsertCoefficient(A, row, row, 1);
						m_Solver.InsertCoefficient(b, row, row, m_VPlume);
					}
					else if (k >= m_zScreenBeginNode && k <= m_zScreenEndNode && i>=m_rScreenBeginNode) {
						// Screen Grid:
						m_Solver.InsertCoefficient(A, row, row, 1);
						m_Solver.InsertCoefficient(b, row, row, m_VScreen);
					}
					else if (k >= m_zAccelBeginNode && k <= m_zAccelEndNode && i >= m_rAccelBeginNode) {
						// Acceleration Grid:
						m_Solver.InsertCoefficient(A, row, row, 1);
						m_Solver.InsertCoefficient(b, row, row, m_VAccel);
					}
					else if (i == 0) {
						// r = 0 surface:
						m_Solver.InsertCoefficient(A, row, row+step_i, 1);
						m_Solver.InsertCoefficient(A, row, row, -1);
					}
					else if (i == m_RadialNodeCount - 1) {
						// r = m_RadialLength surface:
						m_Solver.InsertCoefficient(A, row, row - step_i, -1);
						m_Solver.InsertCoefficient(A, row, row, 1);
					}
					else if (j == 0) {
						// theta = 0 surface:
						m_Solver.InsertCoefficient(A, row, row, -1);
						m_Solver.InsertCoefficient(A, row, row + step_j, 1);
					}
					else if (j == m_ThetaNodeCount - 1) {
						// theta = m_ThetaLength surface:
						m_Solver.InsertCoefficient(A, row, row - step_j, -1);
						m_Solver.InsertCoefficient(A, row, row, 1);
					}
					else{
						// Interior region:
						double r = m_RadialBegin + i * m_dr;

						double jm1 = 1.0 / (r * r * m_dtheta * m_dtheta);												// U_{i,	j-1,	k}
						double im1 = 1.0 / (m_dr * m_dr) - 1.0 / (r * 2.0 * m_dr);										// U_{i-1,	j,		k}
						double km1 = 1.0 / (m_dz * m_dz);																// U_{i,	j,		k-1}
						double ijk = -2.0 / (m_dz * m_dz) - 2.0 / (m_dr * m_dr) - 2.0 / (r * r * m_dtheta * m_dtheta);	// U_{i,	j,		k}
						double kp1 = km1;																				// U_{i,	j,		k+1}
						double ip1 = 1.0 / (m_dr * m_dr) + 1.0 / (r * 2.0 * m_dr);										// U_{i+1,	j,		k}
						double jp1 = jm1;																				// U_{i,	j+1,	k}

						double InvMagnitude = 1/std::sqrt(jm1 * jm1 + im1 * im1 + km1 * km1 + ijk * ijk + kp1 * kp1 + ip1 * ip1 + jp1 * jp1);

						jm1 = jm1*InvMagnitude;
						im1 = im1*InvMagnitude;
						km1 = km1*InvMagnitude;
						ijk = ijk*InvMagnitude;
						kp1 = kp1*InvMagnitude;
						ip1 = ip1*InvMagnitude;
						jp1 = jp1*InvMagnitude;

						int col = row - step_j;
						m_Solver.InsertCoefficient(A, row, col, jm1);

						col = row - step_i;
						m_Solver.InsertCoefficient(A, row, col, im1);

						col = row - 1;
						m_Solver.InsertCoefficient(A, row, col, km1);

						m_Solver.InsertCoefficient(A, row, row, ijk);

						col = row + 1;
						m_Solver.InsertCoefficient(A, row, col, kp1);

						col = row + step_i;
						m_Solver.InsertCoefficient(A, row, col, ip1);

						col = row + step_j;
						m_Solver.InsertCoefficient(A, row, col, jp1);
					}
				}
			}
		}
	}
	void PoissonSolver::SolvePoisson() {
		LOG_INFO("Configuring Finite Difference Method's linear system.");
		Timer timer;
		ConfigureFDSystem();
		LOG_INFO("Configured in {:.5} ms.", timer.ElapsedMillis());

		m_Solver.Solve();

		m_Solver.OutputSolution("Output");
	}

}