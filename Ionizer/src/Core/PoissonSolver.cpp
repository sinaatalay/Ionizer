//C++ Libraries (Precompiled):
#include "Ionpch.h"

#include "PoissonSolver.h"
#include "Log.h"
#include "Timer.h"

namespace Ionizer {

	PoissonSolver::PoissonSolver() {
		// Setting up the Ion Thruster Geometry:
		double pi = 3.141592653589793;
		m_dtheta = pi / 36;
		m_dr = 0.00002;
		m_dz = 0.00002;

		// Physical lengths of the domain:
		m_ThetaLength = pi / 6;	   // [rad] Whole azimuthal length
		m_RadialLength = 0.002;	   // [m] Whole axial length.
		m_AxialLength = 0.005;	   // [m] Whole radial length.

		// Physical lengths of the thruster:
		m_wScreen = 0.0004;						 // [m]
		m_rScreen = 0.001;						 // [m]
		m_wAccel = 0.0008;						 // [m]
		m_rAccel = 0.0006;						 // [m]
		m_AxialDischargeLength = 0.001;			 // [m]
		m_AxialDistanceBetweenGrids = 0.0012;	 // [m]

		// Potentials:
		m_VDischarge = 2266;	// [V] Bulk plasma potential
		m_VPlume = 0;			// [V] Plume plasma potential
		m_VAccel = -400;		// [V] Acceleration grid (second grid) 
		m_VScreen = 2241;		// [V] Screen grid (first grid) potential

		GeometryCalculation();

		m_n = m_TotalNodeCount;
		m_Solver.SetOrder(m_n);
	}
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

		m_AxialDischargeLength = geometry.GetzDischarge();
		m_AxialDistanceBetweenGrids = geometry.GetzDistance();

		m_VDischarge = geometry.GetVDischarge();
		m_VPlume = geometry.GetVPlume();

		GeometryCalculation();

		m_n = m_TotalNodeCount;
		m_Solver.SetOrder(m_n);
	}
	PoissonSolver::~PoissonSolver() {}

	void PoissonSolver::ConfigureFDSystem() {
		enum MatrixSelection : uint8_t {
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
					else if (k >= m_zScreenBeginNode && k <= m_zScreenEndNode && i >= m_rScreenBeginNode) {
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
						m_Solver.InsertCoefficient(A, row, row + step_i, 1);
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
					else {
						// Interior region:
						double r = m_RadialBegin + i * m_dr;

						double jm1 = 1.0 / (r * r * m_dtheta * m_dtheta);												// U_{i,	j-1,	k}
						double im1 = 1.0 / (m_dr * m_dr) - 1.0 / (r * 2.0 * m_dr);										// U_{i-1,	j,		k}
						double km1 = 1.0 / (m_dz * m_dz);																// U_{i,	j,		k-1}
						double ijk = -2.0 / (m_dz * m_dz) - 2.0 / (m_dr * m_dr) - 2.0 / (r * r * m_dtheta * m_dtheta);	// U_{i,	j,		k}
						double kp1 = km1;																				// U_{i,	j,		k+1}
						double ip1 = 1.0 / (m_dr * m_dr) + 1.0 / (r * 2.0 * m_dr);										// U_{i+1,	j,		k}
						double jp1 = jm1;																				// U_{i,	j+1,	k}

						double InvMagnitude = 1 / std::sqrt(jm1 * jm1 + im1 * im1 + km1 * km1 + ijk * ijk + kp1 * kp1 + ip1 * ip1 + jp1 * jp1);

						jm1 = jm1 * InvMagnitude;
						im1 = im1 * InvMagnitude;
						km1 = km1 * InvMagnitude;
						ijk = ijk * InvMagnitude;
						kp1 = kp1 * InvMagnitude;
						ip1 = ip1 * InvMagnitude;
						jp1 = jp1 * InvMagnitude;

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

		timer.Reset();
		m_Solver.Solve();
	}

	void PoissonSolver::OutputSolution(const std::string& FileName) {
		m_Solver.OutputSolution(FileName);
	}

	std::vector<uint32_t> PoissonSolver::GetImage(uint32_t ViewportWidth, uint32_t ViewportHeight) {
		std::vector<double> solution = m_Solver.GetSolution();
		int step_j = m_RadialNodeCount * m_AxialNodeCount;

		double theta = 0;
		int ThetaNode = (theta / m_dtheta + 0.5) + 1;
		solution.erase(solution.begin(), solution.begin() + step_j * (ThetaNode - 1));
		solution.erase(solution.end() - step_j * (m_ThetaNodeCount - ThetaNode), solution.end());


		std::vector<uint32_t> Image;
		Image.resize(step_j);

		int red;

		for (int height = 0; height < m_RadialNodeCount; height++) {
			for (int width = 0; width < m_AxialNodeCount; width++) {
				int heightSol = m_RadialNodeCount - height - 1;
				int widthSol = width;
				double a = solution[heightSol * m_AxialNodeCount + widthSol];
				double b = m_VAccel;
				double epsilon = std::pow(10, -6);
				bool test1 = fabs(a - b) <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
				b = m_VScreen;
				bool test2 = fabs(a - b) <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
				if (test1 || test2) {
					Image[width + height * m_AxialNodeCount] = 0xff575757;
				}
				else {
					red = (a + 405) / (m_VDischarge - m_VAccel + 10) * 255;
					Image[width + height * m_AxialNodeCount] = (red << 24) + ((0 & 0xff) << 16) + ((0 & 0xff) << 8) + (255 & 0xff);
				}
			}
		}

		std::vector<uint32_t> ImageFinal;
		double ratio = ViewportWidth / (double)m_AxialNodeCount;
		uint32_t HeightFinal= m_RadialNodeCount * ratio;
		ImageFinal.resize(ViewportWidth * ViewportHeight);

		for (int height = 0; height < ViewportHeight; height++) {
			for (int width = 0; width < ViewportWidth; width++) {
				if (width / ratio + height / ratio * m_AxialNodeCount > 25350) {
					int a = 5;
				}
				if (height < HeightFinal) {
					ImageFinal[width + height * ViewportWidth] = Image[ (int)(width / ratio) + (int)(height / ratio) * m_AxialNodeCount];
				}
				else {
					ImageFinal[width + height * ViewportWidth] = 0x00000000;
				}

			}
		}

		return ImageFinal;
	}
}