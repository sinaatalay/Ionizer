#pragma once
//C++ Libraries (Precompiled):
#include "Ionpch.h"

#include "IonThrusterGeometry.h"
#include "SparseSystemSolver.h"


namespace Ionizer {

	class PoissonSolver : public IonThrusterGeometry {
	public:
		PoissonSolver();
		explicit PoissonSolver(const IonThrusterGeometry& geometry);
		~PoissonSolver();

		// ConfigureFDSystem() calculates the coefficients of Poisson's equation's Finite Difference Method system, Ax=b, in cylindrical coordinates:
		void ConfigureFDSystem();

		// SolvePoisson() solves the linear system Ax=b created with the Finite Difference Method by using SparseSystemSolver class:
		void SolvePoisson();

		// OutputSolution() outputs the solution vector x to a *.txt file.
		void OutputSolution(const std::string& FileName);

		// GetImage() returns an image of the solution as a std::vector<double>.
		std::vector<uint32_t> GetImage(uint32_t ViewportWidth, uint32_t ViewportHeight, float angle);
	private:
		int m_n;						// The order of the linear system
		SparseSystemSolver m_Solver;	// The sparse linear system solver
	};
}