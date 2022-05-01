#pragma once
#include "IonThrusterGeometry.h"
#include "SparseSystemSolver.h"


namespace picdsmc {

	class PoissonSolver : public IonThrusterGeometry {
	public:
		PoissonSolver();
		explicit PoissonSolver(const IonThrusterGeometry& geometry);
		~PoissonSolver();

		void ConfigureFDSystem();

		void SolvePoisson();
	private:
		int m_n;						// The order of the linear system
		SparseSystemSolver m_Solver;	// The sparse linear system solver
	};
}