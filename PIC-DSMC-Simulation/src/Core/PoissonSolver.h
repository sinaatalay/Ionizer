#pragma once
#include "Geometry.h"

class PoissonSolver : public Geometry {
public:
	explicit PoissonSolver();
	~PoissonSolver();

	void ConfigureCoefficients();

private:
	double* m_A;
};

