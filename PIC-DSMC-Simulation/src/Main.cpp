constexpr double PI = 3.14159265;

//===C++ Libraries===============
#include <iostream>
#include <cmath>

//===Vendor======================
#include "Core/Log.h" // Logger from https://github.com/gabime/spdlog

//===Core========================
#include "Core/Geometry.h"
#include "Core/PoissonSolver.h"

int main() {
	Log::Init(LOG_LEVEL_ALL);
	LOG_INFO("Welcome to PIC-DSMC Simulation!");

	//==============================================================================
	//==============================================================================
	//==============================================================================
	//==============================================================================

	// Try to use multiples of the relevant geometry dimensions.
	// Otherwise the program will change the geometry accordingly.
	double dz = 0.00002;
	double dr = 0.00002;
	double dtheta = PI / 12;

	// Physical lengths of the domain.
	double ThetaLength = PI / 6;	// [rad] Whole azimuthal length.
	double AxialLength = 0.005;		// [m] Whole axial length.
	double RadialLength = 0.002;	// [m] Whole radial length.

	// Physical lengths of the thruster.
	double ScreenGridWidth = 0.0004;			// [m]
	double ScreenGridHoleRadius = 0.001;		// [m]

	double AccelerationGridWidth = 0.0008;		// [m]
	double AccelerationGridHoleRadius = 0.0006;	// [m]

	double DischageRegionAxialLength = 0.001;	// [m]
	double DistanceBetweenGrids = 0.0012;		// [m]

	// Potentials
	double V_Discharge = 2266;		 // Bulk plasma potential 
	double V_Plume = 0;				 // Plume plasma potential
	double V_Accel = -400;			 // Acceleration grid (second grid) potential 
	double V_Screen = 2241;			 // Screen grid (first grid) potential

	// Poisson solver settings
	int MaxIterations = 100;			// Maximum number of outer GMRES iterations to take
	int KrylovSubspaceDimension = 850;	// Maximum number of inner GMRES iterations to take (can't be bigger than the order of the linear system)
	double AbsoluteTolerance = 0.001;	// Absolute tolerance for GMRES

	//==============================================================================
	//==============================================================================
	//==============================================================================
	//==============================================================================

	Geometry geometry = Geometry();

	geometry.Setdz(dz);
	geometry.Setdr(dr);
	geometry.Setdtheta(dtheta);

	geometry.SetThetaLength(ThetaLength);
	geometry.SetAxialLength(AxialLength);
	geometry.SetRadialLength(RadialLength);

	geometry.SetScreenGridWidth(ScreenGridWidth);
	geometry.SetScreenGridRadius(ScreenGridHoleRadius);
	geometry.SetScreenGridVoltage(V_Screen);

	geometry.SetAccelGridWidth(AccelerationGridWidth);
	geometry.SetAccelGridRadius(AccelerationGridHoleRadius);
	geometry.SetAccelGridVoltage(V_Accel);

	geometry.SetAxialDischargeLength(DischageRegionAxialLength);
	geometry.SetDistanceBetweenGrids(DistanceBetweenGrids);

	geometry.SetVDischarge(V_Discharge);
	geometry.SetVPlume(V_Plume);

	//==============================================================================
	//==============================================================================
	//==============================================================================
	//==============================================================================

	PoissonSolver Poisson = PoissonSolver();
	Poisson.SetGeometry(geometry);
	Poisson.LogGeometry();

	Poisson.AllocateMemory();

	Poisson.ConfigureMatrixA();

	Poisson.SolvePoisson(MaxIterations, KrylovSubspaceDimension, AbsoluteTolerance);

	LOG_INFO("Have a nice day!");
	std::cin.get();
}