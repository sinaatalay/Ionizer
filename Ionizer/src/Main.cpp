//C++ Libraries (Precompiled):
#include "Ionpch.h"

//Core:
#include "Core/IonThrusterGeometry.h"
#include "Core/Log.h"
#include "Core/PoissonSolver.h"

using namespace Ionizer;

int main() {
	LOG_INIT(LOG_LEVEL_ALL);
	LOG_INFO("Welcome to PIC-DSMC Simulation!");

	constexpr double pi = 3.141592653589793; //Compile time constant.

	// Setting up the Ion Thruster Geometry:

	// Try to use multiples of the relevant geometry dimensions.
	// Otherwise the program will change the geometry accordingly.
	double dz = 0.00002;
	double dr = 0.00002;
	double dtheta = pi / 36;

	// Physical lengths of the domain:
	double ThetaLength = pi / 6;	// [rad] Whole azimuthal length
	double AxialLength = 0.005;		// [m] Whole axial length.
	double RadialLength = 0.002;	// [m] Whole radial length.

	// Physical lengths of the thruster:
	double ScreenGridWidth = 0.0004;			// [m]
	double ScreenGridHoleRadius = 0.001;		// [m]

	double AccelerationGridWidth = 0.0008;		// [m]
	double AccelerationGridHoleRadius = 0.0006;	// [m]

	double DischageRegionAxialLength = 0.001;	// [m]
	double DistanceBetweenGrids = 0.0012;		// [m]

	// Potentials:
	double V_Discharge = 2266;		 // [V] Bulk plasma potential
	double V_Plume = 0;				 // [V] Plume plasma potential
	double V_Accel = -400;			 // [V] Acceleration grid (second grid) potential
	double V_Screen = 2241;			 // [V] Screen grid (first grid) potential

	IonThrusterGeometry geometry;

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

	// Setting up the poisson solver:

	PoissonSolver poisson(geometry);

	poisson.LogGeometry();

	poisson.SolvePoisson();
	poisson.OutputSolution("PoissonSolution1.txt");

	LOG_INFO("Have a nice day!");

	std::cin.get();
}