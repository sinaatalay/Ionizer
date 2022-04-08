#include "Geometry.h"
#include "Log.h"

Geometry::Geometry() {
	m_ThetaBegin = 0.0;
	m_AxialBegin = 0.0;
	m_RadialBegin = 0.0;
}
Geometry::~Geometry() {}

//===============================================
void Geometry::Setdtheta(double dtheta) { m_dtheta = dtheta; }
void Geometry::Setdr(double dr) { m_dr = dr; }
void Geometry::Setdz(double dz) { m_dz = dz; }
//===============================================
void Geometry::SetThetaLength(double length) { m_ThetaLength = length; }
void Geometry::SetAxialLength(double length) { m_AxialLength = length; }
void Geometry::SetRadialLength(double length) { m_RadialLength = length; }
//===============================================
void Geometry::SetScreenGridWidth(double length) { m_wScreen = length; }
void Geometry::SetScreenGridRadius(double length) { m_rScreen = length; }
//===============================================
void Geometry::SetAccelGridWidth(double length) { m_wAccel = length; }
void Geometry::SetAccelGridRadius(double length) { m_rAccel = length; }
//===============================================
void Geometry::SetAxialDischargeLength(double length) { m_zDischarge = length; }
void Geometry::SetDistanceBetweenGrids(double length) { m_zDistance = length; };
//===============================================
//===============================================
void Geometry::Calculation() {
	m_ThetaEnd = m_ThetaBegin + m_ThetaLength;
	m_RadialEnd = m_RadialBegin + m_RadialLength;
	m_AxialEnd = m_AxialBegin + m_AxialLength;
	//===============================================
	//===============================================
	//===============================================
	// FIXING TOTAL LENGTH PROBLEMS:
	double epsilon = pow(10, -6);

	m_RadialNodeCount = (int)(m_RadialLength / m_dr + 0.5 + 1);
	double a = (m_RadialNodeCount-1) * m_dr;
	if (fabs(a - m_RadialLength) > epsilon) {
		LOG_WARN("(NODE COUNT PROBLEM) Whole radial length has changed from {:.8} m to {:.8} m", m_RadialLength, a);
		m_RadialLength = a;
	}

	m_AxialNodeCount = (int)(m_AxialLength / m_dz + 0.5 + 1);
	a = (m_AxialNodeCount - 1) * m_dz;
	if (fabs(a - m_AxialLength) > epsilon) {
		LOG_WARN("(NODE COUNT PROBLEM) Whole axial length has changed from {:.8} m to {:.8} m", m_AxialLength, a);
		m_AxialLength = a;
	}

	m_ThetaNodeCount = (int)(m_ThetaLength / m_dtheta + 0.5 + 1);
	a = (m_ThetaNodeCount - 1) * m_dtheta;
	if (fabs(a - m_ThetaLength) > epsilon) {
		LOG_WARN("(NODE COUNT PROBLEM) Whole azimuthal length has changed from {:.8} rad to {:.8} rad", m_ThetaLength, a);
		m_ThetaLength = a;
	}

	m_TotalNodeCount = m_AxialNodeCount + m_RadialNodeCount + m_ThetaNodeCount;
	//===============================================
	//===============================================
	//===============================================
	// FIXING RADIAL PROBLEMS:
	m_rScreenBeginNode = (int)(m_rScreen / m_dr + 0.5);
	a = m_rScreenBeginNode * m_dr;
	if (fabs(a - m_rScreen) > epsilon) {
		LOG_WARN("(NODE COUNT PROBLEM) Radius of the screen grid has changed from {:.8} m to {:.8} m", m_rScreen, a);
		m_rScreen = a;
	}

	m_rAccelBeginNode = (int)(m_rAccel / m_dr + 0.5);
	a = m_rAccelBeginNode * m_dr;
	if (fabs(a - m_rAccel) > epsilon) {
		LOG_WARN("(NODE COUNT PROBLEM) Radius of the acceleration grid has changed from {:.8} m to {:.8} m", m_rAccel, a);
		m_rAccel = a;
	}
	//===============================================
	//===============================================
	//===============================================
	// FIXING AXIAL PROBLEMS:
	m_zScreenBeginNode = (int)((m_zDischarge) / m_dz + 0.5);
	a = m_zScreenBeginNode * m_dz;
	if (fabs(a - m_zDischarge) > epsilon) {
		LOG_WARN("(NODE COUNT PROBLEM) Axial length of discharge region has changed from {:.8} m to {:.8} m", m_zDischarge, a);
		m_zDischarge = a;
	}

	m_zScreenEndNode = (int)((m_zDischarge + m_wScreen) / m_dz + 0.5);
	a = (m_zScreenEndNode - m_zScreenBeginNode) * m_dz;
	if (fabs(a - m_wScreen) > epsilon) {
		LOG_WARN("(NODE COUNT PROBLEM) Width of the screen grid has changed from {:.8} m to {:.8} m", m_wScreen, a);
		m_wScreen = a;
	}

	m_zAccelBeginNode = (int)((m_zDischarge + m_wScreen + m_zDistance) / m_dz + 0.5);
	a = (m_zAccelBeginNode - m_zScreenEndNode) * m_dz;
	if ((a - m_zDistance) > epsilon) {
		LOG_WARN("(NODE COUNT PROBLEM) The distance between the grids has changed from {:.8} m to {:.8} m", m_zDistance, a);
		m_zDistance = a;
	}

	m_zAccelEndNode = (int)((m_zDischarge + m_wScreen + m_zDistance + m_wAccel) / m_dz + 0.5);
	a = (m_zAccelEndNode - m_zAccelBeginNode) * m_dz;
	if ((a - m_wAccel) > epsilon) {
		LOG_WARN("(NODE COUNT PROBLEM) Width of the acceleration grid has changed from {:.8} m to {:.8} m", m_wAccel, a);
		m_wAccel = a;
	}

	m_zPlume = m_AxialLength - (m_zDischarge + m_zDistance + m_wAccel + m_wScreen);
}

void Geometry::SetGeometry(const Geometry& geometry) {
	m_dtheta = geometry.Getdtheta();
	m_dr = geometry.Getdr();
	m_dz = geometry.Getdz();

	m_ThetaLength = geometry.GetThetaLength();
	m_RadialLength = geometry.GetRadialLength();
	m_AxialLength = geometry.GetAxialLength();

	m_wScreen = geometry.GetwScreen();
	m_rScreen = geometry.GetrScreen();

	m_wAccel = geometry.GetwAccel();
	m_rAccel = geometry.GetrAccel();

	m_zDischarge = geometry.GetzDischarge();
	m_zDistance = geometry.GetzDistance();
	Calculation();
}

void Geometry::LogGeometry() const {
	LOG_INFO("Whole azimuthal length: {:.5} rad", m_ThetaLength);	// [rad] Whole azimuthal length.
	LOG_INFO("Whole axial length: {:.5} m", m_AxialLength);	// [m] Whole axial length.
	LOG_INFO("Whole radial length: {:.5} m", m_RadialLength);  // [m] Whole radial length.
	//===============================================
	LOG_INFO("Theta node count: {}", m_ThetaNodeCount);
	LOG_INFO("Axial node count: {}", m_AxialNodeCount);
	LOG_INFO("Radial node count: {}", m_RadialNodeCount);
	LOG_INFO("Total node count: {}", m_TotalNodeCount);
	//===============================================
	LOG_INFO("dr: {:.5} m", m_dr);
	LOG_INFO("dz: {:.5} m", m_dz);
	LOG_INFO("dtheta: {:.5} rad", m_dtheta);
	//===============================================
	LOG_INFO("Beginning of the azimuthal domain: {:.5} rad", m_ThetaBegin);
	LOG_INFO("End of the azimuthal domain: {:.5} rad", m_ThetaEnd);
	LOG_INFO("Beginning of the radial domain: {:.5} m", m_RadialBegin);
	LOG_INFO("End of the radial domain: {:.5} m", m_RadialEnd);
	LOG_INFO("Beginning of the axial domain: {:.5} m", m_AxialBegin);
	LOG_INFO("End of the axial domain: {:.5} m", m_AxialEnd);
	//===============================================
	LOG_INFO("Portion of the domain up to the screen grid: {:.5} m", m_zDischarge);
	LOG_INFO("Distance betweeen the screen and accel grid: {:.5} m", m_zDistance);
	LOG_INFO("Portion of the domain that remains in the plume region: {:.5} m", m_zPlume);
	//===============================================
	LOG_INFO("Width of the screen grid: {:.5} m", m_wScreen);
	LOG_INFO("Radius of the screen grid: {:.5} m", m_rScreen);
	//===============================================
	LOG_INFO("Width of the acceleration grid: {:.5} m", m_wAccel);
	LOG_INFO("Radius of the acceleration grid: {:.5} m", m_rAccel);
	//===============================================
	LOG_INFO("Start node of the acceleration grid in radial direction: {}", m_rAccelBeginNode);
	LOG_INFO("Start node of the acceleration grid in axial direction: {}", m_zAccelBeginNode);
	LOG_INFO("End node of the acceleration grid in axial direction: {}", m_zAccelEndNode);
	//===============================================
	LOG_INFO("Start node of the screen grid in radial direction: {}", m_rScreenBeginNode);
	LOG_INFO("Start node of the screen grid in axial direction: {}", m_zScreenBeginNode);
	LOG_INFO("End node of the screen grid in axial direction: {}", m_zScreenEndNode);
}