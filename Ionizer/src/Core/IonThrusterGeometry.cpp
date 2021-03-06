//C++ Libraries (Precompiled):
#include "Ionpch.h"

#include "IonThrusterGeometry.h"
#include "Log.h"

namespace Ionizer {

	IonThrusterGeometry::IonThrusterGeometry() {}
	IonThrusterGeometry::IonThrusterGeometry(const IonThrusterGeometry& geometry) {
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
	}
	IonThrusterGeometry::~IonThrusterGeometry() {}

	void IonThrusterGeometry::Setdtheta(double dtheta) { m_dtheta = dtheta; }
	void IonThrusterGeometry::Setdr(double dr) { m_dr = dr; }
	void IonThrusterGeometry::Setdz(double dz) { m_dz = dz; }

	void IonThrusterGeometry::SetThetaLength(double length) { m_ThetaLength = length; }
	void IonThrusterGeometry::SetAxialLength(double length) { m_AxialLength = length; }
	void IonThrusterGeometry::SetRadialLength(double length) { m_RadialLength = length; }

	void IonThrusterGeometry::SetScreenGridWidth(double length) { m_wScreen = length; }
	void IonThrusterGeometry::SetScreenGridRadius(double length) { m_rScreen = length; }
	void IonThrusterGeometry::SetScreenGridVoltage(double v) { m_VScreen = v; }

	void IonThrusterGeometry::SetAccelGridWidth(double length) { m_wAccel = length; }
	void IonThrusterGeometry::SetAccelGridRadius(double length) { m_rAccel = length; }
	void IonThrusterGeometry::SetAccelGridVoltage(double v) { m_VAccel = v; }

	void IonThrusterGeometry::SetAxialDischargeLength(double length) { m_AxialDischargeLength = length; }
	void IonThrusterGeometry::SetDistanceBetweenGrids(double length) { m_AxialDistanceBetweenGrids = length; };

	void IonThrusterGeometry::SetVDischarge(double v) { m_VDischarge = v; }
	void IonThrusterGeometry::SetVPlume(double v) { m_VPlume = v; }

	void IonThrusterGeometry::GeometryCalculation() {
		m_ThetaEnd = m_ThetaBegin + m_ThetaLength;
		m_RadialEnd = m_RadialBegin + m_RadialLength;
		m_AxialEnd = m_AxialBegin + m_AxialLength;
		//===============================================
		//===============================================
		//===============================================
		// FIXING TOTAL LENGTH PROBLEMS:
		double epsilon = std::pow(10, -6);

		m_RadialNodeCount = (int)(m_RadialLength / m_dr + 0.5 + 1);
		double a = (m_RadialNodeCount - 1) * m_dr;
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

		m_TotalNodeCount = m_AxialNodeCount * m_RadialNodeCount * m_ThetaNodeCount;
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
		m_zScreenBeginNode = (int)((m_AxialDischargeLength) / m_dz + 0.5);
		a = m_zScreenBeginNode * m_dz;
		if (fabs(a - m_AxialDischargeLength) > epsilon) {
			LOG_WARN("(NODE COUNT PROBLEM) Axial length of discharge region has changed from {:.8} m to {:.8} m", m_AxialDischargeLength, a);
			m_AxialDischargeLength = a;
		}

		m_zScreenEndNode = (int)((m_AxialDischargeLength + m_wScreen) / m_dz + 0.5);
		a = (m_zScreenEndNode - m_zScreenBeginNode) * m_dz;
		if (fabs(a - m_wScreen) > epsilon) {
			LOG_WARN("(NODE COUNT PROBLEM) Width of the screen grid has changed from {:.8} m to {:.8} m", m_wScreen, a);
			m_wScreen = a;
		}

		m_zAccelBeginNode = (int)((m_AxialDischargeLength + m_wScreen + m_AxialDistanceBetweenGrids) / m_dz + 0.5);
		a = (m_zAccelBeginNode - m_zScreenEndNode) * m_dz;
		if ((a - m_AxialDistanceBetweenGrids) > epsilon) {
			LOG_WARN("(NODE COUNT PROBLEM) The distance between the grids has changed from {:.8} m to {:.8} m", m_AxialDistanceBetweenGrids, a);
			m_AxialDistanceBetweenGrids = a;
		}

		m_zAccelEndNode = (int)((m_AxialDischargeLength + m_wScreen + m_AxialDistanceBetweenGrids + m_wAccel) / m_dz + 0.5);
		a = (m_zAccelEndNode - m_zAccelBeginNode) * m_dz;
		if ((a - m_wAccel) > epsilon) {
			LOG_WARN("(NODE COUNT PROBLEM) Width of the acceleration grid has changed from {:.8} m to {:.8} m", m_wAccel, a);
			m_wAccel = a;
		}

		m_zPlume = m_AxialLength - (m_AxialDischargeLength + m_AxialDistanceBetweenGrids + m_wAccel + m_wScreen);
	}
	void IonThrusterGeometry::LogGeometry() const {
		LOG_INFO("Whole azimuthal length: {:.5} rad", m_ThetaLength);
		LOG_INFO("Whole axial length: {:.5} m", m_AxialLength);
		LOG_INFO("Whole radial length: {:.5} m", m_RadialLength);
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
		LOG_INFO("Portion of the domain up to the screen grid: {:.5} m", m_AxialDischargeLength);
		LOG_INFO("Voltage of the discharge: {:.5} V", m_VDischarge);
		LOG_INFO("Distance betweeen the screen and accel grid: {:.5} m", m_AxialDistanceBetweenGrids);
		LOG_INFO("Portion of the domain that remains in the plume region: {:.5} m", m_zPlume);
		LOG_INFO("Voltage of the plume: {:.5} V", m_VPlume);
		//===============================================
		LOG_INFO("Width of the screen grid: {:.5} m", m_wScreen);
		LOG_INFO("Radius of the screen grid: {:.5} m", m_rScreen);
		LOG_INFO("Voltage of the screen grid: {:.5} V", m_VScreen);
		//===============================================
		LOG_INFO("Width of the acceleration grid: {:.5} m", m_wAccel);
		LOG_INFO("Radius of the acceleration grid: {:.5} m", m_rAccel);
		LOG_INFO("Voltage of the acceleration grid: {:.5} V", m_VAccel);
		//===============================================
		LOG_INFO("Start node of the acceleration grid in radial direction: {}", m_rAccelBeginNode);
		LOG_INFO("Start node of the acceleration grid in axial direction: {}", m_zAccelBeginNode);
		LOG_INFO("End node of the acceleration grid in axial direction: {}", m_zAccelEndNode);
		//===============================================
		LOG_INFO("Start node of the screen grid in radial direction: {}", m_rScreenBeginNode);
		LOG_INFO("Start node of the screen grid in axial direction: {}", m_zScreenBeginNode);
		LOG_INFO("End node of the screen grid in axial direction: {}", m_zScreenEndNode);
	}

}