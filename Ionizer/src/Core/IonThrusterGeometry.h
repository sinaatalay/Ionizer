#pragma once

namespace picdsmc {

	class IonThrusterGeometry {
	public:
		IonThrusterGeometry();
		explicit IonThrusterGeometry(const IonThrusterGeometry& geometry);
		virtual ~IonThrusterGeometry();

		void Setdtheta(double dtheta);
		void Setdr(double dr);
		void Setdz(double dz);

		void SetThetaLength(double length);
		void SetRadialLength(double length);
		void SetAxialLength(double length);

		void SetScreenGridWidth(double length);
		void SetScreenGridRadius(double length);
		void SetScreenGridVoltage(double v);

		void SetAccelGridWidth(double length);
		void SetAccelGridRadius(double length);
		void SetAccelGridVoltage(double v);

		void SetAxialDischargeLength(double length);
		void SetDistanceBetweenGrids(double length);

		void SetVDischarge(double v);
		void SetVPlume(double v);

		double GetThetaLength() const { return m_ThetaLength; }
		double GetAxialLength() const { return m_AxialLength; }
		double GetRadialLength() const { return m_RadialLength; }

		int GetThetaNodeCount() const { return m_ThetaNodeCount; }
		int GetAxialNodeCount() const { return m_AxialNodeCount; }
		int GetRadialNodeCount() const { return m_RadialNodeCount; }
		int GetTotalNodeCount() const { return m_TotalNodeCount; }

		double Getdr() const { return m_dr; }
		double Getdz() const { return m_dz; }
		double Getdtheta() const { return m_dtheta; }

		double GetThetaBegin() const { return m_ThetaBegin; }
		double GetThetaEnd() const { return m_ThetaEnd; }
		double GetRadialBegin() const { return m_RadialBegin; }
		double GetRadialEnd() const { return m_RadialEnd; }
		double GetAxialBegin() const { return m_AxialBegin; }
		double GetAxialEnd() const { return m_AxialEnd; }

		double GetzDischarge() const { return m_zDischarge; }
		double GetzDistance() const { return m_zDistance; }
		double GetzPlume() const { return m_zPlume; }

		double GetwScreen() const { return m_wScreen; }
		double GetrScreen() const { return m_rScreen; }
		double GetVScreen() const { return m_VScreen; }

		double GetwAccel() const { return m_wAccel; }
		double GetrAccel() const { return m_rAccel; }
		double GetVAccel() const { return m_VAccel; }

		int GetrAccelBeginNode() const { return m_rAccelBeginNode; }
		int GetzAccelBeginNode() const { return m_zAccelBeginNode; }
		int GetzAccelEndNode() const { return m_zAccelEndNode; }

		int GetrScreenBeginNode() const { return m_rScreenBeginNode; }
		int GetzScreenBeginNode() const { return m_zScreenBeginNode; }
		int GetzScreenEndNode() const { return m_zScreenEndNode; }

		double GetVDischarge() const { return m_VDischarge; }
		double GetVPlume() const { return m_VPlume; }

		// GeometryCalculation() computes the other geometry dimensions derived from the set inputs and round them according to the node lengths:
		void GeometryCalculation();

		// LogGeometry() logs all the members to the console:
		void LogGeometry() const;
	protected:
		double m_ThetaLength;	// [rad] Whole azimuthal length
		double m_AxialLength;	// [m] Whole axial length
		double m_RadialLength;  // [m] Whole radial length

		int m_ThetaNodeCount;
		int m_AxialNodeCount;
		int m_RadialNodeCount;
		int m_TotalNodeCount;

		double m_dr, m_dz, m_dtheta; // [m] Differential lengths

		double m_ThetaBegin;	// [m] Beginning of azimuthal domain
		double m_ThetaEnd;		// [m] End of azimuthal domain
		double m_RadialBegin;	// [m] Beginning of radial domain
		double m_RadialEnd;		// [m] End of radial domain
		double m_AxialBegin;	// [m] Beginning of the axial domain
		double m_AxialEnd;		// [m] End of the axial domain

		double m_zDischarge;	// [m] Portion of the domain up to the screen grid
		double m_VDischarge;	// [V] Voltage of the discharge
		double m_zDistance;		// [m] Distance betweeen the screen and accel grid
		double m_zPlume;		// [m] Portion of the domain that remains in the plume region, space region
		double m_VPlume;		// [V] Voltage of the plume

		double m_wScreen;		// [m] Width of the screen grid
		double m_rScreen;		// [m] Radius of the screen grid hole
		double m_VScreen;		// [V] Voltage of the screen grid

		double m_wAccel;		// [m] Width of the acceleration grid
		double m_rAccel;		// [m] Radius of the acceleration grid hole
		double m_VAccel;		// [V] Voltage of the acceleration grid

		int m_rAccelBeginNode;	// Start node of the acceleration grid in radial direction
		int m_zAccelBeginNode;	// Start node of the acceleration grid in axial direction
		int m_zAccelEndNode;	// End node of the acceleration grid in axial direction

		int m_rScreenBeginNode; // Start node of the screen grid in radial direction
		int m_zScreenBeginNode; // Start node of the screen grid in axial direction
		int m_zScreenEndNode;	// End node of the screen grid in axial direction
	};

}