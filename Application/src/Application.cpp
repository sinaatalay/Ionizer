#include "Ionizer.h"
#include "Walnut/Application.h"
#include "Walnut/EntryPoint.h"
#include "Walnut/Image.h"

static void NoteMarker(const char* desc)
{
	ImGui::TextDisabled("(*)");
	if (ImGui::IsItemHovered())
	{
		ImGui::BeginTooltip();
		ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
		ImGui::TextUnformatted(desc);
		ImGui::PopTextWrapPos();
		ImGui::EndTooltip();
	}
}

class IonizerLayer : public Walnut::Layer {
public:
	virtual void OnUIRender() override {
		ImGui::Begin("Menu");
		ImGui::Text("Welcome to Ionizer!");

		static float angle = 0.0f;
		ImGui::SliderFloat(" ", &angle, 0.0f, 30.0f, "Solution at theta = %.3f");
		if (ImGui::Button("Solve Poisson")) {
			Render(angle);
		}

		ImGui::Separator();

		if (ImGui::CollapsingHeader("Configuration"))
		{
			ImGui::LabelText("Label", "Value");
			ImGui::Text("Physical lengths of the nodes:");
			static double dz = 0.00002;
			ImGui::InputDouble("dz [m]", &dz, 0.000002f, 0.0005f, "%.4f");
			ImGui::SameLine;
			NoteMarker("Try to use multiples of the relevant geometry dimensions.\n"
				"Otherwise the program will change the geometry accordingly");

			static double dr = 0.00002;
			ImGui::InputDouble("dr [m]", &dr, 0.000002f, 0.0005f, "%.4f");
			ImGui::SameLine;
			NoteMarker("Try to use multiples of the relevant geometry dimensions.\n"
				"Otherwise the program will change the geometry accordingly");

			static double dtheta = 5;
			ImGui::InputDouble("dtheta [deg]", &dtheta, 0.1, 10, "%.4f");
			ImGui::SameLine;
			NoteMarker("Try to use multiples of the relevant geometry dimensions.\n"
				"Otherwise the program will change the geometry accordingly");
			dtheta = 3.141592653589793 * 2 * dtheta / 360;

			ImGui::Text("Physical lengths of the domain:");
			static double ThetaLength = 30;
			ImGui::InputDouble("theta length [deg]", &ThetaLength, 10.0f, 60.0f, "%.4f");
			ThetaLength = 3.141592653589793 * 2 * ThetaLength / 360;

			static double RadialLength = 0.002;
			ImGui::InputDouble("radial length [m]", &RadialLength, 0.001f, 0.005f, "%.4f");

			static double AxialLength = 0.002;
			ImGui::InputDouble("axial length [m]", &AxialLength, 0.001f, 0.005f, "%.4f");

			ImGui::Text("Physical lengths of the thruster:");
			static double ScreenGridWidth = 0.0004;
			ImGui::InputDouble("screen grid width [m]", &ScreenGridWidth, 0.0001f, 0.001f, "%.4f");

			static double ScreenGridHoleRadius = 0.001;
			ImGui::InputDouble("screen grid hole radius [m]", &ScreenGridHoleRadius, 0.0001f, 0.004f, "%.4f");

			static double AccelerationGridWidth = 0.0008;
			ImGui::InputDouble("acceleration grid width [m]", &AccelerationGridWidth, 0.0001f, 0.0015f, "%.4f");

			static double AccelerationGridHoleRadius = 0.0006;
			ImGui::InputDouble("acceleration grid hole radius [m]", &AccelerationGridHoleRadius, 0.0001f, 0.001f, "%.4f");

			static double DischageRegionAxialLength = 0.001;
			ImGui::InputDouble("dischage region axial length [m]", &DischageRegionAxialLength, 0.0001f, 0.002f, "%.4f");

			static double DistanceBetweenGrids = 0.0012;
			ImGui::InputDouble("distance between grids [m]", &DistanceBetweenGrids, 0.0003f, 0.002f, "%.4f");

			ImGui::Text("Potentials:");
			static double V_Discharge = 2266;
			ImGui::InputDouble("bulk plasma potential [V]", &V_Discharge, 1.0f, 4000.0f, "%.4f");

			static double V_Plume = 0;
			ImGui::InputDouble("plume plasma potential [V]", &V_Plume, -1000.0f, 1000.0f, "%.4f");

			static double V_Accel = -400;
			ImGui::InputDouble("Acceleration grid (second grid) potential [V]", &V_Accel, -1000.0f, 4000.0f, "%.4f");

			static double V_Screen = 2241;
			ImGui::InputDouble("Screen grid (first grid) potential [m]", &V_Screen, 1.0f, 4000.0f, "%.4f");

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
		}

		ImGui::End();

		//ImGui::ShowDemoWindow();

		ImGui::Begin("Viewport");

		m_ViewportWidth = ImGui::GetContentRegionAvail().x;
		m_ViewportHeight = ImGui::GetContentRegionAvail().y;

		if (m_Image) {
			ImGui::Image(m_Image->GetDescriptorSet(), { (float)m_Image->GetWidth(), (float)m_Image->GetHeight() });
		}

		ImGui::End();
	}

	void Render(float angle) {

		if (m_Image == nullptr || m_ViewportWidth != m_Image->GetWidth()) {
			m_Image = std::make_shared<Walnut::Image>(m_ViewportWidth, m_ViewportHeight, Walnut::ImageFormat::RGBA);
			delete[] m_ImageData;
			m_ImageData = new uint32_t[m_ViewportWidth * m_ViewportHeight];
		}

		Ionizer::PoissonSolver poisson(geometry);
		poisson.LogGeometry();
		poisson.SolvePoisson();

		std::vector<uint32_t> hop = poisson.GetImage(m_ViewportWidth, m_ViewportHeight, angle);

		for (uint32_t i = 0; i < m_ViewportWidth*m_ViewportHeight; i++) {
			m_ImageData[i] = hop[i];
		}

		m_Image->SetData(m_ImageData);
	}
private:
	std::shared_ptr<Walnut::Image> m_Image;
	uint32_t* m_ImageData = nullptr;
	uint32_t m_ViewportWidth;
	uint32_t m_ViewportHeight;
	Ionizer::IonThrusterGeometry geometry;
};

Walnut::Application* Walnut::CreateApplication(int argc, char** argv) {
	Walnut::ApplicationSpecification spec;
	spec.Name = "Ionizer";

	LOG_INIT(LOG_LEVEL_ALL);
	LOG_INFO("Welcome to Ionizer!");

	Walnut::Application* app = new Walnut::Application(spec);
	app->PushLayer<IonizerLayer>();
	app->SetMenubarCallback([app]() {
		if (ImGui::BeginMenu("File")) {
			if (ImGui::MenuItem("Exit")) {
				app->Close();
			}
			ImGui::EndMenu();
		}
	});
	return app;
}