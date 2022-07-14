#include "Ionizer.h"
#include "Walnut/Application.h"
#include "Walnut/EntryPoint.h"
#include "Walnut/Image.h"
#include "Walnut/Timer.h"

#include "Renderer.h"

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

		if (ImGui::CollapsingHeader("Configuration")) {
            if (ImGui::TreeNode("Step sizes")) {
				ImGui::LabelText("Label", "Value");
				static double dz = 0.00002;
				ImGui::InputDouble("dz", &dz, 0.00002, 0.0001, "%.5f");
				static double dr = 0.00002;
				ImGui::InputDouble("dr", &dr, 0.00002, 0.0001, "%.5f");
				static double dtheta = 5;
				ImGui::InputDouble("dtheta [deg]", &dtheta, 1, 5, "%.5f");
                ImGui::TreePop();
            }

			if (ImGui::TreeNode("Lengths of the domain")){
				ImGui::Text("Try to use multiples of the step sizes. Otherwise the program will change the geometry accordingly.");
				ImGui::LabelText("Label", "Value");
				static double ThetaLength = 30;
				ImGui::InputDouble("theta length [deg]", &ThetaLength, 1, 5, "%.4f");

				static double RadialLength = 0.002;
				ImGui::InputDouble("radial length", &RadialLength, 0.002, 0.01, "%.4f");

				static double AxialLength = 0.002;
				ImGui::InputDouble("axial length", &AxialLength, 0.002, 0.01, "%.4f");
				ImGui::TreePop();
			}

			if (ImGui::TreeNode("Lengths of the thruster")){
				ImGui::LabelText("Label", "Value");
				static double ScreenGridWidth = 0.0004;
				ImGui::InputDouble("screen grid width [m]", &ScreenGridWidth, 0.0001, 0.002, "%.4f");

				static double ScreenGridHoleRadius = 0.001;
				ImGui::InputDouble("screen grid hole radius [m]", &ScreenGridHoleRadius, 0.001, 0.05, "%.4f");

				static double AccelerationGridWidth = 0.0008;
				ImGui::InputDouble("acceleration grid width [m]", &AccelerationGridWidth, 0.001, 0.004, "%.4f");

				static double AccelerationGridHoleRadius = 0.0006;
				ImGui::InputDouble("acceleration grid hole radius [m]", &AccelerationGridHoleRadius, 0.0001, 0.004, "%.4f");

				static double DischargeRegionAxialLength = 0.001;
				ImGui::InputDouble("discharge region axial length [m]", &DischargeRegionAxialLength, 0.001, 0.005, "%.4f");

				static double DistanceBetweenGrids = 0.0012;
				ImGui::InputDouble("distance between grids [m]", &DistanceBetweenGrids, 0.0001, 0.0005, "%.4f");
				ImGui::TreePop();
			}

			if (ImGui::TreeNode("Potentials")) {
				ImGui::LabelText("Label", "Value");
				static double V_Discharge = 2266;
				ImGui::InputDouble("bulk plasma potential [V]", &V_Discharge, 10, 50, "%.4f");

				static double V_Plume = 0;
				ImGui::InputDouble("plume plasma potential [V]", &V_Plume, 10, 50, "%.4f");

				static double V_Accel = -400;
				ImGui::InputDouble("Acceleration grid (second grid) potential [V]", &V_Accel, 10, 50, "%.4f");

				static double V_Screen = 2241;
				ImGui::InputDouble("Screen grid (first grid) potential [m]", &V_Screen, 10, 50, "%.4f");
				ImGui::TreePop();
			}
		}

		ImGui::End();

		// ImGui::ShowDemoWindow();

		ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
		ImGui::Begin("Viewport");

		m_ViewportWidth = ImGui::GetContentRegionAvail().x;
		m_ViewportHeight = ImGui::GetContentRegionAvail().y;


		auto image = m_Renderer.GetFinalImage();
		if (image) {
			ImGui::Image(image->GetDescriptorSet(), { (float)image->GetWidth(), (float)image->GetHeight() });
		}

		ImGui::End();
		ImGui::PopStyleVar();
	}

	void Render(float angle) {

		Walnut::Timer timer;

		m_Renderer.OnResize(m_ViewportWidth, m_ViewportHeight);
		m_Renderer.Render();

		m_LastRenderTime = timer.ElapsedMillis();

	}
private:
	IonizerApp::Renderer m_Renderer;
	uint32_t m_ViewportWidth=0, m_ViewportHeight = 0;

	float m_LastRenderTime = 0.0f;
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