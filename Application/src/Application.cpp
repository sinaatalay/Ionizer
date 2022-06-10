#include "Ionizer.h"
#include "Walnut/Application.h"
#include "Walnut/EntryPoint.h"
#include "Walnut/Image.h"

class ExampleLayer : public Walnut::Layer {
public:
	virtual void OnUIRender() override {
		ImGui::Begin("Settings");
		if (ImGui::Button("Solve Poisson")) {
			Render();
		}
		ImGui::End();

		// ImGui::ShowDemoWindow();
		ImGui::Begin("Viewport");

		m_ViewportWidth = ImGui::GetContentRegionAvail().x;
		m_ViewportHeight = ImGui::GetContentRegionAvail().y;

		if (m_Image) {
			ImGui::Image(m_Image->GetDescriptorSet(), { (float)m_Image->GetWidth(), (float)m_Image->GetHeight() });
		}

		ImGui::End();
	}

	void Render() {

		if (m_Image == nullptr || m_ViewportWidth != m_Image->GetWidth() || m_ViewportHeight != m_Image->GetHeight()) {
			m_Image = std::make_shared<Walnut::Image>(m_ViewportWidth, m_ViewportHeight, Walnut::ImageFormat::RGBA);
			delete[] m_ImageData;
			m_ImageData = new uint32_t[m_ViewportWidth * m_ViewportHeight];
		}

#if 1
		Ionizer::PoissonSolver poisson;
		poisson.SolvePoisson();

		std::vector<uint32_t> hop = poisson.GetImage(m_ViewportWidth, m_ViewportHeight);

		for (uint32_t i = 0; i < m_ViewportWidth * m_ViewportHeight; i++) {
			m_ImageData[i] = hop[i];
		}
#else
		for (uint32_t i = 0; i < m_ViewportWidth * m_ViewportHeight; i++) {
			if(i<20*m_ViewportWidth){
				m_ImageData[i] = 0xff0000ff;
			}
			else {
				m_ImageData[i] = 0xff000000;
			}
		}
#endif
		m_Image->SetData(m_ImageData);
	}
private:
	std::shared_ptr<Walnut::Image> m_Image;
	uint32_t* m_ImageData = nullptr;
	uint32_t m_ViewportWidth;
	uint32_t m_ViewportHeight;
};

Walnut::Application* Walnut::CreateApplication(int argc, char** argv) {
	Walnut::ApplicationSpecification spec;
	spec.Name = "Ionizer";

	LOG_INIT(LOG_LEVEL_ALL);
	LOG_INFO("Welcome to PIC-DSMC Simulation!");

	Walnut::Application* app = new Walnut::Application(spec);
	app->PushLayer<ExampleLayer>();
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