#pragma once

#include "Walnut/Image.h"
#include <memory>
#include <glm/glm.hpp>

namespace IonizerApp {
	class Renderer {
	public:
		Renderer() = default;

		void Render();

		void OnResize(uint32_t width, uint32_t height);

		std::shared_ptr<Walnut::Image> GetFinalImage() const { return m_FinalImage; }
	private:
		std::shared_ptr<Walnut::Image> m_FinalImage;
		uint32_t* m_ImageData = nullptr;
	};
}

