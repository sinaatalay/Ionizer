#pragma once
//C++ Libraries (Precompiled):
#include "Ionpch.h"

namespace Ionizer {

	class Color {
	public:
		Color() : r(0), g(0), b(0) {}
		Color(int R, int G, int B) : r(R), g(G), b(B) {}

		void SetHSV(float fH, float fS, float fV);
		void SetRGB(int R, int G, int B);
		uint32_t const GetHEX(){ return ((255 & 0xff) << 24) + ((b & 0xff) << 16) + ((g & 0xff) << 8) + (r & 0xff); };
	private:
		unsigned char r, g, b;
	};

}

