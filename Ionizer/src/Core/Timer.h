#pragma once

//C++ Libraries (Precompiled):
#include "Ionpch.h"

namespace Ionizer {

	class Timer {
	public:
		Timer();

		void Reset();
		float Elapsed();
		float ElapsedMillis();
	private:
		std::chrono::time_point<std::chrono::high_resolution_clock> m_Start;
	};

}

