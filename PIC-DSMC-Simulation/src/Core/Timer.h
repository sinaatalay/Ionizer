#pragma once

#include <iostream>
#include <chrono>

namespace picdsmc {

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

