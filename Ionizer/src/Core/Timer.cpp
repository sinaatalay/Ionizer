#pragma once

#include <iostream>
#include <chrono>

#include "Timer.h"
#include "Log.h"

namespace Ionizer {

		Timer::Timer(){
			Reset();
		}

		void Timer::Reset(){
			m_Start = std::chrono::high_resolution_clock::now();
		}
		float Timer::Elapsed(){
			return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - m_Start).count() * 0.001f * 0.001f * 0.001f;
		}
		float Timer::ElapsedMillis(){
			return Elapsed() * 1000.0f;
		}

}

