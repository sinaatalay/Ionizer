#pragma once

#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>

namespace picdsmc{

	class Log {
	public:
		static void Init(int level);
		inline static std::shared_ptr<spdlog::logger>& GetLogger() { return s_Logger; }
	private:
		static std::shared_ptr<spdlog::logger> s_Logger;
	};

}

// Log macros:
#define LOG_INIT(...)			::picdsmc::Log::Init(__VA_ARGS__)

#define LOG_DEBUG(...)			::picdsmc::Log::GetLogger()->debug(__VA_ARGS__)
#define LOG_INFO(...)			::picdsmc::Log::GetLogger()->info(__VA_ARGS__)
#define LOG_WARN(...)			::picdsmc::Log::GetLogger()->warn(__VA_ARGS__)
#define LOG_ERROR(...)			::picdsmc::Log::GetLogger()->error(__VA_ARGS__)
#define LOG_CRITICAL(...)		::picdsmc::Log::GetLogger()->critical(__VA_ARGS__)


// Log level macros:
#define LOG_LEVEL_ALL 0
#define LOG_LEVEL_DEBUG 1
#define LOG_LEVEL_INFO 2
#define LOG_LEVEL_WARN 3
#define LOG_LEVEL_ERROR 4
#define LOG_LEVEL_CRITICAL 5
#define LOG_LEVEL_NONE 6