#include "Log.h"

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>

std::shared_ptr<spdlog::logger> Log::s_Logger;

void Log::Init(int level) {
	std::vector<spdlog::sink_ptr> logSinks;
	logSinks.emplace_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
	logSinks.emplace_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>("Simulation.log", true));

	// Formatting: https://github.com/gabime/spdlog/wiki/3.-Custom-formatting
	logSinks[0]->set_pattern("%^[%T] [%l]: %v%$"); 
	logSinks[1]->set_pattern("[%T] [%l]: %v");

	s_Logger = std::make_shared<spdlog::logger>("PIC-DSMC Simulation", begin(logSinks), end(logSinks));
	spdlog::register_logger(s_Logger);
	s_Logger->set_level((spdlog::level::level_enum)level);
	s_Logger->flush_on((spdlog::level::level_enum)level);
}

