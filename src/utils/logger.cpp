#include "logger.hpp"
#include <chrono>
#include <iomanip>

namespace bonesim {
    void Logger::log(const std::string& level, const std::string& msg) {
        auto now = std::chrono::system_clock::now();
        auto time_t = std::chrono::system_clock::to_time_t(now);
        std::cout << "[" << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S")
                  << "] [" << level << "] " << msg << std::endl;
    }
    void Logger::info(const std::string& msg) { log("INFO", msg); }
    void Logger::warning(const std::string& msg) { log("WARNING", msg); }
    void Logger::error(const std::string& msg) { log("ERROR", msg); }
}