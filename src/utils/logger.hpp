#pragma once
#include <string>
#include <iostream>

namespace bonesim {
    class Logger {
    public:
        static void info(const std::string& msg);
        static void warning(const std::string& msg);
        static void error(const std::string& msg);
    private:
        static void log(const std::string& level, const std::string& msg);
    };
}