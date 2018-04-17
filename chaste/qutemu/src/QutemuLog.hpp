#pragma once

#include <string>
#include <sstream>
#include <iostream>

#define COUT(msg) if (PetscTools::AmMaster()) std::cout << msg << std::endl << std::flush
#define LOG(msg) { \
    std::stringstream __ss; \
    __ss << msg; \
    QutemuLog::Log(__ss.str()); \
}

class QutemuLog
{
public:
    static void Log(const std::string &s);
    static std::string GetLog();
};