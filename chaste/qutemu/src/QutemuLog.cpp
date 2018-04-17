#include "PetscTools.hpp"
#include "QutemuLog.hpp"

std::stringstream log_stream;
void QutemuLog::Log(const std::string &s)
{
    if (PetscTools::AmMaster()) {
        log_stream << s << std::endl;
        COUT(s);
    }
}

std::string QutemuLog::GetLog() {
    return log_stream.str();
}