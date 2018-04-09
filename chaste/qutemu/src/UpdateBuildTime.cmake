message(STATUS "updating buildtime...")
execute_process(COMMAND ${timekeeper_exe})

file(READ build_timestamp build_time)
file(WRITE QutemuVersion.cpp
"#include \"QutemuVersion.hpp\"

const char* QutemuVersion::GetBuildTime() {
    return \"${build_time}\";
}"
)
