file(READ ${Chaste_BINARY_DIR}/global/build_timestamp build_time)
file(WRITE QutemuVersion.cpp
"#include \"QutemuVersion.hpp\"

const char* QutemuVersion::GetBuildTime() {
    return \"${build_time}\";
}"
)
