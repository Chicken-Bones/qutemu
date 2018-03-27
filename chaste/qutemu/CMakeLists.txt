find_package(Chaste COMPONENTS cell_based heart)

#have to pull a bit of the chaste_do_project and chaste_do_common back into here to support cellml files
message("Configuring project qutemu")
set(component project_qutemu)

set(Chaste_${component}_SOURCES "")
file(GLOB cellml_files src/cellml/*.cellml)
foreach(cellml_file ${cellml_files})
    chaste_do_cellml(Chaste_${component}_SOURCES ${cellml_file} OFF  "--output-dir" "${CMAKE_CURRENT_BINARY_DIR}/src/cellml")
endforeach()
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src/cellml)

#build timestamps
list(APPEND Chaste_${component}_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/src/QutemuVersion.cpp)
add_custom_command(OUTPUT
        dummy
        ${CMAKE_CURRENT_BINARY_DIR}/src/QutemuVersion.cpp
        COMMAND ${CMAKE_COMMAND}
            "-DChaste_BINARY_DIR=${Chaste_BINARY_DIR}"
            -P ${CMAKE_CURRENT_SOURCE_DIR}/src/UpdateBuildTime.cmake
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/src"
        COMMENT "Generating Qutemu Build Timestamp"
        VERBATIM
        )

#fix error: ‘template<class> class std::auto_ptr’ is deprecated [-Werror=deprecated-declarations]
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(DEPRECATION_FLAG "-wr1478")
else()
    set(DEPRECATION_FLAG "-Wno-deprecated-declarations")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DEPRECATION_FLAG}")

#because ${ARGN} is not passed through the chaste_do_project macro like it is through the chaste_do_component
#chaste_do_project(qutemu ${Chaste_${component}_SOURCES})


set(Chaste_${component}_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
header_dirs(${Chaste_${component}_SOURCE_DIR} Chaste_${component}_INCLUDE_DIRS)

set(Chaste_${component}_INCLUDE_DIRS ${Chaste_${component}_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR}/src/cellml)

Chaste_DO_COMMON(${component} ${Chaste_${component}_SOURCES})
