cmake_minimum_required(VERSION 3.28)

project(VHLLE VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GSL REQUIRED)

# User provides path to vhlle_params repository
set(VHLLE_PARAMS_DIR "" CACHE PATH "Path to vhlle_params repository")

if(VHLLE_PARAMS_DIR)
    message(STATUS "vhlle_params found at ${VHLLE_PARAMS_DIR} " 
                   "- parameter tables will be symlinked")
else()
    message(WARNING "VHLLE_PARAMS_DIR not set. Parameter tables (eos/, ic/, params/) "
                    "will not be symlinked.\n"
                    "Set it with -DVHLLE_PARAMS_DIR=/path/to/vhlle_params")
endif()

# Subdirectories
add_subdirectory(src)

# Install the exported targets for find_package(vhlle)
install(EXPORT vhlle-targets
    FILE vhlleTargets.cmake
    NAMESPACE vhlle::
    DESTINATION lib/cmake/vhlle
)

# Package config
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/vhlleConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)
install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/vhlleConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/vhlleConfigVersion.cmake"
    DESTINATION lib/cmake/vhlle
)