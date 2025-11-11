cmake_minimum_required(VERSION 3.10)
project(MatrixLib VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Generate compile_commands.json for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Include directories for the library
include_directories(include)

# Library target (header-only library, but we can create an interface target)
add_library(matrixlib INTERFACE)
target_include_directories(matrixlib INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Example executable
add_executable(matrixlib_example examples/app.cpp)
target_link_libraries(matrixlib_example PRIVATE matrixlib)

# For Windows compatibility (if needed)
if(WIN32)
    target_compile_definitions(matrixlib INTERFACE _CRT_SECURE_NO_WARNINGS)
    target_compile_definitions(matrixlib_example PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

# Installation (optional)
install(TARGETS matrixlib EXPORT MatrixLibTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(DIRECTORY include/ DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)

