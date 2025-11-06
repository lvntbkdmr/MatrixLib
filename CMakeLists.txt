cmake_minimum_required(VERSION 3.12)
project(MatrixLib VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard (C++14 minimum for constexpr features)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler warnings for avionics compliance
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
    # Additional warnings for safety-critical code
    add_compile_options(-Wconversion -Wsign-conversion)
endif()

# Header-only library
add_library(MatrixLib INTERFACE)
target_include_directories(MatrixLib INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Install headers
install(DIRECTORY include/ DESTINATION include)

# Example executable
add_executable(matrix_example examples/example.cpp)
target_link_libraries(matrix_example PRIVATE MatrixLib)

# Test executable
add_executable(matrix_tests tests/test_matrix.cpp)
target_link_libraries(matrix_tests PRIVATE MatrixLib)

# Optional: Add tests if Google Test is available
option(BUILD_TESTS "Build unit tests" OFF)
if(BUILD_TESTS)
    enable_testing()
    # Add test executable here when test framework is available
    add_test(NAME MatrixTests COMMAND matrix_tests)
endif()

