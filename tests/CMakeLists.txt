# Test configuration for MatrixLib

# Fetch Google Test framework
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.14.0
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Enable Google Test integration with CTest
include(GoogleTest)

# VectorCls unit tests
add_executable(test_VectorCls test_VectorCls.cpp)
target_link_libraries(test_VectorCls PRIVATE matrixlib gtest gtest_main)
target_include_directories(test_VectorCls PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Discover tests for CTest
gtest_discover_tests(test_VectorCls)

# Add test to CTest
add_test(NAME VectorClsTests COMMAND test_VectorCls)
