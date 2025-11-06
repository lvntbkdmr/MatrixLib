/**
 * @file test_matrix.cpp
 * @brief Unit tests for the Matrix library
 * 
 * Note: This is a basic test framework. For production use,
 * consider integrating with Google Test or similar framework.
 */

#include "matrix.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace avionics::math;

// Simple test framework
static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define ASSERT_EQ(a, b) \
    do { \
        tests_run++; \
        if ((a) == (b)) { \
            tests_passed++; \
        } else { \
            tests_failed++; \
            std::cout << "FAIL: " << __FILE__ << ":" << __LINE__ \
                      << " - Expected " << (a) << " == " << (b) << "\n"; \
        } \
    } while(0)

#define ASSERT_NEAR(a, b, epsilon) \
    do { \
        tests_run++; \
        if (std::abs((a) - (b)) < (epsilon)) { \
            tests_passed++; \
        } else { \
            tests_failed++; \
            std::cout << "FAIL: " << __FILE__ << ":" << __LINE__ \
                      << " - Expected " << (a) << " ≈ " << (b) \
                      << " (epsilon: " << (epsilon) << ")\n"; \
        } \
    } while(0)

#define ASSERT_TRUE(condition) \
    do { \
        tests_run++; \
        if (condition) { \
            tests_passed++; \
        } else { \
            tests_failed++; \
            std::cout << "FAIL: " << __FILE__ << ":" << __LINE__ \
                      << " - Condition is false\n"; \
        } \
    } while(0)

void test_construction() {
    std::cout << "Testing construction...\n";
    
    // Default construction (zero matrix)
    MatrixCls<3, 3, float> mat1;
    ASSERT_EQ(mat1(0, 0), 0.0f);
    ASSERT_EQ(mat1(2, 2), 0.0f);
    
    // Construction with initializer array
    const float init[9] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f};
    MatrixCls<3, 3, float> mat2(init);
    ASSERT_EQ(mat2(0, 0), 1.0f);
    ASSERT_EQ(mat2(0, 1), 2.0f);
    ASSERT_EQ(mat2(2, 2), 9.0f);
    
    // Copy construction
    MatrixCls<3, 3, float> mat3(mat2);
    ASSERT_EQ(mat3(0, 0), 1.0f);
    ASSERT_EQ(mat3(2, 2), 9.0f);
    
    std::cout << "  Construction tests passed\n\n";
}

void test_accessors() {
    std::cout << "Testing accessors...\n";
    
    MatrixCls<3, 3, float> mat;
    mat(0, 0) = 1.0f;
    mat(1, 1) = 2.0f;
    mat(2, 2) = 3.0f;
    
    ASSERT_EQ(mat(0, 0), 1.0f);
    ASSERT_EQ(mat(1, 1), 2.0f);
    ASSERT_EQ(mat(2, 2), 3.0f);
    
    // Test static accessors
    ASSERT_EQ(MatrixCls<3, 3, float>::rows(), 3);
    ASSERT_EQ(MatrixCls<3, 3, float>::cols(), 3);
    ASSERT_EQ(MatrixCls<3, 3, float>::size(), 9);
    
    std::cout << "  Accessor tests passed\n\n";
}

void test_bounds_checking() {
    std::cout << "Testing bounds checking...\n";
    
    MatrixCls<3, 3, float> mat;
    float value;
    
    // Valid access
    ASSERT_EQ(mat.get(1, 1, value), MatrixStatus::SUCCESS);
    ASSERT_EQ(mat.set(1, 1, 5.0f), MatrixStatus::SUCCESS);
    ASSERT_EQ(mat.get(1, 1, value), MatrixStatus::SUCCESS);
    ASSERT_EQ(value, 5.0f);
    
    // Invalid row
    ASSERT_EQ(mat.get(5, 1, value), MatrixStatus::OUT_OF_BOUNDS);
    ASSERT_EQ(mat.set(5, 1, 5.0f), MatrixStatus::OUT_OF_BOUNDS);
    
    // Invalid column
    ASSERT_EQ(mat.get(1, 5, value), MatrixStatus::OUT_OF_BOUNDS);
    ASSERT_EQ(mat.set(1, 5, 5.0f), MatrixStatus::OUT_OF_BOUNDS);
    
    std::cout << "  Bounds checking tests passed\n\n";
}

void test_identity() {
    std::cout << "Testing identity matrix...\n";
    
    MatrixCls<3, 3, float> mat;
    ASSERT_EQ(mat.identity(), MatrixStatus::SUCCESS);
    
    ASSERT_EQ(mat(0, 0), 1.0f);
    ASSERT_EQ(mat(0, 1), 0.0f);
    ASSERT_EQ(mat(1, 1), 1.0f);
    ASSERT_EQ(mat(2, 2), 1.0f);
    
    // Non-square matrix should fail
    MatrixCls<2, 3, float> non_square;
    ASSERT_EQ(non_square.identity(), MatrixStatus::INVALID_DIMENSIONS);
    
    std::cout << "  Identity tests passed\n\n";
}

void test_arithmetic_operators() {
    std::cout << "Testing arithmetic operators...\n";
    
    const float init1[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    const float init2[4] = {5.0f, 6.0f, 7.0f, 8.0f};
    
    MatrixCls<2, 2, float> mat1(init1);
    MatrixCls<2, 2, float> mat2(init2);
    
    // Addition
    auto sum = mat1 + mat2;
    ASSERT_EQ(sum(0, 0), 6.0f);
    ASSERT_EQ(sum(0, 1), 8.0f);
    ASSERT_EQ(sum(1, 0), 10.0f);
    ASSERT_EQ(sum(1, 1), 12.0f);
    
    // Subtraction
    auto diff = mat2 - mat1;
    ASSERT_EQ(diff(0, 0), 4.0f);
    ASSERT_EQ(diff(0, 1), 4.0f);
    ASSERT_EQ(diff(1, 0), 4.0f);
    ASSERT_EQ(diff(1, 1), 4.0f);
    
    // Scalar multiplication
    auto scaled = mat1 * 2.0f;
    ASSERT_EQ(scaled(0, 0), 2.0f);
    ASSERT_EQ(scaled(0, 1), 4.0f);
    
    // Compound assignment
    mat1 += mat2;
    ASSERT_EQ(mat1(0, 0), 6.0f);
    
    std::cout << "  Arithmetic operator tests passed\n\n";
}

void test_matrix_multiplication() {
    std::cout << "Testing matrix multiplication...\n";
    
    MatrixCls<2, 3, float> matA;
    matA(0, 0) = 1.0f; matA(0, 1) = 2.0f; matA(0, 2) = 3.0f;
    matA(1, 0) = 4.0f; matA(1, 1) = 5.0f; matA(1, 2) = 6.0f;
    
    MatrixCls<3, 2, float> matB;
    matB(0, 0) = 7.0f; matB(0, 1) = 8.0f;
    matB(1, 0) = 9.0f; matB(1, 1) = 10.0f;
    matB(2, 0) = 11.0f; matB(2, 1) = 12.0f;
    
    auto result = matA * matB;  // Should be 2x2
    
    ASSERT_EQ(result.rows(), 2);
    ASSERT_EQ(result.cols(), 2);
    ASSERT_NEAR(result(0, 0), 58.0f, 0.001f);   // 1*7 + 2*9 + 3*11
    ASSERT_NEAR(result(0, 1), 64.0f, 0.001f);   // 1*8 + 2*10 + 3*12
    ASSERT_NEAR(result(1, 0), 139.0f, 0.001f);  // 4*7 + 5*9 + 6*11
    ASSERT_NEAR(result(1, 1), 154.0f, 0.001f);  // 4*8 + 5*10 + 6*12
    
    std::cout << "  Matrix multiplication tests passed\n\n";
}

void test_transpose() {
    std::cout << "Testing transpose...\n";
    
    MatrixCls<2, 3, float> mat;
    mat(0, 0) = 1.0f; mat(0, 1) = 2.0f; mat(0, 2) = 3.0f;
    mat(1, 0) = 4.0f; mat(1, 1) = 5.0f; mat(1, 2) = 6.0f;
    
    auto trans = transpose(mat);
    
    ASSERT_EQ(trans.rows(), 3);
    ASSERT_EQ(trans.cols(), 2);
    ASSERT_EQ(trans(0, 0), 1.0f);
    ASSERT_EQ(trans(0, 1), 4.0f);
    ASSERT_EQ(trans(1, 0), 2.0f);
    ASSERT_EQ(trans(1, 1), 5.0f);
    ASSERT_EQ(trans(2, 0), 3.0f);
    ASSERT_EQ(trans(2, 1), 6.0f);
    
    std::cout << "  Transpose tests passed\n\n";
}

void test_determinant() {
    std::cout << "Testing determinant...\n";
    
    MatrixCls<2, 2, float> mat2x2;
    mat2x2(0, 0) = 1.0f; mat2x2(0, 1) = 2.0f;
    mat2x2(1, 0) = 3.0f; mat2x2(1, 1) = 4.0f;
    
    float det = determinant(mat2x2);
    ASSERT_NEAR(det, -2.0f, 0.001f);  // 1*4 - 2*3 = -2
    
    MatrixCls<3, 3, float> mat3x3;
    mat3x3(0, 0) = 1.0f; mat3x3(0, 1) = 2.0f; mat3x3(0, 2) = 3.0f;
    mat3x3(1, 0) = 4.0f; mat3x3(1, 1) = 5.0f; mat3x3(1, 2) = 6.0f;
    mat3x3(2, 0) = 7.0f; mat3x3(2, 1) = 8.0f; mat3x3(2, 2) = 9.0f;
    
    float det3 = determinant(mat3x3);
    ASSERT_NEAR(det3, 0.0f, 0.001f);  // This matrix is singular
    
    std::cout << "  Determinant tests passed\n\n";
}

void test_inverse() {
    std::cout << "Testing inverse...\n";
    
    MatrixCls<2, 2, float> mat;
    mat(0, 0) = 1.0f; mat(0, 1) = 2.0f;
    mat(1, 0) = 3.0f; mat(1, 1) = 4.0f;
    
    MatrixCls<2, 2, float> inv;
    ASSERT_EQ(inverse(mat, inv), MatrixStatus::SUCCESS);
    
    // Verify: mat * inv should be identity
    auto product = mat * inv;
    ASSERT_NEAR(product(0, 0), 1.0f, 0.001f);
    ASSERT_NEAR(product(0, 1), 0.0f, 0.001f);
    ASSERT_NEAR(product(1, 0), 0.0f, 0.001f);
    ASSERT_NEAR(product(1, 1), 1.0f, 0.001f);
    
    // Test singular matrix
    MatrixCls<2, 2, float> singular;
    singular(0, 0) = 1.0f; singular(0, 1) = 2.0f;
    singular(1, 0) = 2.0f; singular(1, 1) = 4.0f;  // det = 0
    
    MatrixCls<2, 2, float> inv_singular;
    ASSERT_EQ(inverse(singular, inv_singular), MatrixStatus::SINGULAR_MATRIX);
    
    std::cout << "  Inverse tests passed\n\n";
}

void test_lu_decomp() {
    std::cout << "Testing LU decomposition...\n";
    
    // Test with a known matrix
    MatrixCls<3, 3, float> mat;
    mat(0, 0) = 2.0f; mat(0, 1) = 1.0f; mat(0, 2) = 1.0f;
    mat(1, 0) = 4.0f; mat(1, 1) = 3.0f; mat(1, 2) = 3.0f;
    mat(2, 0) = 2.0f; mat(2, 1) = 3.0f; mat(2, 2) = 4.0f;
    
    MatrixCls<3, 3, float> L, U;
    ASSERT_EQ(mat.LuDecomp(L, U), MatrixStatus::SUCCESS);
    
    // Verify: L * U should equal original matrix
    auto reconstructed = L * U;
    ASSERT_NEAR(reconstructed(0, 0), mat(0, 0), 0.001f);
    ASSERT_NEAR(reconstructed(0, 1), mat(0, 1), 0.001f);
    ASSERT_NEAR(reconstructed(0, 2), mat(0, 2), 0.001f);
    ASSERT_NEAR(reconstructed(1, 0), mat(1, 0), 0.001f);
    ASSERT_NEAR(reconstructed(1, 1), mat(1, 1), 0.001f);
    ASSERT_NEAR(reconstructed(1, 2), mat(1, 2), 0.001f);
    ASSERT_NEAR(reconstructed(2, 0), mat(2, 0), 0.001f);
    ASSERT_NEAR(reconstructed(2, 1), mat(2, 1), 0.001f);
    ASSERT_NEAR(reconstructed(2, 2), mat(2, 2), 0.001f);
    
    // Test non-square matrix (should fail)
    MatrixCls<2, 3, float> non_square;
    MatrixCls<2, 3, float> L_ns, U_ns;
    ASSERT_EQ(non_square.LuDecomp(L_ns, U_ns), MatrixStatus::INVALID_DIMENSIONS);
    
    // Test singular matrix (should fail)
    MatrixCls<2, 2, float> singular;
    singular(0, 0) = 1.0f; singular(0, 1) = 2.0f;
    singular(1, 0) = 2.0f; singular(1, 1) = 4.0f;  // det = 0
    MatrixCls<2, 2, float> L_sing, U_sing;
    ASSERT_EQ(singular.LuDecomp(L_sing, U_sing), MatrixStatus::SINGULAR_MATRIX);
    
    // Test LU decomposition with size parameter
    MatrixCls<6, 6, float> mat6x6;
    mat6x6(0, 0) = 2.0f; mat6x6(0, 1) = 1.0f; mat6x6(0, 2) = 1.0f;
    mat6x6(1, 0) = 4.0f; mat6x6(1, 1) = 3.0f; mat6x6(1, 2) = 3.0f;
    mat6x6(2, 0) = 2.0f; mat6x6(2, 1) = 3.0f; mat6x6(2, 2) = 4.0f;
    
    MatrixCls<6, 6, float> L_size, U_size;
    ASSERT_EQ(mat6x6.LuDecomp(3, L_size, U_size), MatrixStatus::SUCCESS);
    
    // Verify: L * U should equal original top-left 3x3
    auto reconstructed_size = L_size * U_size;
    ASSERT_NEAR(reconstructed_size(0, 0), mat6x6(0, 0), 0.001f);
    ASSERT_NEAR(reconstructed_size(0, 1), mat6x6(0, 1), 0.001f);
    ASSERT_NEAR(reconstructed_size(0, 2), mat6x6(0, 2), 0.001f);
    ASSERT_NEAR(reconstructed_size(1, 0), mat6x6(1, 0), 0.001f);
    ASSERT_NEAR(reconstructed_size(1, 1), mat6x6(1, 1), 0.001f);
    ASSERT_NEAR(reconstructed_size(1, 2), mat6x6(1, 2), 0.001f);
    ASSERT_NEAR(reconstructed_size(2, 0), mat6x6(2, 0), 0.001f);
    ASSERT_NEAR(reconstructed_size(2, 1), mat6x6(2, 1), 0.001f);
    ASSERT_NEAR(reconstructed_size(2, 2), mat6x6(2, 2), 0.001f);
    
    std::cout << "  LU decomposition tests passed\n\n";
}

void test_backslash() {
    std::cout << "Testing backslash operation...\n";
    
    // Test with a known system: A * x = B
    MatrixCls<3, 3, float> A;
    A(0, 0) = 2.0f; A(0, 1) = 1.0f; A(0, 2) = 1.0f;
    A(1, 0) = 4.0f; A(1, 1) = 3.0f; A(1, 2) = 3.0f;
    A(2, 0) = 2.0f; A(2, 1) = 3.0f; A(2, 2) = 4.0f;
    
    MatrixCls<3, 1, float> B;
    B(0, 0) = 1.0f;
    B(1, 0) = 2.0f;
    B(2, 0) = 3.0f;
    
    MatrixCls<3, 1, float> x;
    ASSERT_EQ(backslash(A, B, x), MatrixStatus::SUCCESS);
    
    // Verify: A * x should equal B
    auto verify = A * x;
    ASSERT_NEAR(verify(0, 0), B(0, 0), 0.001f);
    ASSERT_NEAR(verify(1, 0), B(1, 0), 0.001f);
    ASSERT_NEAR(verify(2, 0), B(2, 0), 0.001f);
    
    // Test backslash with size parameter
    MatrixCls<6, 6, float> A_size, B_size;
    A_size(0, 0) = 2.0f; A_size(0, 1) = 1.0f; A_size(0, 2) = 1.0f;
    A_size(1, 0) = 4.0f; A_size(1, 1) = 3.0f; A_size(1, 2) = 3.0f;
    A_size(2, 0) = 2.0f; A_size(2, 1) = 3.0f; A_size(2, 2) = 4.0f;
    
    B_size(0, 0) = 1.0f;
    B_size(1, 0) = 2.0f;
    B_size(2, 0) = 3.0f;
    
    MatrixCls<6, 1, float> x_size;
    ASSERT_EQ(backslash(A_size, B_size, 3, x_size), MatrixStatus::SUCCESS);
    
    // Verify: A * x should equal B (only top-left 3x3 portion)
    auto verify_size = A_size * x_size;
    ASSERT_NEAR(verify_size(0, 0), B_size(0, 0), 0.001f);
    ASSERT_NEAR(verify_size(1, 0), B_size(1, 0), 0.001f);
    ASSERT_NEAR(verify_size(2, 0), B_size(2, 0), 0.001f);
    
    // Test singular matrix (should fail)
    MatrixCls<2, 2, float> A_sing;
    A_sing(0, 0) = 1.0f; A_sing(0, 1) = 2.0f;
    A_sing(1, 0) = 2.0f; A_sing(1, 1) = 4.0f;  // det = 0
    
    MatrixCls<2, 1, float> B_sing, x_sing;
    B_sing(0, 0) = 1.0f;
    B_sing(1, 0) = 2.0f;
    
    ASSERT_EQ(backslash(A_sing, B_sing, x_sing), MatrixStatus::SINGULAR_MATRIX);
    
    std::cout << "  Backslash operation tests passed\n\n";
}

void test_comparison() {
    std::cout << "Testing comparison operators...\n";
    
    const float init[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    MatrixCls<2, 2, float> mat1(init);
    MatrixCls<2, 2, float> mat2(init);
    MatrixCls<2, 2, float> mat3;
    
    ASSERT_TRUE(mat1 == mat2);
    ASSERT_TRUE(mat1 != mat3);
    
    std::cout << "  Comparison tests passed\n\n";
}

void test_type_aliases() {
    std::cout << "Testing type aliases...\n";
    
    Matrix2x2<float> mat2x2;
    Matrix3x3<float> mat3x3;
    Matrix4x4<float> mat4x4;
    Vector2<float> vec2;
    Vector3<float> vec3;
    Vector4<float> vec4;
    
    ASSERT_EQ(mat2x2.rows(), 2);
    ASSERT_EQ(mat2x2.cols(), 2);
    ASSERT_EQ(mat3x3.rows(), 3);
    ASSERT_EQ(vec3.rows(), 3);
    ASSERT_EQ(vec3.cols(), 1);
    
    std::cout << "  Type alias tests passed\n\n";
}

int main() {
    std::cout << "=== Matrix Library Unit Tests ===\n\n";
    
    test_construction();
    test_accessors();
    test_bounds_checking();
    test_identity();
    test_arithmetic_operators();
    test_matrix_multiplication();
    test_transpose();
    test_determinant();
    test_inverse();
    test_lu_decomp();
    test_backslash();
    test_comparison();
    test_type_aliases();
    
    std::cout << "=== Test Results ===\n";
    std::cout << "Tests run: " << tests_run << "\n";
    std::cout << "Tests passed: " << tests_passed << "\n";
    std::cout << "Tests failed: " << tests_failed << "\n";
    
    if (tests_failed == 0) {
        std::cout << "\nAll tests PASSED!\n";
        return 0;
    } else {
        std::cout << "\nSome tests FAILED!\n";
        return 1;
    }
}
