/**
 * @file example.cpp
 * @brief Example usage of the Matrix library
 */

#include "matrix.hpp"
#include <iostream>
#include <iomanip>

using namespace avionics::math;

void printMatrix(const MatrixCls<3, 3, float>& mat, const char* name) {
    std::cout << name << ":\n";
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(2) << mat(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

int main() {
    std::cout << "=== Matrix Library Example ===\n\n";

    // Create a 3x3 identity matrix
    MatrixCls<3, 3, float> mat1;
    mat1.identity();
    printMatrix(mat1, "Identity matrix (3x3)");

    // Create another matrix with initial values
    const float init_data[9] = {1.0f, 2.0f, 3.0f,
                                 4.0f, 5.0f, 6.0f,
                                 7.0f, 8.0f, 9.0f};
    MatrixCls<3, 3, float> mat2(init_data);
    printMatrix(mat2, "Matrix 2 (3x3)");

    // Matrix operations using function overloading
    auto mat3 = mat1 + mat2;  // Addition
    printMatrix(mat3, "Matrix 1 + Matrix 2");

    auto mat4 = mat2 * 2.0f;  // Scalar multiplication
    printMatrix(mat4, "Matrix 2 * 2.0");

    auto mat5 = transpose(mat2);  // Transpose
    std::cout << "Transpose of Matrix 2:\n";
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(2) << mat5(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    // Access elements
    float value = mat2(1, 2);  // Gets element at row 1, col 2
    std::cout << "Element at (1, 2): " << value << "\n\n";

    // Safe access with bounds checking
    float safe_value;
    if (mat2.get(1, 2, safe_value) == MatrixStatus::SUCCESS) {
        std::cout << "Safe access - Element at (1, 2): " << safe_value << "\n";
    }
    std::cout << "\n";

    // Test bounds checking
    if (mat2.get(5, 2, safe_value) == MatrixStatus::OUT_OF_BOUNDS) {
        std::cout << "Bounds checking works: OUT_OF_BOUNDS returned for invalid index\n";
    }
    std::cout << "\n";

    // Matrix multiplication example
    MatrixCls<2, 3, float> matA;
    matA(0, 0) = 1.0f; matA(0, 1) = 2.0f; matA(0, 2) = 3.0f;
    matA(1, 0) = 4.0f; matA(1, 1) = 5.0f; matA(1, 2) = 6.0f;

    MatrixCls<3, 2, float> matB;
    matB(0, 0) = 7.0f; matB(0, 1) = 8.0f;
    matB(1, 0) = 9.0f; matB(1, 1) = 10.0f;
    matB(2, 0) = 11.0f; matB(2, 1) = 12.0f;

    auto matC = matA * matB;  // Result is 2x2
    std::cout << "Matrix multiplication (2x3 * 3x2 = 2x2):\n";
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(2) << matC(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    // Determinant example (2x2)
    MatrixCls<2, 2, float> mat2x2;
    mat2x2(0, 0) = 1.0f; mat2x2(0, 1) = 2.0f;
    mat2x2(1, 0) = 3.0f; mat2x2(1, 1) = 4.0f;
    
    float det = determinant(mat2x2);
    std::cout << "Determinant of 2x2 matrix: " << det << "\n\n";

    // Inverse example (2x2)
    MatrixCls<2, 2, float> inv_result;
    if (inverse(mat2x2, inv_result) == MatrixStatus::SUCCESS) {
        std::cout << "Inverse of 2x2 matrix:\n";
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                std::cout << std::setw(8) << std::fixed << std::setprecision(2) << inv_result(i, j) << " ";
            }
            std::cout << "\n";
        }
    }
    std::cout << "\n";

    // LU Decomposition example
    MatrixCls<3, 3, float> mat_lu;
    mat_lu(0, 0) = 2.0f; mat_lu(0, 1) = 1.0f; mat_lu(0, 2) = 1.0f;
    mat_lu(1, 0) = 4.0f; mat_lu(1, 1) = 3.0f; mat_lu(1, 2) = 3.0f;
    mat_lu(2, 0) = 2.0f; mat_lu(2, 1) = 3.0f; mat_lu(2, 2) = 4.0f;

    MatrixCls<3, 3, float> L, U;
    if (mat_lu.LuDecomp(L, U) == MatrixStatus::SUCCESS) {
        std::cout << "LU Decomposition:\n";
        std::cout << "Original matrix:\n";
        printMatrix(mat_lu, "");
        
        std::cout << "Lower triangular matrix (L):\n";
        printMatrix(L, "");
        
        std::cout << "Upper triangular matrix (U):\n";
        printMatrix(U, "");
        
        // Verify: L * U should equal original matrix
        auto reconstructed = L * U;
        std::cout << "L * U (should equal original):\n";
        printMatrix(reconstructed, "");
    } else {
        std::cout << "LU Decomposition failed\n\n";
    }

    // Backslash operation example (A \ B)
    std::cout << "Backslash operation (A \\ B) example:\n";
    MatrixCls<3, 3, float> A;
    A(0, 0) = 2.0f; A(0, 1) = 1.0f; A(0, 2) = 1.0f;
    A(1, 0) = 4.0f; A(1, 1) = 3.0f; A(1, 2) = 3.0f;
    A(2, 0) = 2.0f; A(2, 1) = 3.0f; A(2, 2) = 4.0f;

    MatrixCls<3, 1, float> B;
    B(0, 0) = 1.0f;
    B(1, 0) = 2.0f;
    B(2, 0) = 3.0f;

    // Solve A * x = B
    MatrixCls<3, 1, float> x;
    if (backslash(A, B, x) == MatrixStatus::SUCCESS) {
        std::cout << "Solution x to A*x = B:\n";
        for (int i = 0; i < 3; ++i) {
            std::cout << "x(" << i << ") = " << x(i, 0) << "\n";
        }
        
        // Verify: A * x should equal B
        auto verify = A * x;
        std::cout << "\nVerification (A * x):\n";
        for (int i = 0; i < 3; ++i) {
            std::cout << "B(" << i << ") = " << verify(i, 0) << " (expected: " << B(i, 0) << ")\n";
        }
    }
    std::cout << "\n";

    // Type aliases example
    Matrix3x3<float> mat3x3;
    mat3x3.identity();
    Vector3<float> vec3;
    vec3(0, 0) = 1.0f;
    vec3(1, 0) = 2.0f;
    vec3(2, 0) = 3.0f;

    std::cout << "Using type aliases - Vector3:\n";
    for (int i = 0; i < 3; ++i) {
        std::cout << vec3(i, 0) << "\n";
    }

    std::cout << "\n=== Example Complete ===\n";
    return 0;
}
