/**
 * @file submatrix_example.cpp
 * @brief Example of runtime submatrix extraction
 * 
 * This example demonstrates how to work with submatrices when the size
 * is determined at runtime.
 */

#include "matrix.hpp"
#include <iostream>
#include <iomanip>

using namespace avionics::math;

void printSubmatrix(const MatrixCls<6, 6, float>& mat, int size, const char* name) {
    std::cout << name << " (top-left " << size << "x" << size << "):\n";
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(2) << mat(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

int main() {
    std::cout << "=== Runtime Submatrix Extraction Example ===\n\n";

    // Create a 6x6 matrix
    MatrixCls<6, 6, float> mat;
    
    // Fill with some values
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            mat(i, j) = static_cast<float>(i * 6 + j + 1);
        }
    }

    std::cout << "Original 6x6 matrix:\n";
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            std::cout << std::setw(4) << mat(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    // Runtime variable x (can be 1-6)
    int x = 3;  // This can change at runtime

    // Method 1: Extract into a fixed-size result matrix (RECOMMENDED)
    // Use the maximum size (6x6) and work with the top-left x×x portion
    MatrixCls<6, 6, float> submat;
    if (mat.extractTopLeft(x, submat) == MatrixStatus::SUCCESS) {
        std::cout << "Method 1: Extracted " << x << "x" << x << " submatrix:\n";
        printSubmatrix(submat, x, "Submatrix");
        
        // Now you can do operations on submat
        // Only use elements (0:x-1, 0:x-1)
        auto submat_transpose = transpose(submat);
        std::cout << "Transpose of submatrix:\n";
        printSubmatrix(submat_transpose, x, "Transposed");
        
        // Matrix multiplication example (only using top-left x×x)
        MatrixCls<6, 6, float> identity;
        identity.identity();
        auto product = submat * identity;
        std::cout << "Submatrix * Identity:\n";
        printSubmatrix(product, x, "Product");
    }

    // Method 2: Extract into a specific size matrix (if you know the size at compile time)
    // This is more memory efficient but requires compile-time knowledge
    if (x == 3) {
        MatrixCls<3, 3, float> submat3x3;
        if (mat.extractTopLeft(3, submat3x3) == MatrixStatus::SUCCESS) {
            std::cout << "Method 2: Extracted into 3x3 matrix:\n";
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    std::cout << std::setw(8) << std::fixed << std::setprecision(2) 
                              << submat3x3(i, j) << " ";
                }
                std::cout << "\n";
            }
            std::cout << "\n";
            
            // Now you can do operations on the 3x3 matrix
            auto submat3x3_transpose = transpose(submat3x3);
            std::cout << "Transpose of 3x3 submatrix:\n";
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    std::cout << std::setw(8) << std::fixed << std::setprecision(2) 
                              << submat3x3_transpose(i, j) << " ";
                }
                std::cout << "\n";
            }
            std::cout << "\n";
        }
    }

    // Method 3: Work directly with the original matrix (no extraction needed)
    // Just use mat(0:x-1, 0:x-1) in your operations
    std::cout << "Method 3: Working directly with original matrix (x=" << x << "):\n";
    std::cout << "Top-left " << x << "x" << x << " elements:\n";
    for (int i = 0; i < x; ++i) {
        for (int j = 0; j < x; ++j) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(2) 
                      << mat(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    // Example: LU decomposition on submatrix
    if (x == 3) {
        MatrixCls<3, 3, float> submat3x3;
        mat.extractTopLeft(3, submat3x3);
        
        MatrixCls<3, 3, float> L, U;
        if (submat3x3.LuDecomp(L, U) == MatrixStatus::SUCCESS) {
            std::cout << "LU Decomposition of 3x3 submatrix:\n";
            std::cout << "L matrix:\n";
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    std::cout << std::setw(8) << std::fixed << std::setprecision(2) 
                              << L(i, j) << " ";
                }
                std::cout << "\n";
            }
            std::cout << "\nU matrix:\n";
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    std::cout << std::setw(8) << std::fixed << std::setprecision(2) 
                              << U(i, j) << " ";
                }
                std::cout << "\n";
            }
        }
    }

    // Example: Backslash operation with runtime size
    std::cout << "Backslash operation with runtime size (x=" << x << "):\n";
    MatrixCls<6, 6, float> A;
    
    // Fill A with some values
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            A(i, j) = static_cast<float>((i == j) ? (i + 1) * 2 : (i + j) * 0.5f);
        }
    }
    
    // Create B as a vector
    MatrixCls<6, 1, float> B;
    for (int i = 0; i < 6; ++i) {
        B(i, 0) = static_cast<float>(i + 1);
    }
    
    // Solve: A * result = B (only using top-left x×x portion of A, top-left x×1 of B)
    // Note: backslash() automatically extracts the submatrices internally, so no need
    // to call extractTopLeft() beforehand!
    MatrixCls<6, 1, float> result;
    if (backslash(A, B, x, result) == MatrixStatus::SUCCESS) {
        std::cout << "Solution (top-left " << x << "x" << x << " portion):\n";
        for (int i = 0; i < x; ++i) {
            std::cout << "result(" << i << ") = " << result(i, 0) << "\n";
        }
        
        // Verify: A * result should equal B (only top-left x×x portion)
        auto verify = A * result;
        std::cout << "\nVerification (A * result, top-left " << x << " elements):\n";
        for (int i = 0; i < x; ++i) {
            std::cout << "B(" << i << ") = " << verify(i, 0) << " (expected: " << B(i, 0) << ")\n";
        }
    } else {
        std::cout << "Backslash operation failed\n";
    }
    std::cout << "\n";

    std::cout << "\n=== Example Complete ===\n";
    return 0;
}

