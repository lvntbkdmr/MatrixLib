#include "MatrixCls.h"
#include "VectorCls.h"
#include "SubMatPtrCls.h"
#include "SubVecPtrCls.h"
#include "MatrixPkg.h"
#include <iostream>
#include <iomanip>

int main()
{
    std::cout << "=== Forward Substitution with Runtime-Sized Subvectors ===\n\n";
    
    // Create compile-time sized matrix and vectors
    MatrixCls<6, 6, double> L;
    VectorCls<6, double> b;
    VectorCls<6, double> x;
    
    // Initialize L as a lower triangular matrix
    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = 0; j < 6; ++j) {
            if (i >= j) {
                L(i, j) = static_cast<double>(i + 1);  // Lower triangular
            } else {
                L(i, j) = 0.0;
            }
        }
    }
    
    // Initialize b
    for (size_t i = 0; i < 6; ++i) {
        b[i] = static_cast<double>(i + 1);
    }
    
    std::cout << "Full 6x6 matrix L:\n";
    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = 0; j < 6; ++j) {
            std::cout << std::setw(6) << L(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    std::cout << "Vector b = [";
    for (size_t i = 0; i < 6; ++i) {
        std::cout << b[i];
        if (i < 5) std::cout << ", ";
    }
    std::cout << "]\n\n";
    
    // Runtime decision: use a 3x3 submatrix
    size_t size = 3;  // Runtime decision
    
    // Extract runtime-sized submatrix and subvectors
    SubMatPtrCls<double> L_sub = L.submatrix(0, 0, size, size);
    SubVecPtrCls<double> b_sub = b.subvector(0, size);
    SubVecPtrCls<double> x_sub = x.subvector(0, size);
    
    std::cout << "Extracted " << size << "x" << size << " submatrix:\n";
    for (size_t i = 0; i < L_sub.rows(); ++i) {
        for (size_t j = 0; j < L_sub.cols(); ++j) {
            std::cout << std::setw(6) << L_sub(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    std::cout << "Subvector b_sub = [";
    for (size_t i = 0; i < b_sub.size(); ++i) {
        std::cout << b_sub[i];
        if (i < b_sub.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n\n";
    
    // Method 1: Write directly to x_sub (recommended for runtime-sized vectors)
    std::cout << "Method 1: Writing directly to x_sub\n";
    L_sub.forwardSubstitution(b_sub, x_sub);
    
    std::cout << "Solution x_sub = [";
    for (size_t i = 0; i < x_sub.size(); ++i) {
        std::cout << x_sub[i];
        if (i < x_sub.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n\n";
    
    // Method 2: Return VectorCls and assign (requires compile-time size)
    std::cout << "Method 2: Returning VectorCls and assigning\n";
    VectorCls<3, double> x_result = L_sub.forwardSubstitution<3>(b_sub);
    x_sub = x_result;  // Assignment works because sizes match at runtime
    
    std::cout << "Solution x_sub (after assignment) = [";
    for (size_t i = 0; i < x_sub.size(); ++i) {
        std::cout << x_sub[i];
        if (i < x_sub.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n\n";
    
    // Verify: L_sub * x_sub should equal b_sub
    std::cout << "Verification (L_sub * x_sub):\n";
    for (size_t i = 0; i < L_sub.rows(); ++i) {
        double sum = 0.0;
        for (size_t j = 0; j <= i; ++j) {
            sum += L_sub(i, j) * x_sub[j];
        }
        std::cout << "Row " << i << ": " << sum << " (expected: " << b_sub[i] << ")\n";
    }
    std::cout << "\n";
    
    // Example with different size
    std::cout << "=== Example with size = 4 ===\n\n";
    size = 4;
    
    SubMatPtrCls<double> L_sub2 = L.submatrix(0, 0, size, size);
    SubVecPtrCls<double> b_sub2 = b.subvector(0, size);
    SubVecPtrCls<double> x_sub2 = x.subvector(0, size);
    
    L_sub2.forwardSubstitution(b_sub2, x_sub2);
    
    std::cout << "Solution x_sub2 (size=" << size << ") = [";
    for (size_t i = 0; i < x_sub2.size(); ++i) {
        std::cout << x_sub2[i];
        if (i < x_sub2.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n\n";
    
    std::cout << "=== All Examples Completed Successfully! ===\n";
    
    return 0;
}

