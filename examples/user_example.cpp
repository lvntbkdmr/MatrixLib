#include "MatrixCls.h"
#include "VectorCls.h"
#include "SubMatPtrCls.h"
#include "SubVecPtrCls.h"
#include "MatrixPkg.h"
#include <iostream>
#include <iomanip>

int main()
{
    std::cout << "=== User's Desired Usage Pattern ===\n\n";
    
    // User's exact usage pattern
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
    
    size_t size = 3;  // Runtime decision
    
    // User's exact syntax
    SubMatPtrCls<double> L_sub = L.submatrix(0, 0, size, size);
    SubVecPtrCls<double> b_sub = b.subvector(0, size);
    SubVecPtrCls<double> x_sub = x.subvector(0, size);
    
    // User's desired assignment syntax - NOW WORKS!
    x_sub = L_sub.forwardSubstitution(b_sub);
    
    std::cout << "Matrix L_sub (" << L_sub.rows() << "x" << L_sub.cols() << "):\n";
    for (size_t i = 0; i < L_sub.rows(); ++i) {
        for (size_t j = 0; j < L_sub.cols(); ++j) {
            std::cout << std::setw(6) << L_sub(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    std::cout << "Vector b_sub = [";
    for (size_t i = 0; i < b_sub.size(); ++i) {
        std::cout << b_sub[i];
        if (i < b_sub.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n\n";
    
    std::cout << "Solution x_sub = [";
    for (size_t i = 0; i < x_sub.size(); ++i) {
        std::cout << x_sub[i];
        if (i < x_sub.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n\n";
    
    // Verify: L_sub * x_sub should equal b_sub
    std::cout << "Verification (L_sub * x_sub):\n";
    bool all_match = true;
    for (size_t i = 0; i < L_sub.rows(); ++i) {
        double sum = 0.0;
        for (size_t j = 0; j <= i; ++j) {
            sum += L_sub(i, j) * x_sub[j];
        }
        bool match = std::abs(sum - b_sub[i]) < 1e-10;
        std::cout << "Row " << i << ": " << sum << " (expected: " << b_sub[i] << ") " 
                  << (match ? "✓" : "✗") << "\n";
        if (!match) all_match = false;
    }
    std::cout << "\n";
    
    if (all_match) {
        std::cout << "✓ Verification passed! The solution is correct.\n\n";
    } else {
        std::cout << "✗ Verification failed!\n\n";
    }
    
    // Test with different size
    std::cout << "=== Testing with size = 4 ===\n\n";
    size = 4;
    
    SubMatPtrCls<double> L_sub2 = L.submatrix(0, 0, size, size);
    SubVecPtrCls<double> b_sub2 = b.subvector(0, size);
    SubVecPtrCls<double> x_sub2 = x.subvector(0, size);
    
    // Same syntax works!
    x_sub2 = L_sub2.forwardSubstitution(b_sub2);
    
    std::cout << "Solution x_sub2 (size=" << size << ") = [";
    for (size_t i = 0; i < x_sub2.size(); ++i) {
        std::cout << x_sub2[i];
        if (i < x_sub2.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n\n";
    
    std::cout << "=== All Examples Completed Successfully! ===\n";
    
    return 0;
}

