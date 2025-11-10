#include "MatrixCls.h"
#include "VectorCls.h"
#include "SubMatPtrCls.h"
#include "SubVecPtrCls.h"
#include "MatrixPkg.h"
#include <iostream>
#include <iomanip>

int main()
{
    std::cout << "=== Forward Substitution with Runtime-Sized Submatrices ===\n\n";
    
    // Example 1: Fixed-size matrix with runtime submatrix selection
    // Create a 6x6 lower triangular matrix (compile-time size)
    MatrixCls<6, 6, double> L_full;
    
    // Initialize as a lower triangular matrix
    // L[i][j] = 0 for i < j, L[i][j] = i+1 for i >= j (example values)
    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = 0; j < 6; ++j) {
            if (i >= j) {
                L_full(i, j) = static_cast<double>(i + 1);  // Lower triangular
            } else {
                L_full(i, j) = 0.0;
            }
        }
    }
    
    std::cout << "Full 6x6 matrix L:\n";
    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = 0; j < 6; ++j) {
            std::cout << std::setw(6) << L_full(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    // Example 1a: Extract a 3x3 top-left submatrix at runtime
    size_t submatrix_size = 3;  // Runtime decision
    SubMatPtrCls<double> L_3x3 = L_full.submatrix(0, 0, submatrix_size, submatrix_size);
    
    std::cout << "Extracted 3x3 submatrix (top-left):\n";
    for (size_t i = 0; i < L_3x3.rows(); ++i) {
        for (size_t j = 0; j < L_3x3.cols(); ++j) {
            std::cout << std::setw(6) << L_3x3(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    // Create a 3-element vector b
    VectorCls<3, double> b1;
    b1[0] = 1.0;
    b1[1] = 2.0;
    b1[2] = 3.0;
    
    std::cout << "Vector b1 = [";
    for (size_t i = 0; i < 3; ++i) {
        std::cout << b1[i];
        if (i < 2) std::cout << ", ";
    }
    std::cout << "]\n\n";
    
    // Solve L_3x3 * x = b1 using forward substitution
    // Using member function with VectorCls
    VectorCls<3, double> x1 = L_3x3.forwardSubstitution<3>(b1);
    
    std::cout << "Solution x1 (using member function with VectorCls):\n";
    std::cout << "x1 = [";
    for (size_t i = 0; i < 3; ++i) {
        std::cout << x1[i];
        if (i < 2) std::cout << ", ";
    }
    std::cout << "]\n\n";
    
    // Verify: L_3x3 * x1 should equal b1
    std::cout << "Verification (L_3x3 * x1):\n";
    VectorCls<3, double> verify1;
    for (size_t i = 0; i < 3; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j <= i; ++j) {
            sum += L_3x3(i, j) * x1[j];
        }
        verify1[i] = sum;
        std::cout << "Row " << i << ": " << verify1[i] << " (expected: " << b1[i] << ")\n";
    }
    std::cout << "\n";
    
    // Example 1b: Using SubVecPtrCls with the submatrix
    // Create a subvector view from a larger vector
    VectorCls<6, double> b_full;
    b_full[0] = 1.0;
    b_full[1] = 2.0;
    b_full[2] = 3.0;
    b_full[3] = 4.0;
    b_full[4] = 5.0;
    b_full[5] = 6.0;
    
    // Create a subvector view of the first 3 elements
    SubVecPtrCls<3, double> b_sub(&b_full.Elements[0]);
    
    std::cout << "Subvector b_sub (first 3 elements of b_full):\n";
    std::cout << "b_sub = [";
    for (size_t i = 0; i < 3; ++i) {
        std::cout << b_sub[i];
        if (i < 2) std::cout << ", ";
    }
    std::cout << "]\n\n";
    
    // Solve using member function with SubVecPtrCls
    VectorCls<3, double> x2 = L_3x3.forwardSubstitution<3>(b_sub);
    
    std::cout << "Solution x2 (using member function with SubVecPtrCls):\n";
    std::cout << "x2 = [";
    for (size_t i = 0; i < 3; ++i) {
        std::cout << x2[i];
        if (i < 2) std::cout << ", ";
    }
    std::cout << "]\n\n";
    
    // Example 1c: Using the generic function from MatrixPkg
    VectorCls<3, double> x3 = forwardSubstitution(L_3x3, b1);
    
    std::cout << "Solution x3 (using generic function from MatrixPkg):\n";
    std::cout << "x3 = [";
    for (size_t i = 0; i < 3; ++i) {
        std::cout << x3[i];
        if (i < 2) std::cout << ", ";
    }
    std::cout << "]\n\n";
    
    // Example 2: Dynamic submatrix size based on user input or runtime condition
    std::cout << "=== Example 2: Dynamic Submatrix Size ===\n\n";
    
    // Simulate a runtime decision (e.g., based on problem size)
    size_t problem_size = 4;  // Could come from user input, config file, etc.
    
    // Extract a problem_size x problem_size submatrix
    SubMatPtrCls<double> L_dynamic = L_full.submatrix(0, 0, problem_size, problem_size);
    
    std::cout << "Extracted " << problem_size << "x" << problem_size 
              << " submatrix (runtime size):\n";
    for (size_t i = 0; i < L_dynamic.rows(); ++i) {
        for (size_t j = 0; j < L_dynamic.cols(); ++j) {
            std::cout << std::setw(6) << L_dynamic(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    // Create corresponding vector
    VectorCls<4, double> b2;
    for (size_t i = 0; i < 4; ++i) {
        b2[i] = static_cast<double>(i + 1);
    }
    
    std::cout << "Vector b2 = [";
    for (size_t i = 0; i < 4; ++i) {
        std::cout << b2[i];
        if (i < 3) std::cout << ", ";
    }
    std::cout << "]\n\n";
    
    // Solve - note: template parameter must match the vector size
    VectorCls<4, double> x4 = L_dynamic.forwardSubstitution<4>(b2);
    
    std::cout << "Solution x4:\n";
    std::cout << "x4 = [";
    for (size_t i = 0; i < 4; ++i) {
        std::cout << x4[i];
        if (i < 3) std::cout << ", ";
    }
    std::cout << "]\n\n";
    
    // Example 3: Using the convenience operator() for submatrix extraction
    std::cout << "=== Example 3: Using operator() for Submatrix ===\n\n";
    
    // Extract submatrix from (0,0) to (2,2) inclusive using operator()
    SubMatPtrCls<double> L_op = L_full(0, 0, 2, 2);
    
    std::cout << "Submatrix extracted using operator()(0, 0, 2, 2):\n";
    for (size_t i = 0; i < L_op.rows(); ++i) {
        for (size_t j = 0; j < L_op.cols(); ++j) {
            std::cout << std::setw(6) << L_op(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    VectorCls<3, double> x5 = L_op.forwardSubstitution<3>(b1);
    
    std::cout << "Solution x5:\n";
    std::cout << "x5 = [";
    for (size_t i = 0; i < 3; ++i) {
        std::cout << x5[i];
        if (i < 2) std::cout << ", ";
    }
    std::cout << "]\n\n";
    
    // Example 4: Non-top-left submatrix
    std::cout << "=== Example 4: Non-Top-Left Submatrix ===\n\n";
    
    // Extract a 2x2 submatrix starting at (2,2)
    SubMatPtrCls<double> L_center = L_full.submatrix(2, 2, 2, 2);
    
    std::cout << "2x2 submatrix starting at (2,2):\n";
    for (size_t i = 0; i < L_center.rows(); ++i) {
        for (size_t j = 0; j < L_center.cols(); ++j) {
            std::cout << std::setw(6) << L_center(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    VectorCls<2, double> b3;
    b3[0] = 1.0;
    b3[1] = 2.0;
    
    VectorCls<2, double> x6 = L_center.forwardSubstitution<2>(b3);
    
    std::cout << "Solution x6:\n";
    std::cout << "x6 = [";
    for (size_t i = 0; i < 2; ++i) {
        std::cout << x6[i];
        if (i < 1) std::cout << ", ";
    }
    std::cout << "]\n\n";
    
    std::cout << "=== All Examples Completed Successfully! ===\n";
    
    return 0;
}

