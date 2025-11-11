/**
 * @file app.cpp
 * @brief Test application for the Matrix Library LeftDivide functionality
 * 
 * This application tests the LeftDivide method (solving Ax = b) for various
 * matrix sizes and types. It verifies the correctness of the solution by
 * computing the residual ||Ax - b||.
 * 
 * The tests include:
 * - Small matrices (2x2, 3x3, 4x4) with known solutions
 * - Larger matrices (5x5, 6x6) with tridiagonal structure
 * - Verification using residual computation
 * 
 * @note This is a test/demo application. For production use, integrate
 * the matrix library into your application code.
 */

#include "MatrixPkg.h"
#include <iostream>
#include <iomanip>
#include <cmath>

/**
 * @brief Test LeftDivide function for a specific matrix size
 * @param size Matrix size (size x size)
 * @param testName Name of the test for output purposes
 * @return true if test passed (residual < 1e-10), false otherwise
 * 
 * This function:
 * 1. Creates a test matrix A and vector b
 * 2. Solves Ax = b using LeftDivide
 * 3. Computes the residual ||Ax - b||
 * 4. Verifies the solution is correct
 * 
 * For small matrices (size <= 4), it also compares against expected solutions.
 */
bool testLeftDivide(UINT32 size, const char* testName)
{
    const UINT32 MAX_N = 6;  // Maximum matrix size for stack allocation
    
    // Create matrices and vectors with maximum size MAX_N
    // The actual problem size will be 'size' (set via setSubMatrixSize)
    MatrixCls<MAX_N, MAX_N, REAL64> A;
    VectorCls<MAX_N, REAL64> b;
    VectorCls<MAX_N, REAL64> x;
    VectorCls<MAX_N, REAL64> x_expected;  // Expected solution (for verification)

    // Set the runtime dimensions to the test size
    // This allows using only a portion of the allocated matrix
    A.setSubMatrixSize(size, size);

    // ========================================================================
    // Test Case Setup: Well-conditioned matrices with known solutions
    // ========================================================================
    // For size 2: A = [4, 1; 2, 3], b = [5; 11], solution: x = [0.4; 3.4]
    // For size 3: A = [2, 1, 0; 1, 2, 1; 0, 1, 2], b = [3; 4; 3], solution: x = [1; 1; 1]
    // For size 4: A = [4, 1, 0, 0; 1, 4, 1, 0; 0, 1, 4, 1; 0, 0, 1, 4], b = [5; 6; 6; 5], solution: x = [1; 1; 1; 1]
    // For size 5, 6: Tridiagonal matrices with diagonal dominance
    
    if (size == 2)
    {
        // A = [4, 1; 2, 3], b = [5; 11]
        // Solution: x = [0.4; 3.4]
        // Verification: 4*0.4 + 1*3.4 = 1.6 + 3.4 = 5.0 ✓
        //               2*0.4 + 3*3.4 = 0.8 + 10.2 = 11.0 ✓
        A(0, 0) = 4.0; A(0, 1) = 1.0;
        A(1, 0) = 2.0; A(1, 1) = 3.0;
        b(0) = 5.0;
        b(1) = 11.0;
        x_expected(0) = 0.4;
        x_expected(1) = 3.4;
    }
    else if (size == 3)
    {
        // A = [2, 1, 0; 1, 2, 1; 0, 1, 2], b = [3; 4; 3]
        A(0, 0) = 2.0; A(0, 1) = 1.0; A(0, 2) = 0.0;
        A(1, 0) = 1.0; A(1, 1) = 2.0; A(1, 2) = 1.0;
        A(2, 0) = 0.0; A(2, 1) = 1.0; A(2, 2) = 2.0;
        b(0) = 3.0;
        b(1) = 4.0;
        b(2) = 3.0;
        x_expected(0) = 1.0;
        x_expected(1) = 1.0;
        x_expected(2) = 1.0;
    }
    else if (size == 4)
    {
        // A = [4, 1, 0, 0; 1, 4, 1, 0; 0, 1, 4, 1; 0, 0, 1, 4]
        // b = [5; 6; 6; 5]
        A(0, 0) = 4.0; A(0, 1) = 1.0; A(0, 2) = 0.0; A(0, 3) = 0.0;
        A(1, 0) = 1.0; A(1, 1) = 4.0; A(1, 2) = 1.0; A(1, 3) = 0.0;
        A(2, 0) = 0.0; A(2, 1) = 1.0; A(2, 2) = 4.0; A(2, 3) = 1.0;
        A(3, 0) = 0.0; A(3, 1) = 0.0; A(3, 2) = 1.0; A(3, 3) = 4.0;
        b(0) = 5.0;
        b(1) = 6.0;
        b(2) = 6.0;
        b(3) = 5.0;
        x_expected(0) = 1.0;
        x_expected(1) = 1.0;
        x_expected(2) = 1.0;
        x_expected(3) = 1.0;
    }
    else if (size == 5)
    {
        // A = Identity-like with small perturbations, b = [1; 2; 3; 4; 5]
        for (UINT32 i = 0; i < size; i++)
        {
            for (UINT32 j = 0; j < size; j++)
            {
                if (i == j)
                    A(i, j) = 5.0 + i * 0.1;  // Diagonal dominance
                else if (i == j + 1 || i == j - 1)
                    A(i, j) = 0.5;  // Tridiagonal
                else
                    A(i, j) = 0.0;
            }
            b(i) = static_cast<REAL64>(i + 1);
        }
        // Expected solution computed manually or left as verification
        // For this test, we'll just verify it solves correctly
    }
    else if (size == 6)
    {
        // A = Tridiagonal matrix with diagonal dominance
        for (UINT32 i = 0; i < size; i++)
        {
            for (UINT32 j = 0; j < size; j++)
            {
                if (i == j)
                    A(i, j) = 6.0;  // Strong diagonal
                else if (i == j + 1 || i == j - 1)
                    A(i, j) = 1.0;  // Tridiagonal
                else
                    A(i, j) = 0.0;
            }
            b(i) = static_cast<REAL64>(i + 1) * 2.0;
        }
    }

    // ========================================================================
    // Solve the linear system Ax = b
    // ========================================================================
    
    // Create a submatrix view for verification (to compute residual)
    SubMatPtrCls<MAX_N, MAX_N, REAL64> A_sub = A(0, 0, size, size);

    // Copy submatrix to temporary for LeftDivide
    // Note: In production code, you could use the submatrix view directly
    // if LeftDivide supported SubMatPtrCls (currently requires MatrixCls)
    MatrixCls<MAX_N, MAX_N, REAL64> A_temp;
    A_temp.setSubMatrixSize(size, size);
    for (UINT32 i = 0; i < size; i++)
    {
        for (UINT32 j = 0; j < size; j++)
        {
            A_temp(i, j) = A_sub(i, j);
        }
    }

    // Create vectors for the sub-problem (copy from full vectors)
    // Set the runtime length to match the problem size
    VectorCls<MAX_N, REAL64> b_sub, x_sub;
    b_sub.setLength(size);
    x_sub.setLength(size);
    for (UINT32 i = 0; i < size; i++)
    {
        b_sub(i) = b(i);
        x_sub(i) = x(i);
    }

    // Solve using LeftDivide (instance method)
    // This automatically selects the best algorithm based on matrix properties
    bool success = A_temp.LeftDivide(b_sub, x_sub);

    // Copy result back to full vector
    for (UINT32 i = 0; i < size; i++)
    {
        x(i) = x_sub(i);
    }

    if (!success)
    {
        std::cout << testName << " FAILED: LeftDivide returned false (singular matrix)" << std::endl;
        return false;
    }

    // Print results
    std::cout << "\n" << testName << " (Size " << size << "x" << size << "):" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Solution x:" << std::endl;
    for (UINT32 i = 0; i < size; i++)
    {
        std::cout << "  x[" << i << "] = " << std::setw(10) << x(i);
        if (size <= 4 && x_expected.getLength() > i)
        {
            REAL64 error = std::abs(x(i) - x_expected(i));
            std::cout << " (expected: " << std::setw(10) << x_expected(i) 
                      << ", error: " << std::setw(10) << error << ")";
        }
        std::cout << std::endl;
    }

    // ========================================================================
    // Verify solution by computing residual: ||Ax - b||
    // ========================================================================
    // A good solution should have a very small residual (close to zero)
    REAL64 maxResidual = 0.0;
    for (UINT32 i = 0; i < size; i++)
    {
        // Compute Ax for row i
        REAL64 residual = 0.0;
        for (UINT32 j = 0; j < size; j++)
        {
            residual += A_sub(i, j) * x_sub(j);  // Use x_sub, not x
        }
        // Compute (Ax - b) for row i
        residual -= b_sub(i);
        // Track the maximum absolute residual
        REAL64 absResidual = std::abs(residual);
        if (absResidual > maxResidual)
            maxResidual = absResidual;
    }
    
    std::cout << "  Max residual ||Ax - b||: " << maxResidual << std::endl;
    
    // Consider test passed if residual is very small (machine precision)
    // This indicates the solution is numerically accurate
    bool passed = (maxResidual < 1e-10);
    if (passed)
        std::cout << "  ✓ PASSED" << std::endl;
    else
        std::cout << "  ✗ FAILED (residual too large)" << std::endl;
    
    return passed;
}

/**
 * @brief Main function - runs all LeftDivide tests
 * @return 0 if all tests pass, 1 otherwise
 * 
 * This function runs a series of tests for the LeftDivide method
 * with different matrix sizes (2x2 through 6x6). Each test:
 * - Creates a test matrix and vector
 * - Solves the linear system
 * - Verifies the solution accuracy
 */
int main()
{
    std::cout << "=== Testing LeftDivide Function ===" << std::endl;
    std::cout << "This test suite verifies the LeftDivide method (solving Ax = b)" << std::endl;
    std::cout << "for various matrix sizes and types." << std::endl << std::endl;
    
    bool allPassed = true;
    
    // Test different matrix sizes
    // Each test uses a well-conditioned matrix with known properties
    allPassed &= testLeftDivide(2, "Test 1: 2x2 matrix");
    allPassed &= testLeftDivide(3, "Test 2: 3x3 tridiagonal matrix");
    allPassed &= testLeftDivide(4, "Test 3: 4x4 tridiagonal matrix");
    allPassed &= testLeftDivide(5, "Test 4: 5x5 tridiagonal matrix");
    allPassed &= testLeftDivide(6, "Test 5: 6x6 tridiagonal matrix");
    
    // Print summary
    std::cout << "\n=== Summary ===" << std::endl;
    if (allPassed)
    {
        std::cout << "All tests PASSED! ✓" << std::endl;
        std::cout << "The LeftDivide function is working correctly." << std::endl;
        return 0;
    }
    else
    {
        std::cout << "Some tests FAILED! ✗" << std::endl;
        std::cout << "Please review the output above for details." << std::endl;
        return 1;
    }
}